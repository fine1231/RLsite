import torch.nn as nn
import torch
import torch.nn.functional as F
from torch.autograd import Variable

def conv3x3(in_channels, out_channels, stride=1):
    return nn.Conv3d(in_channels, out_channels, kernel_size=3, stride=stride, padding=1, bias=True)
def conv1x1(in_channels, out_channels, stride=1):
    return nn.Conv3d(in_channels, out_channels, kernel_size=1, stride=stride, padding=0, bias=True)

class ResidualBlock(nn.Module):
    def __init__(self, in_channels, out_channels, stride=1, downsample=None):
        super(ResidualBlock, self).__init__()
        self.conv1 = conv3x3(in_channels, out_channels, stride)
        self.bn1 = nn.InstanceNorm3d(out_channels, affine=True)
        self.elu = nn.ReLU(inplace=True)
        self.conv2 = conv3x3(out_channels, out_channels)
        self.bn2 = nn.InstanceNorm3d(out_channels, affine=True)
        self.downsample = downsample

    def forward(self, x):
        residual = x
        out = self.conv1(x)
        out = self.bn1(out)
        out = self.elu(out)
        out = self.conv2(out)
        out = self.bn2(out)

        if self.downsample is not None:
            residual = self.downsample(x)
        out += residual
        out = self.elu(out)
        return out

def init_weights(m):
    if isinstance(m, nn.Linear):
        nn.init.kaiming_uniform_(m.weight, nonlinearity='relu')
        if m.bias is not None:
            nn.init.zeros_(m.bias)
    elif isinstance(m, nn.Conv3d):
        nn.init.kaiming_normal_(m.weight, mode='fan_out', nonlinearity='relu')
        if m.bias is not None:
            nn.init.zeros_(m.bias)
    elif isinstance(m, nn.LSTM):
        for name, param in m.named_parameters():
            if 'weight' in name:
                nn.init.kaiming_normal_(param, nonlinearity='relu')
            elif 'bias' in name:
                nn.init.zeros_(param)

class RLsite(nn.Module):
    def __init__(self, args, block):
        super(RLsite, self).__init__()
        self.batch_size = args['batch_size']
        self.in_channels = args['in_channels']
        self.type = args['task_type']
        self.r = args['r']
        self.d = args['d']
        self.lstm_hid_dim = args['lstm_hid_dim']
        
        # CNN Layers
        self.conv = conv3x3(1, self.in_channels)
        self.bn = nn.InstanceNorm3d(self.in_channels, affine=True)
        self.elu = nn.ReLU(inplace=False)
        self.layer1 = self.make_layer(block, args['cnn_channels'], args['cnn_layers'])
        self.layer2 = self.make_layer(block, args['cnn_channels'], args['cnn_layers'])
        self.linear_first_seq = nn.Linear(args['cnn_channels'], args['d_a'])
        self.linear_second_seq = nn.Linear(args['d_a'], self.r)
        self.linear_final_step = nn.Linear(self.lstm_hid_dim*2+args['d_a'], args['dense_hid'])
        self.linear_final = nn.Linear(args['dense_hid'], args['n_classes'])
        self.layer_norm_2d = nn.LayerNorm(normalized_shape=(10, 18))

        # LSTM Layers
        self.lstm = torch.nn.LSTM(args['emb_dim'],args['lstm_hid_dim'],2,batch_first=True,bidirectional=True,dropout=args['dropout'])
        self.linear_first = nn.Linear(2*self.lstm_hid_dim,args['d_a'])
        self.linear_second = nn.Linear(args['d_a'],args['r'])
        self.layer_norm = torch.nn.LayerNorm(args['emb_dim'])
        self.apply(init_weights)

    def softmax(self, input, axis=1):
        input_size = input.size()
        trans_input = input.transpose(axis, len(input_size)-1)
        trans_size = trans_input.size()
        input_2d = trans_input.contiguous().view(-1, trans_size[-1])
        soft_max_2d = F.softmax(input_2d, dim=-1)
        soft_max_nd = soft_max_2d.view(*trans_size)
        return soft_max_nd.transpose(axis, len(input_size)-1)

    def init_hidden(self, batch_size):
        return (Variable(torch.zeros(4, batch_size, self.lstm_hid_dim).cuda()),
                Variable(torch.zeros(4, batch_size, self.lstm_hid_dim).cuda()))

    def make_layer(self, block, out_channels, blocks, stride=1):
        downsample = None
        if stride != 1 or self.in_channels != out_channels:
            downsample = nn.Sequential(
                conv1x1(self.in_channels, out_channels, stride=stride),
                nn.InstanceNorm3d(out_channels, affine=True))
        layers = [block(self.in_channels, out_channels, stride, downsample)]
        self.in_channels = out_channels
        layers.extend([block(out_channels, out_channels) for _ in range(1, blocks)])
        return nn.Sequential(*layers)

    def forward(self, x1, x2):
        x1 = x1.unsqueeze(1)
        mean = x1.mean(dim=(0, 2, 3, 4), keepdim=True)
        std = x1.std(dim=(0, 2, 3, 4), keepdim=True)
        x_normalized = (x1 - mean) / (std + 1e-7)
        pic = self.conv(x_normalized)
        pic = self.bn(pic)
        pic = self.elu(pic)
        pic = self.layer1(pic)
        pic = self.layer2(pic)
        pic_emb = torch.mean(pic, 3).permute(0, 2, 3, 1)
        pic_emb = torch.min(pic_emb, dim=1)[0]
        pic_emb = pic_emb.squeeze(1) 
        seq_att = F.tanh(self.linear_first_seq(pic_emb))
        seq_att = self.linear_second_seq(seq_att)
        seq_att = self.softmax(seq_att, 1).transpose(1, 2)
        seq_embed = torch.matmul(seq_att, pic_emb)
        avg_seq_embed = torch.sum(seq_embed, 1) / self.r


        x2 = x2.float().cuda()
        x_normalized = self.layer_norm(x2)
        batch_size = x2.size(0)
        self.hidden_state = self.init_hidden(batch_size)
        outputs, self.hidden_state = self.lstm(x_normalized, self.hidden_state)
        sentence_att = F.tanh(self.linear_first(outputs))
        sentence_att = self.linear_second(sentence_att)
        sentence_att = self.softmax(sentence_att, 1).transpose(1, 2)

        sentence_embed = torch.matmul(sentence_att, outputs)
        avg_sentence_embed = torch.sum(sentence_embed, 1) / self.r
        sscomplex = torch.cat([avg_sentence_embed, avg_seq_embed], dim=1)
        sscomplex = F.relu(self.linear_final_step(sscomplex))
        output = torch.sigmoid(self.linear_final(sscomplex))
        return output, seq_att  # not 0
