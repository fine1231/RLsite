import torch
from subfunction import *
from dataPre import *
from sklearn import metrics
import numpy as np
import torch
from torch.autograd import Variable

def load_model(pkl_path, attention_model):
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    state_dict = torch.load(pkl_path, map_location=device)
    new_state_dict = {}
    for key, value in state_dict.items():
        new_key = key.replace("module.", "") 
        new_state_dict[new_key] = value
    attention_model.load_state_dict(new_state_dict, strict=False)
    attention_model.to(device)
    attention_model.eval()
    return attention_model

def testPerNucleic(testArgs):
    result = {}
    pkl_file_path = testArgs['pkl_file_path'] 
    attention_model = modelArgs['model']
    attention_model = load_model(pkl_file_path, attention_model)
    for x in testArgs['test_nucleics']:
        print('\nCurrent test nucleic: ', x.strip())
        x = x.strip()
        data = testArgs['testDataDict'][x]
        engFeatureDict1 = testArgs['engFeatureDict']
        test_dataset = ProDataset(dataSet=data, engFeatureDict=testArgs['engFeatureDict'])
        test_loader = DataLoader(dataset=test_dataset, batch_size=1, shuffle=False, drop_last=True)
        testArgs['test_loader'] = test_loader
        testArgs['model'] = attention_model
        testRecall, testPrecision, testAuc, all_pred, all_target = test(testArgs)
        result[x] = [testRecall, testPrecision, testAuc, all_pred, all_target]
    return result

def test(testArgs):
    test_loader = testArgs['test_loader']
    attention_model = testArgs['model']
    n_batches = 0
    correct = 0
    threshold = 0.505
    all_pred = np.array([])
    all_target = np.array([])

    with torch.no_grad():
        for batch_idx, (Name, featureMap, properties) in enumerate(test_loader):
            y = properties
            attention_model.hidden_state = attention_model.init_hidden(batch_size=1)
            featureMap = featureMap.cuda()
            featureSeq = featureMap[:, :, -13:, 0]
            featureSeq_filit = featureSeq[:, :, :]
            featureSeq = featureSeq_filit
            featureEng = featureMap[:, :, :10, :]
            y_pred, att = attention_model(featureEng, featureSeq)
            round_y_pred = y_pred.round(decimals=3)
#            print("Name =", Name, "y =", properties, "y_pred =", round_y_pred)
            if not bool(attention_model.type):
                pred = torch.round(y_pred.type(torch.DoubleTensor).squeeze(1))
                correct += torch.eq(torch.round(y_pred.type(torch.DoubleTensor).squeeze(1)), y.type(torch.DoubleTensor)).data.sum()
                all_pred = np.concatenate((all_pred, y_pred.data.cpu().squeeze(1).numpy()), axis=0)
                all_target = np.concatenate((all_target, y.data.cpu().numpy()), axis=0)
    testSize = round(len(test_loader.dataset), 3)
    testRecall = round(metrics.recall_score(all_target, (all_pred > threshold).astype(int)), 3)
    testPrecision = round(metrics.precision_score(all_target, (all_pred > threshold).astype(int)), 3)
    testAuc = round(metrics.roc_auc_score(all_target, all_pred), 3)
    testMcc = round(metrics.matthews_corrcoef(all_target, (all_pred > threshold).astype(int)), 3)
    print(testArgs['name'],"=",testSize,"  test recall =",testRecall,"  test precision =",testPrecision,"  test auc =",testAuc,"  test Mcc = ",testMcc)
    return testRecall,testPrecision,testAuc,all_pred,all_target
