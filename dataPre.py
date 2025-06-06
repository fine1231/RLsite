from model import *
from subfunction import *
from torch.utils.data import Dataset, DataLoader
import numpy as np

class ProDataset(Dataset):
    def __init__(self, dataSet, engFeatureDict):
        self.dataSet = dataSet 
        self.len = len(dataSet)
        self.dict = engFeatureDict 
        self.properties = [int (x[1]) for x in dataSet]
        self.properties_list = list(sorted(set(self.properties)))
    def __getitem__(self, index):
        Name,label = self.dataSet[index]
        featureMap = self.dict[Name]
        return Name, featureMap, int(label)
    def __len__(self):
        return self.len
    def get_properties(self):
        return self.property_list
    def get_property(self, id):
        return self.property_list[id]
    def get_property_id(self, property):
        return self.property_list.index(property)

print('load path datas...')
testFoldPath = './data/dataPre/Test/nativeTE18'
valListPath = './data/dataPre/Test/TE9' 
JLListPath = './data/dataPre/Test/JL10' 
TLListPath = './data/dataPre/Test/TL12' 
featurePath = './data/engFeature'
featureDictPath = './data/dataPre/FeatureList/list_ph'

print('get feature datas...')
engFeatureDict = getFeatureDict(featurePath, featureDictPath, num_prev=5, num_next=5)

print('get TE18_nucleic dic...')
testNucleicList = getTestNucleicList(testFoldPath)
dataDict = getDataDict(testNucleicList)

print('get TE9_nucleic dic...')
valList = getTestNucleicList(valListPath)
valDataDict = getValDataDict(valList)

print('get JL10_nucleic dic...')
JLList = getTestNucleicList(JLListPath)
JLDataDict = getValDataDict(JLList)

print('get TL12_nucleic dic...')
TLList = getTestNucleicList(TLListPath)
TLDataDict = getValDataDict(TLList) 

print('model args...')
modelArgs = {}
modelArgs['batch_size'] = 1
modelArgs['d_a'] = 32
modelArgs['r'] = 10
modelArgs['d'] = 11
modelArgs['in_channels'] = 8
modelArgs['cnn_channels'] = 32
modelArgs['cnn_layers'] = 4 
modelArgs['lstm_hid_dim'] = 32 
modelArgs['pkl'] = './pkl/model.pkl'
modelArgs['dense_hid'] = 64
modelArgs['task_type'] = 0
modelArgs['n_classes'] = 1
modelArgs['emb_dim'] = 13 
modelArgs['dropout'] = 0.2
modelArgs['model'] =  RLsite(modelArgs,block = ResidualBlock).cuda()
modelArgs['engFeatureDict'] = engFeatureDict
modelArgs['doTE18'] = False
modelArgs['doRB9'] = False
modelArgs['doJL10'] = True
modelArgs['doTL12'] = False
modelArgs['test_nucleics'] = testNucleicList
modelArgs['testDataDict'] = dataDict
modelArgs['val_nucleics'] = valList
modelArgs['valDataDict'] = valDataDict
modelArgs['JL_nucleics'] = JLList
modelArgs['JLDataDict'] = JLDataDict
modelArgs['TL_nucleics'] = TLList
modelArgs['TLDataDict'] = TLDataDict
