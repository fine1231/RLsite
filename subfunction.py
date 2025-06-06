import numpy as np
import re
import torch
import warnings
from torch.autograd import Variable

def create_variable(tensor):
    if torch.cuda.is_available():
        return Variable(tensor.cuda())
    else:
        return Variable(tensor)

def readLinesStrip(lines):
    for i in range(len(lines)):
        lines[i] = lines[i].rstrip('\n').replace('\t',' ').replace('        ',' ').replace('       ',' ').replace('      ',' ').replace('     ',' ').replace('   ',' ').replace('  ',' ')
    return lines

def getNucleic(path, featureMapName, featureMap = True):
    nucleics = open(path+"/"+featureMapName).readlines()
    nucleics = readLinesStrip(nucleics)
    Name = featureMapName 
    if(featureMap):
        featureMap = []
        for i in range(1,len(nucleics)):
            featureMap.append(nucleics[i])
        return Name, featureMap

def read_index_file(index_file_path):
    merge_instructions = []
    with open(index_file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            target_index = int(parts[0]) - 1 
            merge_indexes = [int(x) for x in parts[1:]] 
            merge_instructions.append((target_index, merge_indexes))
    return merge_instructions

def getFeatureDict(featurePath, featureDictPath, num_prev, num_next):
    featureDict = open(featureDictPath).readlines()
    engFeatureDict = {}
    engFeatureDict_SN = {}
    engFeatureDict_SS = {}
    
    for data in featureDict:
        nucleicName = data.strip()
        eachFeaturePath = "./data/dataPre/STRUC/" + nucleicName
        indexFilePath = "./data/dataPre/SN/" + nucleicName + "_sn"
        eachMap = open(eachFeaturePath, 'r').readlines()
        eachFeatureMap = [x.rstrip('\n').strip() for x in eachMap]
        merge_instructions = read_index_file(indexFilePath)
        eachFeature2D = {} 
        for i, data in enumerate(eachFeatureMap):
            featureMapName = data.strip().split(' ')[0]
            Name, featureMap = getNucleic(featurePath, featureMapName)
            featuremap_np = [list(map(float, x.strip(' ').split(' '))) for x in featureMap]
            feature2D = np.expand_dims(featuremap_np, axis=0)
            feature2D = torch.FloatTensor(feature2D)
            eachFeature2D[Name] = feature2D
        for target_index, merge_indexes in merge_instructions:
            if not (0 <= target_index < len(eachFeatureMap)):
                raise IndexError(f"Target index {target_index + 1} is out of range for feature map length {len(eachFeatureMap)}.")
            target_featureMapName = eachFeatureMap[target_index].strip().split(' ')[0]
            target_feature2D = eachFeature2D.get(target_featureMapName, torch.zeros((1, feature2D.shape[1], feature2D.shape[2])))
            merge_feature2Ds = [target_feature2D]
            for idx in merge_indexes:
                if not (0 <= idx < len(eachFeatureMap)):
                    merge_feature2D = torch.zeros_like(target_feature2D)
                else:
                    merge_featureMapName = eachFeatureMap[idx].strip().split(' ')[0]
                    merge_feature2D = eachFeature2D.get(merge_featureMapName, torch.zeros_like(target_feature2D))               
                merge_feature2Ds.append(merge_feature2D)
            combined_feature2D = torch.cat(merge_feature2Ds, dim=0)
            engFeatureDict_SN[target_featureMapName] = combined_feature2D[:, :10, :]
        for i, data in enumerate(eachFeatureMap):
            featureMapName = data.strip().split(' ')[0]
            feature2D = eachFeature2D[featureMapName]
            prev_feature2Ds = []
            for j in range(num_prev):
                prev_idx = num_prev-1-j
                prev_feature_idx = i - prev_idx - 1
                if prev_feature_idx >= 0:
                    prev_featureMapName = eachFeatureMap[prev_feature_idx].strip().split(' ')[0] 
                    prev_feature2D = eachFeature2D.get(prev_featureMapName, None)
                    if prev_feature2D is None:
                        prev_feature2D = torch.zeros_like(feature2D)
                else:
                    prev_feature2D = torch.zeros_like(feature2D)
                prev_feature2Ds.append(prev_feature2D)
            next_feature2Ds = []
            for j in range(num_next):
                next_idx = j
                next_feature_idx = i + next_idx + 1
                if next_feature_idx < len(eachFeatureMap):
                    next_featureMapName = eachFeatureMap[next_feature_idx].strip().split(' ')[0]
                    next_feature2D = eachFeature2D.get(next_featureMapName, None)
                    if next_feature2D is None:
                        next_feature2D = torch.zeros_like(feature2D)
                else:
                    next_feature2D = torch.zeros_like(feature2D)
                next_feature2Ds.append(next_feature2D)
            combined_feature2D = torch.cat(prev_feature2Ds + [feature2D] + next_feature2Ds, dim=0)
            engFeatureDict_SS[featureMapName] = combined_feature2D[:, -13:, :]
        for i, data in enumerate(eachFeatureMap):
            featureMapName = data.strip().split(' ')[0]
            data_SN = engFeatureDict_SN[featureMapName]
            data_SS = engFeatureDict_SS[featureMapName]
            engFeatureDict[featureMapName] = np.concatenate((data_SN, data_SS),axis=1)  
    return engFeatureDict

def getTestNucleicList(testFoldPath):
    with open(testFoldPath, 'r') as f:
        testNucleicList_cpi = f.read().strip().split('\n') 
    testNucleicList = [cpi.strip().split()[0] for cpi in testNucleicList_cpi]
    return testNucleicList

def getDataDict(testNucleicList):
    dataDict = {}
    for x in testNucleicList:
        xData = []
        nucleicName = x.strip(' ')
        nucleicPath = "./data/dataPre/ACDE/"+nucleicName
        nucleicTest = open(nucleicPath,'r').readlines()
        nucleic = [x.strip(' ').split(' ')[0] for x in nucleicTest]
        Label = [x.rstrip('\n').strip(' ').split('   ')[1] for x in nucleicTest]
        for i in range(len(nucleic)):
            xData.append([nucleic[i], Label[i]])
        dataDict[nucleicName] = xData
    return dataDict

def getValDataDict(valNucleicList):
    valDataDict = {}
    for x in valNucleicList:
        xData = []
        nucleicName = x.strip(' ')
        nucleicPath = "./data/dataPre/ValDB/"+nucleicName
        nucleicVal = open(nucleicPath, 'r').readlines()
        nucleic = [x.strip(' ').split(' ')[0] for x in nucleicVal]
        vLabel = [x.rstrip('\n').strip(' ').split('   ')[1] for x in nucleicVal]
        for i in range(len(nucleic)):
            xData.append([nucleic[i], vLabel[i]])
        valDataDict[nucleicName] = xData
    return valDataDict
