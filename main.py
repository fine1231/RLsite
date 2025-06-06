import torch
from sklearn import metrics
from predict import *
import warnings
warnings.filterwarnings("ignore")
torch.cuda.set_device(0)
print('cuda size == 1')

def main():
    """
    Parsing command line parameters, reading data, fitting and scoring model.
    """

    if(modelArgs['doTE18']):
        testArgs = {
            'pkl_file_path': modelArgs['pkl'],
            'model': modelArgs['model'], 
            'name': 'TE18 size',
            'test_nucleics': modelArgs['test_nucleics'],
            'testDataDict': modelArgs['testDataDict'],
            'engFeatureDict': modelArgs['engFeatureDict'],
        }
        result = testPerNucleic(testArgs)

    if(modelArgs['doRB9']):
        testArgs = {
            'pkl_file_path': modelArgs['pkl'],
            'model': modelArgs['model'], 
            'name': 'RB9 size',
            'test_nucleics': modelArgs['val_nucleics'],
            'testDataDict': modelArgs['valDataDict'],
            'engFeatureDict': modelArgs['engFeatureDict'],
        }
        result = testPerNucleic(testArgs)

    if(modelArgs['doJL10']):
        testArgs = {
            'pkl_file_path': modelArgs['pkl'],
            'model': modelArgs['model'], 
            'name': 'JL10 size',
            'test_nucleics': modelArgs['JL_nucleics'],
            'testDataDict': modelArgs['JLDataDict'],
            'engFeatureDict': modelArgs['engFeatureDict'],
        }
        result = testPerNucleic(testArgs)

    if(modelArgs['doTL12']):
        testArgs = {
            'pkl_file_path': modelArgs['pkl'],
            'model': modelArgs['model'],
            'name': 'TL12 size',
            'test_nucleics': modelArgs['TL_nucleics'],
            'testDataDict': modelArgs['TLDataDict'],
            'engFeatureDict': modelArgs['engFeatureDict'],
        }

        result = testPerNucleic(testArgs)

if __name__ == "__main__":
    main()

