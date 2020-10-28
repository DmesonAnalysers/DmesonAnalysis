'''
python script to convert ML models from .pkl to .model
'''

import os
import sys
import argparse

from hipe4ml.model_handler import ModelHandler

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFilePkl', metavar='text', default='model.pkl', help='input pickle file to be converted')
args = parser.parse_args()

ModelPath = os.path.expanduser(args.inFilePkl)
print(f'Loaded saved model: {ModelPath}')
ModelHandl = ModelHandler()
ModelHandl.load_model_handler(ModelPath)

if '.pickle' in ModelPath:
    outFileName = ModelPath.replace('.pickle', '.model')
elif '.pkl' in ModelPath:
    outFileName = ModelPath.replace('.pkl', '.model')
else:
    print(f'ERROR: invalid input file {ModelHandl}, please check it! Exit')
    sys.exit()

ModelHandl.dump_original_model(outFileName, True)
print(f'Saved model: {outFileName}')
