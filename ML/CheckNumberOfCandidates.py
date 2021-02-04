'''
python script to check the number of candidates in the skimmed data for ML
run: python CheckNumberOfCandidates.py config_training_FileName.yml
'''
import sys
import argparse
import yaml
import pandas as pd

sys.path.append('..')
from utils.DfUtils import LoadDfFromRootOrParquet

# inputs
parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='config_training_FileName.yml',
                    help='config file used for the training')
args = parser.parse_args()

# load configfiles
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

# load dataframes
print('Reading input files')
bkg = LoadDfFromRootOrParquet(
    inputCfg['input']['data'], inTreeNames=inputCfg['input']['treename'])
bkg = bkg.query(inputCfg['data_prep']['filt_bkg_mass'])
prompt = LoadDfFromRootOrParquet(
    inputCfg['input']['prompt'], inTreeNames=inputCfg['input']['treename'])
if inputCfg['input']['FD']:
    FD = LoadDfFromRootOrParquet(
        inputCfg['input']['FD'], inTreeNames=inputCfg['input']['treename'])

# loop over training pt bins
for ptMin, ptMax, bkg_mult in zip(inputCfg['pt_ranges']['min'],
                                  inputCfg['pt_ranges']['max'],
                                  inputCfg['data_prep']['bkg_mult']):
    print(f'\nPt bin {ptMin}-{ptMax} GeV/c, available candidates:')
    numBkg = len(bkg.query(f'{ptMin} < pt_cand < {ptMax}'))
    numPrompt = len(prompt.query(f'{ptMin} < pt_cand < {ptMax}'))
    print(f'  - bkg -> {numBkg}\n  - prompt -> {numPrompt}')
    numFD = 0
    if inputCfg['input']['FD']:
        numFD = len(FD.query(f'{ptMin} < pt_cand < {ptMax}'))
        print(f'  - non-prompt -> {numFD}')
    print(f'  - fraction of bkg used -> {100 * (numPrompt + numFD) * bkg_mult / numBkg :.2f}%\n')
