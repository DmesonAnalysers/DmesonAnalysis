'''
python script to check the number of candidates in the skimmed data for ML
run: python CheckNumberOfCandidates.py config_training_FileName.yml
'''
import argparse
import yaml
import pandas as pd

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
bkg = pd.read_parquet(inputCfg['input']['data']).query(inputCfg['data_prep']['filt_bkg_mass'])
prompt = pd.read_parquet(inputCfg['input']['prompt'])
if inputCfg['input']['FD']:
    FD = pd.read_parquet(inputCfg['input']['FD'])

# loop over training pt bins
for ptMin, ptMax in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max']):
    print(f'\nPt bin {ptMin}-{ptMax} GeV/c, available candidates:')
    numBkg = len(bkg.query(f'{ptMin} < pt_cand < {ptMax}'))
    numPrompt = len(prompt.query(f'{ptMin} < pt_cand < {ptMax}'))
    print(f'  - bkg -> {numBkg}\n  - prompt -> {numPrompt}')
    if inputCfg['input']['FD']:
        numFD = len(FD.query(f'{ptMin} < pt_cand < {ptMax}'))
        print(f'  - non-prompt -> {numFD}\n')
