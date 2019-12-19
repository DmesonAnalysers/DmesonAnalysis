'''
python script to filter tree from task output and save output trees in root (default) or parquet files for ML studies
run: python FilterTree4ML.py cfgFileName.yml
with option --parquet the output files are saved into parquet files instead of root ones
'''

import sys
import argparse
import yaml
import uproot
sys.path.append('..')
from utils.DfUtils import WriteTree, FilterBitDf #pylint: disable=wrong-import-position,import-error,no-name-in-module


bitSignal = 0
bitBkg = 1
bitPrompt = 2
bitFD = 3
bitRefl = 4

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='cfgFileName.yml',
                    help='input config yaml file name')
parser.add_argument('--parquet', default=False, action='store_true',
                    help='flag to save output files into parquet files')

args = parser.parse_args()

with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inFileNames = cfg['infile']['filename']
if not isinstance(inFileNames, list):
    inFileNames = [inFileNames]
inDirName = cfg['infile']['dirname']
inTreeName = cfg['infile']['treename']
isMC = cfg['infile']['isMC']
outTreeName = cfg['outfile']['treename']
outSuffix = cfg['outfile']['suffix']
outDirName = cfg['outfile']['dirpath']

preSelections = cfg['skimming']['preselections']
colsToKeep = cfg['skimming']['colstokeep']

if colsToKeep and 'inv_mass' not in colsToKeep:
    print('Warning: invariant mass branch (inv_mass) disabled. Are you sure you don\'t want to keep it?')
if colsToKeep and 'pt_cand' not in colsToKeep:
    print('Warning: pt branch (pt_cand) disabled. Are you sure you don\'t want to keep it?')

PtMin = cfg['skimming']['pt']['min']
PtMax = cfg['skimming']['pt']['max']

for iFile, inFile in enumerate(inFileNames):
    inTree = uproot.open(inFile)['{0}/{1}'.format(inDirName, inTreeName)]
    if iFile == 0:
        dataFrame = inTree.pandas.df()
    else:
        dataFrame.append(inTree.pandas.df())

if not colsToKeep:
    colsToKeep = list(dataFrame.columns)
    colsToKeep.remove('cand_type')

print('Applying pT selection')
dataFramePtCut = dataFrame.query('pt_cand > {0} & pt_cand < {1}'.format(PtMin, PtMax))
del dataFrame
dataFramePtCutSel = dataFramePtCut.query(preSelections)
del dataFramePtCut

if isMC:
    print('Getting bkg dataframe')
    dataFramePtCutSelBkg = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitBkg])
    print('Getting prompt dataframe')
    dataFramePtCutSelPrompt = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSignal, bitPrompt], 'and')
    print('Getting FD dataframe')
    dataFramePtCutSelFD = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSignal, bitFD], 'and')
    print('Getting reflected prompt dataframe')
    dataFramePtCutSelPromptRefl = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSignal, bitPrompt, bitRefl], 'and')
    print('Getting reflected signal dataframe')
    dataFramePtCutSelFDRefl = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSignal, bitFD, bitRefl], 'and')
    del dataFramePtCutSel

    if not args.parquet:
        if not dataFramePtCutSelBkg.empty:
            print('Saving bkg tree')
            WriteTree(dataFramePtCutSelBkg, colsToKeep, outTreeName, \
                '{0}/Bkg{1}_pT_{2:.0f}_{3:.0f}.root'.format(outDirName, outSuffix, PtMin, PtMax))
        if not dataFramePtCutSelPrompt.empty:
            print('Saving prompt tree')
            WriteTree(dataFramePtCutSelPrompt, colsToKeep, outTreeName, \
                '{0}/Prompt{1}_pT_{2:.0f}_{3:.0f}.root'.format(outDirName, outSuffix, PtMin, PtMax))
        if not dataFramePtCutSelFD.empty:
            print('Saving FD tree')
            WriteTree(dataFramePtCutSelFD, colsToKeep, outTreeName, \
                '{0}/FD{1}_pT_{2:.0f}_{3:.0f}.root'.format(outDirName, outSuffix, PtMin, PtMax))
        if not dataFramePtCutSelPromptRefl.empty:
            print('Saving prompt refl tree')
            WriteTree(dataFramePtCutSelPromptRefl, colsToKeep, outTreeName, \
                '{0}/PromptRefl{1}_pT_{2:.0f}_{3:.0f}.root'.format(outDirName, outSuffix, PtMin, PtMax))
        if not dataFramePtCutSelFDRefl.empty:
            print('Saving FD refl tree')
            WriteTree(dataFramePtCutSelFDRefl, colsToKeep, outTreeName, \
                '{0}/FDRefl{1}_pT_{2:.0f}_{3:.0f}.root'.format(outDirName, outSuffix, PtMin, PtMax))
    else:
        if not dataFramePtCutSelBkg.empty:
            print('Saving bkg parquet')
            dataFramePtCutSelBkg[colsToKeep].to_parquet('{0}/Bkg{1}_pT_{2:.0f}_{3:.0f}.parquet.gzip'.format(\
                outDirName, outSuffix, PtMin, PtMax), compression='gzip')
        if not dataFramePtCutSelPrompt.empty:
            print('Saving prompt parquet')
            dataFramePtCutSelPrompt[colsToKeep].to_parquet('{0}/Prompt{1}_pT_{2:.0f}_{3:.0f}.parquet.gzip'.format(\
                outDirName, outSuffix, PtMin, PtMax), compression='gzip')
        if not dataFramePtCutSelFD.empty:
            print('Saving FD parquet')
            dataFramePtCutSelFD[colsToKeep].to_parquet('{0}/FD{1}_pT_{2:.0f}_{3:.0f}.parquet.gzip'.format(\
                outDirName, outSuffix, PtMin, PtMax), compression='gzip')
        if not dataFramePtCutSelPromptRefl.empty:
            print('Saving prompt refl parquet')
            dataFramePtCutSelPromptRefl[colsToKeep].to_parquet(\
                '{0}/PromptRefl{1}_pT_{2:.0f}_{3:.0f}.parquet.gzip'.format(\
                    outDirName, outSuffix, PtMin, PtMax), compression='gzip')
        if not dataFramePtCutSelFDRefl.empty:
            print('Saving FD refl parquet')
            dataFramePtCutSelFDRefl[colsToKeep].to_parquet('{0}/FDRefl{1}_pT_{2:.0f}_{3:.0f}.parquet.gzip'.format(\
                outDirName, outSuffix, PtMin, PtMax), compression='gzip')
else:
    if not args.parquet:
        print('Saving data tree')
        WriteTree(dataFramePtCutSel, colsToKeep, outTreeName, \
            '{0}/Data{1}_pT_{2:.0f}_{3:.0f}.root'.format(outDirName, outSuffix, PtMin, PtMax))
    else:
        print('Saving data parquet')
        dataFramePtCutSel[colsToKeep].to_parquet('{0}/Data{1}_pT_{2:.0f}_{3:.0f}.parquet.gzip'.format(\
            outDirName, outSuffix, PtMin, PtMax), compression='gzip')
