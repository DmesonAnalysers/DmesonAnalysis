'''
python script to filter tree from task output and save output trees in root (default) or parquet files for ML studies
run: python FilterTree4ML.py cfgFileName.yml
with option --root the output files are saved into root files instead of parquet ones
'''

import sys
import argparse
import numpy as np
import yaml
sys.path.append('..')
from utils.DfUtils import WriteTree, FilterBitDf, LoadDfFromRootOrParquet, GetMind0 #pylint: disable=wrong-import-position,import-error,no-name-in-module


bitSignal = 0
bitBkg = 1
bitPrompt = 2
bitFD = 3
bitRefl = 4
bitSecPeak = 9

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='cfgFileName.yml',
                    help='input config yaml file name')
parser.add_argument('--root', default=False, action='store_true',
                    help='flag to save output files into root files instead of parquet')

args = parser.parse_args()
print('Opening input file')
with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inFileNames = cfg['infile']['filename']
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

dataFrame = LoadDfFromRootOrParquet(inFileNames, inDirName, inTreeName)

if not colsToKeep:
    colsToKeep = list(dataFrame.columns)
    colsToKeep.remove('cand_type')

print('Applying selections')
dataFramePtCut = dataFrame.query(f'pt_cand > {PtMin} & pt_cand < {PtMax}')
del dataFrame
if preSelections:
    dataFramePtCutSel = dataFramePtCut.astype(float).query(preSelections)
    del dataFramePtCut
else:
    dataFramePtCutSel = dataFramePtCut

if cfg['missingvalues']['enable']:
    dataFramePtCutSel = dataFramePtCutSel.replace(cfg['missingvalues']['toreplace'], value=np.nan)

if cfg['singletrackvars'] and cfg['singletrackvars']['addAODfiltervars']:
    # this assumes that we are analysing a 3 prong!
    if set(['pt_prong0', 'pt_prong1', 'pt_prong2']).issubset(dataFramePtCutSel.columns):
        dataFramePtCutSel['pt_prong_min'] = dataFramePtCutSel[['pt_prong0', 'pt_prong1', 'pt_prong2']].min(axis=1)
        colsToKeep.append('pt_prong_min')
        if set(['imp_par_prong0', 'imp_par_prong1', 'imp_par_prong2']).issubset(dataFramePtCutSel.columns):
            dataFramePtCutSel['imp_par_min_ptgtr2'] = dataFramePtCutSel.apply(lambda x: GetMind0(
                [x['pt_prong0'], x['pt_prong1'], x['pt_prong2']],
                [x['imp_par_prong0'], x['imp_par_prong1'], x['imp_par_prong2']], 2), axis=1)
            colsToKeep.append('imp_par_min_ptgtr2')

if isMC:
    print('Getting bkg dataframe')
    dataFramePtCutSelBkg = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitBkg])
    print('Getting prompt dataframe')
    dataFramePtCutSelPrompt = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSignal, bitPrompt], 'and')
    dataFramePtCutSelPrompt = FilterBitDf(dataFramePtCutSelPrompt, 'cand_type', [bitRefl], 'not')
    print('Getting FD dataframe')
    dataFramePtCutSelFD = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSignal, bitFD], 'and')
    dataFramePtCutSelFD = FilterBitDf(dataFramePtCutSelFD, 'cand_type', [bitRefl], 'not')
    print('Getting reflected prompt dataframe')
    dataFramePtCutSelPromptRefl = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSignal, bitPrompt, bitRefl], 'and')
    print('Getting reflected signal dataframe')
    dataFramePtCutSelFDRefl = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSignal, bitFD, bitRefl], 'and')
    print('Getting second-peak prompt dataframe')
    dataFramePtCutSelSecPeakPrompt = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSecPeak, bitPrompt], 'and')
    dataFramePtCutSelSecPeakPrompt = FilterBitDf(dataFramePtCutSelSecPeakPrompt, 'cand_type', [bitRefl], 'not')
    print('Getting second-peak FD dataframe')
    dataFramePtCutSelSecPeakFD = FilterBitDf(dataFramePtCutSel, 'cand_type', [bitSecPeak, bitFD], 'and')
    dataFramePtCutSelSecPeakFD = FilterBitDf(dataFramePtCutSelSecPeakFD, 'cand_type', [bitRefl], 'not')
    del dataFramePtCutSel

    if args.root:
        if not dataFramePtCutSelBkg.empty:
            print('Saving bkg tree')
            WriteTree(dataFramePtCutSelBkg, colsToKeep, outTreeName,
                      f'{outDirName}/Bkg{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.root')
        if not dataFramePtCutSelPrompt.empty:
            print('Saving prompt tree')
            WriteTree(dataFramePtCutSelPrompt, colsToKeep, outTreeName,
                      f'{outDirName}/Prompt{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.root')
        if not dataFramePtCutSelFD.empty:
            print('Saving FD tree')
            WriteTree(dataFramePtCutSelFD, colsToKeep, outTreeName,
                      f'{outDirName}/FD{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.root')
        if not dataFramePtCutSelPromptRefl.empty:
            print('Saving prompt refl tree')
            WriteTree(dataFramePtCutSelPromptRefl, colsToKeep, outTreeName,
                      f'{outDirName}/PromptRefl{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.root')
        if not dataFramePtCutSelFDRefl.empty:
            print('Saving FD refl tree')
            WriteTree(dataFramePtCutSelFDRefl, colsToKeep, outTreeName,
                      f'{outDirName}/FDRefl{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.root')
        if not dataFramePtCutSelSecPeakPrompt.empty:
            print('Saving prompt tree')
            WriteTree(dataFramePtCutSelSecPeakPrompt, colsToKeep, outTreeName,
                      f'{outDirName}/SecPeakPrompt{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.root')
        if not dataFramePtCutSelSecPeakFD.empty:
            print('Saving FD tree')
            WriteTree(dataFramePtCutSelSecPeakFD, colsToKeep, outTreeName,
                      f'{outDirName}/SecPeakFD{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.root')
    else:
        if not dataFramePtCutSelBkg.empty:
            print('Saving bkg parquet')
            dataFramePtCutSelBkg[colsToKeep].to_parquet(\
                f'{outDirName}/Bkg{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip', compression='gzip')
        if not dataFramePtCutSelPrompt.empty:
            print('Saving prompt parquet')
            dataFramePtCutSelPrompt[colsToKeep].to_parquet(\
                f'{outDirName}/Prompt{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip', compression='gzip')
        if not dataFramePtCutSelFD.empty:
            print('Saving FD parquet')
            dataFramePtCutSelFD[colsToKeep].to_parquet(\
                f'{outDirName}/FD{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip', compression='gzip')
        if not dataFramePtCutSelPromptRefl.empty:
            print('Saving prompt refl parquet')
            dataFramePtCutSelPromptRefl[colsToKeep].to_parquet(\
                f'{outDirName}/PromptRefl{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip', compression='gzip')
        if not dataFramePtCutSelFDRefl.empty:
            print('Saving FD refl parquet')
            dataFramePtCutSelFDRefl[colsToKeep].to_parquet(\
                f'{outDirName}/FDRefl{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip', compression='gzip')
        if not dataFramePtCutSelSecPeakPrompt.empty:
            print('Saving prompt parquet')
            dataFramePtCutSelSecPeakPrompt[colsToKeep].to_parquet(\
                f'{outDirName}/SecPeakPrompt{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip', compression='gzip')
        if not dataFramePtCutSelSecPeakFD.empty:
            print('Saving FD parquet')
            dataFramePtCutSelSecPeakFD[colsToKeep].to_parquet(\
                f'{outDirName}/SecPeakFD{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip', compression='gzip')
else:
    if args.root:
        print('Saving data tree')
        WriteTree(dataFramePtCutSel, colsToKeep, outTreeName,
                  f'{outDirName}/Data{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.root')
    else:
        print('Saving data parquet')
        dataFramePtCutSel[colsToKeep].to_parquet(\
            f'{outDirName}/Data{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip', compression='gzip')
