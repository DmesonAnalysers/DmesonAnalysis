'''
python script to filter tree from task output and save output trees in parquet files for ML studies
run: python FilterTree4ML.py cfgFileName.yml
'''

import sys
import argparse
import numpy as np
import yaml
sys.path.append('..')
from utils.DfUtils import FilterBitDf, LoadDfFromRootOrParquet, GetMind0 #pylint: disable=wrong-import-position,import-error

# common bits
bitSignal = 0
bitBkg = 1
bitPrompt = 2
bitFD = 3
bitRefl = 4
# channel specific bits
# Ds
bitSecPeakDs = 9
# LctopK0s
bitLctopK0s = 9
# LctopiL
bitLctopiL = 10
# LctopKpi
bitLcNonRes = 9
bitLcLambda1520 = 10
bitLcKStar = 11
bitLcDelta = 12

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('configfile', metavar='text', default='cfgFileName.yml',
                    help='input config yaml file name')

args = parser.parse_args()
print('Opening input file')
with open(args.configfile, 'r') as ymlCfgFile:
    cfg = yaml.load(ymlCfgFile, yaml.FullLoader)

channel = cfg['channel']
if channel not in ['Ds', 'D0', 'Dplus', 'Dstar', 'LctopKpi', 'LctopK0s', 'LctopiL']:
    print('Error: only Ds, D0, Dplus, Dstar, LctopKpi, LctopK0s, and LctopiL channels are implemented! Exit')
    sys.exit()

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

bitsForSel, labelsContr = ({} for _ in range(2))
labelsContr = {'bkg': 'Bkg', 'prompt_sig': 'Prompt', 'FD_sig': 'FD',
               'prompt_sig_refl': 'PromptRefl', 'FD_sig_refl': 'FDRefl',
               'prompt_sec_peak': 'SecPeakPrompt', 'FD_sec_peak': 'SecPeakFD',
               'prompt_sig_nonreso': 'PromptNonRes', 'FD_sig_nonreso': 'FDNonRes',
               'prompt_sig_Lambda1520': 'PromptLambda1520', 'FD_sig_Lambda1520': 'FDLambda1520',
               'prompt_sig_KStar': 'PromptKStar', 'FD_sig_KStar': 'FDKStar',
               'prompt_sig_Delta': 'PromptDelta', 'FD_sig_Delta': 'FDDelta'}

if isMC:
    if 'Ds' in channel:
        bitsForSel = {'bkg': [bitBkg],
                      'prompt_sig': [bitSignal, bitPrompt], 'FD_sig': [bitSignal, bitFD],
                      'prompt_sig_refl': [bitSignal, bitPrompt, bitRefl], 'FD_sig_refl': [bitSignal, bitFD, bitRefl],
                      'prompt_sec_peak': [bitSecPeakDs, bitPrompt], 'FD_sec_peak': [bitSecPeakDs, bitFD]}
    elif 'Dplus' in channel:
        bitsForSel = {'bkg': [bitBkg], 'prompt_sig': [bitSignal, bitPrompt], 'FD_sig': [bitSignal, bitFD]}
    elif 'Dstar' in channel:
        bitsForSel = {'bkg': [bitBkg], 'prompt_sig': [bitSignal, bitPrompt], 'FD_sig': [bitSignal, bitFD]}
    elif 'D0' in channel:
        bitsForSel = {'bkg': [bitBkg], 'prompt_sig': [bitSignal, bitPrompt], 'FD_sig': [bitSignal, bitFD], 'prompt_sig_refl': [bitSignal, bitPrompt, bitRefl], 'FD_sig_refl': [bitSignal, bitFD, bitRefl]}
    elif 'LctopKpi' in channel:
        bitsForSel = {'bkg': [bitBkg],
                      'prompt_sig_nonreso': [bitSignal, bitLcNonRes, bitPrompt],
                      'FD_sig_nonreso': [bitSignal, bitLcNonRes, bitFD],
                      'prompt_sig_Lambda1520': [bitSignal, bitLcLambda1520, bitPrompt],
                      'FD_sig_Lambda1520': [bitSignal, bitLcLambda1520, bitFD],
                      'prompt_sig_KStar': [bitSignal, bitLcKStar, bitPrompt],
                      'FD_sig_KStar': [bitSignal, bitLcKStar, bitFD],
                      'prompt_sig_Delta': [bitSignal, bitLcDelta, bitPrompt],
                      'FD_sig_Delta': [bitSignal, bitLcDelta, bitFD],
                      'prompt_sig_refl': [bitSignal, bitPrompt, bitRefl], 'FD_sig': [bitSignal, bitFD, bitRefl]}
    elif 'LctopK0s' in channel:
        bitsForSel = {'bkg': [bitBkg],
                      'prompt_sig': [bitSignal, bitLctopK0s, bitPrompt], 'FD_sig': [bitSignal, bitLctopK0s, bitFD]}
    elif 'LctopiL' in channel:
        bitsForSel = {'bkg': [bitBkg],
                      'prompt_sig': [bitSignal, bitLctopiL, bitPrompt], 'FD_sig': [bitSignal, bitLctopiL, bitFD]}

    for contr in bitsForSel:
        print(f'Getting {labelsContr[contr]} dataframe')
        dataFramePtCutSelContr = FilterBitDf(dataFramePtCutSel, 'cand_type', bitsForSel[contr], 'and')
        # always check that it is not reflected, unless is the reflection contribution
        if 'refl' not in contr:
            dataFramePtCutSelContr = FilterBitDf(dataFramePtCutSelContr, 'cand_type', [bitRefl], 'not')

        if not dataFramePtCutSelContr.empty:
            print(f'Saving {labelsContr[contr]} parquet')
            dataFramePtCutSelContr[colsToKeep].to_parquet(
                f'{outDirName}/{labelsContr[contr]}{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip',
                compression='gzip')
else:
    print('Saving data to parquet')
    dataFramePtCutSel[colsToKeep].to_parquet(f'{outDirName}/Data{outSuffix}_pT_{PtMin:.0f}_{PtMax:.0f}.parquet.gzip',
                                             compression='gzip')
