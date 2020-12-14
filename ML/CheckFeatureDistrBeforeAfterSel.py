'''
python script to check variable distribution using the hipe4ml package
run: python ChecckDistr.py cfgFileName.yml cutSetFileName.yml
'''
import sys
import argparse
import yaml
import matplotlib.pyplot as plt
from hipe4ml import plot_utils

sys.path.append('..')
from utils.DfUtils import LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error,no-name-in-module

# inputs
parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with path of input dataframes for check')
parser.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                    help='cut set file name')
parser.add_argument('--outputDir', metavar='text', default='.',
                    help='output directory for plots')
args = parser.parse_args()

# input dataframes
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
isMC = inputCfg['isMC']
if isMC:
    dfPrompt = LoadDfFromRootOrParquet(inputCfg['tree']['filenamePrompt'],
                                       inputCfg['tree']['dirname'], inputCfg['tree']['treename'])
    dfFD = LoadDfFromRootOrParquet(inputCfg['tree']['filenameFD'],
                                   inputCfg['tree']['dirname'], inputCfg['tree']['treename'])
else:
    dfAll = LoadDfFromRootOrParquet(inputCfg['tree']['filenameAll'],
                                    inputCfg['tree']['dirname'], inputCfg['tree']['treename'])

# selections to be applied
with open(args.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']
selToApply = []
for iPt, _ in enumerate(cutVars['Pt']['min']):
    selToApply.append('')
    for varName in cutVars:
        if varName == 'InvMass':
            continue
        if selToApply[iPt] != '':
            selToApply[iPt] += ' & '
        selToApply[iPt] += f"{cutVars[varName]['min'][iPt]}<{cutVars[varName]['name']}<{cutVars[varName]['max'][iPt]}"

LegLabels = ['before selection', 'after selection']
varsToRemove = ['pt_B'] # HARD CODED

for (cuts, ptMin, ptMax) in zip(selToApply, cutVars['Pt']['min'], cutVars['Pt']['max']):
    print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
    if isMC:
        dfPromptList = [dfPrompt.query(f'{ptMin} < pt_cand < {ptMax}'), dfPrompt.astype(float).query(cuts)]
        dfFDList = [dfFD.query(f'{ptMin} < pt_cand < {ptMax}'), dfFD.astype(float).query(cuts)]

        varsToDraw = list(dfPromptList[0].columns)
        for varToRemove in varsToRemove:
            if varToRemove in varsToDraw:
                varsToDraw.remove(varToRemove)


        plot_utils.plot_distr(dfPromptList, varsToDraw, 100, LegLabels, figsize=(12, 7), density=True)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(f'{args.outputDir}/PromptDistrCompBeforeAfterSel_pT_{ptMin}_{ptMax}.pdf')
        plt.close('all')
        del dfPromptList

        plot_utils.plot_distr(dfFDList, varsToDraw, 100, LegLabels, figsize=(12, 7), density=True)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(f'{args.outputDir}/FDDistrCompBeforeAfterSel_pT_{ptMin}_{ptMax}.pdf')
        plt.close('all')
        del dfFDList

    else:
        dfAllList = [dfAll.query(f'{ptMin} < pt_cand < {ptMax}'), dfAll.astype(float).query(cuts)]

        varsToDraw = list(dfAllList[0].columns)
        for varToRemove in varsToRemove:
            if varToRemove in varsToDraw:
                varsToDraw.remove(varToRemove)

        plot_utils.plot_distr(dfAllList, varsToDraw, 100, LegLabels, figsize=(12, 7), density=True)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(f'{args.outputDir}/DataDistrCompBeforeAfterSel_pT_{ptMin}_{ptMax}.pdf')
        plt.close('all')
        del dfAllList
