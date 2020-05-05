'''
python script to check ML distribution with and without improver using the hipe4ml package
run: python CheckMLProbImprover.py cfgFileName.yml
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
args = parser.parse_args()

# input dataframes
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inputDict = inputCfg['input']

dfPrompt = LoadDfFromRootOrParquet(inputDict['filenamePrompt'], inputDict['dirname'], inputDict['treename'])
dfFD = LoadDfFromRootOrParquet(inputDict['filenameFD'], inputDict['dirname'], inputDict['treename'])
dfPromptNoImp = LoadDfFromRootOrParquet(inputDict['filenamePromptNoImp'], inputDict['dirname'], inputDict['treename'])
dfFDNoImp = LoadDfFromRootOrParquet(inputDict['filenameFDNoImp'], inputDict['dirname'], inputDict['treename'])

LegLabels = inputCfg['output']['leg_labels']
varsToRemove = inputCfg['var_to_remove']

for (ptMin, ptMax) in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max']):
    print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
    dfPromptList = [dfPrompt.query(f'{ptMin} < pt_cand < {ptMax}'), dfPromptNoImp.query(f'{ptMin} < pt_cand < {ptMax}')]
    dfFDList = [dfFD.query(f'{ptMin} < pt_cand < {ptMax}'), dfFDNoImp.query(f'{ptMin} < pt_cand < {ptMax}')]

    varsToDraw = list(dfPromptList[0].columns)
    for varToRemove in varsToRemove:
        if varToRemove in varsToDraw:
            varsToDraw.remove(varToRemove)

    outputDir = inputCfg['output']['dir']
    plot_utils.plot_distr(dfPromptList, varsToDraw, (12, 10), 100, True, LegLabels)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig(f'{outputDir}/PromptDistrCompImprover_pT_{ptMin}_{ptMax}.pdf')
    plt.close('all')
    del dfPromptList

    plot_utils.plot_distr(dfFDList, varsToDraw, (12, 10), 100, True, LegLabels)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig(f'{outputDir}/FDDistrCompImprover_pT_{ptMin}_{ptMax}.pdf')
    plt.close('all')
    del dfFDList
