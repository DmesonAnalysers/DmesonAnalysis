'''
python script to check variable distribution using the hipe4ml package
run: python ChecckDistr.py cfgFileNameCheck.yml
'''
import argparse
import yaml
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hipe4ml import plot_utils

sys.path.insert(0, '..')
from utils.DfUtils import LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error

def main():
    # read config file
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileNameCheck.yml',
                        help='config file name for check')
    args = parser.parse_args()

    print('Loading check configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading check configuration: Done!')

    print('Loading data files: ...', end='\r')
    DfList = []
    inTreeName = inputCfg['input']['treename']
    for filePath in inputCfg['input']['files']:
        DfList.append(LoadDfFromRootOrParquet(filePath, "", inTreeName))

    print('Loading data files: Done!')

    print('Appling simple pre-filtering: ...', end='\r')
    DfListSel = []
    for df, query in zip(DfList, inputCfg['queries']):
        DfListSel.append(df.query(query))
    print('Pre-filtering: Done!')
    del DfList

    VarsToDraw = inputCfg['plotting_columns']
    LegLabels = inputCfg['output']['leg_labels']
    Colors = inputCfg['output']['colors']
    OutPutDir = inputCfg['output']['dir']

    for PtMin, PtMax, LimMin, LimMax in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'],
                                            inputCfg['plot_lim_min'], inputCfg['plot_lim_max']):
        print(f'Plot variable distributions --- {PtMin} < pT < {PtMax} GeV/c')
        DfListPt = []
        for df in DfListSel:
            DfListPt.append(df.query(f'{PtMin} < pt_cand < {PtMax}'))
        #print(len(DfListPt), len(Colors))
        DistrPlot = plot_utils.plot_distr(DfListPt, VarsToDraw, 1000, LegLabels, figsize=(6, 6), density=True,
                                          histtype='stepfilled', grid=False, log=True, colors=Colors, alpha=0.3)
        plt.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.95, hspace=0.4)
        if not isinstance(DistrPlot, np.ndarray):
            DistrPlot = np.array([DistrPlot])
        print (len(DistrPlot),len (LimMin), len(LimMax), len(inputCfg['xaxes_label']))
        for ax, minVar, maxVar, xLabel in zip(DistrPlot, LimMin, LimMax, inputCfg['xaxes_label']):
             ax.set_xlim(minVar, maxVar)

             ax.set_xlabel(xLabel, fontsize=10, ha='right', position=(1, 20))
             ax.set_ylabel('Counts (arb. units)', fontsize=10, ha='right', position=(20, 1))
             plt.legend(frameon=False, fontsize=10, loc='best')

             ax.set_title('')
             '''
             textstr = r'pp, $\sqrt{s}$ = 5.02 TeV'
             textstr2 = r'$3 < p_{\mathrm{T}} < 4~\mathrm{GeV}/c$'

             ax.text(0.56, 0.75, textstr, transform=ax.transAxes, fontsize=15,
                    verticalalignment='top')
             ax.text(0.56, 0.69, textstr2, transform=ax.transAxes, fontsize=15,
                    verticalalignment='top')
            '''
             plt.tight_layout()
        plt.savefig(f'{OutPutDir}/NsigzoomDistrComp_pT_{PtMin}_{PtMax}.pdf')
        plt.close('all')
        del DfListPt

    del DfListSel

main()
