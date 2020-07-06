'''
python script to check variable distribution using the hipe4ml package
run: python ChecckDistr.py cfgFileNameCheck.yml
'''
import argparse
import yaml
import pandas as pd
import matplotlib.pyplot as plt
from hipe4ml import plot_utils

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
    for filePath in inputCfg['input']['files']:
        DfList.append(pd.read_parquet(filePath))
    print('Loading data files: Done!')

    VarsToDraw = inputCfg['plotting_columns']
    LegLabels = inputCfg['output']['leg_labels']
    OutPutDir = inputCfg['output']['dir']

    for PtMin, PtMax, LimMin, LimMax in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'],
                                            inputCfg['plot_lim_min'], inputCfg['plot_lim_max']):
        print(f'Plot variable distributions --- {PtMin} < pT < {PtMax} GeV/c')
        DfListPt = []
        for df in DfList:
            DfListPt.append(df.query(f'{PtMin} < pt_cand < {PtMax}'))
        DistrPlot = plot_utils.plot_distr(DfListPt, VarsToDraw, 100, LegLabels, figsize=(12, 12), density=True,
                                          histtype='step', grid=False, log=True)
        plt.subplots_adjust(left=0.1, bottom=0.05, right=0.95, top=0.95, hspace=0.4)
        for ax, minVar, maxVar, xLabel in zip(DistrPlot, LimMin, LimMax, inputCfg['xaxes_label']):
            ax.set_xlim(minVar, maxVar)    
            ax.set_xlabel(xLabel)
            ax.set_title('')
        plt.savefig(f'{OutPutDir}/DistrComp_pT_{PtMin}_{PtMax}.pdf')
        plt.close('all')
        del DfListPt

    del DfList


main()
