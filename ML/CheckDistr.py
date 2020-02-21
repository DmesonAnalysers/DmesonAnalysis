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

    for (PtMin, PtMax) in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max']):
        print(f'Plot variable distributions --- {PtMin} < pT < {PtMax} GeV/c')
        DfListPt = []
        for df in DfList:
            DfListPt.append(df.query(f'{PtMin} < pt_cand < {PtMax}'))
        VarsToDraw = inputCfg['plotting_columns']
        LegLabels = inputCfg['output']['leg_labels']
        OutPutDir = inputCfg['output']['dir']
        plot_utils.plot_distr(DfListPt, VarsToDraw, (12, 7), 100, True, LegLabels)
        plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
        plt.savefig(f'{OutPutDir}/DistrComp_pT_{PtMin}_{PtMax}.pdf')
        plt.close('all')
        del DfListPt

    del DfList


main()
