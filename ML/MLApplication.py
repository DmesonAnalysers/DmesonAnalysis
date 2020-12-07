'''
python standalone script to apply trained models to data using the hipe4ml package
run: python MLApplication.py cfgFileNameML.yml
Use the same config file as that used in the classification
'''
import os
import sys
import argparse
import yaml

from hipe4ml.model_handler import ModelHandler
from hipe4ml.tree_handler import TreeHandler

def main(): #pylint: disable=too-many-statements, too-many-branches
    # read config file
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileNameML.yml', help='config file name for ml')
    args = parser.parse_args()

    print('Loading analysis configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading analysis configuration: Done!')

    PtBins = [[a, b] for a, b in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'])]
    OutputLabels = inputCfg['output']['out_labels']
    ColumnsToSave = inputCfg['appl']['column_to_save_list']
    ModelList = inputCfg['ml']['saved_models']
    ModelHandls = []
    for iBin in range(len(PtBins)):
        ModelPath = ModelList[iBin]
        if not isinstance(ModelPath, str):
            print('\033[91mERROR: path to model not correctly defined!\033[0m')
            sys.exit()
        ModelPath = os.path.expanduser(ModelPath)
        print(f'Loaded saved model: {ModelPath}')
        ModelHandl = ModelHandler()
        ModelHandl.load_model_handler(ModelPath)
        ModelHandls.append(ModelHandl)

    for inputFile, outName in zip(inputCfg['standalone_appl']['inputs'], inputCfg['standalone_appl']['output_names']):
        print(f'Loading and preparing data file {inputFile}: ...', end='\r')
        DataHandler = TreeHandler(inputFile)
        DataHandler.slice_data_frame('pt_cand', PtBins, True)
        print(f'Loading and preparing data files {inputFile}: Done!')

        print('Applying ML model to dataframes: ...', end='\r')
        for iBin, PtBin in enumerate(PtBins):
            OutPutDirPt = os.path.join(os.path.expanduser(inputCfg['standalone_appl']['output_dir']),
                                                          f'pt{PtBin[0]}_{PtBin[1]}')
            if os.path.isdir(OutPutDirPt):
                print((f'\033[93mWARNING: Output directory \'{OutPutDirPt}\' already exists,'
                     ' overwrites possibly ongoing!\033[0m'))
            else:
                os.makedirs(OutPutDirPt)
            DataDfPtSel = DataHandler.get_slice(iBin)
            yPred = ModelHandls[iBin].predict(DataDfPtSel, inputCfg['ml']['raw_output'])
            ColumnsToSaveFinal = ColumnsToSave
            if not isinstance(ColumnsToSaveFinal, list):
                print('\033[91mERROR: column_to_save_list must be defined!\033[0m')
                sys.exit()
            if 'inv_mass' not in ColumnsToSaveFinal:
                print('\033[93mWARNING: inv_mass is not going to be saved in the output dataframe!\033[0m')
            if 'pt_cand' not in ColumnsToSaveFinal:
                print('\033[93mWARNING: pt_cand is not going to be saved in the output dataframe!\033[0m')
            if 'pt_B' in ColumnsToSaveFinal and 'pt_B' not in DataDfPtSel.columns:
                ColumnsToSaveFinal.remove('pt_B') # only in MC
            DataDfPtSel = DataDfPtSel.loc[:, ColumnsToSaveFinal]
            if ModelHandls[iBin].get_n_classes() < 3:
                DataDfPtSel['ML_output'] = yPred
            else:
                for Pred, Lab in enumerate(OutputLabels):
                    DataDfPtSel[f'ML_output_{Lab}'] = yPred[:, Pred]
            DataDfPtSel.to_parquet(f'{OutPutDirPt}/{outName}_pT_{PtBin[0]}_{PtBin[1]}_ModelApplied.parquet.gzip')
            del DataDfPtSel
        print('Applying ML model to dataframes: Done!')

main()
