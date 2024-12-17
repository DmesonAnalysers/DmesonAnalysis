'''
python standalone script to apply trained models to data using the hipe4ml package
run: python MLApplication.py cfgFileNameML.yml
Use the same config file as that used in the classification
'''
import os
import sys
import argparse
import yaml
import matplotlib.pyplot as plt

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

    if inputCfg.get('savecfg'):
        # Save the YAML file to the folder
        if not os.path.isdir(os.path.expanduser(inputCfg['standalone_appl']['output_dir'])):
            os.makedirs(os.path.expanduser(inputCfg['standalone_appl']['output_dir']))
        with open(f'{os.path.expanduser(inputCfg["standalone_appl"]["output_dir"])}/cfg.yml', 'w') as ymlOutFile:
            yaml.dump(inputCfg, ymlOutFile, default_flow_style=False)

    PtBins = [[a, b] for a, b in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'])]
    OutputLabels = [inputCfg['output']['out_labels']['Bkg'],
                    inputCfg['output']['out_labels']['Prompt']]
    if inputCfg['output']['out_labels']['FD'] is not None:
        OutputLabels.append(inputCfg['output']['out_labels']['FD'])
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

    treename = inputCfg['standalone_appl']['treename']
    for inputFile, outName in zip(inputCfg['standalone_appl']['inputs'], inputCfg['standalone_appl']['output_names']):
        print(f'Loading and preparing data file {inputFile}: ...', end='\r')
        if treename is None:
            DataHandler = TreeHandler(inputFile)
        else:
            DataHandler = TreeHandler(inputFile, treename)

        DataHandler.slice_data_frame('fPt', PtBins, True)
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
            if 'fM' not in ColumnsToSaveFinal:
                print('\033[93mWARNING: fM is not going to be saved in the output dataframe!\033[0m')
            if 'fPt' not in ColumnsToSaveFinal:
                print('\033[93mWARNING: fPt is not going to be saved in the output dataframe!\033[0m')
            if 'pt_B' in ColumnsToSaveFinal and 'pt_B' not in DataDfPtSel.columns:
                ColumnsToSaveFinal.remove('pt_B') # only in MC
            DataDfPtSel = DataDfPtSel.loc[:, ColumnsToSaveFinal]
            if ModelHandls[iBin].get_n_classes() < 3:
                DataDfPtSel['ML_output'] = yPred
            else:
                for Pred, Lab in enumerate(OutputLabels):
                    DataDfPtSel[f'ML_output_{Lab}'] = yPred[:, Pred]
            DataDfPtSel.to_parquet(f'{OutPutDirPt}/{outName}_pT_{PtBin[0]}_{PtBin[1]}_ModelApplied.parquet.gzip')
            
            if inputCfg.get('savedistrs'):
                plt.figure(figsize=(10, 6))
                for col in DataDfPtSel.columns:
                    if 'ML_output' in col:
                        plt.hist(DataDfPtSel[col], bins=100, alpha=0.5, label=col, log=True)
                plt.title(f'Distributions of ML Outputs for {outName}')
                plt.xlabel('Score')
                plt.ylabel('Frequency (log scale)')
                plt.legend()
                plt.savefig(f"{OutPutDirPt}/{outName}Distrs.pdf", format="pdf", bbox_inches="tight")
            
            del DataDfPtSel
        print('Applying ML model to dataframes: Done!')

main()
