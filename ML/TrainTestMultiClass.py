'''
python script to run basic training and application using the hipe4ml package
run: python TrainTestMulticlass.py cfgFileNameML.yml [--train, --apply]
--train -> to perform only the training and save the models in pkl
--apply -> to perform only the application loading saved models
'''
import os
import sys
import argparse
import yaml
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler

def data_prep(inputCfg, iBin, PtMin, PtMax, OutPutDirPt, DataDf, PromptDf, FDDf): #pylint: disable=too-many-statements
    '''
    function for data preparation
    '''
    DataDfPtSel = DataDf.query(f'{PtMin} < pt_cand < {PtMax}')
    BkgDfPtSel = DataDfPtSel.query(inputCfg['data_prep']['filt_bkg_mass'])
    PromptDfPtSel = PromptDf.query(f'{PtMin} < pt_cand < {PtMax}')
    FDDfPtSel = FDDf.query(f'{PtMin} < pt_cand < {PtMax}')

    nPrompt = len(PromptDfPtSel)
    nFD = len(FDDfPtSel)
    nBkg = len(BkgDfPtSel)
    print((f'Number of available candidates in {PtMin} < pT < {PtMax} GeV/c:\n     Prompt: {nPrompt}'
           f'\n     FD: {nFD}\n     Bkg: {nBkg}'))

    dataset_opt = inputCfg['data_prep']['dataset_opt']
    seed_split = inputCfg['data_prep']['seed_split']
    test_f = inputCfg['data_prep']['test_fraction']

    if dataset_opt == 'equal':

        nCandToKeep = min([nPrompt, nFD, nBkg])
        print(('Keep same number of prompt, FD, and background (minimum) for training and '
               f'testing ({1 - test_f}-{test_f}): {nCandToKeep}'))
        print(f'Fraction of real data candidates used for ML: {nCandToKeep/nBkg:.5f}')

        if nPrompt > nCandToKeep:
            print((f'Remaining prompt candidates ({nPrompt - nCandToKeep})'
                   'will be used for the efficiency together with test set'))
        if nFD > nCandToKeep:
            print((f'Remaining FD candidates ({nFD - nCandToKeep}) will be used for the '
                   'efficiency together with test set'))

        TotDfPtSel = pd.concat([BkgDfPtSel.iloc[:nCandToKeep], PromptDfPtSel.iloc[:nCandToKeep],
                                FDDfPtSel.iloc[:nCandToKeep]], sort=True)
        LabelsArray = [0] * nCandToKeep + [1] * nCandToKeep + [2] * nCandToKeep
        TrainSet, TestSet, yTrain, yTest = train_test_split(TotDfPtSel, LabelsArray, test_size=test_f,
                                                            random_state=seed_split)
        TrainTestData = [TrainSet, yTrain, TestSet, yTest]
        CandTypeFlags = pd.Series(yTest)
        PromptDfPtSelForEff = pd.concat([PromptDfPtSel.iloc[nCandToKeep:], TestSet[CandTypeFlags.values == 1]],
                                        sort=False)
        FDDfPtSelForEff = pd.concat([FDDfPtSel.iloc[nCandToKeep:], TestSet[CandTypeFlags.values == 2]], sort=False)
        del TotDfPtSel

    elif dataset_opt == 'max_signal':

        nCandBkg = round(inputCfg['ml']['bkg_mult'][iBin] * (nPrompt + nFD))
        print((f'Keep all prompt and FD and use {nCandBkg} bkg candidates for training and '
               f'testing ({1 - test_f}-{test_f})'))
        if nCandBkg >= nBkg:
            nCandBkg = nBkg
            print('WARNING: using all bkg available, not good!')
        print(f'Fraction of real data candidates used for ML: {nCandBkg/nBkg:.5f}')

        TotDfPtSel = pd.concat([BkgDfPtSel.iloc[:nCandBkg], PromptDfPtSel, FDDfPtSel], sort=True)
        LabelsArray = [0] * nCandBkg + [1] * nPrompt + [2] * nFD
        TrainSet, TestSet, yTrain, yTest = train_test_split(TotDfPtSel, LabelsArray, test_size=test_f,
                                                            random_state=seed_split)
        TrainTestData = [TrainSet, yTrain, TestSet, yTest]
        CandTypeFlags = pd.Series(yTest)
        PromptDfPtSelForEff = TestSet[CandTypeFlags.values == 1]
        FDDfPtSelForEff = TestSet[CandTypeFlags.values == 2]
        del TotDfPtSel

    else:
        print(f'ERROR: {dataset_opt} is not a valid option!')
        sys.exit()

    # plots
    VarsToDraw = inputCfg['ml']['plotting_columns']
    LegLabels = inputCfg['output']['leg_labels']
    OutputLabels = inputCfg['output']['out_labels']
    #_____________________________________________
    plot_utils.plot_distr([BkgDfPtSel, PromptDfPtSel, FDDfPtSel], VarsToDraw, (12, 7), 100, True, LegLabels, 0.3)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig(f'{OutPutDirPt}/DistributionsAll_pT_{PtMin}_{PtMax}.pdf')
    plt.close('all')
    #_____________________________________________
    CorrMatrixFig = plot_utils.plot_corr([BkgDfPtSel, PromptDfPtSel, FDDfPtSel], VarsToDraw, LegLabels)
    for Fig, Lab in zip(CorrMatrixFig, OutputLabels):
        plt.figure(Fig.number)
        plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
        Fig.savefig(f'{OutPutDirPt}/CorrMatrix{Lab}_pT_{PtMin}_{PtMax}.pdf')

    del BkgDfPtSel, PromptDfPtSel, FDDfPtSel
    return TrainTestData, DataDfPtSel, PromptDfPtSelForEff, FDDfPtSelForEff


def train_test(inputCfg, PtMin, PtMax, OutPutDirPt, TrainTestData):
    '''
    function for model training and testing
    '''
    modelClf = xgb.XGBClassifier()
    TrainCols = inputCfg['ml']['training_columns']
    HyperPars = inputCfg['ml']['hyper_par']
    if not isinstance(TrainCols, list):
        print('ERROR: training columns must be defined!')
        sys.exit()
    if not isinstance(HyperPars, dict):
        print('ERROR: hyper-parameters must be defined or be an empty dict!')
        sys.exit()
    ModelHandl = ModelHandler(modelClf, TrainCols, HyperPars)

    # hyperparams optimization --> not working with multi-class classification at the moment
    #HypRanges = {
    #    # # defines the maximum depth of a single tree (regularization)
    #    'max_depth': (1, 30),
    #    'learning_rate': (0.01, 0.3),  # learning rate
    #    'n_estimators': (50, 1000)  # number of boosting trees
    #}
    #ModelHandl.optimize_params_bayes(TrainTestData, HypRanges, None)

    # train and test the model with the updated hyperparameters
    ModelHandl.train_test_model(TrainTestData)
    yPredTest = ModelHandl.predict(TrainTestData[2], inputCfg['ml']['raw_output'], True)

    # save model handler in pickle
    ModelHandl.dump_model_handler(f'{OutPutDirPt}/ModelHandler_pT_{PtMin}_{PtMax}.pickle')

    #plots
    LegLabels = inputCfg['output']['leg_labels']
    OutputLabels = inputCfg['output']['out_labels']
    #_____________________________________________
    plt.rcParams["figure.figsize"] = (10, 7)
    MLOutputFig = plot_utils.plot_output_train_test(ModelHandl, TrainTestData, 80, inputCfg['ml']['raw_output'], 
                                                    LegLabels, True, inputCfg['plots']['train_test_log'], density=True)
    for Fig, Lab in zip(MLOutputFig, OutputLabels):
        Fig.savefig(f'{OutPutDirPt}/MLOutputDistr{Lab}_pT_{PtMin}_{PtMax}.pdf')
    #_____________________________________________
    plt.rcParams["figure.figsize"] = (8, 7)
    ROCCurveFig = plot_utils.plot_roc(TrainTestData[3], yPredTest, LegLabels)
    ROCCurveFig.savefig(f'{OutPutDirPt}/ROCCurveAll_pT_{PtMin}_{PtMax}.pdf')
    #_____________________________________________
    PrecisionRecallFig = plot_utils.plot_precision_recall(TrainTestData[3], yPredTest, LegLabels)
    PrecisionRecallFig.savefig(f'{OutPutDirPt}/PrecisionRecallAll_pT_{PtMin}_{PtMax}.pdf')
    #_____________________________________________
    plt.rcParams["figure.figsize"] = (12, 7)
    FeaturesImportanceFig = plot_utils.plot_feature_imp(TrainTestData[2][TrainCols], TrainTestData[3], ModelHandl)
    for iFig, Fig in enumerate(FeaturesImportanceFig):
        if iFig < 3:
            Fig.savefig(f'{OutPutDirPt}/FeatureImportance{OutputLabels[iFig]}_pT_{PtMin}_{PtMax}.pdf')
        else:
            Fig.savefig(f'{OutPutDirPt}/FeatureImportanceAll_pT_{PtMin}_{PtMax}.pdf')

    return ModelHandl


def appl(inputCfg, PtMin, PtMax, OutPutDirPt, ModelHandl, DataDfPtSel, PromptDfPtSelForEff, FDDfPtSelForEff):
    OutputLabels = inputCfg['output']['out_labels']
    print('Applying ML model to prompt dataframe: ...', end='\r')
    yPredPromptEff = ModelHandl.predict(PromptDfPtSelForEff, inputCfg['ml']['raw_output'], True)
    PromptDfPtSelForEff = PromptDfPtSelForEff.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        PromptDfPtSelForEff[f'ML_output_{Lab}'] = yPredPromptEff[:, Pred]
    PromptDfPtSelForEff.to_parquet(f'{OutPutDirPt}/Prompt_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to prompt dataframe: Done!')

    print('Applying ML model to FD dataframe: ...', end='\r')
    yPredFDEff = ModelHandl.predict(FDDfPtSelForEff, inputCfg['ml']['raw_output'], True)
    FDDfPtSelForEff = FDDfPtSelForEff.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        FDDfPtSelForEff[f'ML_output_{Lab}'] = yPredFDEff[:, Pred]
    FDDfPtSelForEff.to_parquet(f'{OutPutDirPt}/FD_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to FD dataframe: Done!')

    print('Applying ML model to data dataframe: ...', end='\r')
    yPredData = ModelHandl.predict(DataDfPtSel, inputCfg['ml']['raw_output'], True)
    DataDfPtSel = DataDfPtSel.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        DataDfPtSel[f'ML_output_{Lab}'] = yPredData[:, Pred]
    DataDfPtSel.to_parquet(f'{OutPutDirPt}/Data_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to data dataframe: Done!')


def main():
    # read config file
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileNameML.yml', help='config file name for ml')
    parser.add_argument("--train", help="perform only training and testing", action="store_true")
    parser.add_argument("--apply", help="perform only application", action="store_true")
    args = parser.parse_args()

    print('Loading analysis configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading analysis configuration: Done!')

    print('Loading data files: ...', end='\r')
    PromptDf = pd.read_parquet(inputCfg['input']['prompt'])
    FDDf = pd.read_parquet(inputCfg['input']['FD'])
    DataDf = pd.read_parquet(inputCfg['input']['data'])
    print('Loading data files: Done!')

    for iBin, (PtMin, PtMax) in enumerate(zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max'])):

        print(f'\n\033[94mStarting ML analysis --- {PtMin} < pT < {PtMax} GeV/c\033[0m')

        OutPutDirPt = os.path.join(inputCfg['output']['dir'], f'pt{PtMin}_{PtMax}')
        if os.path.isdir(OutPutDirPt):
            print('Output directory already exists, overwrites possibly ongoing!')
        else:
            os.mkdir(OutPutDirPt)

        # data preparation
        #_____________________________________________
        TrainTestData, DataDfPtSel, PromptDfPtSelForEff, FDDfPtSelForEff = data_prep( \
            inputCfg, iBin, PtMin, PtMax, OutPutDirPt, DataDf, PromptDf, FDDf)

        # training, testing
        #_____________________________________________
        if not args.apply:
            ModelHandl = train_test(inputCfg, PtMin, PtMax, OutPutDirPt, TrainTestData)
        else:
            ModelList = inputCfg['ml']['saved_models']
            ModelPath = ModelList[iBin]
            if not isinstance(ModelPath, str):
                print(f'ERROR: path to model not correctly defined!')
                sys.exit()
            print(f'Loaded saved model: {ModelPath}')
            ModelHandl = ModelHandler()
            ModelHandl.load_model_handler(ModelPath)

        # model application
        #_____________________________________________
        if not args.train:
            appl(inputCfg, PtMin, PtMax, OutPutDirPt, ModelHandl, DataDfPtSel, PromptDfPtSelForEff, FDDfPtSelForEff)

        # delete dataframes to release memory
        for data in TrainTestData:
            del data
        del DataDfPtSel, PromptDfPtSelForEff, FDDfPtSelForEff


main()
