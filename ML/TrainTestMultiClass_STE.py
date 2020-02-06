'''
python script to run basic training and application using the hipe4ml package
run: python TrainTestMulticlass.py cfgFileNameML.yml
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

# read config file
parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileNameML.yml',
                    help='config file name for ml')
args = parser.parse_args()

print('Loading analysis configuration')
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

PromptDf = pd.read_parquet(inputCfg['input']['prompt'])
FDDf = pd.read_parquet(inputCfg['input']['FD'])
BkgDf = pd.read_parquet(inputCfg['input']['Bkg'])

LegLabels = inputCfg['output']['leg_labels']
OutputLabels = inputCfg['output']['out_labels']
OutPutDir = inputCfg['output']['dir']

# data preparation
#_____________________________________________
for (PtMin, PtMax) in zip(inputCfg['pt_ranges']['min'], inputCfg['pt_ranges']['max']):

    print(f'\n\nStarting analysis --- {PtMin} < pT < {PtMax} GeV/c ')

    OutPutDirPt = os.path.join(OutPutDir, f'pt{PtMin}_{PtMax}')

    if os.path.isdir(OutPutDirPt):
        print('\nOutput directory already exists, overwrites possibly ongoing!')
    if not os.path.exists(OutPutDirPt):
        os.makedirs(OutPutDirPt)

    BkgDfPtSel = BkgDf.query(f'{PtMin} < pt_cand < {PtMax}')
    PromptDfPtSel = PromptDf.query(f'{PtMin} < pt_cand < {PtMax}')
    FDDfPtSel = FDDf.query(f'{PtMin} < pt_cand < {PtMax}')

    '''
    #hiding nPr=nFD=nBkg

    nPrompt = len(PromptDfPtSel)
    nFD = len(FDDfPtSel)
    nBkg = len(BkgDfPtSel)
    print(f'\nNumber of available candidates in {PtMin} < pT < {PtMax} GeV/c:')
    print('Prompt:', nPrompt)
    print('FD:', nFD)
    print('Bkg:', nBkg)

    nCandToKeep = min([nPrompt, nFD, nBkg])
    print('\nKeep same number of prompt, FD, and background (minimum) for training and testing (50-50):',
          nCandToKeep)
    
    print('Fraction of real data candidates used for ML: ', nCandToKeep / nBkg)

    if nPrompt > nFD:
        print(f'Remaining prompt candidates ({nPrompt - nCandToKeep})' \
              'will be used for the efficiency together with test set')
    elif nFD > nPrompt:
        print(f'Remaining FD candidates ({nFD - nCandToKeep}) will be used for the efficiency together with test set\n')
    
    #hiding Set implementation
    TotDfPtSel = pd.concat([BkgDfPtSel.iloc[:nCandToKeep].copy(), PromptDfPtSel.iloc[:nCandToKeep].copy(),
                            FDDfPtSel.iloc[:nCandToKeep]PtMin.copy()], sort=True)

    LabelsArray = [0] * nCandToKeep + [1] * nCandToKeep + [2] * nCandToKeep

    TrainSet, TestSet, yTrain, yTest = train_test_split(TotDfPtSel, LabelsArray, test_size=0.5, random_state=42)
    TrainTestData = [TrainSet, yTrain, TestSet, yTest]

    CandTypeFlags = pd.Series(yTest)
    PromptDfPtSelForEff = pd.concat([PromptDfPtSel.iloc[nCandToKeep:].copy(), TestSet[CandTypeFlags.values == 1]],
                                    sort=False)
    FDDfPtSelForEff = pd.concat([FDDfPtSel.iloc[nCandToKeep:].copy(), TestSet[CandTypeFlags.values == 2]], sort=False)
    
    excluding training and testing 
    # training, testing, and model application
    #_____________________________________________
    modelClf = xgb.XGBClassifier()
    TrainCols = inputCfg['ml']['training_columns']
    if not TrainCols:
        print('ERROR: training columns must be defined!')
        sys.exit()
    ModelHandl = ModelHandler(modelClf, TrainCols)

    # hyperparams optimization --> not working with multi-class classification at the moment
    #HypRanges = {
    #    # # defines the maximum depth of a single tree (regularization)
    #    'max_depth': (1, 30),
    #    'learning_rate': (0.01, 0.3),  # learning rate
    #    'n_estimators': (50, 1000)  # number of boosting trees
    #}
    #ModelHandl.optimize_params_bayes(TrainTestData, HypRanges, None)

    ModelHandl.set_model_params(inputCfg['ml']['hyper_par'])

    # train and test the model with the updated hyperparameters
    ModelHandl.train_test_model(TrainTestData)
    yPredTest = ModelHandl.predict(TrainTestData[2])

    # apply model
    print('Applying ML model to prompt dataframe: ...', end='\r')
    yPredPromptEff = ModelHandl.predict(PromptDfPtSelForEff)
    PromptDfPtSelForEff = PromptDfPtSelForEff.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        PromptDfPtSelForEff[f'ML_output_{Lab}'] = yPredPromptEff[:, Pred]
    PromptDfPtSelForEff.to_parquet(f'{OutPutDirPt}/Prompt_Dpluspp5TeV_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to prompt dataframe: Done!')

    print('Applying ML model to FD dataframe: ...', end='\r')
    yPredFDEff = ModelHandl.predict(FDDfPtSelForEff)
    FDDfPtSelForEff = FDDfPtSelForEff.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        FDDfPtSelForEff[f'ML_output_{Lab}'] = yPredFDEff[:, Pred]
    FDDfPtSelForEff.to_parquet(f'{OutPutDirPt}/FD_Dpluspp5TeV_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to FD dataframe: Done!')

    print('Applying ML model to data dataframe: ...', end='\r')
    yPredData = ModelHandl.predict(DataDfPtSel)
    DataDfPtSel = DataDfPtSel.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        DataDfPtSel[f'ML_output_{Lab}'] = yPredData[:, Pred]
    DataDfPtSel.to_parquet(f'{OutPutDirPt}/Data_Dpluspp5TeV_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to data dataframe: Done!')

    # save model handler in pickle
    #_____________________________________________
    ModelHandl.dump_model_handler(f'{OutPutDirPt}/ModelHandler_pT_{PtMin}_{PtMax}.pickle')
    
    '''
    # plots
    #_____________________________________________
    VarsToDraw = inputCfg['ml']['plotting_columns']
    #_____________________________________________
    FeatDistrFig = plot_utils.plot_distr([BkgDfPtSel, PromptDfPtSel, FDDfPtSel], VarsToDraw, (12, 7),
                                         100, True, LegLabels)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig(f'{OutPutDirPt}/DistributionsAll_pT_{PtMin}_{PtMax}.pdf')
    plt.close('all')

    #_____________________________________________
    CorrMatrixFig = plot_utils.plot_corr([BkgDfPtSel, PromptDfPtSel, FDDfPtSel], VarsToDraw, LegLabels)
    for Fig, Lab in zip(CorrMatrixFig, OutputLabels):
        plt.figure(Fig.number)
        plt.subplots_adjust(left=0.2, bottom=0.1, right=0.95, top=0.9)
        Fig.savefig(f'{OutPutDirPt}/CorrMatrix{Lab}_pT_{PtMin}_{PtMax}.pdf')
    
    '''
    excludng training plot
    #_____________________________________________
    plt.rcParams["figure.figsize"] = (10, 7)
    MLOutputFig = plot_utils.plot_output_train_test(ModelHandl, TrainTestData, 80, True, LegLabels)
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
    FeaturesImportanceFig = plot_utils.plot_feature_imp(TrainSet[TrainCols], yTrain, ModelHandl)
    for iFig, Fig in enumerate(FeaturesImportanceFig):
        if iFig < 3:
            Fig.savefig(f'{OutPutDirPt}/FeatureImportance{OutputLabels[iFig]}_pT_{PtMin}_{PtMax}.pdf')
        else:
            Fig.savefig(f'{OutPutDirPt}/FeatureImportanceAll_pT_{PtMin}_{PtMax}.pdf')
    
    # delete dataframes to release memory
    del DataDfPtSel, BkgDfPtSel, PromptDfPtSel, FDDfPtSel, TrainSet, TestSet, yTrain, yTest
    '''
