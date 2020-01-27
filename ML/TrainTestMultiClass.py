'''
python script to run basic training and application using the hipe4ml package
run: python TrainTestMulticlass.py cfgFileName.yml
'''

import pandas as pd
import xgboost as xgb
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt

from hipe4ml import plot_utils
from hipe4ml.model_handler import ModelHandler

# data preparation
#_____________________________________________

# TODO: change hard coded parameters with parameters set ina config yaml file
PtMins = [1, 2, 4, 6, 12]
PtMaxs = [2, 4, 6, 12, 50]

PromptDf = pd.read_parquet('../../AnalysisNonPromptDpp2017/Dplus/MC/LHC18a4a2/Prompt_Dpluspp5TeV_pT_1_50.parquet.gzip')
FDDf = pd.read_parquet('../../AnalysisNonPromptDpp2017/Dplus/MC/LHC18a4a2/FD_Dpluspp5TeV_pT_1_50.parquet.gzip')
DataDf = pd.read_parquet('../../AnalysisNonPromptDpp2017/Dplus/Data/LHC17pq/Data_Dpluspp5TeV_pT_1_50.parquet.gzip')

LegLabels = ['Background', r'Prompt D$^+$', r'Feed-down D$^+$']
OutputLabels = ['Bkg', 'Prompt', 'FD']

OutPutDir = '../../AnalysisNonPromptDpp2017/Dplus/MLoutput'

for (PtMin, PtMax) in zip(PtMins, PtMaxs):

    print('\nStarting ML analysis')

    DataDfPtSel = DataDf.query(f'{PtMin} < pt_cand < {PtMax}')
    BkgDfPtSel = DataDfPtSel.query('inv_mass < 1.82 or 1.92 < inv_mass < 2.00')
    PromptDfPtSel = PromptDf.query(f'{PtMin} < pt_cand < {PtMax}')
    FDDfPtSel = FDDf.query(f'{PtMin} < pt_cand < {PtMax}')

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

    TotDfPtSel = pd.concat([BkgDfPtSel.iloc[:nCandToKeep].copy(), PromptDfPtSel.iloc[:nCandToKeep].copy(),
                           FDDfPtSel.iloc[:nCandToKeep].copy()], sort=True)

    LabelsArray = [0 for _ in range(nCandToKeep)] + [1 for _ in range(nCandToKeep)] + [2 for _ in range(nCandToKeep)]

    TrainSet, TestSet, yTrain, yTest = train_test_split(TotDfPtSel, LabelsArray, test_size=0.5, random_state=42)
    TrainTestData = [TrainSet, yTrain, TestSet, yTest]

    CandTypeFlags = pd.Series(yTest)
    PromptDfPtSelForEff = pd.concat([PromptDfPtSel.iloc[nCandToKeep:].copy(), TestSet[CandTypeFlags.values == 1]],
                                    sort=False)
    FDDfPtSelForEff = pd.concat([FDDfPtSel.iloc[nCandToKeep:].copy(), TestSet[CandTypeFlags.values == 2]], sort=False)

    #remove pT and inv-mass from training columns --> training columns to be set in a config yaml file
    TrainCols = list(TotDfPtSel.columns)
    TrainCols.remove('inv_mass')
    TrainCols.remove('pt_cand')

    # training, testing, and model application
    #_____________________________________________
    modelClf = xgb.XGBClassifier()
    ModelHandl = ModelHandler(modelClf, TrainCols)

    # hyperparams optimization --> not working with multi-class classification at the moment
    #HypRanges = {
    #    # # defines the maximum depth of a single tree (regularization)
    #    'max_depth': (1, 30),
    #    'learning_rate': (0.01, 0.3),  # learning rate
    #    'n_estimators': (50, 1000)  # number of boosting trees
    #}
    #ModelHandl.optimize_params_bayes(TrainTestData, HypRanges, None)

    # TODO: change hard coded hyperparameters!
    HypPars = {'max_depth':3, 'learning_rate':0.12, 'n_estimators':250, 'min_child_weight':5, 'colsample':0.9}
    ModelHandl.set_model_params(HypPars)

    # train and test the model with the updated hyperparameters
    ModelHandl.train_test_model(TrainTestData)
    yPredTest = ModelHandl.predict(TrainTestData[2])

    # apply model
    print('Applying ML model to prompt dataframe: ...', end='\r')
    yPredPromptEff = ModelHandl.predict(PromptDfPtSelForEff)
    PromptDfPtSelForEff = PromptDfPtSelForEff.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        PromptDfPtSelForEff[f'ML_output_{Lab}'] = yPredPromptEff[:, Pred]
    PromptDfPtSelForEff.to_parquet(f'{OutPutDir}/Prompt_Dpluspp5TeV_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to prompt dataframe: Done!')

    print('Applying ML model to FD dataframe: ...', end='\r')
    yPredFDEff = ModelHandl.predict(FDDfPtSelForEff)
    FDDfPtSelForEff = FDDfPtSelForEff.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        FDDfPtSelForEff[f'ML_output_{Lab}'] = yPredFDEff[:, Pred]
    FDDfPtSelForEff.to_parquet(f'{OutPutDir}/FD_Dpluspp5TeV_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to FD dataframe: Done!')

    print('Applying ML model to data dataframe: ...', end='\r')
    yPredData = ModelHandl.predict(DataDfPtSel)
    DataDfPtSel = DataDfPtSel.loc[:, ['inv_mass', 'pt_cand']]
    for Pred, Lab in enumerate(OutputLabels):
        DataDfPtSel[f'ML_output_{Lab}'] = yPredData[:, Pred]
    DataDfPtSel.to_parquet(f'{OutPutDir}/Data_Dpluspp5TeV_pT_{PtMin}_{PtMax}_ModelApplied.parquet.gzip')
    print('Applying ML model to data dataframe: Done!')


    # TODO: add possibility to save models instead or with the model application

    # plots
    #_____________________________________________
    VarsToDraw = list(TotDfPtSel.columns)

    #_____________________________________________
    FeatDistrFig = plot_utils.plot_distr([BkgDfPtSel, PromptDfPtSel, FDDfPtSel], VarsToDraw, (12, 7), 
                                         100, True, LegLabels)
    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.99, top=0.96, hspace=0.55, wspace=0.55)
    plt.savefig(f'{OutPutDir}/DistributionsAll_pT_{PtMin}_{PtMax}.pdf')
    plt.close('all')

    #_____________________________________________
    CorrMatrixFig = plot_utils.plot_corr([BkgDfPtSel, PromptDfPtSel, FDDfPtSel], VarsToDraw, LegLabels)
    for Fig, Lab in zip(CorrMatrixFig, OutputLabels):
        plt.figure(Fig.number)
        plt.subplots_adjust(left=0.2, bottom=0.25, right=0.95, top=0.9)
        Fig.savefig(f'{OutPutDir}/CorrMatrix{Lab}_pT_{PtMin}_{PtMax}.pdf')

    #_____________________________________________
    plt.rcParams["figure.figsize"] = (10, 7)
    MLOutputFig = plot_utils.plot_output_train_test(ModelHandl, TrainTestData, 80, True, LegLabels)
    for Fig, Lab in zip(MLOutputFig, OutputLabels):
        Fig.savefig(f'{OutPutDir}/MLOutputDistr{Lab}_pT_{PtMin}_{PtMax}.pdf')

    #_____________________________________________
    plt.rcParams["figure.figsize"] = (8, 7)
    ROCCurveFig = plot_utils.plot_roc(TrainTestData[3], yPredTest, LegLabels)
    ROCCurveFig.savefig(f'{OutPutDir}/ROCCurveAll_pT_{PtMin}_{PtMax}.pdf')

    #_____________________________________________
    PrecisionRecallFig = plot_utils.plot_precision_recall(TrainTestData[3], yPredTest, LegLabels)
    PrecisionRecallFig.savefig(f'{OutPutDir}/PrecisionRecallAll_pT_{PtMin}_{PtMax}.pdf')

    #_____________________________________________
    plt.rcParams["figure.figsize"] = (12, 7)
    FeaturesImportanceFig = plot_utils.plot_feature_imp(TrainSet[TrainCols], yTrain, ModelHandl)
    for iFig, Fig in enumerate(FeaturesImportanceFig):
        if iFig < 3:
            Fig.savefig(f'{OutPutDir}/FeatureImportance{OutputLabels[iFig]}_pT_{PtMin}_{PtMax}.pdf')
        else:
            Fig.savefig(f'{OutPutDir}/FeatureImportanceAll_pT_{PtMin}_{PtMax}.pdf')

    # delete dataframes to release memory
    del DataDfPtSel, BkgDfPtSel, PromptDfPtSel, FDDfPtSel, TrainSet, TestSet, yTrain, yTest
