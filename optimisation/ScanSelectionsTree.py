'''
python script for the optimisation of the working point using TTrees or pandas dataframes as input
run: python ScanSelectionTree.py cfgFileName.yml outFileName.root
'''

import sys
import argparse
import time
import itertools
import numpy as np
import yaml
from root_numpy import fill_hist
from ROOT import TFile, TH1F, TH2F, TF1, TCanvas, TNtuple, TDirectoryFile  # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, kRainBow, kBlack, kFullCircle  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.AnalysisUtils import ComputeEfficiency, GetPromptFDFractionFc, GetExpectedBkgFromSideBands  #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import  GetExpectedBkgFromSideBandsImp, GetExpectedBkgFromMC, GetExpectedSignal  #pylint: disable=wrong-import-position,import-error
from utils.FitUtils import SingleGaus #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle  #pylint: disable=wrong-import-position,import-error
from utils.DfUtils import LoadDfFromRootOrParquet  #pylint: disable=wrong-import-position,import-error
from utils.ReadModel import ReadTAMU, ReadPHSD, ReadMCatsHQ, ReadCatania  #pylint: disable=wrong-import-position,import-error

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('outFileName', metavar='text', default='outFile.root',
                    help='output root file name')
parser.add_argument("--batch", help="suppress video output",
                    action="store_true")
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

# load dataframes from input files
dfPrompt = LoadDfFromRootOrParquet(inputCfg['infiles']['signal']['prompt']['filename'],
                                   inputCfg['infiles']['signal']['prompt']['dirname'],
                                   inputCfg['infiles']['signal']['prompt']['treename'])
dfFD = LoadDfFromRootOrParquet(inputCfg['infiles']['signal']['feeddown']['filename'],
                               inputCfg['infiles']['signal']['feeddown']['dirname'],
                               inputCfg['infiles']['signal']['feeddown']['treename'])
dfBkg = LoadDfFromRootOrParquet(inputCfg['infiles']['background']['filename'],
                                inputCfg['infiles']['background']['dirname'],
                                inputCfg['infiles']['background']['treename'])
if inputCfg['infiles']['secpeak']['prompt']['filename']:
    dfSecPeakPrompt = LoadDfFromRootOrParquet(inputCfg['infiles']['secpeak']['prompt']['filename'],
                                              inputCfg['infiles']['secpeak']['prompt']['dirname'],
                                              inputCfg['infiles']['secpeak']['prompt']['treename'])
else:
    dfSecPeakPrompt = None
if inputCfg['infiles']['secpeak']['feeddown']['filename']:
    dfSecPeakFD = LoadDfFromRootOrParquet(inputCfg['infiles']['secpeak']['feeddown']['filename'],
                                          inputCfg['infiles']['secpeak']['feeddown']['dirname'],
                                          inputCfg['infiles']['secpeak']['feeddown']['treename'])
else:
    dfSecPeakFD = None

# reshuffle bkg and take only a fraction of it, seed fixed for reproducibility
dfBkg = dfBkg.sample(frac=inputCfg['infiles']['background']['fractiontokeep'], random_state=42).reset_index(drop=True)

# load cut values to scan
ptMins = inputCfg['ptmin']
ptMaxs = inputCfg['ptmax']
if not isinstance(ptMins, list):
    ptMins = [ptMins]
if not isinstance(ptMaxs, list):
    ptMaxs = [ptMaxs]
ParCutsName = inputCfg['dfparametercuts']['name']
EnableParCuts = inputCfg['dfparametercuts']['enable']
if EnableParCuts:
    ParCutMins = inputCfg['dfparametercuts']['min']
    ParCutMaxs = inputCfg['dfparametercuts']['max']
    if not isinstance(ParCutMins, list):
        CutParMins = [ParCutMins]
    if not isinstance(ParCutMaxs, list):
        ParCutMaxs = [ParCutMaxs]
else:
    ParCutsName = 'Integral'
    ParCutMins, ParCutMaxs = [], []
    ParCutMins.append(-1.e10)
    ParCutMaxs.append(1.e10)

cutVars = inputCfg['cutvars']
cutRanges, upperLowerCuts, varNames = [], [], []
for var in cutVars:
    cutRanges.append(np.arange(cutVars[var]['min'], cutVars[var]['max'] +
                               cutVars[var]['step'] / 10, cutVars[var]['step']).tolist())
    if cutVars[var]['upperlowercut'] == 'Upper':
        upperLowerCuts.append('<')
    else:
        upperLowerCuts.append('>')
    varNames.append(var)

# load preselection efficiency
if inputCfg['infiles']['preseleff']['filename']:
    infilePreselEff = TFile.Open(inputCfg['infiles']['preseleff']['filename'])
    hPreselEffPrompt = infilePreselEff.Get(inputCfg['infiles']['preseleff']['prompthistoname'])
    hPreselEffFD = infilePreselEff.Get(inputCfg['infiles']['preseleff']['feeddownhistoname'])

# load acceptance
infileAcc = TFile.Open(inputCfg['infiles']['acceptance'])
hPtGenAcc = infileAcc.Get('hPtGenAcc')
hPtGenLimAcc = infileAcc.Get('hPtGenLimAcc')

# load cross sections
inFileCrossSec = TFile.Open(inputCfg['predictions']['crosssec']['filename'])
hCrossSecPrompt = inFileCrossSec.Get(inputCfg['predictions']['crosssec']['histonames']['prompt'])
hCrossSecFD = inFileCrossSec.Get(inputCfg['predictions']['crosssec']['histonames']['feeddown'])

# load RAA
RaaPrompt_config = inputCfg['predictions']['Raa']['prompt']
if not isinstance(RaaPrompt_config, float) and not isinstance(RaaPrompt_config, int):
    if not isinstance(RaaPrompt_config, str):
        print('ERROR: RAA must be at least a string or a number. Exit')
        sys.exit()
    else:
        Raa_model_name = inputCfg['predictions']['Raa']['model']
        if Raa_model_name not in ['phsd', 'Catania', 'tamu', 'MCatsHQ']:
            print('ERROR: wrong model name, please check the list of avaliable models. Exit')
            sys.exit()
        else:
            if Raa_model_name == 'phsd':
                RaaPromptSpline, _, ptMinRaaPrompt, ptMaxRaaPrompt = ReadPHSD(RaaPrompt_config)
            elif Raa_model_name == 'Catania':
                RaaPromptSpline, _, ptMinRaaPrompt, ptMaxRaaPrompt = ReadCatania(RaaPrompt_config)
            elif Raa_model_name == 'MCatsHQ':
                RaaPromptSpline, _, ptMinRaaPrompt, ptMaxRaaPrompt = ReadMCatsHQ(RaaPrompt_config)
            elif Raa_model_name == 'tamu':
                RaaPromptSpline, _, ptMinRaaPrompt, ptMaxRaaPrompt = ReadTAMU(RaaPrompt_config)
else:
    RaaPrompt = RaaPrompt_config

RaaFD_config = inputCfg['predictions']['Raa']['feeddown']
if not isinstance(RaaFD_config, float) and not isinstance(RaaFD_config, int):
    if not isinstance(RaaFD_config, str):
        print('ERROR: RAA must be at least a string or a number. Exit')
        sys.exit()
    else:
        Raa_model_name = inputCfg['predictions']['Raa']['model']
        if Raa_model_name not in ['phsd', 'Catania', 'tamu', 'MCatsHQ']:
            print('ERROR: wrong model name, please check the list of avaliable models. Exit')
            sys.exit()
        else:
            if Raa_model_name == 'phsd':
                RaaFDSpline, _, ptMinRaaFD, ptMaxRaaFD = ReadPHSD(RaaFD_config)
            elif Raa_model_name == 'Catania':
                RaaFDSpline, _, ptMinRaaFD, ptMaxRaaFD = ReadCatania(RaaFD_config)
            elif Raa_model_name == 'MCatsHQ':
                RaaFDSpline, _, ptMinRaaFD, ptMaxRaaFD = ReadMCatsHQ(RaaFD_config)
            elif Raa_model_name == 'tamu':
                RaaFDSpline, _, ptMinRaaFD, ptMaxRaaFD = ReadTAMU(RaaFD_config)
else:
    RaaFD = RaaFD_config

# load constant terms
nExpEv = inputCfg['nExpectedEvents']
Taa = inputCfg['Taa']
sigmaMB = inputCfg['sigmaMB']

# set batch mode if enabled
if args.batch:
    gROOT.SetBatch(True)
    gROOT.ProcessLine("gErrorIgnoreLevel = kFatal;")

# output file with TNtuple
outFile = TFile(args.outFileName, 'recreate')
outDirFitSB = TDirectoryFile('SBfits', 'SBfits')
outDirFitSB.Write()
outDirFitSBPt = []

outFile.cd()
outDirPlots = TDirectoryFile('plots', 'plots')
outDirPlots.Write()
outDirPlotsPt = []

estNames = {'Signif': 'expected significance', 'SoverB': 'S/B', 'S': 'expected signal', 'B': 'expected background',
            'EffAccPrompt': '(Acc#times#font[152]{e})_{prompt}', 'EffAccFD': '(Acc#times#font[152]{e})_{FD}',
            'fPrompt': '#it{f}_{ prompt}^{ fc}', 'fFD': '#it{f}_{ FD}^{ fc}'}

varsName4Tuple = (':'.join(cutVars) + ':PtMin:PtMax:ParCutMin:ParCutMax:EffAccPromptError:EffAccFDError:SError:BError'
                  ':SignifError:SoverBError:' + ':'.join(estNames.keys()))
tSignif = TNtuple('tSignif', 'tSignif', varsName4Tuple)

totSets = 1
for cutRange in cutRanges:
    totSets *= len(cutRange)
print(f'Total number of sets per pT bin: {totSets}')

SetGlobalStyle(padleftmargin=0.12, padrightmargin=0.2, padbottommargin=0.15, padtopmargin=0.075,
               titleoffset=1., palette=kRainBow, titlesize=0.06, labelsize=0.055, maxdigits=4)

cSignifVsRest, hSignifVsRest, cEstimVsCut, hEstimVsCut = [], [], [], []
counter = 0
for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    outDirFitSB.cd()
    outDirFitSBPt.append(TDirectoryFile(f'pT{ptMin}-{ptMax}', f'pT{ptMin}-{ptMax}'))
    outDirFitSBPt[iPt].Write()
    outDirPlots.cd()
    outDirPlotsPt.append(TDirectoryFile(f'pT{ptMin}-{ptMax}', f'pT{ptMin}-{ptMax}'))
    outDirPlotsPt[iPt].Write()
    dfPromptPt = dfPrompt.query(f'{ptMin} < pt_cand < {ptMax}')
    dfFDPt = dfFD.query(f'{ptMin} < pt_cand < {ptMax}')
    dfBkgPt = dfBkg.query(f'{ptMin} < pt_cand < {ptMax}')

    # Raa
    ptCent = (ptMax + ptMin) / 2.
    if isinstance(RaaPrompt_config, str):
        if ptMinRaaPrompt < ptCent < ptMaxRaaPrompt:
            RaaPrompt = RaaPromptSpline['yCent'](ptCent)
        elif ptCent > ptMaxRaaPrompt:
            RaaPrompt = RaaPromptSpline['yCent'](ptMaxRaaPrompt)
        else:
            RaaPrompt = RaaPromptSpline['yCent'](ptMinRaaPrompt)
    if isinstance(RaaFD_config, str):
        if ptMinRaaFD < ptCent < ptMaxRaaFD:
            RaaFD = RaaFDSpline['yCent'](ptCent)
        elif ptCent > ptMaxRaaFD:
            RaaFD = RaaFDSpline['yCent'](ptMaxRaaFD)
        else:
            RaaFD = RaaFDSpline['yCent'](ptMinRaaFD)

    # denominator for efficiency
    nTotPrompt = len(dfPromptPt)
    nTotFD = len(dfFDPt)

    # preselection efficiency (if input provided)
    if inputCfg['infiles']['preseleff']['filename']:
        ptBinPreselEff = hPreselEffPrompt.GetXaxis().FindBin(ptMin*1.0001)
        preselEffPrompt = hPreselEffPrompt.GetBinContent(ptBinPreselEff)
        preselEffFD = hPreselEffFD.GetBinContent(ptBinPreselEff)
        preselEffPromptUnc = hPreselEffPrompt.GetBinError(ptBinPreselEff)
        preselEffFDUnc = hPreselEffFD.GetBinError(ptBinPreselEff)
    else:
        preselEffPrompt = 1.
        preselEffFD = 1.

    # acceptance
    ptBinAccMin = hPtGenAcc.GetXaxis().FindBin(ptMin*1.0001)
    ptBinAccMax = hPtGenAcc.GetXaxis().FindBin(ptMax*0.9999)
    numAcc = hPtGenAcc.Integral(ptBinAccMin, ptBinAccMax)
    denAcc = hPtGenLimAcc.Integral(ptBinAccMin, ptBinAccMax)
    acc, accUnc = ComputeEfficiency(numAcc, denAcc, np.sqrt(numAcc), np.sqrt(denAcc))

    # cross section from theory
    ptBinCrossSecMin = hCrossSecPrompt.GetXaxis().FindBin(ptMin*1.0001)
    ptBinCrossSecMax = hCrossSecPrompt.GetXaxis().FindBin(ptMax*0.9999)
    crossSecPrompt = hCrossSecPrompt.Integral(ptBinCrossSecMin, ptBinCrossSecMax, 'width') / (ptMax - ptMin)
    crossSecFD = hCrossSecFD.Integral(ptBinCrossSecMin, ptBinCrossSecMax, 'width') / (ptMax - ptMin)

    # signal histograms
    hMassSignal = TH1F(f'hMassSignal_pT{ptMin}-{ptMax}', ';#it{M} (GeV/#it{c});Counts', 400,
                       min(dfPromptPt['inv_mass']), max(dfPromptPt['inv_mass']))
    fill_hist(hMassSignal, np.concatenate((dfPromptPt['inv_mass'].values, dfFDPt['inv_mass'].values)))
    funcSignal = TF1('funcSignal', SingleGaus, 1.6, 2.2, 3)
    funcSignal.SetParameters(hMassSignal.Integral('width'), hMassSignal.GetMean(), hMassSignal.GetRMS())
    hMassSignal.Fit('funcSignal', 'Q0')
    mean = funcSignal.GetParameter(1)
    sigma = funcSignal.GetParameter(2)
    # SecPeak
    meanSecPeak = inputCfg['infiles']['secpeak']['mean']
    sigmaSecPeak = inputCfg['infiles']['secpeak']['sigma']
    if dfSecPeakPrompt and dfSecPeakFD:
        hMassSecPeak = TH1F(f'hMassSignal_pT{ptMin}-{ptMax}', ';#it{M} (GeV/#it{c});Counts', 400,
                            min(dfSecPeakPrompt['inv_mass']), max(dfSecPeakPrompt['inv_mass']))
        fill_hist(hMassSecPeak, np.concatenate((dfSecPeakPrompt['inv_mass'].values, dfSecPeakFD['inv_mass'].values)))
        funcSignal.SetParameters(hMassSecPeak.Integral('width'), hMassSecPeak.GetMean(), hMassSecPeak.GetRMS())
        hMassSecPeak.Fit('funcSignal', 'Q0')
        meanSecPeak = funcSignal.GetParameter(1)
        sigmaSecPeak = funcSignal.GetParameter(2)

    # output histos
    hSignifVsRest.append(dict())
    hEstimVsCut.append(dict())
    for iParCut, (ParCutMin, ParCutMax) in enumerate(zip(ParCutMins, ParCutMaxs)):
        if len(varNames) == 1:
            for est in estNames:
                minVar = cutVars[varNames[0]]['min'] - cutVars[varNames[0]]['step'] / 2
                maxVar = cutVars[varNames[0]]['max'] + cutVars[varNames[0]]['step'] / 2
                nBinsVar = int((maxVar - minVar) / cutVars[varNames[0]]['step'])
                hEstimVsCut[iPt][est] = TH1F(f'h{est}VsCut_pT{ptMin}-{ptMax}_{ParCutsName}{ParCutMin}-{ParCutMax}',
                                             f';{varNames[0]};{estNames[est]}', nBinsVar, minVar, maxVar)
                SetObjectStyle(hEstimVsCut[iPt][est], color=kBlack, marker=kFullCircle, linewidth=1)
        elif len(varNames) == 2:
            for est in estNames:
                minVar0 = cutVars[varNames[0]]['min'] - cutVars[varNames[0]]['step'] / 2
                minVar1 = cutVars[varNames[1]]['min'] - cutVars[varNames[1]]['step'] / 2
                maxVar0 = cutVars[varNames[0]]['max'] + cutVars[varNames[0]]['step'] / 2
                maxVar1 = cutVars[varNames[1]]['max'] + cutVars[varNames[1]]['step'] / 2
                nBinsVar0 = int((maxVar0 - minVar0) / cutVars[varNames[0]]['step'])
                nBinsVar1 = int((maxVar1 - minVar1) / cutVars[varNames[1]]['step'])
                hEstimVsCut[iPt][est] = TH2F(f'h{est}VsCut_pT{ptMin}-{ptMax}_{ParCutsName}{ParCutMin}-{ParCutMax}',
                                             f';{varNames[0]};{varNames[1]};{estNames[est]}',
                                             nBinsVar0, minVar0, maxVar0, nBinsVar1, minVar1, maxVar1)
        startTime = time.time()
        for iSet, cutSet in enumerate(itertools.product(*cutRanges)):
            for iCut, (cut, upperLower, varName) in enumerate(zip(cutSet, upperLowerCuts, varNames)):
                if iCut == 0:
                    selToApply = f'{varName}{upperLower}{cut}'
                else:
                    selToApply += f' & {varName}{upperLower}{cut}'
            if ParCutsName and EnableParCuts:
                selToApply += f' & {ParCutMin} < {ParCutsName} < {ParCutMax}'

            if (iSet+1) % 100 == 0:
                print(f'Time elapsed to test up to cut set number {iSet+1}: {time.time()-startTime:.2f}s', end='\r')

            dfPromptPtSel = dfPromptPt.query(selToApply)
            dfFDPtSel = dfFDPt.query(selToApply)
            dfBkgPtSel = dfBkgPt.query(selToApply)
            effPrompt, effPromptUnc = ComputeEfficiency(len(dfPromptPtSel), nTotPrompt,
                                                        np.sqrt(len(dfPromptPtSel)), np.sqrt(nTotPrompt))
            effFD, effFDUnc = ComputeEfficiency(len(dfFDPtSel), nTotFD, np.sqrt(len(dfFDPtSel)), np.sqrt(nTotFD))
            effTimesAccPrompt = effPrompt * preselEffPrompt * acc
            effTimesAccFD = effFD * preselEffFD * acc
            fPrompt, fFD = GetPromptFDFractionFc(effTimesAccPrompt, effTimesAccFD,
                                                 crossSecPrompt, crossSecFD, RaaPrompt, RaaFD)
            hMassBkg = TH1F(f'hMassBkg_pT{ptMin}-{ptMax}_cutSet{iSet}', ';#it{M} (GeV/#it{c});Counts', 200,
                            min(dfBkgPtSel['inv_mass']), max(dfBkgPtSel['inv_mass']))
            fill_hist(hMassBkg, dfBkgPtSel['inv_mass'].values)

            # expected signal, BR already included in cross section
            if inputCfg['expectedSignalFrom'] == 'prompt':
                expSignal = GetExpectedSignal(crossSecPrompt, ptMax-ptMin, 1., effTimesAccPrompt,
                                              fPrompt[0], 1., 1., nExpEv, sigmaMB, Taa, RaaPrompt)
            elif inputCfg['expectedSignalFrom'] == 'feeddown':
                expSignal = GetExpectedSignal(crossSecFD, ptMax-ptMin, 1., effTimesAccFD,
                                              fFD[0], 1., 1., nExpEv, sigmaMB, Taa, RaaFD)

            # expected background
            bkgConfig = inputCfg['infiles']['background']
            outDirFitSBPt[iPt].cd()
            expBkg = 0.
            errExpBkg = 0.
            if bkgConfig['isMC']:
                expBkg, errExpBkg, hMassBkg = GetExpectedBkgFromMC(hMassBkg, mean, sigma)
            elif bkgConfig['impFit']:
                expBkg, errExpBkg, hMassBkg = GetExpectedBkgFromSideBandsImp(hMassBkg, bkgConfig['fitFunc'],
                                                                             bkgConfig['nSigma'], mean, sigma,
                                                                             meanSecPeak, sigmaSecPeak)
            else:
                expBkg, errExpBkg, hMassBkg = GetExpectedBkgFromSideBands(hMassBkg, bkgConfig['fitFunc'],
                                                                          bkgConfig['nSigma'], mean, sigma,
                                                                          meanSecPeak, sigmaSecPeak)
            hMassBkg.Write()
            expBkg *= nExpEv / bkgConfig['nEvents'] / bkgConfig['fractiontokeep']
            errExpBkg *= nExpEv / bkgConfig['nEvents'] / bkgConfig['fractiontokeep']

            if inputCfg['infiles']['background']['corrfactor']['filename']:
                inFile = TFile.Open(inputCfg['infiles']['background']['corrfactor']['filename'])
                hBkgCorrFactor = inFile.Get(inputCfg['infiles']['background']['corrfactor']['histoname'])
                expBkg *= hBkgCorrFactor.GetBinContent(hBkgCorrFactor.FindBin(ptCent))
                errExpBkg *= hBkgCorrFactor.GetBinContent(hBkgCorrFactor.FindBin(ptCent))

            # S/B and significance
            expSoverB = 0.
            expSignif = 0.
            errS = 0. # TODO: think how to define a meaningful error on the estimated signal and propagate it
            errSoverB = 0.
            errSignif = 0.
            if expBkg > 0:
                expSoverB = expSignal / expBkg
                expSignif = expSignal / np.sqrt(expSignal + expBkg)
                errSoverB = expSoverB * errExpBkg / expBkg
                errSignif = expSignif * 0.5 * errExpBkg / (expSignal + expBkg)

            # Efficiency
            EffAccFDError = np.sqrt((effFDUnc/effFD)**2
                                    + (preselEffFDUnc/preselEffFD)**2
                                    + (accUnc/acc)**2)*effTimesAccFD
            EffAccPromptError = np.sqrt((effPromptUnc/effPrompt)**2
                                        + (preselEffPromptUnc/preselEffPrompt)**2
                                        + (accUnc/acc)**2)*effTimesAccPrompt

            tupleForNtuple = cutSet + (ptMin, ptMax, ParCutMin, ParCutMax, EffAccPromptError, EffAccFDError,
                                       errS, errExpBkg, errSignif, errSoverB, expSignif, expSoverB, expSignal, expBkg,
                                       effTimesAccPrompt, effTimesAccFD, fPrompt[0], fFD[0])
            tSignif.Fill(np.array(tupleForNtuple, 'f'))
            estValues = {'Signif': expSignif, 'SoverB': expSoverB, 'S': expSignal, 'B': expBkg,
                         'EffAccPrompt': effTimesAccPrompt, 'EffAccFD': effTimesAccFD,
                         'fPrompt': fPrompt[0], 'fFD': fFD[0]}
            estValuesErr = {'SignifError': errSignif, 'SoverBError': errSoverB, 'SError': errS, 'BError': errExpBkg,
                            'EffAccPromptError': EffAccPromptError, 'EffAccFDError': EffAccFDError}
            if len(varNames) == 1:
                binVar = hEstimVsCut[iPt]['Signif'].GetXaxis().FindBin(cutSet[0])
                for est in estValues:
                    hEstimVsCut[iPt][est].SetBinContent(binVar, estValues[est])
                    if f'{est}Error' in estValuesErr:
                        hEstimVsCut[iPt][est].SetBinError(binVar, estValuesErr[f'{est}Error'])
            if len(varNames) == 2:
                binVar0 = hEstimVsCut[iPt]['Signif'].GetXaxis().FindBin(cutSet[0])
                binVar1 = hEstimVsCut[iPt]['Signif'].GetYaxis().FindBin(cutSet[1])
                for est in estValues:
                    hEstimVsCut[iPt][est].SetBinContent(binVar0, binVar1, estValues[est])

        if ParCutsName != 'Integral':
            print(f'Time elapsed to test cut sets for pT bin {ptMin}-{ptMax} '
                  f'and {ParCutsName} bin {ParCutMin}-{ParCutMax}: {time.time()-startTime:.2f}s')
        else:
            print(f'Time elapsed to test cut sets for pT bin {ptMin}-{ptMax}: {time.time()-startTime:.2f}s')
        # plots
        outDirPlotsPt[iPt].mkdir(f'{ParCutsName}{ParCutMin}-{ParCutMax}')
        cSignifVsRest.append(TCanvas(f'cSignifVsRest_pT{ptMin}-{ptMax}_{ParCutsName}{ParCutMin}-{ParCutMax}',
                                     '', 800, 1000))
        cSignifVsRest[counter].Divide(2, 4)
        for iPad, est in enumerate(estNames):
            if est != 'Signif':
                hFrame = cSignifVsRest[counter].cd(iPad).DrawFrame(tSignif.GetMinimum(est)*0.8,
                                                                   tSignif.GetMinimum('Signif')*0.8,
                                                                   tSignif.GetMaximum(est)*1.2,
                                                                   tSignif.GetMaximum('Signif')*1.2,
                                                                   f";{estNames[est]};{estNames['Signif']}")
                hFrame.GetXaxis().SetDecimals()
                hFrame.GetYaxis().SetDecimals()
                hSignifVsRest[iPt][est] = (TH2F((f'hSignifVs{est}_pT{ptMin}-{ptMax}_{ParCutsName}'
                                                 f'{ParCutMin}-{ParCutMax}'),
                                                f";{estNames[est]};{estNames['Signif']}", 50,
                                                tSignif.GetMinimum(est)*0.8, tSignif.GetMaximum(est)*1.2, 50,
                                                tSignif.GetMinimum('Signif')*0.8, tSignif.GetMaximum('Signif')*1.))
                tSignif.Draw(f'Signif:{est}>>hSignifVs{est}_pT{ptMin}-{ptMax}_{ParCutsName}{ParCutMin}-{ParCutMax}',
                             f'PtMin == {ptMin} && PtMax == {ptMax}', 'colz same')
                cSignifVsRest[counter].Update()
                cSignifVsRest[counter].Modified()
                outDirPlotsPt[iPt].cd(f'{ParCutsName}{ParCutMin}-{ParCutMax}')
                hSignifVsRest[iPt][est].Write()
        outDirPlotsPt[iPt].cd(f'{ParCutsName}{ParCutMin}-{ParCutMax}')
        cSignifVsRest[counter].Write()
        if 1 <= len(varNames) <= 2:
            if len(varNames) == 1:
                cEstimVsCut.append(TCanvas(
                    f'cEstimVsCut_pT{ptMin}-{ptMax}_{ParCutsName}{ParCutMin}-{ParCutMax}', '', 800, 1000))
                cEstimVsCut[counter].Divide(2, 4)
                for iPad, est in enumerate(hEstimVsCut[iPt]):
                    hFrame = cEstimVsCut[counter].cd(iPad+1).DrawFrame(minVar, tSignif.GetMinimum(est)*0.8,
                                                                       maxVar, tSignif.GetMaximum(est)*1.2,
                                                                       f';{varNames[0]};{estNames[est]}')
                    if 'Eff' in est:
                        cEstimVsCut[counter].cd(iPad+1).SetLogy()
                        hFrame.GetYaxis().SetMoreLogLabels()
                    hFrame.GetXaxis().SetNdivisions(505)
                    hFrame.GetXaxis().SetDecimals()
                    hFrame.GetYaxis().SetDecimals()
                    hEstimVsCut[iPt][est].DrawCopy('psame')
                    outDirPlotsPt[iPt].cd(f'{ParCutsName}{ParCutMin}-{ParCutMax}')
                    hEstimVsCut[iPt][est].Write()
            elif len(varNames) == 2:
                cEstimVsCut.append(TCanvas(f'cEstimVsCut_pT{ptMin}-{ptMax}_{ParCutsName}{ParCutMin}-{ParCutMax}',
                                           '', 800, 1000))
                cEstimVsCut[counter].Divide(2, 4)
                for iPad, est in enumerate(hEstimVsCut[iPt]):
                    minVar0 = cutVars[varNames[0]]['min'] - cutVars[varNames[0]]['step'] / 2
                    minVar1 = cutVars[varNames[1]]['min'] - cutVars[varNames[1]]['step'] / 2
                    maxVar0 = cutVars[varNames[0]]['max'] + cutVars[varNames[0]]['step'] / 2
                    maxVar1 = cutVars[varNames[1]]['max'] + cutVars[varNames[1]]['step'] / 2
                    hFrame = cEstimVsCut[counter].cd(iPad+1).DrawFrame(minVar0, minVar1, maxVar0, maxVar1,
                                                                       f';{varNames[0]};{varNames[1]};{estNames[est]}')
                    if 'Eff' in est:
                        cEstimVsCut[counter].cd(iPad+1).SetLogz()
                        hFrame.GetZaxis().SetMoreLogLabels()
                    hFrame.GetXaxis().SetNdivisions(505)
                    hFrame.GetYaxis().SetNdivisions(505)
                    hFrame.GetXaxis().SetDecimals()
                    hFrame.GetYaxis().SetDecimals()
                    hEstimVsCut[iPt][est].DrawCopy('colzsame')
                    outDirPlotsPt[iPt].cd(f'{ParCutsName}{ParCutMin}-{ParCutMax}')
                    hEstimVsCut[iPt][est].Write()
            cEstimVsCut[counter].Update()
            cEstimVsCut[counter].Modified()
            outDirPlotsPt[iPt].cd(f'{ParCutsName}{ParCutMin}-{ParCutMax}')
            cEstimVsCut[counter].Write()
        counter += 1
outFile.cd()
tSignif.Write()
outFile.Close()

if not args.batch:
    input('Press enter to exit')
