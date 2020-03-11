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
from ROOT import TFile, TH1F, TH2F, TCanvas, TNtuple, TDirectoryFile, TAxis, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, kRainBow, kBlack, kBlue, kRed, kFullCircle  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
# pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.AnalysisUtils import ComputeEfficiency, GetPromptFDFractionFc, GetExpectedBkgFromSideBands, \
    GetExpectedBkgFromMC, GetExpectedSignal
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.DfUtils import LoadDfFromRootOrParquet
from utils.ReadModel import InterpolateModel, ReadFONLL, ReadGMVFNS, ReadKtFact, ReadTAMU, ReadPHSD, \
    ReadMCatsHQ, ReadCatania

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('outFileName', metavar='text', default='outFile.root',
                    help='output root file name')
parser.add_argument("--batch", help="suppress video output", action="store_true")
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
    dfSecPeakPrompt = LoadDfFromRootOrParquet(inputCfg['infiles']['signal']['prompt']['filename'],
                                              inputCfg['infiles']['signal']['prompt']['dirname'],
                                              inputCfg['infiles']['signal']['prompt']['treename'])
else:
    dfSecPeakPrompt = None
if inputCfg['infiles']['secpeak']['feeddown']['filename']:
    dfSecPeakFD = LoadDfFromRootOrParquet(inputCfg['infiles']['signal']['feeddown']['filename'],
                                          inputCfg['infiles']['signal']['feeddown']['dirname'],
                                          inputCfg['infiles']['signal']['feeddown']['treename'])
else:
    dfSecPeakFD = None

# reshuffle bkg and take only a fraction of it
dfBkg = dfBkg.sample(frac=inputCfg['infiles']['background']['fractiontokeep']).reset_index(drop=True)

# load cut values to scan
ptMins = inputCfg['ptmin']
ptMaxs = inputCfg['ptmax']
if not isinstance(ptMins, list):
    ptMins = [ptMins]
if not isinstance(ptMaxs, list):
    ptMaxs = [ptMaxs]
cutVars = inputCfg['cutvars']
cutRanges, upperLowerCuts, varNames = [], [], []
for _, var in enumerate(cutVars):
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
# TODO: add possibility to pass prediction from file instead of fixed value
RaaPrompt_config = inputCfg['predictions']['Raa']['prompt']
if not isinstance(RaaPrompt_config, float) and not isinstance(RaaPrompt_config, int):
    if not isinstance(RaaPrompt_config, str):
        print('ERROR: RAA must be at least a string or a number. Exit')
        exit()
    else:
        Raa_model_name = inputCfg['predictions']['Raa']['model']
        if Raa_model_name not in ['phsd', 'KtFact', 'GMVFNS', 'Catania', 'fonll', 'tamu']:
            print('ERROR: wrong model name, please check the list of avaliable models. Exit')
            exit()
        else:
            if Raa_model_name == 'phsd':
                RaaPromptSpline, _ = ReadPHSD(RaaPrompt_config)
            elif Raa_model_name == 'KtFact':
                RaaPromptSpline, _ = ReadKtFact(RaaPrompt_config)
            elif Raa_model_name == 'GMVFNS':
                RaaPromptSpline, _ = ReadGMVFNS(RaaPrompt_config)
            elif Raa_model_name == 'Catania':
                RaaPromptSpline, _ = ReadCatania(RaaPrompt_config)
            elif Raa_model_name == 'fonll':
                RaaPromptSpline, _ = ReadFONLL(RaaPrompt_config)
            elif Raa_model_name == 'tamu':
                RaaPromptSpline, _ = ReadTAMU(RaaPrompt_config)

else:
    RaaPrompt = RaaPrompt_config

RaaFD_config = inputCfg['predictions']['Raa']['feeddown']
if not isinstance(RaaFD_config, float) and not isinstance(RaaFD_config, int):
    if not isinstance(RaaFD_config, str):
        print('ERROR: RAA must be at least a string or a number. Exit')
        exit()
    else:
        Raa_model_name = inputCfg['predictions']['Raa']['model']
        if Raa_model_name not in ['phsd', 'KtFact', 'GMVFNS', 'Catania', 'fonll', 'tamu', 'MCatsHQ']:
            print('ERROR: wrong model name, please check the list of avaliable models. Exit')
            exit()
        else:
            if Raa_model_name == 'phsd':
                RaaFDSpline, _ = ReadPHSD(RaaFD_config)
            elif Raa_model_name == 'KtFact':
                RaaFDSpline, _ = ReadKtFact(RaaFD_config)
            elif Raa_model_name == 'GMVFNS':
                RaaFDSpline, _ = ReadGMVFNS(RaaFD_config)
            elif Raa_model_name == 'Catania':
                RaaFDSpline, _ = ReadCatania(RaaFD_config)
            elif Raa_model_name == 'fonll':
                RaaFDSpline, _ = ReadFONLL(RaaFD_config)
            elif Raa_model_name == 'MCatsHQ':
                RaaFDSpline, _ = ReadMCatsHQ(RaaFD_config)
            elif Raa_model_name == 'tamu':
                RaaFDSpline, _ = ReadTAMU(RaaFD_config)

else:
    RaaFD = RaaFD_config

# load constant terms
nExpEv = inputCfg['nExpectedEvents']
Taa = inputCfg['Taa']
BR = inputCfg['BR']
sigmaMB = inputCfg['sigmaMB']

# set batch mode if enabled
if args.batch:
    gROOT.SetBatch(True)
    gROOT.ProcessLine("gErrorIgnoreLevel = kFatal;")

# output file with TNtuple
outFile = TFile(args.outFileName, 'recreate')
if not inputCfg['infiles']['background']['isMC']:
    outDirFitSB = TDirectoryFile('SBfits', 'SBfits')
    outDirFitSB.Write()
    outDirFitSBPt = []

outFile.cd()
outDirPlots = TDirectoryFile('plots', 'plots')
outDirPlots.Write()
outDirPlotsPt = []

estNames = {'Signif':'expected significance', 'SoverB':'S/B', 'EffAccPrompt':'(Acc#times#font[152]{e})_{prompt}',
            'EffAccFD':'(Acc#times#font[152]{e})_{FD}', 'fPrompt':'#it{f}_{ prompt}^{ fc}', 'fFD':'#it{f}_{ FD}^{ fc}'}

varsName4Tuple = ':'.join(cutVars) + ':PtMin:PtMax:S:B:' + ':'.join(estNames.keys())
tSignif = TNtuple('tSignif', 'tSignif', varsName4Tuple)

totSets = 1
for cutRange in cutRanges:
    totSets *= len(cutRange)
print(f'Total number of sets per pT bin: {totSets}')

SetGlobalStyle(padleftmargin=0.2, padrightmargin=0.2, padbottommargin=0.15,
               titleoffset=1.4, titleoffsety=1.8, palette=kRainBow)

cSignifVsRest, hSignifVsRest, cEstimVsCut, hEstimVsCut = [], [], [], []
for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    if not inputCfg['infiles']['background']['isMC']:
        outDirFitSB.cd()
        outDirFitSBPt.append(TDirectoryFile(f'pT{ptMin}-{ptMax}', f'pT{ptMin}-{ptMax}'))
        outDirFitSBPt[iPt].Write()
    outDirPlots.cd()
    outDirPlotsPt.append(TDirectoryFile(f'pT{ptMin}-{ptMax}', f'pT{ptMin}-{ptMax}'))
    outDirPlotsPt[iPt].Write()
    dfPromptPt = dfPrompt.query(f'{ptMin} < pt_cand < {ptMax}')
    dfFDPt = dfFD.query(f'{ptMin} < pt_cand < {ptMax}')
    dfBkgPt = dfBkg.query(f'{ptMin} < pt_cand < {ptMax}')
    ptCent = (ptMax-ptMin)/2.
    RaaPrompt = RaaPromptSpline['yCent'](ptCent)
    RaaFD = RaaFDSpline['yCent'](ptCent)
    # denominator for efficiency
    nTotPrompt = len(dfPromptPt)
    nTotFD = len(dfFDPt)
    # preselection efficiency (if input provided)
    if inputCfg['infiles']['preseleff']['filename']:
        ptBinPreselEff = hPreselEffPrompt.GetXaxis().FindBin(ptMin*1.0001)
        preselEffPrompt = hPreselEffPrompt.GetBinContent(ptBinPreselEff)
        preselEffFD = hPreselEffFD.GetBinContent(ptBinPreselEff)
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
    # output histos
    hSignifVsRest.append(dict())
    hEstimVsCut.append(dict())
    if len(varNames) == 1:
        for est in estNames:
            minVar = cutVars[varNames[0]]['min'] - cutVars[varNames[0]]['step'] / 2
            maxVar = cutVars[varNames[0]]['max'] + cutVars[varNames[0]]['step'] / 2
            nBinsVar = int((maxVar - minVar) / cutVars[varNames[0]]['step'])
            hEstimVsCut[iPt][est] = TH1F(f'h{est}VsCut_pT{ptMin}-{ptMax}', f';{varNames[0]};{estNames[est]}',
                                         nBinsVar, minVar, maxVar)
            SetObjectStyle(hEstimVsCut[iPt][est], color=kBlack, marker=kFullCircle)
    elif len(varNames) == 2:
        for est in estNames:
            minVar0 = cutVars[varNames[0]]['min'] - cutVars[varNames[0]]['step'] / 2
            minVar1 = cutVars[varNames[1]]['min'] - cutVars[varNames[0]]['step'] / 2
            maxVar0 = cutVars[varNames[0]]['max'] + cutVars[varNames[0]]['step'] / 2
            maxVar1 = cutVars[varNames[1]]['max'] + cutVars[varNames[0]]['step'] / 2
            nBinsVar0 = int((maxVar0 - minVar0) / cutVars[varNames[0]]['step'])
            nBinsVar1 = int((maxVar1 - minVar1) / cutVars[varNames[1]]['step'])
            hEstimVsCut[iPt][est] = TH2F(f'h{est}VsCut_pT{ptMin}-{ptMax}',
                                         f';{varNames[0]};{varNames[1]};{estNames[est]}',
                                         nBinsVar0, minVar0, maxVar0, nBinsVar1, minVar1, maxVar1)
    startTime = time.time()

    for iSet, cutSet in enumerate(itertools.product(*cutRanges)):
        selToApply = ''
        for iCut, (cut, upperLower, varName) in enumerate(zip(cutSet, upperLowerCuts, varNames)):
            if iCut == 0:
                selToApply = f'{varName}{upperLower}{cut}'
            else:
                selToApply += f' & {varName}{upperLower}{cut}'
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
        fPrompt, fFD = GetPromptFDFractionFc(
            effTimesAccPrompt, effTimesAccFD, crossSecPrompt, crossSecFD, RaaPrompt, RaaFD)
        hMassBkg = TH1F(f'hMassBkg_pT{ptMin}-{ptMax}_cutSet{iSet}', ';#it{M} (GeV/#it{c});Counts',
                        400, min(dfBkgPtSel['inv_mass']), max(dfBkgPtSel['inv_mass']))
        hMassSignal = TH1F(f'hMassSignal_pT{ptMin}-{ptMax}_cutSet{iSet}', ';#it{M} (GeV/#it{c});Counts', 400, min(
            dfPromptPtSel['inv_mass']), max(dfPromptPtSel['inv_mass']))

        fill_hist(hMassBkg, dfBkgPtSel['inv_mass'].values)
        fill_hist(hMassSignal, list(dfPromptPtSel['inv_mass'].values)+list(dfFDPtSel['inv_mass'].values))

        # SecPeak
        if dfSecPeakPrompt and dfSecPeakFD: 
            hMassSeaPeak = TH1F(f'hMassSignal_pT{ptMin}-{ptMax}_cutSet{iSet}', ';#it{M} (GeV/#it{c});Counts', 400, min(
                dfSecPeakPrompt['inv_mass']), max(dfSecPeakPrompt['inv_mass']))
            fill_hist(hMassSeaPeak, list(dfSecPeakPrompt['inv_mass'].values)+list(dfSecPeakFD['inv_mass'].values))
        else:
            hMassSeaPeak = None
        # expected signal
        expSignal = GetExpectedSignal(crossSecPrompt, ptMax-ptMin, 1., effTimesAccPrompt,
                                      fPrompt[0], 1., 1., nExpEv, sigmaMB, Taa, spline((ptMax-ptMin)/2))
        # BR already included in cross section

        # expected background
        if inputCfg['infiles']['background']['isMC']:
            expBkg = GetExpectedBkgFromMC(hMassBkg)
        else:
            if hMassSeaPeak:
                expBkg, hMassSB = GetExpectedBkgFromSideBands(hMassBkg, 'pol2', 4, hMassSignal, -1., -1., hMassSeaPeak)
            else:
                expBkg, hMassSB = GetExpectedBkgFromSideBands(hMassBkg, 'pol2', 4, hMassSignal)
            outDirFitSBPt[iPt].cd()
            hMassSB.Write()
        expBkg *= nExpEv / inputCfg['infiles']['background']['nEvents'] / \
            inputCfg['infiles']['background']['fractiontokeep']
        # S/B and significance
        if expBkg > 0:
            expSoverB = expSignal / expBkg
        else:
            expSoverB = -1.
        if expSignal + expBkg > 0:
            expSignif = expSignal / (np.sqrt(expSignal + expBkg))
        else:
            expSignif = -1.
        tupleForNtuple = cutSet + (ptMin, ptMax, expSignal, expBkg, expSignif,
                                   expSoverB, effTimesAccPrompt, effTimesAccFD, fPrompt[0], fFD[0])
        tSignif.Fill(np.array(tupleForNtuple, 'f'))
        estValues = {'Signif': expSignif, 'SoverB': expSoverB, 'EffAccPrompt': effTimesAccPrompt,
                     'EffAccFD': effTimesAccFD, 'fPrompt': fPrompt[0], 'fFD': fFD[0]}
        if len(varNames) == 1:
            binVar = hEstimVsCut[iPt]['Signif'].GetXaxis().FindBin(cutSet[0])
            for est in estValues:
                hEstimVsCut[iPt][est].SetBinContent(binVar, estValues[est])
        if len(varNames) == 2:
            binVar0 = hEstimVsCut[iPt]['Signif'].GetXaxis().FindBin(cutSet[0])
            binVar1 = hEstimVsCut[iPt]['Signif'].GetYaxis().FindBin(cutSet[1])
            for est in estValues:
                hEstimVsCut[iPt][est].SetBinContent(binVar0, binVar1, estValues[est])
    print(f'Time elapsed to test cut sets for pT bin {ptMin}-{ptMax}: {time.time()-startTime:.2f}')
    # plots
    cSignifVsRest.append(TCanvas(f'cSignifVsRest_pT{ptMin}-{ptMax}', '', 1920, 1080))
    cSignifVsRest[iPt].Divide(3, 2)
    for iPad, est in enumerate(estNames):
        if est != 'Signif':
            hFrame = cSignifVsRest[iPt].cd(iPad).DrawFrame(tSignif.GetMinimum(est)*0.8,
                                                           tSignif.GetMinimum('Signif')*0.8,
                                                           tSignif.GetMaximum(est)*1.2,
                                                           tSignif.GetMaximum('Signif')*1.2,
                                                           f";{estNames[est]};{estNames['Signif']}")
            hFrame.GetXaxis().SetDecimals()
            hFrame.GetYaxis().SetDecimals()
            hSignifVsRest[iPt][est] = (TH2F(f'hSignifVs{est}_pT{ptMin}-{ptMax}',
                                            f";{estNames[est]};{estNames['Signif']}", 50,
                                            tSignif.GetMinimum(est)*0.8, tSignif.GetMaximum(est)*1.2, 50,
                                            tSignif.GetMinimum('Signif')*0.8, tSignif.GetMaximum('Signif')*1.))
            tSignif.Draw(f'Signif:{est}>>hSignifVs{est}_pT{ptMin}-{ptMax}', f'PtMin == {ptMin} && PtMax == {ptMax}',
                         'colz same')
            cSignifVsRest[iPt].Update()
            cSignifVsRest[iPt].Modified()
            outDirPlotsPt[iPt].cd()
            hSignifVsRest[iPt][est].Write()
    outDirPlotsPt[iPt].cd()
    cSignifVsRest[iPt].Write()
    if 1 <= len(varNames) <= 2:
        if len(varNames) == 1:
            cEstimVsCut.append(TCanvas(f'cEstimVsCut_pT{ptMin}-{ptMax}', '', 1920, 1080))
            cEstimVsCut[iPt].Divide(3, 2)
            for iPad, est in enumerate(hEstimVsCut[iPt]):
                hFrame = cEstimVsCut[iPt].cd(iPad+1).DrawFrame(minVar, tSignif.GetMinimum(est)*0.8, 
                                                               maxVar, tSignif.GetMaximum(est)*1.2,
                                                               f';{varNames[0]};{estNames[est]}')
                if 'Eff' in est:
                    cEstimVsCut[iPt].cd(iPad+1).SetLogy()
                    hFrame.GetYaxis().SetMoreLogLabels()
                hFrame.GetXaxis().SetNdivisions(505)
                hFrame.GetXaxis().SetDecimals()
                hFrame.GetYaxis().SetDecimals()
                hEstimVsCut[iPt][est].DrawCopy('psame')
                outDirPlotsPt[iPt].cd()
                hEstimVsCut[iPt][est].Write()
        elif len(varNames) == 2:
            cEstimVsCut.append(TCanvas(f'cEstimVsCut_pT{ptMin}-{ptMax}', '', 1920, 1080))
            cEstimVsCut[iPt].Divide(3, 2)
            for iPad, est in enumerate(hEstimVsCut[iPt]):
                minVar0 = cutVars[varNames[0]]['min'] - cutVars[varNames[0]]['step'] / 2
                minVar1 = cutVars[varNames[1]]['min'] - cutVars[varNames[0]]['step'] / 2
                maxVar0 = cutVars[varNames[0]]['max'] + cutVars[varNames[0]]['step'] / 2
                maxVar1 = cutVars[varNames[1]]['max'] + cutVars[varNames[0]]['step'] / 2
                hFrame = cEstimVsCut[iPt].cd(iPad+1).DrawFrame(minVar0, minVar1, maxVar0, maxVar1,
                                                               f';{varNames[0]};{varNames[1]};{estNames[est]}')
                if 'Eff' in est:
                    cEstimVsCut[iPt].cd(iPad+1).SetLogz()
                    hFrame.GetZaxis().SetMoreLogLabels()
                hFrame.GetXaxis().SetNdivisions(505)
                hFrame.GetYaxis().SetNdivisions(505)
                hFrame.GetXaxis().SetDecimals()
                hFrame.GetYaxis().SetDecimals()
                hEstimVsCut[iPt][est].DrawCopy('colz same')
                outDirPlotsPt[iPt].cd()
                hEstimVsCut[iPt][est].Write()
        cEstimVsCut[iPt].Update()
        cEstimVsCut[iPt].Modified()
        outDirPlotsPt[iPt].cd()
        cEstimVsCut[iPt].Write()

outFile.cd()
tSignif.Write()
outFile.Close()

if not args.batch:
    input('Press enter to exit')
