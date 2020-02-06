'''
python script for the optimisation of the working point using TTrees or pandas dataframes as input
run: python ScanSelectionTree.py cfgFileName.yml outFileName.root
'''

import sys
import argparse
import itertools
import numpy as np
import yaml
from root_numpy import fill_hist
from ROOT import TFile, TH1F, TH2F, TCanvas, TLegend, TNtuple, TDirectoryFile # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, kRainBow  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
#pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.AnalysisUtils import ComputeEfficiency, GetPromptFDFractionFc, GetExpectedBkgFromSideBands, \
    GetExpectedBkgFromMC, GetExpectedSignal
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.DfUtils import LoadDfFromRootOrParquet

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('outFileName', metavar='text', default='outFile.root',
                    help='output root file name')
parser.add_argument('--batch', action='store_true', default=False,
                    help='flag to enable batch mode')
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

#reshuffle bkg and take only a fraction of it
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
    hPreselEffFD = infilePreselEff.Get(inputCfg['infiles']['preseleff']['FDhistoname'])

# load acceptance
infileAcc = TFile.Open(inputCfg['infiles']['acceptance'])
hPtGenAcc = infileAcc.Get('hPtGenAcc')
hPtGenLimAcc = infileAcc.Get('hPtGenLimAcc')

# load cross sections
inFileCrossSec = TFile.Open(inputCfg['predictions']['crosssec']['filename'])
hCrossSecPrompt = inFileCrossSec.Get(
    f"{inputCfg['predictions']['crosssec']['histonames']['prompt']}_central")
hCrossSecFD = inFileCrossSec.Get(
    f"{inputCfg['predictions']['crosssec']['histonames']['feeddown']}_central_corr")

# load RAA
# TODO: add possibility to pass prediction from file instead of fixed value
RaaPrompt = inputCfg['predictions']['Raa']['prompt']
if not isinstance(RaaPrompt, float) and not isinstance(RaaPrompt, int):
    print('ERROR: currently only fixed values of RAA are supported. Exit')
    exit()
RaaFD = inputCfg['predictions']['Raa']['feeddown']
if not isinstance(RaaFD, float) and not isinstance(RaaFD, int):
    print('ERROR: currently only fixed values of RAA are supported. Exit')
    exit()

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
    outDirPt = []

varsName4Tuple = ':'.join(cutVars) + ':PtMin:PtMax:S:B:Signif:SoverB:EffPrompt:EffFD:fprompt:fFD'
tSignif = TNtuple('tSignif', 'tSignif', varsName4Tuple)

totSets = 1
for cutRange in cutRanges:
    totSets *= len(cutRange)
print(f'Total number of sets per pT bin:{totSets}')

for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):

    if not inputCfg['infiles']['background']['isMC']:
        outDirFitSB.cd()
        outDirPt.append(TDirectoryFile(f'pT{ptMin}-{ptMax}', f'pT{ptMin}-{ptMax}'))
        outDirPt[iPt].Write()

    dfPromptPt = dfPrompt.query(f'{ptMin} < pt_cand < {ptMax}')
    dfFDPt = dfFD.query(f'{ptMin} < pt_cand < {ptMax}')
    dfBkgPt = dfBkg.query(f'{ptMin} < pt_cand < {ptMax}')

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
    acc, accUnc = ComputeEfficiency(numAcc, np.sqrt(numAcc), denAcc, np.sqrt(denAcc))

    # cross section from theory
    ptBinCrossSecMin = hCrossSecPrompt.GetXaxis().FindBin(ptMin*1.0001)
    ptBinCrossSecMax = hCrossSecPrompt.GetXaxis().FindBin(ptMax*0.9999)
    crossSecPrompt = hCrossSecPrompt.Integral(ptBinCrossSecMin, ptBinCrossSecMax) / (ptMax - ptMin)
    crossSecFD = hCrossSecFD.Integral(ptBinCrossSecMin, ptBinCrossSecMax) / (ptMax - ptMin)

    for iSet, cutSet in enumerate(itertools.product(*cutRanges)):
        selToApply = ''
        for iCut, (cut, upperLower, varName) in enumerate(zip(cutSet, upperLowerCuts, varNames)):
            if iCut == 0:
                selToApply = f'{varName}{upperLower}{cut}'
            else:
                selToApply += f' & {varName}{upperLower}{cut}'

        if (iSet+1) % 100 == 0:
            print(f'Testing cut set number {iSet+1}: {selToApply}') 

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
        fill_hist(hMassSignal, dfPromptPtSel['inv_mass'].values)#  + dfFDPtSel['inv_mass'].values

        # expected signal
        expSignal = GetExpectedSignal(crossSecPrompt, ptMax-ptMin, 1., effTimesAccPrompt,
                                      fPrompt[0], BR, 1., nExpEv, sigmaMB, Taa, RaaPrompt)

        # expected background
        if inputCfg['infiles']['background']['isMC']:
            expBkg = GetExpectedBkgFromMC(hMassBkg)
        else:
            expBkg, hMassSB = GetExpectedBkgFromSideBands(hMassBkg, 'pol2', 4, hMassSignal)
            outDirPt[iPt].cd()
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


# SetGlobalStyle(padleftmargin=0.15, padtopmargin=0.08, titleoffsetx=1., titleoffsety=1.4, opttitle=1, palette=kRainBow)
# TODO: add plots

outFile.cd()
tSignif.Write()
outFile.Close()
