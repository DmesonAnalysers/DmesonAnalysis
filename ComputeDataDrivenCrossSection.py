'''
python script for the computation of the production cross section of prompt or feed-down D
run: python ComputeDataDrivenCrossSection.py rawYieldFile.root effAccFile.root fracFile.root outFile.root
                                             [--prompt] [--FD] [--Dplus] [--Ds] [--system] [--energy] [--batch]
prompt or FD and Dplus or Ds must be specified
'''

import sys
import argparse
import numpy as np
from ROOT import TFile, TCanvas, TLegend, TGraphErrors, gROOT  # pylint: disable=import-error,no-name-in-module
from ROOT import AliHFSystErr  # pylint: disable=import-error,no-name-in-module
from utils.AnalysisUtils import ComputeCrossSection, GetPromptFDFractionCutSet
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('rawYieldFileName', metavar='text', default='rawYieldFile.root', help='root file with raw yields')
parser.add_argument('effAccFileName', metavar='text', default='effAccFile.root',
                    help='root file with efficiency and acceptance')
parser.add_argument('fracFileName', metavar='text', default='fracFile.root',
                    help='root file with prompt (FD) fraction')
parser.add_argument('outFileName', metavar='text', default='outFile.root', help='root output file name')
parser.add_argument('--system', metavar='text', default='pp', help='collision system (pp, pPb, PbPb)')
parser.add_argument('--energy', metavar='text', default='5.02', help='energy (5.02)')
parser.add_argument('--centrality', metavar='text', default='010', help='centrality (010, 3050)')
parser.add_argument("--prompt", action='store_true', help='flag to compute prompt cross section', default=False)
parser.add_argument("--FD", action='store_true', help='flag to compute FD cross section', default=False)
parser.add_argument("--Dplus", action='store_true', help='flag to compute D+ cross section', default=False)
parser.add_argument("--Ds", action='store_true', help='flag to compute Ds cross section', default=False)
parser.add_argument("--batch", action='store_true', help='suppress video output', default=False)
args = parser.parse_args()

# Define arguments (systematic uncertainties from AliHFSystErr)
systErr = AliHFSystErr()
systErr.SetIsDataDrivenFDAnalysis(True)

if args.system == 'pp':
    axisTitle = ';#it{p}_{T} (GeV/#it{c}); d#sigma/d#it{p}_{T} #times BR (#mub GeV^{-1} #it{c})'
    histoName = 'CrossSection'
    systErr.SetCollisionType(0)
    systErr.SetRunNumber(17)
    if args.energy == '5.02':
        if args.Dplus:
            systErr.SetIs5TeVAnalysis(True)
        sigmaMB = 50.87e+3 # ub
        lumiUnc = 0.021
    else:
        print(f'Energy {args.energy} not implemented! Exit')
        sys.exit()
elif args.system == 'PbPb':
    axisTitle = ';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T} #times BR (GeV^{-1} #it{c})'
    histoName = 'CorrYield'
    if args.centrality in ['010', '3050']:
        systErr.SetCentrality(args.centrality)
    else:
        print('ERROR: only 0-10 and 30-50 centrality classes implemented! Exit')
        sys.exit()
    systErr.SetCollisionType(1)
    systErr.SetRunNumber(18)
    sigmaMB = 1. # yields in case of PbPb

if args.Dplus and args.Ds:
    print('ERROR: not possible to select more than one meson at the same time! Exit')
    sys.exit()
elif args.Dplus:
    systErr.Init(2)
elif args.Ds:
    systErr.Init(4)
else:
    print('ERROR: Dplus or Ds must be specified! Exit')
    sys.exit()

if args.prompt and args.FD:
    print('ERROR: not possible to select prompt and FD at the same time! Exit')
    sys.exit()

# load input file
rawYieldFile = TFile.Open(args.rawYieldFileName)
hRawYields = rawYieldFile.Get('hRawYields')
hEvForNorm = rawYieldFile.Get('hEvForNorm')
nEv = hEvForNorm.GetBinContent(1)

effAccFile = TFile.Open(args.effAccFileName)
hEffAccPrompt = effAccFile.Get('hAccEffPrompt')
hEffAccFD = effAccFile.Get('hAccEffFD')

fracFile = TFile.Open(args.fracFileName)
hCorrYieldPrompt = fracFile.Get('hCorrYieldPrompt')
hCorrYieldFD = fracFile.Get('hCorrYieldFD')
hCovPromptPrompt = fracFile.Get('hCovPromptPrompt')
hCovPromptFD = fracFile.Get('hCovPromptFD')
hCovFDFD = fracFile.Get('hCovFDFD')

# TODO: improve protection checking the limits of the bins besides the number
if hRawYields.GetNbinsX() != hEffAccPrompt.GetNbinsX() or hRawYields.GetNbinsX() != hCorrYieldPrompt.GetNbinsX():
    print('ERROR: inconsistent number of bins in input objects! Exit')
    sys.exit()

hCrossSection = hRawYields.Clone(f'h{histoName}')
hCrossSection.SetTitle(axisTitle)

# systematic uncertainties --> total taken as sum in quadrature (uncorrelated)
systGetter = {'YieldExtr': 'GetRawYieldErr', 'SelEff': 'GetCutsEffErr', 'TrEff': 'GetTrackingEffErr',
              'PIDEff': 'GetPIDEffErr', 'PtShape': 'GetMCPtShapeErr', 'FD': 'GetDataDrivenFDErr',
              'Tot': 'GetTotalSystErr'}
systColors = {'YieldExtr': 'kAzure+4', 'SelEff': 'kRed', 'TrEff': 'kBlue+1',
              'PIDEff': 'kOrange+7', 'PtShape': 'kRed', 'FD': 'kAzure+4',
              'Tot': 'kBlack'}
gCrossSectionSyst = {}
for systSource in systGetter:
    gCrossSectionSyst[systSource] = TGraphErrors(0)
    gCrossSectionSyst[systSource].SetName(f'g{histoName}Syst{systSource}')
    gCrossSectionSyst[systSource].SetTitle(axisTitle)
    SetObjectStyle(gCrossSectionSyst[systSource], color=GetROOTColor(systColors[systSource]), fillstyle=0)

hPromptFrac = hEffAccPrompt.Clone('hPromptFrac')
hFDFrac = hEffAccFD.Clone('hFDFrac')
hPromptFrac.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{f}_{prompt}')
hFDFrac.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{f}_{FD}')

for iPt in range(hCrossSection.GetNbinsX()):
    ptMin = hRawYields.GetBinLowEdge(iPt+1)
    ptMax = ptMin+hRawYields.GetBinWidth(iPt+1)
    ptCent = hRawYields.GetBinCenter(iPt+1)
    rawYield = hRawYields.GetBinContent(iPt+1)
    rawYieldUnc = hRawYields.GetBinError(iPt+1)
    effAccPrompt = hEffAccPrompt.GetBinContent(iPt+1)
    effAccFD = hEffAccFD.GetBinContent(iPt+1)
    effAccPromptUnc = hEffAccPrompt.GetBinError(iPt+1)
    effAccFDUnc = hEffAccFD.GetBinError(iPt+1)

    # ingredients for (prompt or FD) fraction computation
    corrYieldPrompt = hCorrYieldPrompt.GetBinContent(iPt+1)
    corrYieldFD = hCorrYieldFD.GetBinContent(iPt+1)
    covPromptPrompt = hCovPromptPrompt.GetBinContent(iPt+1)
    covPromptFD = hCovPromptFD.GetBinContent(iPt+1)
    covFDFD = hCovFDFD.GetBinContent(iPt+1)

    # prompt and FD and fractions
    fracPromptFD, uncFracPromptFD = GetPromptFDFractionCutSet(effAccPrompt, effAccFD, corrYieldPrompt, corrYieldFD,
                                                              covPromptPrompt, covFDFD, covPromptFD)

    hPromptFrac.SetBinContent(iPt+1, fracPromptFD[0])
    hPromptFrac.SetBinError(iPt+1, uncFracPromptFD[0])
    hFDFrac.SetBinContent(iPt+1, fracPromptFD[1])
    hFDFrac.SetBinError(iPt+1, uncFracPromptFD[1])

    if args.prompt:
        effAcc = effAccPrompt
        uncEffAcc = effAccPromptUnc
        frac = fracPromptFD[0]
        uncFrac = uncFracPromptFD[0]
    else:
        effAcc = effAccFD
        uncEffAcc = effAccFDUnc
        frac = fracPromptFD[1]
        uncFrac = uncFracPromptFD[1]

    if args.FD:
        crossSec, crossSecUnc = ComputeCrossSection(rawYield, rawYieldUnc, frac, uncFrac, effAcc,
                                                    ptMax - ptMin, 1., sigmaMB, nEv, 1., 'corr')
    else:
        # TODO: check if uncorrelated is the right option or anti-correlated is better
        crossSec, crossSecUnc = ComputeCrossSection(rawYield, rawYieldUnc, frac, uncFrac, effAcc,
                                                    ptMax - ptMin, 1., sigmaMB, nEv, 1., 'uncorr')

    hCrossSection.SetBinContent(iPt+1, crossSec)
    hCrossSection.SetBinError(iPt+1, crossSecUnc)

    # systematic uncertainties (statistical uncertainty on eff included)
    for systSource in systGetter:
        gCrossSectionSyst[systSource].SetPoint(iPt, ptCent, crossSec)
        if 'SelEff' not in systSource:
            relSyst = getattr(systErr, systGetter[systSource])(ptCent)
        else:
            relSyst = np.sqrt(getattr(systErr, systGetter[systSource])(ptCent)**2 + (uncEffAcc/effAcc)**2)
        gCrossSectionSyst[systSource].SetPointError(iPt, 0.4, relSyst * crossSec)

if args.system == 'pp':
    gCrossSectionSystLumi = TGraphErrors(0)
    gCrossSectionSystLumi.SetName('gCrossSectionSystLumi')
    gCrossSectionSystLumi.SetTitle('Luminosity syst. unc.;;')
    gCrossSectionSystLumi.SetPoint(0, 1., 1.)
    gCrossSectionSystLumi.SetPointError(0, 0.4, lumiUnc)
    SetObjectStyle(gCrossSectionSystLumi, color=GetROOTColor('kBlue'), fillstyle=0)

gROOT.SetBatch(args.batch)
SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14)

cCrossSec = TCanvas(f'c{histoName}', '', 700, 800)
cCrossSec.SetLogy()
hCrossSection.Draw()
gCrossSectionSyst['Tot'].Draw('2')

legFrac = TLegend(0.2, 0.84, 0.4, 0.94)
legFrac.SetBorderSize(0)
legFrac.SetFillStyle(0)
legFrac.SetTextSize(0.045)
legFrac.AddEntry(hPromptFrac, 'Prompt', 'p')
legFrac.AddEntry(hFDFrac, 'Non-prompt', 'p')

legEff = legFrac.Clone('legEff')
legEff.SetY1(0.2)
legEff.SetY2(0.4)
cCrossSec.Update()

cFrac = TCanvas('cFrac', '', 800, 800)
cFrac.DrawFrame(hPromptFrac.GetBinLowEdge(1), 0., ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); fraction')
hPromptFrac.Draw('same')
hFDFrac.Draw('same')
legFrac.Draw()
cFrac.Update()

cEff = TCanvas('cEff', '', 800, 800)
cEff.DrawFrame(hPromptFrac.GetBinLowEdge(1), 1.e-4, ptMax, 1., ';#it{p}_{T} (GeV/#it{c}); (Acc#times#font[152]{e})')
cEff.SetLogy()
hEffAccPrompt.Draw('same')
hEffAccFD.Draw('same')
legEff.Draw()
cEff.Update()

outFile = TFile(args.outFileName, 'recreate')
hCrossSection.Write()
for systSource in systGetter:
    gCrossSectionSyst[systSource].Write()
if args.system == 'pp':
    gCrossSectionSystLumi.Write()
hRawYields.Write()
hEffAccPrompt.Write()
hEffAccFD.Write()
hPromptFrac.Write()
hFDFrac.Write()
hEvForNorm.Write()
cCrossSec.Write()
cFrac.Write()
cEff.Write()
systErr.Write()
outFile.Close()

if not args.batch:
    input('Press enter to exit')
