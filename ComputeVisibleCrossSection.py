'''
python script to compute visible cross section
run: python ComputeVisibleCrossSection.py inFilePtDiffCrossSec.root outFile.root [--Dplus] [--Ds]
'''

import sys
import argparse
import ctypes
import numpy as np
from ROOT import TFile, TH1F, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kFullCircle # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.AnalysisUtils import ScaleGraph

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFileNamePtDiffCrossSec', metavar='text', default='inFilePtDiffCrossSec.root',
                    help='root file with pT-differential cross section')
parser.add_argument('outFileName', metavar='text', default='outFile.root', help='output root file')
parser.add_argument('--Dplus', action='store_true', default=False, help='enable calculation for D+')
parser.add_argument('--Ds', action='store_true', default=False, help='enable calculation for Ds')
args = parser.parse_args()

if not args.Dplus and not args.Ds:
    print('ERROR: you should enable the comparison for either D+ or Ds! Exit')
    sys.exit()
elif args.Dplus and args.Ds:
    print('ERROR: you cannot enable the comparison for both D+ and Ds! Exit')
    sys.exit()

if args.Dplus:
    BR = 0.0898
    mesonName = 'D^{+}'
elif args.Ds:
    BR = 0.0227
    mesonName = 'D_{s}^{+}'

SetGlobalStyle()

# load inputs
inFileCross = TFile.Open(args.inFileNamePtDiffCrossSec)
hCrossSection = inFileCross.Get('hCrossSection')
isDataDriven = True
if not hCrossSection:
    hCrossSection = inFileCross.Get('histoSigmaCorr')
    hCrossSection.Scale(1.e-6 / BR)
    hCrossSection.SetStats(0)
    gCrossSectionFDSyst = inFileCross.Get('gSigmaCorrConservative')
    gCrossSectionFDSyst.RemovePoint(0)
    ScaleGraph(gCrossSectionFDSyst, 1.e-6 / BR)
    isDataDriven = False
hCrossSection.SetDirectory(0)
systErr = inFileCross.Get('AliHFSystErr')
inFileCross.Close()

# compute visible cross section
visCrossSec, visCrossSecStatUnc, visCrossSecUncorrSystUnc, \
    visCrossSecCorrSystUncLow, visCrossSecCorrSystUncHigh = (0 for _ in range(5))

for iPt in range(hCrossSection.GetNbinsX()):

    binWidth = hCrossSection.GetBinWidth(iPt+1)
    ptCent = hCrossSection.GetBinCenter(iPt+1)
    visCrossSec += hCrossSection.GetBinContent(iPt+1) * binWidth
    visCrossSecStatUnc += (hCrossSection.GetBinError(iPt+1) * binWidth)**2

    # uncorrelated systematic uncertainty (yield extraction)
    visCrossSecUncorrSystUnc += systErr.GetRawYieldErr(ptCent)**2

    # correlated systematic uncertainty (sel eff, PID eff, gen pT shape, tracking, FD)
    totCorrSystUncSqLow = (systErr.GetCutsEffErr(ptCent)**2 + systErr.GetMCPtShapeErr(ptCent)**2 + \
        systErr.GetPIDEffErr(ptCent)**2 + systErr.GetTrackingEffErr(ptCent)**2) \
            * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    totCorrSystUncSqHigh = (systErr.GetCutsEffErr(ptCent)**2 + systErr.GetMCPtShapeErr(ptCent)**2 + \
        systErr.GetPIDEffErr(ptCent)**2 + systErr.GetTrackingEffErr(ptCent)**2) \
            * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    if isDataDriven:
        totCorrSystUncSqLow += systErr.GetDataDrivenFDErr(ptCent)**2 \
            * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
        totCorrSystUncSqHigh += systErr.GetDataDrivenFDErr(ptCent)**2 \
            * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    else:
        for iPtFD in range(gCrossSectionFDSyst.GetN()):
            ptCentFD, sigma = ctypes.c_double(), ctypes.c_double()
            gCrossSectionFDSyst.GetPoint(iPtFD, ptCentFD, sigma)
            if abs(ptCentFD.value-ptCent) < 0.01:
                break
        totCorrSystUncSqLow += gCrossSectionFDSyst.GetErrorYlow(iPtFD)**2 * binWidth**2
        totCorrSystUncSqLow += gCrossSectionFDSyst.GetErrorYhigh(iPtFD)**2 * binWidth**2

    visCrossSecCorrSystUncLow += np.sqrt(totCorrSystUncSqLow)
    visCrossSecCorrSystUncHigh += np.sqrt(totCorrSystUncSqHigh)

visCrossSecStatUnc = np.sqrt(visCrossSecStatUnc)
visCrossSecUncorrSystUnc = np.sqrt(visCrossSecUncorrSystUnc)
visCrossSecTotSystUncLow = np.sqrt(visCrossSecUncorrSystUnc**2 + visCrossSecCorrSystUncLow**2)
visCrossSecTotSystUncHigh = np.sqrt(visCrossSecUncorrSystUnc**2 + visCrossSecCorrSystUncHigh**2)

# fill histos and graphs
ptMin = hCrossSection.GetBinLowEdge(1)
ptMax = hCrossSection.GetXaxis().GetBinUpEdge(hCrossSection.GetNbinsX())
hVisibleCrossSectionStat = TH1F('hVisibleCrossSectionStat',
                                f';;#sigma (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})', 1, 0.5, 1.5)
gVisibleCrossSectionUncorrSyst = TGraphAsymmErrors(0)
gVisibleCrossSectionCorrSyst = TGraphAsymmErrors(0)
gVisibleCrossSectionTotSyst = TGraphAsymmErrors(0)
gVisibleCrossSectionUncorrSyst.SetNameTitle('gVisibleCrossSectionUncorrSyst',
                                            f';;#sigma (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')
gVisibleCrossSectionCorrSyst.SetNameTitle('gVisibleCrossSectionCorrSyst',
                                          f';;#sigma (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')
gVisibleCrossSectionTotSyst.SetNameTitle('gVisibleCrossSectionTotSyst',
                                         f';;#sigma (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')

SetObjectStyle(hVisibleCrossSectionStat, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(gVisibleCrossSectionUncorrSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisibleCrossSectionCorrSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisibleCrossSectionTotSyst, color=kBlack, fillstyle=0)

hVisibleCrossSectionStat.SetBinContent(1, visCrossSec)
hVisibleCrossSectionStat.SetBinError(1, visCrossSecStatUnc)
gVisibleCrossSectionUncorrSyst.SetPoint(0, 1., visCrossSec)
gVisibleCrossSectionCorrSyst.SetPoint(0, 1., visCrossSec)
gVisibleCrossSectionTotSyst.SetPoint(0, 1., visCrossSec)
gVisibleCrossSectionUncorrSyst.SetPointError(0, 0.3, 0.3, visCrossSecUncorrSystUnc, visCrossSecUncorrSystUnc)
gVisibleCrossSectionCorrSyst.SetPointError(0, 0.3, 0.3, visCrossSecCorrSystUncLow, visCrossSecCorrSystUncHigh)
gVisibleCrossSectionTotSyst.SetPointError(0, 0.3, 0.3, visCrossSecTotSystUncLow, visCrossSecTotSystUncHigh)

# otput file
outFile = TFile.Open(args.outFileName, 'recreate')
hVisibleCrossSectionStat.Write()
gVisibleCrossSectionUncorrSyst.Write()
gVisibleCrossSectionCorrSyst.Write()
gVisibleCrossSectionTotSyst.Write()
outFile.Close()

print(f'Saved output file: {args.outFileName}')
