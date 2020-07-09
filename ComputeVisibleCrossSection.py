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
visCrossSec, statUnc, uncorrSystUnc, corrSystUncLow, corrSystUncHigh, \
    FDSystUncLow, FDSystUncHigh, trackSystUnc = (0 for _ in range(8))

for iPt in range(hCrossSection.GetNbinsX()):

    binWidth = hCrossSection.GetBinWidth(iPt+1)
    ptCent = hCrossSection.GetBinCenter(iPt+1)
    visCrossSec += hCrossSection.GetBinContent(iPt+1) * binWidth
    statUnc += (hCrossSection.GetBinError(iPt+1) * binWidth)**2

    # uncorrelated systematic uncertainty (yield extraction)
    uncorrSystUnc += systErr.GetRawYieldErr(ptCent)**2

    # correlated systematic uncertainty (sel eff, PID eff, gen pT shape, tracking, FD)
    trackSystUnc += systErr.GetTrackingEffErr(ptCent) * binWidth
    totCorrSystUncSqLow = (systErr.GetCutsEffErr(ptCent)**2 + systErr.GetMCPtShapeErr(ptCent)**2 + \
        systErr.GetPIDEffErr(ptCent)**2 + systErr.GetTrackingEffErr(ptCent)**2) \
            * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    totCorrSystUncSqHigh = (systErr.GetCutsEffErr(ptCent)**2 + systErr.GetMCPtShapeErr(ptCent)**2 + \
        systErr.GetPIDEffErr(ptCent)**2 + systErr.GetTrackingEffErr(ptCent)**2) \
            * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
    if isDataDriven:
        FDSystUncLow += systErr.GetDataDrivenFDErr(ptCent) * binWidth
        FDSystUncHigh += systErr.GetDataDrivenFDErr(ptCent) * binWidth
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
        FDSystUncLow += gCrossSectionFDSyst.GetErrorYlow(iPtFD) * binWidth
        FDSystUncHigh += gCrossSectionFDSyst.GetErrorYlow(iPtFD) * binWidth
        totCorrSystUncSqLow += gCrossSectionFDSyst.GetErrorYlow(iPtFD)**2 * binWidth**2
        totCorrSystUncSqLow += gCrossSectionFDSyst.GetErrorYhigh(iPtFD)**2 * binWidth**2

    corrSystUncLow += np.sqrt(totCorrSystUncSqLow)
    corrSystUncHigh += np.sqrt(totCorrSystUncSqHigh)

statUnc = np.sqrt(statUnc)
uncorrSystUnc = np.sqrt(uncorrSystUnc)
totSystUncLow = np.sqrt(uncorrSystUnc**2 + corrSystUncLow**2)
totSystUncHigh = np.sqrt(uncorrSystUnc**2 + corrSystUncHigh**2)

# fill histos and graphs
ptMin = hCrossSection.GetBinLowEdge(1)
ptMax = hCrossSection.GetXaxis().GetBinUpEdge(hCrossSection.GetNbinsX())
hVisibleCrossSectionStat = TH1F('hVisibleCrossSectionStat',
                                f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})',
                                1, 0.5, 1.5)

gVisCrossSecUncorrSyst, gVisCrossSecCorrSyst, gVisCrossSecTotSyst, gVisCrossSecSystWoTrAndFD, \
    gVisCrossSecSystTracking, gVisCrossSecSystFD = (TGraphAsymmErrors(0) for _ in range(6))
gVisCrossSecUncorrSyst.SetNameTitle('gVisCrossSecUncorrSyst',
                                    f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')
gVisCrossSecCorrSyst.SetNameTitle('gVisCrossSecCorrSyst',
                                  f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')
gVisCrossSecTotSyst.SetNameTitle('gVisCrossSecTotSyst',
                                 f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')
gVisCrossSecSystWoTrAndFD.SetNameTitle('gVisCrossSecSystWoTrackingAndFD',
                                       f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')
gVisCrossSecSystTracking.SetNameTitle('gVisCrossSecSystTracking',
                                      f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')
gVisCrossSecSystFD.SetNameTitle('gVisCrossSecSystFD',
                                f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})')

SetObjectStyle(hVisibleCrossSectionStat, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(gVisCrossSecUncorrSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecCorrSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecTotSyst, color=kBlack, fillstyle=0)

hVisibleCrossSectionStat.SetBinContent(1, visCrossSec)
hVisibleCrossSectionStat.SetBinError(1, statUnc)
gVisCrossSecUncorrSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecCorrSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecTotSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystWoTrAndFD.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystTracking.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystFD.SetPoint(0, 1., visCrossSec)
gVisCrossSecUncorrSyst.SetPointError(0, 0.3, 0.3, uncorrSystUnc, uncorrSystUnc)
gVisCrossSecCorrSyst.SetPointError(0, 0.3, 0.3, corrSystUncLow, corrSystUncHigh)
gVisCrossSecTotSyst.SetPointError(0, 0.3, 0.3, totSystUncLow, totSystUncHigh)
gVisCrossSecSystWoTrAndFD.SetPointError(0, 0.3, 0.3, np.sqrt(totSystUncLow**2 - trackSystUnc**2 - FDSystUncLow**2),
                                        np.sqrt(totSystUncHigh**2 - trackSystUnc**2 - FDSystUncHigh**2))
gVisCrossSecSystTracking.SetPointError(0, 0.3, 0.3, trackSystUnc, trackSystUnc)
gVisCrossSecSystFD.SetPointError(0, 0.3, 0.3, FDSystUncLow, FDSystUncHigh)

# otput file
outFile = TFile.Open(args.outFileName, 'recreate')
hVisibleCrossSectionStat.Write()
gVisCrossSecUncorrSyst.Write()
gVisCrossSecCorrSyst.Write()
gVisCrossSecTotSyst.Write()
gVisCrossSecSystWoTrAndFD.Write()
gVisCrossSecSystTracking.Write()
gVisCrossSecSystFD.Write()
outFile.Close()

print(f'Saved output file: {args.outFileName}')
