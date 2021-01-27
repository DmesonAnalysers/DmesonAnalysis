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
parser.add_argument('--Dzero', action='store_true', default=False, help='enable calculation for D0')
parser.add_argument('--Dstar', action='store_true', default=False, help='enable calculation for D*')
parser.add_argument('--PbPb', action='store_true', default=True, help='enable calculation for PbPb collisions')
parser.add_argument('--ptmin', default=0., type=float, help='minimum pT')
parser.add_argument('--ptmax', default=1.e10, type=float, help='maximum pT')
args = parser.parse_args()

argChannels = [args.Dplus, args.Ds, args.Dzero, args.Dstar]
if not any(argChannels):
    print('ERROR: you should enable the comparison for either Dplus, Ds, Dzero, or Dstar! Exit')
    sys.exit()
elif len(np.where(argChannels)[0]) > 1:
    print('ERROR: you can enable only one meson at a time! Exit')
    sys.exit()

if not args.PbPb:
    lumiUnc = 0.021 # pp@5.02 TeV
else:
    lumiUnc = 0.

if args.Dplus:
    BR = 0.0938
    BRunc = 0.0016
    mesonName = 'D^{+}'
elif args.Ds:
    BR = 0.0224
    BRunc = 0.0008
    mesonName = 'D_{s}^{+}'
elif args.Dzero:
    BR = 0.03950
    BRunc = 0.00031
    mesonName = 'D^{0}'
elif args.Dstar:
    BR = 0.03950*0.677
    BRunc = np.sqrt(0.00031**2/0.03950**2 + 0.005*2/0.677**2) * BR
    mesonName = 'D*^{+}'

SetGlobalStyle()

# load inputs
HFsystErr = []
inFileCross = TFile.Open(args.inFileNamePtDiffCrossSec)
hCrossSection = inFileCross.Get('hCrossSection')
isDataDriven = True
if not hCrossSection:
    if not args.PbPb:
        hCrossSection = inFileCross.Get('histoSigmaCorr')
        gCrossSectionFDSyst = inFileCross.Get('gSigmaCorrConservative')
        gCrossSectionFDSyst.RemovePoint(0)
        hCrossSection.Scale(1.e-6 / BR)
        ScaleGraph(gCrossSectionFDSyst, 1.e-6 / BR)
    else:
        hCrossSection = inFileCross.Get('hAAC')
        gCrossSectionFDSyst = inFileCross.Get('gaaCsystFD')
    hCrossSection.SetStats(0)
    isDataDriven = False
hCrossSection.SetDirectory(0)
if not args.PbPb:
    if not args.Dzero:
        HFsystErr.append(inFileCross.Get('AliHFSystErr'))
    else:
        HFsystErr.append(inFileCross.Get('AliHFSystErrTopol'))
        HFsystErr.append(inFileCross.Get('AliHFSystErrLowpt'))
else:
    if not args.Dzero:
        HFsystErr.append(inFileCross.Get('AliHFSystErrAA'))
    else:
        if inFileCross.Get('AliHFSystErrAA1'):
            HFsystErr.append(inFileCross.Get('AliHFSystErrAA1'))
            HFsystErr.append(inFileCross.Get('AliHFSystErrAA2'))
        else:
            HFsystErr.append(inFileCross.Get('AliHFSystErrAA'))

inFileCross.Close()

# compute visible cross section
visCrossSec, statUnc, uncorrSystUnc, corrSystUncLow, corrSystUncHigh, \
    FDSystUncLow, FDSystUncHigh, trackSystUnc, systWoTrFDBRLumi = (0 for _ in range(9))

for iPt in range(hCrossSection.GetNbinsX()):

    binWidth = hCrossSection.GetBinWidth(iPt+1)
    ptCent = hCrossSection.GetBinCenter(iPt+1)

    # integrate only between ptmin and ptmax
    if ptCent-binWidth/2 < args.ptmin:
        continue
    if ptCent+binWidth/2 > args.ptmax:
        continue

    visCrossSec += hCrossSection.GetBinContent(iPt+1) * binWidth
    statUnc += (hCrossSection.GetBinError(iPt+1) * binWidth)**2

    # get correct systematic uncertainties
    if ptCent < 1 and args.Dzero:
        systErr = HFsystErr[1] # Dzero low pT
    else:
        systErr = HFsystErr[0] # topological

    # uncorrelated systematic uncertainty (yield extraction)
    uncorrSystUnc += systErr.GetRawYieldErr(ptCent)**2 * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2

    # correlated systematic uncertainty (sel eff, PID eff, gen pT shape, tracking, FD, BR)
    trackSystUnc += systErr.GetTrackingEffErr(ptCent) * hCrossSection.GetBinContent(iPt+1) * binWidth
    systWoTrFDBRLumi += (systErr.GetCutsEffErr(ptCent)**2 + systErr.GetMCPtShapeErr(ptCent)**2 + \
        systErr.GetPIDEffErr(ptCent)**2) * hCrossSection.GetBinContent(iPt+1)**2 * binWidth**2
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
        totCorrSystUncSqHigh += gCrossSectionFDSyst.GetErrorYhigh(iPtFD)**2 * binWidth**2

    corrSystUncLow += np.sqrt(totCorrSystUncSqLow)
    corrSystUncHigh += np.sqrt(totCorrSystUncSqHigh)

statUnc = np.sqrt(statUnc)
uncorrSystUnc = np.sqrt(uncorrSystUnc)
corrSystUncLow = np.sqrt(corrSystUncLow**2 + (lumiUnc**2 + BRunc**2/BR**2)*visCrossSec)
corrSystUncHigh = np.sqrt(corrSystUncHigh**2 + (lumiUnc**2 + BRunc**2/BR**2)*visCrossSec)
totSystUncLow = np.sqrt(uncorrSystUnc**2 + corrSystUncLow**2)
totSystUncHigh = np.sqrt(uncorrSystUnc**2 + corrSystUncHigh**2)

# fill histos and graphs
if args.ptmin <= hCrossSection.GetBinLowEdge(1):
    ptMin = hCrossSection.GetBinLowEdge(1)
else:
    ptMin = args.ptmin
if args.ptmax >= hCrossSection.GetXaxis().GetBinUpEdge(hCrossSection.GetNbinsX()):
    ptMax = hCrossSection.GetXaxis().GetBinUpEdge(hCrossSection.GetNbinsX())
else:
    ptMax = args.ptmax

graphName = 'VisCrossSec'
utilStr = f';;d#sigma/d#it{{y}} (#mub) ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})'
if args.PbPb:
    graphName = 'PtVisYield'
    utilStr = f';;d#it{{N}}/d#it{{y}} ({ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}})'

hVisibleCrossSectionStat = TH1F(f'h{graphName}', utilStr, 1, 0.5, 1.5)

gVisCrossSecUncorrSyst, gVisCrossSecCorrSyst, gVisCrossSecTotSyst, gVisCrossSecSystWoTrFDBRAndLumi, \
    gVisCrossSecSystTracking, gVisCrossSecSystFD, gVisCrossSecSystLumi, gVisCrossSecSystBR = \
        (TGraphAsymmErrors(0) for _ in range(8))
gVisCrossSecUncorrSyst.SetNameTitle(f'g{graphName}UncorrSys', utilStr)
gVisCrossSecCorrSyst.SetNameTitle(f'g{graphName}PtCorrSys', utilStr)
gVisCrossSecTotSyst.SetNameTitle(f'g{graphName}TotSys', utilStr)
gVisCrossSecSystWoTrFDBRAndLumi.SetNameTitle(f'g{graphName}SysWoTrFDBRAndLumi', utilStr)
gVisCrossSecSystTracking.SetNameTitle('gVisCrossSecSysTracking', utilStr)
gVisCrossSecSystFD.SetNameTitle(f'g{graphName}SysFD', utilStr)
gVisCrossSecSystLumi.SetNameTitle(f'g{graphName}SysLumi', utilStr)
gVisCrossSecSystBR.SetNameTitle(f'g{graphName}SysBR', utilStr)

SetObjectStyle(hVisibleCrossSectionStat, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(gVisCrossSecUncorrSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecCorrSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecTotSyst, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystLumi, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystBR, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystFD, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystTracking, color=kBlack, fillstyle=0)
SetObjectStyle(gVisCrossSecSystWoTrFDBRAndLumi, color=kBlack, fillstyle=0)

hVisibleCrossSectionStat.SetBinContent(1, visCrossSec)
hVisibleCrossSectionStat.SetBinError(1, statUnc)
gVisCrossSecUncorrSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecCorrSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecTotSyst.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystWoTrFDBRAndLumi.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystTracking.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystFD.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystLumi.SetPoint(0, 1., visCrossSec)
gVisCrossSecSystBR.SetPoint(0, 1., visCrossSec)
gVisCrossSecUncorrSyst.SetPointError(0, 0.3, 0.3, uncorrSystUnc, uncorrSystUnc)
gVisCrossSecCorrSyst.SetPointError(0, 0.3, 0.3, corrSystUncLow, corrSystUncHigh)
gVisCrossSecTotSyst.SetPointError(0, 0.3, 0.3, totSystUncLow, totSystUncHigh)
gVisCrossSecSystWoTrFDBRAndLumi.SetPointError(0, 0.3, 0.3, systWoTrFDBRLumi, systWoTrFDBRLumi)
gVisCrossSecSystTracking.SetPointError(0, 0.3, 0.3, trackSystUnc, trackSystUnc)
gVisCrossSecSystFD.SetPointError(0, 0.3, 0.3, FDSystUncLow, FDSystUncHigh)
gVisCrossSecSystLumi.SetPointError(0, 0.3, 0.3, visCrossSec*lumiUnc, visCrossSec*lumiUnc)
gVisCrossSecSystBR.SetPointError(0, 0.3, 0.3, visCrossSec*BRunc/BR, visCrossSec*BRunc/BR)

# otput file
outFile = TFile.Open(args.outFileName, 'recreate')
hVisibleCrossSectionStat.Write()
gVisCrossSecUncorrSyst.Write()
gVisCrossSecCorrSyst.Write()
gVisCrossSecTotSyst.Write()
gVisCrossSecSystWoTrFDBRAndLumi.Write()
gVisCrossSecSystTracking.Write()
gVisCrossSecSystFD.Write()
gVisCrossSecSystBR.Write()
if not args.PbPb:
    gVisCrossSecSystLumi.Write()
outFile.Close()

print(f'Saved output file: {args.outFileName}')
