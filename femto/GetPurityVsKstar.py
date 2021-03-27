'''
python script for the evaluation of the D-meson S/B as a function of k*
'''

import sys
import os
import argparse
import ctypes
import numpy as np
import yaml
from ROOT import TFile, TCanvas, TGraphAsymmErrors, TLatex, TF1, TH1F, TGaxis, TLegend, TDirectoryFile, gPad # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kAzure, kOrange, kBlack, kGreen, kOpenCircle # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.AnalysisUtils import ComputeWeightedAverage
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle # pylint: disable=import-error,wrong-import-position

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='yaml config file name')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlFile:
    cfg = yaml.load(ymlFile, yaml.FullLoader)

rawYiedlFileName = cfg['rawyields']['filename']
SoverBHistoName = cfg['rawyields']['histoname']
inFileName = cfg['input']['filename']
suffix = cfg['input']['suffix']
prefix = cfg['input']['prefix']
HFsuffix = cfg['input']['HFsuffix']
evaluateFracFromBeauty = cfg['fractions']['beauty']['enable']
beautyFracFileName = cfg['fractions']['beauty']['filename']
beautyFracHistoName = cfg['fractions']['beauty']['histoname']
evaluateFracFromDstar = cfg['fractions']['Dstar']['enable']
DstarFracFileName = cfg['fractions']['Dstar']['filename']
DstarFracHistoName = cfg['fractions']['Dstar']['histoname']
DstarFracGraphName = cfg['fractions']['Dstar']['graphname']

outDirName = cfg['output']['directory']
outSuffix = cfg['output']['suffix']
addTDirectory = cfg['output']['addTDirectory']

# setting global style
SetGlobalStyle(padbottommargin=0.13, padtopmargin=0.075, titleoffsety=1.2, maxdigits=2)

# load S/B histo and compute purity
inFile = TFile.Open(rawYiedlFileName)
hSoverB = inFile.Get(SoverBHistoName)
SetObjectStyle(hSoverB)
hSignal = inFile.Get(SoverBHistoName.replace('SoverB', 'Signal'))
hBkg = inFile.Get(SoverBHistoName.replace('SoverB', 'Bkg'))
hSoverB.SetDirectory(0)
hSignal.SetDirectory(0)
hBkg.SetDirectory(0)
hSoverB.SetTitle(';#it{p}_{T}(D^{#pm}) (GeV/#it{c});S/B')
hPurity = hSoverB.Clone('hPurity')
hPurity.GetYaxis().SetTitle('S/(S+B)')
hPurity.SetDirectory(0)
for iPt in range(hSoverB.GetNbinsX()):
    bkg = hBkg.GetBinContent(iPt+1)
    bkgUnc = hBkg.GetBinError(iPt+1)
    sgn = hSignal.GetBinContent(iPt+1)
    sgnUnc = hSignal.GetBinError(iPt+1)
    sgnPlusBkg = sgn+bkg
    purity = sgn/sgnPlusBkg
    purityUnc = np.sqrt((bkg/sgnPlusBkg**2)**2 * sgnUnc**2 + (sgn/sgnPlusBkg**2)**2 * bkgUnc**2)
    hPurity.SetBinContent(iPt+1, purity)
    hPurity.SetBinError(iPt+1, purityUnc)
SetObjectStyle(hSoverB, color=kRed+1, fillstyle=0)
SetObjectStyle(hPurity, color=kRed+1, fillstyle=0)
ptMax = hSoverB.GetXaxis().GetBinUpEdge(hSoverB.GetNbinsX())
inFile.Close()
# load systematics if present
systPurity = cfg['rawyields']['systpurity']
if systPurity is not None:
    if len(systPurity) != hPurity.GetNbinsX():
        print('ERROR: binning of S/S+B systematics do not match histogram! Exit')
        sys.exit()

# if enabled, load fraction from beauty
if evaluateFracFromBeauty:
    inFile = TFile.Open(beautyFracFileName)
    hFracFromB = inFile.Get(beautyFracHistoName)
    hFracFromB.SetDirectory(0)
    SetObjectStyle(hFracFromB, color=kAzure+4, fillstyle=0)
    hFracFromC = hFracFromB.Clone('hFracFromC')
    hFracFromC.SetDirectory(0)
    for iPt in range(hFracFromC.GetNbinsX()):
        hFracFromC.SetBinContent(iPt+1, 1-hFracFromC.GetBinContent(iPt+1))
    inFile.Close()
    # load systematics if present
    systBeauty = cfg['fractions']['beauty']['systbeauty']
    if systBeauty is not None:
        if len(systBeauty) != hFracFromB.GetNbinsX():
            print('ERROR: binning of f_nonprompt systematics do not match histogram! Exit')
            sys.exit()

# if enabled, load fraction from Dstar
if evaluateFracFromDstar:
    inFile = TFile.Open(DstarFracFileName)
    inFile.ls()
    hFracFromDstar = inFile.Get(DstarFracHistoName)
    gFracFromDstar = inFile.Get(DstarFracGraphName)
    gFracFromDstar.Draw('ap')
    gFracFromDstar.RemovePoint(0)
    hFracFromDstar.SetDirectory(0)
    SetObjectStyle(hFracFromDstar, color=kGreen+2, fillstyle=0)
    inFile.Close()

# load corr histos
corrNames = ['Particle0_Particle2', 'Particle1_Particle3', 'Particle0_Particle3', 'Particle1_Particle2']
corrTitles = {'Particle0_Particle2': 'p - D^{+}',
              'Particle0_Particle3': 'p - D^{#font[122]{-}}',
              'Particle1_Particle2': '#bar{p} - D^{+}',
              'Particle1_Particle3': '#bar{p} - D^{#font[122]{-}}'}

inFile = TFile.Open(inFileName)
listName = f'{prefix}_CharmFemto{HFsuffix}_ResultQA{suffix}/{prefix}_CharmFemto{HFsuffix}_ResultQA{suffix}'
inList = inFile.Get(listName)

kStarDelta = 0.2
kStarMins = [iK * kStarDelta for iK in range(15)]
kStarMaxs = [(iK+1) * kStarDelta for iK in range(15)]

hSEPairVsPtVsKstar, hSEPairVsPt, hSEPairSignalVsPt = ({} for _ in range(3))
for corrName in corrNames:
    listCorr = inList.FindObject(f'QA_{corrName}')
    hSEPairVsPtVsKstar[corrName], hSEPairVsPt[corrName] = [], []
    for kStarMin, kStarMax in zip(kStarMins, kStarMaxs):
        hSEPairVsPtVsKstar[corrName].append(listCorr.FindObject(f'KstarPtSEPartTwo_{corrName}'))
        hSEPairVsPtVsKstar[corrName][-1].SetDirectory(0)
        kStarMaxBin = hSEPairVsPtVsKstar[corrName][-1].GetXaxis().FindBin(kStarMax)
        # assuming constant pT binning
        ptRebin = round(hSoverB.GetBinWidth(1) / hSEPairVsPtVsKstar[corrName][-1].GetYaxis().GetBinWidth(1))
        hSEPairVsPt[corrName].append(
            hSEPairVsPtVsKstar[corrName][-1].ProjectionY(f'hSEPairVsPt_{corrName}', 1, kStarMaxBin))
        hSEPairVsPt[corrName][-1].SetDirectory(0)
        hSEPairVsPt[corrName][-1].Rebin(ptRebin)
        SetObjectStyle(hSEPairVsPt[corrName][-1], markerstyle=kOpenCircle)
        if 'Particle2' in corrName:
            hSEPairVsPt[corrName][-1].SetTitle(';#it{p}_{T}(D^{+}) (GeV/#it{c});SE pairs')
        else:
            hSEPairVsPt[corrName][-1].SetTitle(';#it{p}_{T}(D^{#font[122]{-}}) (GeV/#it{c});SE pairs')
        hSEPairVsPt[corrName][-1].GetXaxis().SetTitleSize(0.05)
        hSEPairVsPt[corrName][-1].GetXaxis().SetLabelSize(0.05)
        hSEPairVsPt[corrName][-1].GetYaxis().SetTitleSize(0.05)
        hSEPairVsPt[corrName][-1].GetYaxis().SetLabelSize(0.05)
        hSEPairVsPt[corrName][-1].GetYaxis().SetDecimals()
        ptMaxPair = hSEPairVsPt[corrName][-1].GetXaxis().GetBinUpEdge(hSEPairVsPt[corrName][-1].GetNbinsX())
        if ptMaxPair < ptMax:
            ptMax = ptMaxPair
inFile.Close()

for corrName in ['Particle0_Particle2_plus_Particle1_Particle3', 'Particle0_Particle3_plus_Particle1_Particle2']:
    corrNamePart = corrName.split('_plus_')
    hSEPairVsPt[corrName], hSEPairSignalVsPt[corrName] = ([] for _ in range(2))
    for iK, (kStarMin, kStarMax) in enumerate(zip(kStarMins, kStarMaxs)):
        hSEPairVsPt[corrName].append(hSEPairVsPt[corrNamePart[0]][iK].Clone(f'hSEPairVsPt_{corrName}'))
        hSEPairVsPt[corrName][iK].Add(hSEPairVsPt[corrNamePart[1]][iK])
        hSEPairVsPt[corrName][iK].SetDirectory(0)
        hSEPairVsPt[corrName][iK].SetTitle(';#it{p}_{T}(D^{#pm}) (GeV/#it{c});SE pairs')
        if iK == 0:
            corrTitles[corrName] = corrTitles[corrNamePart[0]] + ' #oplus ' + corrTitles[corrNamePart[1]]

        hSEPairSignalVsPt[corrName].append(hSEPairVsPt[corrName][iK].Clone(f'hSEPairSignalVsPt_{corrName}'))
        SetObjectStyle(hSEPairSignalVsPt[corrName][iK], color=kRed+1, fillstyle=0)

# compute average S/B and purities
gAverageSoverB, gSoverBAvPt, gAveragePurity, gAveragePuritySyst, gPurityAvPt = ({} for _ in range(5))
for corrName in ['Particle0_Particle2_plus_Particle1_Particle3', 'Particle0_Particle3_plus_Particle1_Particle2']:
    if corrName == 'Particle0_Particle2_plus_Particle1_Particle3':
        suffix = 'DpPr'
    else:
        suffix = 'DmPr'
    gAverageSoverB[corrName] = TGraphAsymmErrors(0)
    gAverageSoverB[corrName].SetNameTitle(f'gAverageSoverB_{suffix}', ';#it{k}* (GeV/#it{c}); #LT S/B #GT')
    gSoverBAvPt[corrName] = TGraphAsymmErrors(0)
    gSoverBAvPt[corrName].SetNameTitle(f'gSoverBAvPt_{suffix}',
                                       ';#it{k}* (GeV/#it{c}); S/B (#LT #it{p}_{T}(D^{#pm}) #GT)')
    gAveragePurity[corrName] = TGraphAsymmErrors(0)
    gAveragePurity[corrName].SetNameTitle(f'gAveragePurity_{suffix}', ';#it{k}* (GeV/#it{c}); #LT S/(S+B) #GT')
    gAveragePuritySyst[corrName] = TGraphAsymmErrors(0)
    gAveragePuritySyst[corrName].SetNameTitle(f'gAveragePuritySyst_{suffix}', ';#it{k}* (GeV/#it{c}); #LT S/(S+B) #GT')
    gPurityAvPt[corrName] = TGraphAsymmErrors(0)
    gPurityAvPt[corrName].SetNameTitle(f'gPurityAvPt_{suffix}',
                                       ';#it{k}* (GeV/#it{c}); S/(S+B) (#LT #it{p}_{T}(D^{#pm}) #GT)')
    SetObjectStyle(gSoverBAvPt[corrName], color=kRed+1, fillstyle=0)
    SetObjectStyle(gAverageSoverB[corrName], markerstyle=kOpenCircle)
    SetObjectStyle(gPurityAvPt[corrName], color=kRed+1, fillstyle=0)
    SetObjectStyle(gAveragePurity[corrName], markerstyle=kOpenCircle)

    for iK, (kStarMin, kStarMax) in enumerate(zip(kStarMins, kStarMaxs)):
        nPairs, nPairsUnc, SoverB, purity, SoverBUnc, purityUnc = ([] for _ in range(6))
        nBinsToAv = hSoverB.GetXaxis().FindBin(ptMax*0.9999)
        for iPt in range(nBinsToAv):
            ptCent = hSoverB.GetBinCenter(iPt+1)
            ptBinPair = hSEPairVsPt[corrName][iK].GetXaxis().FindBin(ptCent)
            nPairs.append(hSEPairVsPt[corrName][iK].GetBinContent(ptBinPair))
            nPairsUnc.append(hSEPairVsPt[corrName][iK].GetBinError(ptBinPair))
            nPairsRelUnc = nPairsUnc[iPt] / nPairs[iPt] if nPairs[iPt] > 0 else 0.
            SoverB.append(hSoverB.GetBinContent(iPt+1))
            SoverBUnc.append(hSoverB.GetBinError(iPt+1))
            purity.append(hPurity.GetBinContent(iPt+1))
            purityUnc.append(hPurity.GetBinError(iPt+1))
            purityRelUnc = purityUnc[iPt] / purity[iPt] if purity[iPt] > 0 else 0.

            nPairsSignal = nPairs[iPt] * purity[iPt]
            hSEPairSignalVsPt[corrName][iK].SetBinContent(ptBinPair, nPairsSignal)
            hSEPairSignalVsPt[corrName][iK].SetBinError(ptBinPair,
                                                        np.sqrt(purityRelUnc**2 + nPairsRelUnc**2) * nPairsSignal)

        weighAvSoverB, weighAvSoverBUnc = ComputeWeightedAverage(SoverB, nPairs, SoverBUnc, nPairsUnc)
        weighAvPurity, weighAvPurityUnc = ComputeWeightedAverage(purity, nPairs, purityUnc, nPairsUnc)

        # propagate syst unc as pT correlated
        if systPurity is not None:
            weighAvPuritySysLow, _ = ComputeWeightedAverage(np.multiply(purity, (1-np.array(systPurity[:nBinsToAv]))),
                                                            nPairs, purityUnc, nPairsUnc)
            weighAvPuritySysHigh, _ = ComputeWeightedAverage(np.multiply(purity, (1+np.array(systPurity[:nBinsToAv]))),
                                                             nPairs, purityUnc, nPairsUnc)
            systLow = weighAvPurity - weighAvPuritySysLow
            systHigh = weighAvPuritySysHigh - weighAvPurity
            if abs(systHigh-systLow) < 1.e-3 * weighAvPurity:
                systHigh = systLow
            else:
                print('WARNING: asymmetric S/S+B systematics, might be unexpected')

        avPt = hSEPairVsPt[corrName][iK].GetMean()
        avPtUnc = hSEPairVsPt[corrName][iK].GetMeanError()
        # fit data points around mean value
        fSoverB = TF1(f'fSoverB_kstar{iK}', 'pol2', 0., 10.)
        hSoverB.Fit(f'fSoverB_kstar{iK}', 'Q0')
        fPurity = TF1(f'fPurity_kstar{iK}', 'pol4', 0., 10.)
        hPurity.Fit(f'fPurity_kstar{iK}', 'Q0')

        kStarCent = (kStarMax + kStarMin) / 2
        gAverageSoverB[corrName].SetPoint(iK, kStarCent, weighAvSoverB)
        gSoverBAvPt[corrName].SetPoint(iK, kStarCent, fSoverB.Eval(avPt))
        gAverageSoverB[corrName].SetPointError(iK, kStarDelta/2, kStarDelta/2, weighAvSoverBUnc, weighAvSoverBUnc)
        gSoverBAvPt[corrName].SetPointError(iK, kStarDelta/2, kStarDelta/2,
                                            fSoverB.Eval(avPt)-fSoverB.Eval(avPt-avPtUnc),
                                            fSoverB.Eval(avPt)-fSoverB.Eval(avPt+avPtUnc))
        gAveragePurity[corrName].SetPoint(iK, kStarCent, weighAvPurity)
        gPurityAvPt[corrName].SetPoint(iK, kStarCent, fPurity.Eval(avPt))
        gAveragePurity[corrName].SetPointError(iK, kStarDelta/2, kStarDelta/2, weighAvPurityUnc, weighAvPurityUnc)
        gPurityAvPt[corrName].SetPointError(iK, kStarDelta/2, kStarDelta/2,
                                            fPurity.Eval(avPt)-fPurity.Eval(avPt-avPtUnc),
                                            fPurity.Eval(avPt)-fPurity.Eval(avPt+avPtUnc))
        if systPurity is not None:
            gAveragePuritySyst[corrName].SetPoint(iK, kStarCent, weighAvPurity)
            gAveragePuritySyst[corrName].SetPointError(iK, kStarDelta/3, kStarDelta/3, systLow, systHigh)

# compute lambda parameters
hFractions, gFractionsSyst, gBeautyFrac, gDstarFrac = ({} for _ in range(4))
hSEPairVsPtReb, hSEPairSignalVsPtReb, hSEPairBkgVsPt, hSEPairBeautyVsPt, \
    hSEPairDirectVsPt, hSEPairFromDstarVsPt = ({} for _ in range(6))

for corrName in ['Particle0_Particle2_plus_Particle1_Particle3', 'Particle0_Particle3_plus_Particle1_Particle2']:
    if corrName == 'Particle0_Particle2_plus_Particle1_Particle3':
        suffix = 'DpPr'
    else:
        suffix = 'DmPr'

    hFractions[corrName] = TH1F(f'hFractions_{suffix}', ';;fraction', 3, 0.5, 3.5)
    hFractions[corrName].SetDirectory(0)
    hFractions[corrName].GetXaxis().SetBinLabel(1, 'signal / all')
    hFractions[corrName].GetXaxis().SetBinLabel(2, 'b #rightarrow D^{#pm} / signal')
    hFractions[corrName].GetXaxis().SetBinLabel(3, 'c #rightarrow D*^{#pm} #rightarrow D^{#pm} / signal')
    SetObjectStyle(hFractions[corrName], fillstyle=0)

    gFractionsSyst[corrName] = TGraphAsymmErrors(0)
    gFractionsSyst[corrName].SetNameTitle(f'gFractionsSyst_{suffix}', ';;fraction')
    SetObjectStyle(gFractionsSyst[corrName], fillstyle=0)

    gBeautyFrac[corrName], gDstarFrac[corrName] = (TGraphAsymmErrors(0) for _ in range(2))
    SetObjectStyle(gBeautyFrac[corrName], color=kAzure+4, fillstyle=0)
    SetObjectStyle(gDstarFrac[corrName], color=kGreen+2, fillstyle=0)
    gBeautyFrac[corrName].SetName(f'gBeautyFrac_{suffix}')
    gDstarFrac[corrName].SetName(f'gDstarFrac_{suffix}')

    hSEPairVsPtReb[corrName], hSEPairSignalVsPtReb[corrName], hSEPairBkgVsPt[corrName], hSEPairBeautyVsPt[corrName], \
        hSEPairDirectVsPt[corrName], hSEPairFromDstarVsPt[corrName] = ([] for _ in range(6))

    for iK, (kStarMin, kStarMax) in enumerate(zip(kStarMins, kStarMaxs)):
        # rebin SE pair distributions to match fractions
        if evaluateFracFromBeauty:
            ptLimsReb = []
            ptMaxPair = hSEPairVsPt[corrName][-1].GetXaxis().GetBinUpEdge(hSEPairVsPt[corrName][-1].GetNbinsX())
            for iPt in range(hFracFromB.GetNbinsX()):
                ptLowEdge = hFracFromB.GetBinLowEdge(iPt+1)
                if ptLowEdge < ptMaxPair:
                    ptLimsReb.append(ptLowEdge)
            hSEPairVsPtReb[corrName].append(TH1F(f'{hSEPairVsPt[corrName][iK].GetName()}_kstar{iK}_reb',
                                                 '', len(ptLimsReb)-1, np.array(ptLimsReb, 'd')))
            hSEPairSignalVsPtReb[corrName].append(TH1F(f'{hSEPairSignalVsPt[corrName][iK].GetName()}_kstar{iK}_reb',
                                                       '', len(ptLimsReb)-1, np.array(ptLimsReb, 'd')))
            for iPt in range(hSEPairVsPtReb[corrName][iK].GetNbinsX()):
                ptLowEdge = hSEPairVsPtReb[corrName][iK].GetBinLowEdge(iPt+1)
                ptUpEdge = hSEPairVsPtReb[corrName][iK].GetXaxis().GetBinUpEdge(iPt+1)
                ptFineMin = hSEPairVsPt[corrName][iK].GetXaxis().FindBin(ptLowEdge*1.0001)
                ptFineMax = hSEPairVsPt[corrName][iK].GetXaxis().FindBin(ptUpEdge*1.0001)
                contentAll, errorAll, contentSig, errorSig = (0. for _ in range(4))
                for iPtFine in range(ptFineMin, ptFineMax):
                    contentAll += hSEPairVsPt[corrName][iK].GetBinContent(iPtFine)
                    errorAll += hSEPairVsPt[corrName][iK].GetBinError(iPtFine)**2
                    contentSig += hSEPairSignalVsPt[corrName][iK].GetBinContent(iPtFine)
                    errorSig += hSEPairSignalVsPt[corrName][iK].GetBinError(iPtFine)**2
                hSEPairVsPtReb[corrName][iK].SetBinContent(iPt+1, contentAll)
                hSEPairVsPtReb[corrName][iK].SetBinError(iPt+1, np.sqrt(errorAll))
                hSEPairSignalVsPtReb[corrName][iK].SetBinContent(iPt+1, contentSig)
                hSEPairSignalVsPtReb[corrName][iK].SetBinError(iPt+1, np.sqrt(errorSig))
            SetObjectStyle(hSEPairVsPtReb[corrName][iK], markerstyle=kOpenCircle)
            SetObjectStyle(hSEPairSignalVsPtReb[corrName][iK], color=kRed+1, fillstyle=0)
        else:
            hSEPairVsPtReb[corrName].append(hSEPairVsPt[corrName][iK].Clone(
                f'{hSEPairVsPt[corrName][iK].GetName()}_kstar{iK}_reb'))
            hSEPairSignalVsPtReb[corrName].append(hSEPairSignalVsPt[corrName][iK].Clone(
                f'{hSEPairSignalVsPt[corrName][iK].GetName()}_kstar{iK}_reb'))

        hSEPairBkgVsPt[corrName].append(hSEPairVsPtReb[corrName][iK].Clone(
            f'hSEPairBkgVsPt_{corrName}_kstar{iK}_reb'))
        hSEPairBeautyVsPt[corrName].append(hSEPairVsPtReb[corrName][iK].Clone(
            f'hSEPairBeautyVsPt_{corrName}_kstar{iK}_reb'))
        hSEPairDirectVsPt[corrName].append(hSEPairVsPtReb[corrName][iK].Clone(
            f'hSEPairDirectVsPt_{corrName}_kstar{iK}_reb'))
        hSEPairFromDstarVsPt[corrName].append(hSEPairVsPtReb[corrName][iK].Clone(
            f'hSEPairFromDstarVsPt_{corrName}_kstar{iK}_reb'))

        SetObjectStyle(hSEPairBkgVsPt[corrName][iK], color=kOrange+7, fillstyle=0)
        SetObjectStyle(hSEPairBeautyVsPt[corrName][iK], color=kAzure+4, fillstyle=0)
        SetObjectStyle(hSEPairFromDstarVsPt[corrName][iK], color=kGreen+2, fillstyle=0)
        SetObjectStyle(hSEPairDirectVsPt[corrName][iK], color=kRed+1, fillstyle=0)

        allCounts, allCountsUnc, sgnCounts, sgnCountsUnc, bkgCounts, bkgCountsUnc, bCounts, bCountsUnc, \
            DstarCounts, DstarCountsUnc, directCounts, directCountsUnc = ([] for _ in range(12))
        bFrac, bFracUnc, DstarFrac, DstarFracUnc, DstarFracSysLow, DstarFracSysHigh = ([] for _ in range(6))

        for iPt in range(hSEPairVsPtReb[corrName][iK].GetNbinsX()):
            ptCent = hSEPairVsPtReb[corrName][iK].GetBinCenter(iPt+1)
            allCounts.append(hSEPairVsPtReb[corrName][iK].GetBinContent(iPt+1))
            allCountsUnc.append(hSEPairVsPtReb[corrName][iK].GetBinError(iPt+1))
            sgnCounts.append(hSEPairSignalVsPtReb[corrName][iK].GetBinContent(iPt+1))
            sgnCountsUnc.append(hSEPairSignalVsPtReb[corrName][iK].GetBinError(iPt+1))
            bkgCounts.append(allCounts[iPt] - sgnCounts[iPt])
            bkgCountsUnc.append(sgnCountsUnc[iPt])
            hSEPairBkgVsPt[corrName][iK].SetBinContent(iPt+1, bkgCounts[iPt])
            hSEPairBkgVsPt[corrName][iK].SetBinError(iPt+1, bkgCountsUnc[iPt])

            if evaluateFracFromBeauty:
                bFrac.append(hFracFromB.GetBinContent(iPt+1))
                bFracUnc.append(hFracFromB.GetBinError(iPt+1))
                bCounts.append(sgnCounts[iPt] * bFrac[iPt])
                bCountsUnc.append(np.sqrt((
                    bFracUnc[iPt]/bFrac[iPt])**2 + (sgnCountsUnc[iPt]/sgnCounts[iPt])**2) * bCounts[iPt])
                hSEPairBeautyVsPt[corrName][iK].SetBinContent(iPt+1, bCounts[iPt])
                hSEPairBeautyVsPt[corrName][iK].SetBinError(iPt+1, bCountsUnc[iPt])
            else:
                bFrac.append(0.)
                bFracUnc.append(0.)
                bCounts.append(0.)
                bCountsUnc.append(0.)
                hSEPairBeautyVsPt[corrName][iK].SetBinContent(iPt+1, 0.)
                hSEPairBeautyVsPt[corrName][iK].SetBinError(iPt+1, 0.)

            if evaluateFracFromDstar:
                ptBinDstar = hFracFromDstar.GetXaxis().FindBin(ptCent)
                for ptBinDstarGraph in range(gFracFromDstar.GetN()):
                    ptGraph, fracGraph = ctypes.c_double(), ctypes.c_double()
                    gFracFromDstar.GetPoint(ptBinDstarGraph, ptGraph, fracGraph)
                    ptGraph = ptGraph.value
                    if abs(ptGraph - ptCent) < 1.e-3:
                        break
                # Dstar fraction wrt prompt signal
                DstarFrac.append(hFracFromDstar.GetBinContent(ptBinDstar) * (1-bFrac[iPt]))
                DstarFracUnc.append(hFracFromDstar.GetBinError(ptBinDstar) * (1-bFrac[iPt]))
                DstarFracSysLow.append(gFracFromDstar.GetErrorYlow(ptBinDstarGraph) * (1-bFrac[iPt]))
                DstarFracSysHigh.append(gFracFromDstar.GetErrorYhigh(ptBinDstarGraph) * (1-bFrac[iPt]))
                DstarCounts.append(sgnCounts[iPt] * DstarFrac[iPt])
                DstarCountsUnc.append(np.sqrt((
                    DstarFracUnc[iPt]/DstarFrac[iPt])**2 + (sgnCountsUnc[iPt]/sgnCounts[iPt])**2) * DstarCounts[iPt])
                hSEPairFromDstarVsPt[corrName][iK].SetBinContent(iPt+1, DstarCounts[iPt])
                hSEPairFromDstarVsPt[corrName][iK].SetBinError(iPt+1, DstarCountsUnc[iPt])
            else:
                DstarCounts.append(0.)
                DstarCountsUnc.append(0.)
                hSEPairFromDstarVsPt[corrName][iK].SetBinContent(iPt+1, 0.)
                hSEPairFromDstarVsPt[corrName][iK].SetBinError(iPt+1, 0.)

            directCounts.append(allCounts[iPt] - (bkgCounts[iPt] + bCounts[iPt] + DstarCounts[iPt]))
            directCountsUnc.append(np.sqrt(bkgCountsUnc[iPt]**2 + bCountsUnc[iPt]**2 + DstarCountsUnc[iPt]**2))

            hSEPairDirectVsPt[corrName][iK].SetBinContent(iPt+1, directCounts[iPt])
            hSEPairDirectVsPt[corrName][iK].SetBinError(iPt+1, directCountsUnc[iPt])


        kstarPurity, purity = (ctypes.c_double() for _ in range(2))
        gAveragePurity[corrName].GetPoint(iK, kstarPurity, purity)
        purity = purity.value
        purityUnc = gAveragePurity[corrName].GetErrorYlow(iK)
        puritySystUncLow = gAveragePuritySyst[corrName].GetErrorYlow(iK)
        puritySystUncHigh = gAveragePuritySyst[corrName].GetErrorYhigh(iK)

        fracBeauty, fracBeautyUnc = ComputeWeightedAverage(bFrac, sgnCounts, bFracUnc, sgnCountsUnc)
        # propagate syst unc as pT correlated
        if systBeauty is not None:
            fracBeautyAvSysLow, _ = ComputeWeightedAverage(np.multiply(bFrac, (1-np.array(systBeauty[:len(bFrac)]))),
                                                           sgnCounts, bFracUnc, sgnCountsUnc)
            fracBeautyAvSysHigh, _ = ComputeWeightedAverage(np.multiply(bFrac, (1+np.array(systBeauty[:len(bFrac)]))),
                                                            sgnCounts, bFracUnc, sgnCountsUnc)
            fracBeautySysLow = fracBeauty - fracBeautyAvSysLow
            fracBeautySysHigh = fracBeautyAvSysHigh - fracBeauty
            if abs(fracBeautySysHigh-fracBeautySysLow) < 1.e-3 * fracBeauty:
                fracBeautySysHigh = fracBeautySysLow
            else:
                print('WARNING: asymmetric f_nonprompt systematics, might be unexpected')

        fracDstar, fracDstarUnc = ComputeWeightedAverage(DstarFrac, sgnCounts, DstarFracUnc, sgnCountsUnc)
        fracDstarAvSysLow, _ = ComputeWeightedAverage(np.add(DstarFrac, -np.array(DstarFracSysLow)),
                                                      sgnCounts, DstarFracUnc, sgnCountsUnc)
        fracDstarAvSysHigh, _ = ComputeWeightedAverage(np.add(DstarFrac, DstarFracSysHigh),
                                                       sgnCounts, DstarFracUnc, sgnCountsUnc)
        fracDstarAvSysLow = fracDstar - fracDstarAvSysLow
        fracDstarAvSysHigh = fracDstarAvSysHigh - fracDstar

        kStar = (kStarMax+kStarMin) / 2
        kStarDelta = (kStarMax-kStarMin) / 2

        gBeautyFrac[corrName].SetPoint(iK, kStar, fracBeauty)
        gDstarFrac[corrName].SetPoint(iK, kStar, fracDstar)

        gBeautyFrac[corrName].SetPointError(iK, kStarDelta, kStarDelta, fracBeautyUnc, fracBeautyUnc)
        gDstarFrac[corrName].SetPointError(iK, kStarDelta, kStarDelta, fracDstarUnc, fracDstarUnc)

        if iK == 0: # only for k* < 200 MeV/c^2
            hFractions[corrName].SetBinContent(1, purity)
            hFractions[corrName].SetBinContent(2, fracBeauty)
            hFractions[corrName].SetBinContent(3, fracDstar)
            hFractions[corrName].SetBinError(1, purityUnc)
            hFractions[corrName].SetBinError(2, fracBeautyUnc)
            hFractions[corrName].SetBinError(3, fracDstarUnc)
            gFractionsSyst[corrName].SetPoint(0, 1, purity)
            gFractionsSyst[corrName].SetPoint(1, 2, fracBeauty)
            gFractionsSyst[corrName].SetPoint(2, 3, fracDstar)
            if systPurity is not None:
                gFractionsSyst[corrName].SetPointError(0, 0.3, 0.3, puritySystUncLow, puritySystUncHigh)
            if systBeauty is not None:
                gFractionsSyst[corrName].SetPointError(1, 0.3, 0.3, fracBeautySysHigh, fracBeautySysLow)
            gFractionsSyst[corrName].SetPointError(2, 0.3, 0.3, fracDstarAvSysLow, fracDstarAvSysHigh)

# plots
legSoverB = TLegend(0.18, 0.7, 0.5, 0.85)
legSoverB.SetTextSize(0.045)
legSoverB.SetBorderSize(0)
legSoverB.SetFillStyle(0)
legSoverB.AddEntry(gAverageSoverB[corrName], 'weighted average', 'lp')
legSoverB.AddEntry(gSoverBAvPt[corrName], 'S/B (#LT #it{p}_{T}(D^{#pm}) #GT) - pol2 parm', 'lp')

legPurity = TLegend(0.18, 0.2, 0.5, 0.45)
legPurity.SetTextSize(0.045)
legPurity.SetBorderSize(0)
legPurity.SetFillStyle(0)
legPurity.AddEntry(gAverageSoverB[corrName], 'weighted average', 'lp')
legPurity.AddEntry(gSoverBAvPt[corrName], 'S/(S+B) (#LT #it{p}_{T}(D^{#pm}) #GT) - pol4 parm', 'lp')

cSoverBvsKstar = TCanvas('cSoverBvsKstar', '', 1000, 500)
cSoverBvsKstar.Divide(2, 1)
for iPad, corrName in enumerate(gAverageSoverB):
    cSoverBvsKstar.cd(iPad+1).DrawFrame(0., 0., 3., hSoverB.GetMaximum()*1.2, ';#it{k}* (GeV/#it{c});#LT S/B #GT')
    gAverageSoverB[corrName].Draw('p')
    gSoverBAvPt[corrName].Draw('p')
    legSoverB.Draw()
cSoverBvsKstar.Modified()
cSoverBvsKstar.Update()

cPurityvsKstar = TCanvas('cPurityvsKstar', '', 1000, 500)
cPurityvsKstar.Divide(2, 1)
for iPad, corrName in enumerate(gAveragePurity):
    cPurityvsKstar.cd(iPad+1).DrawFrame(0., 0., 3., hPurity.GetMaximum()*1.2, ';#it{k}* (GeV/#it{c});#LT S/(S+B) #GT')
    gAveragePurity[corrName].Draw('p')
    gPurityAvPt[corrName].Draw('p')
    legPurity.Draw()
cPurityvsKstar.Modified()
cPurityvsKstar.Update()

info = TLatex()
info.SetTextColor(kBlack)
info.SetTextSize(0.045)
info.SetTextFont(42)
info.SetNDC()

padOrder = {'Particle0_Particle2': 1,
            'Particle1_Particle3': 2,
            'Particle0_Particle2_plus_Particle1_Particle3': 3,
            'Particle0_Particle3': 4,
            'Particle1_Particle2': 5,
            'Particle0_Particle3_plus_Particle1_Particle2': 6}

cSEDistr = TCanvas('cSEDistr', '', 1200, 900)
cSEDistr.Divide(3, 2)
for corrName in hSEPairVsPt:
    cSEDistr.cd(padOrder[corrName])
    hSEPairVsPt[corrName][0].GetYaxis().SetRangeUser(0., hSEPairVsPt[corrName][0].GetMaximum()*2)
    hSEPairVsPt[corrName][0].DrawCopy('e')
    info.DrawLatex(0.2, 0.86, corrTitles[corrName])
    DName = ''
    if 'Particle2' in corrName and 'Particle3' in corrName:
        DName = 'D^{#pm}'
    elif 'Particle2' in corrName:
        DName = 'D^{+}'
    elif 'Particle3' in corrName:
        DName = 'D^{#font[122]{-}}'

    info.DrawLatex(0.2, 0.80, f'#LT #it{{p}}_{{T}} ({DName}) #GT = '
                   f'{hSEPairVsPt[corrName][0].GetMean():0.2f} #pm {hSEPairVsPt[corrName][0].GetMeanError():0.2f}')
    info.DrawLatex(0.2, 0.74, '#it{k}* < 200 MeV/#it{c}')
cSEDistr.Modified()
cSEDistr.Update()

cSoverBvsDistr = TCanvas('cSoverBvsDistr', '', 1000, 500)
cSoverBvsDistr.Divide(2, 1)
axisPairSoverB, fSoverBToDraw = {}, {}
for iPad, corrName in enumerate(['Particle0_Particle2_plus_Particle1_Particle3',
                                 'Particle0_Particle3_plus_Particle1_Particle2']):

    maxSoverB = hSoverB.GetMaximum()
    maxSEPair = hSEPairVsPt[corrName][0].GetMaximum()

    hFrameSoverB = cSoverBvsDistr.cd(iPad+1).DrawFrame(0., 0., ptMax, maxSEPair,
                                                       ';#it{p}_{T}(D^{#pm}) (GeV/#it{c});SE pairs')
    hFrameSoverB.GetYaxis().SetDecimals()
    hSEPairVsPt[corrName][0].DrawCopy('esame')
    hSoverBToDraw = hSoverB.Clone()
    hSoverBToDraw.Scale(maxSEPair / maxSoverB)
    fSoverBToDraw[corrName] = TF1(f'fSoverB_{corrName}', 'pol2', 0., 10.)
    fSoverBToDraw[corrName].SetLineColor(kRed-7)
    hSoverBToDraw.Fit(f'fSoverB_{corrName}', 'Q0')
    hSoverBToDraw.DrawCopy('same')
    fSoverBToDraw[corrName].Draw('same')
    kStar, avSoverB, SoverBavPt = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    gAverageSoverB[corrName].GetPoint(0, kStar, avSoverB)
    gSoverBAvPt[corrName].GetPoint(0, kStar, SoverBavPt)
    info.DrawLatex(0.2, 0.86, corrTitles[corrName])
    info.DrawLatex(0.2, 0.80, '#it{k}* < 200 MeV/#it{c}')
    info.DrawLatex(0.2, 0.74, f'#LT S/B #GT = '
                   f'{avSoverB.value:0.2f} #pm {gAverageSoverB[corrName].GetErrorYlow(0):0.2f}')
    info.DrawLatex(0.2, 0.68, 'S/B (#LT #it{p}_{T}(D^{#pm}) #GT) = '
                   f'{SoverBavPt.value:0.2f} #pm {gSoverBAvPt[corrName].GetErrorYlow(0):0.2f}')
    cSoverBvsDistr.cd(iPad+1).SetRightMargin(0.14)
    cSoverBvsDistr.cd(iPad+1).SetTicky(0)
    cSoverBvsDistr.cd(iPad+1).Modified()
    cSoverBvsDistr.cd(iPad+1).Update()
    axisPairSoverB[corrName] = TGaxis(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),
                                0., maxSoverB, 510, "+L")
    axisPairSoverB[corrName].SetLineColor(kRed+1)
    axisPairSoverB[corrName].SetLabelColor(kRed+1)
    axisPairSoverB[corrName].SetLabelFont(42)
    axisPairSoverB[corrName].SetLabelSize(0.045)
    axisPairSoverB[corrName].SetTitle('S/B')
    axisPairSoverB[corrName].SetTitleOffset(1.4)
    axisPairSoverB[corrName].SetLabelOffset(0.012)
    axisPairSoverB[corrName].SetTitleColor(kRed+1)
    axisPairSoverB[corrName].SetTitleFont(42)
    axisPairSoverB[corrName].SetTitleSize(0.05)
    axisPairSoverB[corrName].SetMaxDigits(3)
    axisPairSoverB[corrName].Draw()
cSoverBvsDistr.Modified()
cSoverBvsDistr.Update()

cPurityvsDistr = TCanvas('cPurityvsDistr', '', 1000, 500)
cPurityvsDistr.Divide(2, 1)
axisPairPurity, fPurityToDraw = {}, {}
for iPad, corrName in enumerate(['Particle0_Particle2_plus_Particle1_Particle3',
                                 'Particle0_Particle3_plus_Particle1_Particle2']):

    maxPurity = 2.
    maxSEPair = hSEPairVsPt[corrName][0].GetMaximum()

    hFramePurity = cPurityvsDistr.cd(iPad+1).DrawFrame(0., 0., ptMax, maxSEPair,
                                                       ';#it{p}_{T}(D^{#pm}) (GeV/#it{c});SE pairs')
    hFramePurity.GetYaxis().SetDecimals()
    hSEPairVsPt[corrName][0].DrawCopy('esame')
    hPurityToDraw = hPurity.Clone()
    hPurityToDraw.Scale(maxSEPair / maxPurity)
    fPurityToDraw[corrName] = TF1(f'fPurity_{corrName}', 'pol4', 0., 10.)
    fPurityToDraw[corrName].SetLineColor(kRed-7)
    hPurityToDraw.Fit(f'fPurity_{corrName}', 'Q0')
    hPurityToDraw.DrawCopy('esame')
    fPurityToDraw[corrName].Draw('same')
    kStar, avPurity, PurityavPt = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    gAveragePurity[corrName].GetPoint(0, kStar, avPurity)
    gPurityAvPt[corrName].GetPoint(0, kStar, PurityavPt)
    info.DrawLatex(0.2, 0.86, corrTitles[corrName])
    info.DrawLatex(0.2, 0.80, '#it{k}* < 200 MeV/#it{c}')
    info.DrawLatex(0.2, 0.74, f'#LT S/(S+B) #GT = '
                   f'{avPurity.value:0.2f} #pm {gAveragePurity[corrName].GetErrorYlow(0):0.2f}')
    info.DrawLatex(0.2, 0.68, 'S/(S+B) (#LT #it{p}_{T}(D^{#pm}) #GT) = '
                   f'{PurityavPt.value:0.2f} #pm {gPurityAvPt[corrName].GetErrorYlow(0):0.2f}')
    cPurityvsDistr.cd(iPad+1).SetRightMargin(0.14)
    cPurityvsDistr.cd(iPad+1).SetTicky(0)
    cPurityvsDistr.cd(iPad+1).Modified()
    cPurityvsDistr.cd(iPad+1).Update()
    axisPairPurity[corrName] = TGaxis(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),
                                0., maxPurity, 510, "+L")
    axisPairPurity[corrName].SetLineColor(kRed+1)
    axisPairPurity[corrName].SetLabelColor(kRed+1)
    axisPairPurity[corrName].SetLabelFont(42)
    axisPairPurity[corrName].SetLabelSize(0.045)
    axisPairPurity[corrName].SetTitle('S/(S+B)')
    axisPairPurity[corrName].SetTitleOffset(1.4)
    axisPairPurity[corrName].SetLabelOffset(0.012)
    axisPairPurity[corrName].SetTitleColor(kRed+1)
    axisPairPurity[corrName].SetTitleFont(42)
    axisPairPurity[corrName].SetTitleSize(0.05)
    axisPairPurity[corrName].SetMaxDigits(3)
    axisPairPurity[corrName].SetDecimals()
    axisPairPurity[corrName].Draw()
cPurityvsDistr.Modified()
cPurityvsDistr.Update()

cBeautyFracvsDistr = TCanvas('cBeautyFracvsDistr', '', 1000, 500)
cBeautyFracvsDistr.Divide(2, 1)
axisPairBeautyFrac, fBeautyFracToDraw = {}, {}
for iPad, corrName in enumerate(['Particle0_Particle2_plus_Particle1_Particle3',
                                 'Particle0_Particle3_plus_Particle1_Particle2']):

    maxBeautyFrac = 0.5
    maxSEPair = hSEPairSignalVsPtReb[corrName][0].GetMaximum()*2

    hFrameBeautyFrac = cBeautyFracvsDistr.cd(iPad+1).DrawFrame(0., 0., ptMax, maxSEPair,
                                                               ';#it{p}_{T}(D^{#pm}) (GeV/#it{c});SE pairs (signal)')
    hFrameBeautyFrac.GetYaxis().SetDecimals()
    SetObjectStyle(hSEPairSignalVsPtReb[corrName][0], markerstyle=kOpenCircle)
    hSEPairSignalVsPtReb[corrName][0].DrawCopy('esame')
    hBeautyFracToDraw = hFracFromB.Clone()
    hBeautyFracToDraw.Scale(maxSEPair / maxBeautyFrac)
    hBeautyFracToDraw.DrawCopy('esame')
    kStar, avBeautyFrac, BeautyFracavPt = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    gBeautyFrac[corrName].GetPoint(0, kStar, avBeautyFrac)
    gBeautyFrac[corrName].GetPoint(0, kStar, BeautyFracavPt)
    info.DrawLatex(0.2, 0.86, corrTitles[corrName])
    info.DrawLatex(0.2, 0.80, '#it{k}* < 200 MeV/#it{c}')
    info.DrawLatex(0.2, 0.74, f'#LT #it{{f}}_{{non-prompt}} #GT = '
                   f'{avBeautyFrac.value:0.3f} #pm {gBeautyFrac[corrName].GetErrorYlow(0):0.3f}')
    cBeautyFracvsDistr.cd(iPad+1).SetRightMargin(0.14)
    cBeautyFracvsDistr.cd(iPad+1).SetTicky(0)
    cBeautyFracvsDistr.cd(iPad+1).Modified()
    cBeautyFracvsDistr.cd(iPad+1).Update()
    axisPairBeautyFrac[corrName] = TGaxis(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),
                                0., maxBeautyFrac, 510, "+L")
    axisPairBeautyFrac[corrName].SetLineColor(kAzure+4)
    axisPairBeautyFrac[corrName].SetLabelColor(kAzure+4)
    axisPairBeautyFrac[corrName].SetLabelFont(42)
    axisPairBeautyFrac[corrName].SetLabelSize(0.045)
    axisPairBeautyFrac[corrName].SetTitle('#it{f}_{non-prompt}')
    axisPairBeautyFrac[corrName].SetTitleOffset(1.4)
    axisPairBeautyFrac[corrName].SetLabelOffset(0.012)
    axisPairBeautyFrac[corrName].SetTitleColor(kAzure+4)
    axisPairBeautyFrac[corrName].SetTitleFont(42)
    axisPairBeautyFrac[corrName].SetTitleSize(0.05)
    axisPairBeautyFrac[corrName].SetMaxDigits(3)
    axisPairBeautyFrac[corrName].SetDecimals()
    axisPairBeautyFrac[corrName].Draw()
cBeautyFracvsDistr.Modified()
cBeautyFracvsDistr.Update()

cDstarFracvsDistr = TCanvas('cDstarFracvsDistr', '', 1000, 500)
cDstarFracvsDistr.Divide(2, 1)
axisPairDstarFrac, fDstarFracToDraw = {}, {}
for iPad, corrName in enumerate(['Particle0_Particle2_plus_Particle1_Particle3',
                                 'Particle0_Particle3_plus_Particle1_Particle2']):

    maxDstarFrac = 0.5
    maxSEPair = hSEPairSignalVsPtReb[corrName][0].GetMaximum()*2

    hFrameDstarFrac = cDstarFracvsDistr.cd(iPad+1).DrawFrame(0., 0., ptMax, maxSEPair,
                                                               ';#it{p}_{T}(D^{#pm}) (GeV/#it{c});SE pairs (signal)')
    hFrameDstarFrac.GetYaxis().SetDecimals()
    hSEPairSignalVsPtReb[corrName][0].DrawCopy('esame')
    SetObjectStyle(hSEPairSignalVsPtReb[corrName][0], color=kRed+1, fillstyle=0)
    hDstarFracToDraw = hFracFromDstar.Clone()
    for iPt in range(hDstarFracToDraw.GetNbinsX()):
        hDstarFracToDraw.SetBinContent(iPt+1, hFracFromDstar.GetBinContent(iPt+1) * (1-hFracFromB.GetBinContent(iPt+1)))
    hDstarFracToDraw.SetBinContent(1, -1000.)
    SetObjectStyle(hDstarFracToDraw, color=kGreen+2, fillstyle=0)    
    hDstarFracToDraw.Scale(maxSEPair / maxDstarFrac)
    hDstarFracToDraw.DrawCopy('esame')
    kStar, avDstarFrac, DstarFracavPt = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    gDstarFrac[corrName].GetPoint(0, kStar, avDstarFrac)
    gDstarFrac[corrName].GetPoint(0, kStar, DstarFracavPt)
    info.DrawLatex(0.2, 0.86, corrTitles[corrName])
    info.DrawLatex(0.2, 0.80, '#it{k}* < 200 MeV/#it{c}')
    info.DrawLatex(0.2, 0.74, f'#LT D^{{#pm}} #leftarrow D*^{{#pm}} / D^{{#pm}} #GT = '
                   f'{avDstarFrac.value:0.3f} #pm {gDstarFrac[corrName].GetErrorYlow(0):0.3f}')
    cDstarFracvsDistr.cd(iPad+1).SetRightMargin(0.14)
    cDstarFracvsDistr.cd(iPad+1).SetTicky(0)
    cDstarFracvsDistr.cd(iPad+1).Modified()
    cDstarFracvsDistr.cd(iPad+1).Update()
    axisPairDstarFrac[corrName] = TGaxis(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),
                                0., maxDstarFrac, 510, "+L")
    axisPairDstarFrac[corrName].SetLineColor(kGreen+2)
    axisPairDstarFrac[corrName].SetLabelColor(kGreen+2)
    axisPairDstarFrac[corrName].SetLabelFont(42)
    axisPairDstarFrac[corrName].SetLabelSize(0.045)
    axisPairDstarFrac[corrName].SetTitle('D^{#pm} #leftarrow D*^{#pm} / D^{#pm}')
    axisPairDstarFrac[corrName].SetTitleOffset(1.4)
    axisPairDstarFrac[corrName].SetLabelOffset(0.012)
    axisPairDstarFrac[corrName].SetTitleColor(kGreen+2)
    axisPairDstarFrac[corrName].SetTitleFont(42)
    axisPairDstarFrac[corrName].SetTitleSize(0.05)
    axisPairDstarFrac[corrName].SetMaxDigits(3)
    axisPairDstarFrac[corrName].SetDecimals()
    axisPairDstarFrac[corrName].Draw()
cDstarFracvsDistr.Modified()
cDstarFracvsDistr.Update()

legContr = TLegend(0.6, 0.5, 0.9, 0.75)
legContr.SetTextSize(0.04)
legContr.SetBorderSize(0)
legContr.SetFillStyle(0)
legContr.AddEntry(hSEPairVsPtReb[corrName][0], 'measured', 'pl')
legContr.AddEntry(hSEPairBkgVsPt[corrName][0], 'bkg', 'pl')
legContr.AddEntry(hSEPairBeautyVsPt[corrName][0], 'b #rightarrow D^{#pm}', 'pl')
legContr.AddEntry(hSEPairFromDstarVsPt[corrName][0], 'c #rightarrow D*^{#pm} #rightarrow D^{#pm}', 'pl')
legContr.AddEntry(hSEPairDirectVsPt[corrName][0], 'c #rightarrow D^{#pm}', 'pl')

cSEPairsDiffContr, cLambda, cFractions = {}, {}, {}
for corrName in ['Particle0_Particle2_plus_Particle1_Particle3', 'Particle0_Particle3_plus_Particle1_Particle2']:

    if corrName == 'Particle0_Particle2_plus_Particle1_Particle3':
        suffix = 'DpPr'
    else:
        suffix = 'DmPr'

    cSEPairsDiffContr[corrName] = TCanvas(f'cSEPairsDiffContr_{suffix}', '', 500, 500)
    hFramePurity = cSEPairsDiffContr[corrName].cd().DrawFrame(0., 0., ptMax,
                                                              hSEPairVsPtReb[corrName][0].GetMaximum()*1.5,
                                                              ';#it{p}_{T}(D^{#pm}) (GeV/#it{c});SE pairs')
    hFramePurity.GetYaxis().SetDecimals()
    hSEPairVsPtReb[corrName][0].Draw('esame')
    hSEPairBkgVsPt[corrName][0].Draw('esame')
    hSEPairBeautyVsPt[corrName][0].Draw('esame')
    hSEPairFromDstarVsPt[corrName][0].Draw('esame')
    hSEPairDirectVsPt[corrName][0].Draw('esame')
    info.DrawLatex(0.5, 0.86, corrTitles[corrName])
    info.DrawLatex(0.5, 0.80, '#it{k}* < 200 MeV/#it{c}')
    legContr.Draw()
    cSEPairsDiffContr[corrName].Modified()
    cSEPairsDiffContr[corrName].Update()

    cFractions[corrName] = TCanvas(f'cFractions_{suffix}', '', 500, 500)
    hFractions[corrName].GetYaxis().SetRangeUser(0., 1.)
    hFractions[corrName].Draw('e')
    gFractionsSyst[corrName].Draw('2')
    info.DrawLatex(0.5, 0.86, corrTitles[corrName])
    info.DrawLatex(0.5, 0.80, '#it{k}* < 200 MeV/#it{c}')
    cFractions[corrName].Modified()
    cFractions[corrName].Update()

outFileName = 'Purity_vs_kstar' + outSuffix + '.root'
outFileName = os.path.join(outDirName, outFileName)
outFilePurity = TFile(outFileName, 'recreate')
if addTDirectory:
    outDir = TDirectoryFile(outSuffix, outSuffix)
    outDir.Write()
    outDir.cd()
cSEDistr.Write()
cSoverBvsDistr.Write()
cSoverBvsKstar.Write()
cPurityvsDistr.Write()
cPurityvsKstar.Write()
for corrName in gAverageSoverB:
    gAverageSoverB[corrName].Write()
    gSoverBAvPt[corrName].Write()
    gAveragePurity[corrName].Write()
    gPurityAvPt[corrName].Write()
if addTDirectory:
    outDir.Close()
outFilePurity.Close()
print('\nOutput file saved: ', outFileName)

outFileNamePDF = outFileName.replace('.root', '.pdf')
cSEDistr.SaveAs(outFileNamePDF.replace('Purity_vs_kstar', 'PairDistr'))
cSoverBvsDistr.SaveAs(outFileNamePDF.replace('Purity_vs_kstar', 'SoverB_vs_PairDistr_kstar200'))
cPurityvsDistr.SaveAs(outFileNamePDF.replace('Purity_vs_kstar', 'Purity_vs_PairDistr_kstar200'))
cBeautyFracvsDistr.SaveAs(outFileNamePDF.replace('Purity_vs_kstar', 'fBeauty_vs_PairDistr_kstar200'))
cDstarFracvsDistr.SaveAs(outFileNamePDF.replace('Purity_vs_kstar', 'fFromDstar_vs_PairDistr_kstar200'))
cSoverBvsKstar.SaveAs(outFileNamePDF.replace('Purity', 'SoverB'))
cPurityvsKstar.SaveAs(outFileNamePDF)

outFileName = 'Fractions' + outSuffix + '.root'
outFileName = os.path.join(outDirName, outFileName)
outFileFractions = TFile(outFileName, 'recreate')
if addTDirectory:
    outDir = TDirectoryFile(outSuffix, outSuffix)
    outDir.Write()
    outDir.cd()
for corrName in hFractions:
    hFractions[corrName].Write()
    gFractionsSyst[corrName].Write()
    cFractions[corrName].Write()
for corrName in hFractions:
    gAveragePurity[corrName].Write()
    gBeautyFrac[corrName].Write()
    gDstarFrac[corrName].Write()
if addTDirectory:
    outDir.Close()
outFileFractions.Close()
print('\nOutput file saved: ', outFileName)

for corrName in hFractions:
    if corrName == 'Particle0_Particle2_plus_Particle1_Particle3':
        suffix = 'DpPr'
    else:
        suffix = 'DmPr'
    outFileNamePDF = outFileName.replace('.root', f'_{suffix}.pdf')
    cFractions[corrName].SaveAs(outFileNamePDF)
    cSEPairsDiffContr[corrName].SaveAs(outFileNamePDF.replace('Fractions', 'SEparis_diffContr'))

input('Press enter to exit')
