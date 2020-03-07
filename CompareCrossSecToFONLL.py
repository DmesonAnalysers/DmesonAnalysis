'''
python script to compare measured cross sections with FONLL
run: python CompareCrossSecToFONLL.py FONLL.root [--Dplus] [--Ds] [--prompt CrossSecPrompt.root] [--FD CrossSecFD.root] [--logx]
'''

import sys
import argparse
import numpy as np
from ROOT import TFile, TCanvas, TGraphAsymmErrors, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure, kFullCircle, kFullSquare # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('FONLLFileName', metavar='text',
                    default='FONLL.root', help='root file FONLL predictions')
parser.add_argument('--Dplus', action='store_true', default=False,
                    help='enable comparison for D+')
parser.add_argument('--Ds', action='store_true', default=False,
                    help='enable comparison for Ds')
parser.add_argument('--prompt', metavar='text', default=None,
                    help='enable comparison for prompt D and pass Cross section file name')
parser.add_argument('--FD', metavar='text', default=None,
                    help='enable comparison for FD D and pass Cross section file name')
parser.add_argument('--logx', action='store_true', default=False,
                    help='enable log scale for x axis')
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

SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14, titleoffsety=1.6, maxdigits=2, optstat=0)

if args.prompt:
    infilePrompt = TFile.Open(args.prompt)
    hCrossSectionPrompt = infilePrompt.Get('hCrossSection')
    if not hCrossSectionPrompt:
        hCrossSectionPrompt = infilePrompt.Get('histoSigmaCorr')
        hCrossSectionPrompt.Scale(BR/1.e6)
        hCrossSectionPrompt.SetStats(0)
    hCrossSectionPrompt.SetName('hCrossSectionPrompt')
    hCrossSectionPrompt.SetDirectory(0)
    SetObjectStyle(hCrossSectionPrompt, color=kBlack, markerstyle=kFullCircle)
    infilePrompt.Close()

    ptLimitsPrompt = np.array(hCrossSectionPrompt.GetXaxis().GetXbins(), 'd')
    ptMin = list(ptLimitsPrompt)[0]
    ptMax = list(ptLimitsPrompt)[-1]
    sigmaMin = hCrossSectionPrompt.GetMinimum()*0.2
    sigmaMax = hCrossSectionPrompt.GetMaximum()*5

    infileFONLL = TFile.Open(args.FONLLFileName)
    if args.Dplus:
        hFONLLPromptCentral = infileFONLL.Get('hDpluskpipipred_central')
        hFONLLPromptMin = infileFONLL.Get('hDpluskpipipred_min')
        hFONLLPromptMax = infileFONLL.Get('hDpluskpipipred_max')
    elif args.Ds:
        hFONLLPromptCentral = infileFONLL.Get('hDsKkpipred_central')
        hFONLLPromptMin = infileFONLL.Get('hDsKkpipred_min')
        hFONLLPromptMax = infileFONLL.Get('hDsKkpipred_max')
    hFONLLPromptCentral = hFONLLPromptCentral.Rebin(
        hCrossSectionPrompt.GetNbinsX(), 'hFONLLPromptCentral', ptLimitsPrompt)
    hFONLLPromptMin = hFONLLPromptMin.Rebin(hCrossSectionPrompt.GetNbinsX(), 'hFONLLPromptMin', ptLimitsPrompt)
    hFONLLPromptMax = hFONLLPromptMax.Rebin(hCrossSectionPrompt.GetNbinsX(), 'hFONLLPromptMax', ptLimitsPrompt)
    hFONLLPromptCentral.Scale(BR/2.e7, 'width')
    hFONLLPromptMin.Scale(BR/2.e7, 'width')
    hFONLLPromptMax.Scale(BR/2.e7, 'width')
    hFONLLPromptCentral.SetDirectory(0)
    hFONLLPromptMin.SetDirectory(0)
    hFONLLPromptMax.SetDirectory(0)
    hFONLLPromptCentral.SetStats(0)
    hFONLLPromptMin.SetStats(0)
    hFONLLPromptMax.SetStats(0)
    gFONLLPrompt = TGraphAsymmErrors(0)
    for iPt in range(hFONLLPromptCentral.GetNbinsX()):
        gFONLLPrompt.SetPoint(iPt, hFONLLPromptCentral.GetBinCenter(iPt+1), hFONLLPromptCentral.GetBinContent(iPt+1))
        gFONLLPrompt.SetPointError(iPt, hFONLLPromptCentral.GetBinWidth(iPt+1)/2,
                                   hFONLLPromptCentral.GetBinWidth(iPt+1)/2,
                                   hFONLLPromptCentral.GetBinContent(iPt+1)-hFONLLPromptMin.GetBinContent(iPt+1),
                                   hFONLLPromptMax.GetBinContent(iPt+1)-hFONLLPromptCentral.GetBinContent(iPt+1))
    SetObjectStyle(hFONLLPromptCentral, color=kRed+1, markerstyle=0)
    SetObjectStyle(gFONLLPrompt, color=kRed+1, fillalpha=0.2, fillstyle=1000)
    infileFONLL.Close()

if args.FD:
    infileFD = TFile.Open(args.FD)
    hCrossSectionFD = infileFD.Get('hCrossSection')
    if not hCrossSectionFD:
        hCrossSectionFD = infileFD.Get('histoSigmaCorr')
        hCrossSectionFD.SetStats(0)
    hCrossSectionFD.Scale(BR/1.e6)
    hCrossSectionFD.SetName('hCrossSectionFD')
    hCrossSectionFD.SetDirectory(0)
    SetObjectStyle(hCrossSectionFD, color=kBlack, markerstyle=kFullSquare)
    infileFD.Close()

    ptLimitsFD = np.array(hCrossSectionFD.GetXaxis().GetXbins(), 'd')
    if not args.prompt:
        ptMin = list(ptLimitsFD)[0]
        ptMax = list(ptLimitsFD)[-1]
        sigmaMin = hCrossSectionFD.GetBinContent(hCrossSectionFD.GetNbinsX())*0.2
        sigmaMax = hCrossSectionFD.GetBinContent(1)*5
    else:
        if list(ptLimitsFD)[0] < ptMin:
            ptMin = list(ptLimitsFD)[0]
        if list(ptLimitsFD)[-1] > ptMax:
            ptMax = list(ptLimitsFD)[-1]
        if hCrossSectionFD.GetBinContent(hCrossSectionFD.GetNbinsX())*0.2 < sigmaMin:
            sigmaMin = hCrossSectionFD.GetBinContent(hCrossSectionFD.GetNbinsX())*0.2
        if hCrossSectionFD.GetBinContent(1)*5 > sigmaMax:
            sigmaMax = hCrossSectionFD.GetBinContent(1)*5

    infileFONLL = TFile.Open(args.FONLLFileName)
    if args.Dplus:
        hFONLLFDCentral = infileFONLL.Get('hDpluskpipifromBpred_central_corr')
        hFONLLFDMin = infileFONLL.Get('hDpluskpipifromBpred_min_corr')
        hFONLLFDMax = infileFONLL.Get('hDpluskpipifromBpred_max_corr')
    elif args.Ds:
        hFONLLFDCentral = infileFONLL.Get('hDsPhipitoKkpifromBpred_central_corr')
        hFONLLFDMin = infileFONLL.Get('hDsPhipitoKkpifromBpred_min_corr')
        hFONLLFDMax = infileFONLL.Get('hDsPhipitoKkpifromBpred_max_corr')
    hFONLLFDCentral = hFONLLFDCentral.Rebin(hCrossSectionFD.GetNbinsX(), 'hFONLLFDCentral', ptLimitsFD)
    hFONLLFDMin = hFONLLFDMin.Rebin(hCrossSectionFD.GetNbinsX(), 'hFONLLFDMin', ptLimitsFD)
    hFONLLFDMax = hFONLLFDMax.Rebin(hCrossSectionFD.GetNbinsX(), 'hFONLLFDMax', ptLimitsFD)
    hFONLLFDCentral.Scale(BR/2.e7, 'width')
    hFONLLFDMin.Scale(BR/2.e7, 'width')
    hFONLLFDMax.Scale(BR/2.e7, 'width')
    hFONLLFDCentral.SetDirectory(0)
    hFONLLFDMin.SetDirectory(0)
    hFONLLFDMax.SetDirectory(0)
    hFONLLFDCentral.SetStats(0)
    hFONLLFDMin.SetStats(0)
    hFONLLFDMax.SetStats(0)
    gFONLLFD = TGraphAsymmErrors(0)
    for iPt in range(hFONLLFDCentral.GetNbinsX()):
        gFONLLFD.SetPoint(iPt, hFONLLFDCentral.GetBinCenter(iPt+1), hFONLLFDCentral.GetBinContent(iPt+1))
        gFONLLFD.SetPointError(iPt, hFONLLFDCentral.GetBinWidth(iPt+1)/2,
                               hFONLLFDCentral.GetBinWidth(iPt+1)/2,
                               hFONLLFDCentral.GetBinContent(iPt+1)-hFONLLFDMin.GetBinContent(iPt+1),
                               hFONLLFDMax.GetBinContent(iPt+1)-hFONLLFDCentral.GetBinContent(iPt+1))
    SetObjectStyle(hFONLLFDCentral, color=kAzure+4, markerstyle=0)
    SetObjectStyle(gFONLLFD, color=kAzure+4, fillalpha=0.2, fillstyle=1000)
    infileFONLL.Close()

cCrossSec = TCanvas('cCrossSec', '', 700, 800)
hFrame = cCrossSec.DrawFrame(ptMin, sigmaMin, ptMax, sigmaMax,
                             ';#it{p}_{T} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (#mub GeV #it{c}^{-1})')
cCrossSec.SetLogy()
if args.logx:
    cCrossSec.SetLogx()


if args.logx:
    legPrompt = TLegend(0.65, 0.75, 0.85, 0.9)
    legFD = TLegend(0.25, 0.2, 0.45, 0.35)
else:
    legPrompt = TLegend(0.55, 0.75, 0.75, 0.9)
    legFD = TLegend(0.55, 0.55, 0.75, 0.7)
legPrompt.SetTextSize(0.045)
legPrompt.SetBorderSize(0)
legPrompt.SetFillStyle(0)
legFD.SetTextSize(0.045)
legFD.SetBorderSize(0)
legFD.SetFillStyle(0)

if args.prompt:
    gFONLLPrompt.Draw('2')
    hFONLLPromptCentral.Draw('same')
    hCrossSectionPrompt.Draw('esame')
    legPrompt.AddEntry('', f'Prompt {mesonName}', '')
    legPrompt.AddEntry(hCrossSectionPrompt, 'Data', 'p')
    legPrompt.AddEntry(gFONLLPrompt, 'FONLL', 'f')
    legPrompt.Draw()
if args.FD:
    gFONLLFD.Draw('2')
    hFONLLFDCentral.Draw('same')
    hCrossSectionFD.Draw('esame')
    legFD.AddEntry('', f'Non-prompt {mesonName}', '')
    legFD.AddEntry(hCrossSectionFD, 'Data', 'p')
    legFD.AddEntry(gFONLLFD, 'FONLL', 'f')
    legFD.Draw()

input('Press enter to exit')
