'''
Script for the computation of the extra systematic uncertainty
on the theory-driven fprompt due to the multiplicity dependence
run: python ComputeTheoryDrivenFpromptMultSyst.py fracFile.root multDepFile.root
'''

import sys
import argparse
import ctypes
import numpy as np
import uproot
from scipy.interpolate import InterpolatedUnivariateSpline
from ROOT import TCanvas, TFile, TGraphAsymmErrors, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kAzure, kRed, kGreen # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error,no-name-in-module

# load inputs
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('fracFileName', metavar='text', default='fracFile.root',
                    help='root file with prompt fraction (from HFPtSpectum)')
parser.add_argument('multDepFileName', metavar='text', default='multDepFile.root',
                    help='root file containing B/D multiplicity dependence')
parser.add_argument('outFileName', metavar='text', default='outFile.root',
                    help='output root file')
parser.add_argument('--multFactor', type=float, default=1.,
                    help='value of multiplicity wrt to MB')
args = parser.parse_args()

# set global style
SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14, titleoffsety=1.7)

inFileFrac = TFile.Open(args.fracFileName)
gFrac = inFileFrac.Get('gFcConservative')
if not gFrac:
    print('ERROR: graph with prompt fraction not present in input file! Exit')
    sys.exit()
inFileFrac.Close()

sFracRatioVsMult, multFactor = [], []
tunes = [0, 3, 7, 8, 9]
for tune in tunes:
    fRatioVsMult = uproot.open(args.multDepFileName)[f'hfnonprompt_ratioMB_meson_tune{tune}']
    multCent = [(fRatioVsMult.edges[iBin]+fRatioVsMult.edges[iBin+1])/2 for iBin in range(len(fRatioVsMult.edges)-1)]
    sFracRatioVsMult.append(InterpolatedUnivariateSpline(multCent, fRatioVsMult.values))
    multFactor.append(sFracRatioVsMult[-1](args.multFactor))

gFracFONLLUnc, gFracMultUnc, gFracTotUnc = (TGraphAsymmErrors(0) for _ in range(3))
SetObjectStyle(gFracFONLLUnc, color=kAzure+4, fillalpha=0.2)
SetObjectStyle(gFracMultUnc, color=kGreen+2, fillalpha=0.2)
SetObjectStyle(gFracTotUnc, color=kRed+1, fillalpha=0.2)
ptMin, ptMax = -1., 1.
for iPt in range(gFrac.GetN()-1):

    promptFrac, pT = ctypes.c_double(), ctypes.c_double()
    gFrac.GetPoint(iPt+1, pT, promptFrac)
    FONLLUncLow = gFrac.GetErrorYlow(iPt+1)
    FONLLUncHigh = gFrac.GetErrorYhigh(iPt+1)
    ptUncLow = gFrac.GetErrorXlow(iPt+1)
    ptUncHigh = gFrac.GetErrorXhigh(iPt+1)

    nonPromptFrac = 1 - promptFrac.value
    nonPromptFracLow = nonPromptFrac * min(multFactor)
    nonPromptFracHigh = nonPromptFrac * max(multFactor)
    uncMultDepLow = nonPromptFrac-nonPromptFracLow
    uncMultDepHigh = nonPromptFracHigh-nonPromptFrac

    # low for FD is high for prompt and viceversa
    totUncLow = np.sqrt(FONLLUncLow**2 + uncMultDepHigh**2)
    totUncHigh = np.sqrt(FONLLUncHigh**2 + uncMultDepLow**2)

    gFracFONLLUnc.SetPoint(iPt, pT.value, promptFrac.value)
    gFracMultUnc.SetPoint(iPt, pT.value, promptFrac.value)
    gFracTotUnc.SetPoint(iPt, pT.value, promptFrac.value)
    gFracFONLLUnc.SetPointError(iPt, ptUncLow, ptUncHigh, FONLLUncLow, FONLLUncHigh)
    gFracMultUnc.SetPointError(iPt, ptUncLow, ptUncHigh, uncMultDepHigh, uncMultDepLow)
    gFracTotUnc.SetPointError(iPt, ptUncLow, ptUncHigh, totUncLow, totUncHigh)

    if iPt == 0:
        ptMin = pT.value - ptUncLow - 1
    if iPt == gFrac.GetN()-2:
        ptMax = pT.value + ptUncHigh + 1

leg = TLegend(0.4, 0.2, 0.7, 0.4)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextSize(0.045)
leg.AddEntry(gFracFONLLUnc, 'FONLL sys.', 'fp')
leg.AddEntry(gFracMultUnc, 'Mult. dep. sys.', 'fp')
leg.AddEntry(gFracTotUnc, 'Tot. sys.', 'fp')

cFrac = TCanvas('cFrac', '', 800, 800)
cFrac.DrawFrame(ptMin, 0.6, ptMax, 1., ';#it{p}_{T} (GeV/#it{c});#it{f}_{prompt}')
gFracTotUnc.Draw('2')
gFracTotUnc.Draw('p')
gFracFONLLUnc.Draw('2')
gFracFONLLUnc.Draw('p')
gFracMultUnc.Draw('2')
gFracMultUnc.Draw('p')
leg.Draw()

cFrac.SaveAs(args.outFileName.replace('.root', '.pdf'))

outFile = TFile(args.outFileName, 'recreate')
cFrac.Write()
gFracTotUnc.Write()
gFracFONLLUnc.Write()
gFracMultUnc.Write()
outFile.Close()

input('Press enter to exit')
