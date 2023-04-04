'''
Script to fit and combine the correlation functions of the left and right sidebands.
'''

import sys
import argparse
import re
import os
import ctypes
import yaml
import numpy as np
from ROOT import TCanvas, TFile, gStyle, TF1, TGraphErrors, TGraphAsymmErrors, gRandom, TVirtualFitter, kBlue, kOrange, TLegend, TMath, kRed, \
    TSpline3
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text',
                    help='yaml config file name')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlFile:
    cfg = yaml.load(ymlFile, yaml.FullLoader)
spec = cfg['input']['specifier']
outputFileName = cfg['output']
fitMethod = cfg['fit']['method']

if 'PIplusDminus' in spec or 'PIplusDplus' in spec:
    suffixData = '_DPi'
elif 'KplusDminus' in spec or 'KplusDplus' in spec:
    suffixData = '_DK'
else:
    print(
        f'\033[91mError:\033[0m the channel {spec} is not implemented. Exit!')
    sys.exit()

# open files with CF
fSBL = TFile(os.path.join(cfg['input']['sideband_dir'], f'CFOutput_{spec}_SBLeft.root'))
fSBR = TFile(os.path.join(cfg['input']['sideband_dir'], f'CFOutput_{spec}_SBRight.root'))
rebin = cfg['rebin']
gSBL = fSBL.Get(f'Graph_from_hCk_Reweighted_{rebin}MeV')
gSBR = fSBR.Get(f'Graph_from_hCk_Reweighted_{rebin}MeV')
gSBL.SetTitle(';#it{k}* (MeV/c);C(#it{k}*)')
gSBR.SetTitle(';#it{k}* (MeV/c);C(#it{k}*)')

fitRangeMin = cfg['fit']['xrange'][0]
fitRangeMax = cfg['fit']['xrange'][1]
deg = int(re.findall(r'\d{1}', fitMethod)[0])

if 'legendre' in fitMethod:
    # build series of legendre polynomials from degree 0 to l
    fitFunctionStr = ''
    for l in range(deg + 1):
        if l > 0:
            fitFunctionStr += ' + '

        # For unclear reasons, the fit only works if fitRangeMin=0
        shift = (fitRangeMax + fitRangeMin) / 2
        semiRange = (fitRangeMax - fitRangeMin) / 2
        fitFunctionStr += f'[{l}]*ROOT::Math::legendre({l}, (x - {shift})/{semiRange})'
elif 'pol' in fitMethod:
    fitFunctionStr = fitMethod
else:
    print(
        f'\033[91mError:\033[0m the fit mehod {fitMethod} is not implemented. Exit!')
    sys.exit()

fitFuncSBL = TF1('fitFuncSBL', fitFunctionStr, fitRangeMin, fitRangeMax)
fitFuncSBR = TF1('fitFuncSBR', fitFunctionStr, fitRangeMin, fitRangeMax)

for l in range(0, deg + 1):
    fitFuncSBL.SetParLimits(l, cfg['fit']['parmin'][l], cfg['fit']['parmax'][l])
    fitFuncSBL.SetParameter(l, (cfg['fit']['parmin'][l] + cfg['fit']['parmax'][l]) / 2)
    fitFuncSBR.SetParLimits(l, cfg['fit']['parmin'][l], cfg['fit']['parmax'][l])
    fitFuncSBR.SetParameter(l, (cfg['fit']['parmin'][l] + cfg['fit']['parmax'][l]) / 2)

# Create a TGraphErrors to hold the confidence intervals
grintSBL, grintSBR = TGraphErrors(1), TGraphErrors(1)
grintSBL.SetTitle('Fit to CF SB Left with 1#sigma confidence band')
grintSBL.SetName('gSBLFit')
grintSBR.SetTitle('Fit to CF SB Right with 1#sigma confidence band')
grintSBR.SetName('gSBRFit')
nPointsCLSBL = 0
nPointsCLSBR = 0

for i in range(gSBL.GetN()):
    if fitRangeMin <= gSBL.GetX()[i] and gSBL.GetX()[i] <= fitRangeMax:
        grintSBL.SetPoint(nPointsCLSBL, gSBL.GetX()[i], 0)
        nPointsCLSBL += 1

gSBL.Fit('fitFuncSBL', 'MR+', '', fitRangeMin, fitRangeMax)
TVirtualFitter.GetFitter().GetConfidenceIntervals(
    grintSBL, TMath.Erf(1. / TMath.Sqrt(2)))  # 1sigma conf. band

for i in range(gSBR.GetN()):
    if fitRangeMin <= gSBR.GetX()[i] and gSBR.GetX()[i] <= fitRangeMax:
        grintSBR.SetPoint(nPointsCLSBR, gSBR.GetX()[i], 0)
        nPointsCLSBR += 1
gSBR.Fit('fitFuncSBR', 'MR+', '', fitRangeMin, fitRangeMax)
TVirtualFitter.GetFitter().GetConfidenceIntervals(grintSBR, TMath.Erf(1 / TMath.Sqrt(2)))  # 1sigma symmetric conf. band

# Interpolate confidence bands with splines
nSBL = grintSBL.GetN()

yConfBandSBL = (ctypes.c_double * nSBL)(*[grintSBL.GetErrorYhigh(i) for i in range(nSBL)])  # Cast python list to ctypes
splineUncSBL = TSpline3('splUncCFSBL', grintSBL.GetX(), yConfBandSBL, nSBL)

nSBR = grintSBR.GetN()
yConfBandSBR = (ctypes.c_double * nSBR)(*[grintSBR.GetErrorYlow(i) for i in range(nSBL)])
splineUncSBR = TSpline3('splUncCFSBR', grintSBR.GetX(), yConfBandSBR, nSBR)

# Combination of the SBL and SBR CFs
combMethod = cfg['sidebands']['combmethod']
if combMethod == 'kAverage':
    weightSBL = cfg['sidebands']['weightleftsb']
    weightSBR = cfg['sidebands']['weightrightsb']

    sumWeights = weightSBL + weightSBR
    if sumWeights - 1 > 1.e-9:
        print(f'\033[93mWarning:\033[0m The sum of the weights is not 1. Renormalisinsg them.')
        weightSBL /= sumWeights
        weightSBR /= sumWeights

    def modelSBCombFunc(x, par):
        return weightSBL * fitFuncSBL.Eval(x[0]) + weightSBR * fitFuncSBR.Eval(x[0])
    modelSBComb = TF1('funcCFSBComb', modelSBCombFunc, fitRangeMin, fitRangeMax)

    def modelSBCombFuncUnc(x, par):
        return np.sqrt((weightSBL * splineUncSBL.Eval(x[0])) ** 2 + (weightSBR * splineUncSBR.Eval(x[0])) ** 2)
    modelSBUnc = TF1('funcCFSBCombUnc', modelSBCombFuncUnc, fitRangeMin, fitRangeMax)

else:
    print(
        f'\033[91mError:\033[0m the combination mehod {combMethod} is not implemented. Exit!')
    sys.exit()

# Make plots and drawings
oFile = TFile(f'{outputFileName}.root', 'RECREATE')
gStyle.SetOptFit(11111111)
SetGlobalStyle(titlesize=0.04, labelsize=0.035)

# Save Sideband template
modelSBComb.Write()
modelSBUnc.Write()

cSBCombined = TCanvas('cSBCombined', 'Left and right combined', 600, 600)
modelSBComb.SetTitle(gSBL.GetTitle())  # todo: check if works
modelSBComb.Draw()


def modelSBCombFuncUncHigh(x, par):
    return modelSBComb.Eval(x[0]) + modelSBUnc.Eval(x[0])


modelSBCombUncHigh = TF1('funcCFSBCombUncHigh', modelSBCombFuncUncHigh, fitRangeMin, fitRangeMax)


def modelSBCombFuncUncLow(x, par):
    return modelSBComb.Eval(x[0]) - modelSBUnc.Eval(x[0])


modelSBCombUncLow = TF1('funcCFSBCombUncLow', modelSBCombFuncUncLow, fitRangeMin, fitRangeMax)

modelSBCombUncHigh.SetLineStyle(9)
modelSBCombUncHigh.SetLineColor(kBlue)
modelSBCombUncHigh.Draw('same')
modelSBCombUncLow.SetLineStyle(9)
modelSBCombUncLow.SetLineColor(kBlue)
modelSBCombUncLow.Draw('same')

legComb = TLegend(0.3, 0.12, '')
legComb.AddEntry(modelSBComb, 'Combined SB')
legComb.AddEntry(modelSBCombUncHigh, '1#sigma conf. band', 'l')
legComb.Draw('same')

# Left sideband
cSBL = TCanvas('cCFSBL', 'Correlation function SBL', 600, 600)
SetObjectStyle(gSBL)
gSBL.GetYaxis().SetRangeUser(cfg['fit']['yrange'][0], cfg['fit']['yrange'][1])
gSBL.Draw('alpe')
fitFuncSBL.SetLineColor(kRed)
fitFuncSBL.SetMarkerColor(kRed)
grintSBL.SetLineColor(kOrange + 10)
grintSBL.SetMarkerColor(kOrange + 10)
grintSBL.SetFillColorAlpha(kOrange + 6, 0.7)
grintSBL.Draw('3 same')

legSBL = TLegend(0.35, 0.15, '')
legSBL.AddEntry(fitFuncSBL, 'Leg. Pol. fit')
legSBL.AddEntry(grintSBL, '1#sigma CL')
legSBL.AddEntry(gSBL, 'SB left')
legSBL.Draw()
cSBL.Print(f'{outputFileName}_fitSBL.pdf')
fitFuncSBL.Write()
grintSBL.Write()

# Right sideband
cSBR = TCanvas('cCFSBR', 'Correlation function SBR', 600, 600)
SetObjectStyle(gSBR)
gSBR.GetYaxis().SetRangeUser(cfg['fit']['yrange'][0], cfg['fit']['yrange'][1])
gSBR.Draw('alpe')
fitFuncSBR.SetLineColor(kRed)
fitFuncSBR.SetMarkerColor(kRed)
grintSBR.SetLineColor(kOrange + 10)
grintSBR.SetMarkerColor(kOrange + 10)
grintSBR.SetFillColorAlpha(kOrange + 6, 0.7)
grintSBR.Draw('3 same')

legSBR = TLegend(0.35, 0.15, '')
legSBR.AddEntry(fitFuncSBR, 'Leg. Pol. fit')
legSBR.AddEntry(grintSBR, '1#sigma CL')
legSBR.AddEntry(gSBR, 'SB right')
legSBR.Draw()
cSBR.Print(f'{outputFileName}_fitSBR.pdf')
fitFuncSBR.Write()
grintSBR.Write()
oFile.Close()
