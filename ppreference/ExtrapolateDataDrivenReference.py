'''
python script for the pT extrapolation of data-driven pp references
run: python ExtrapolateDataDrivenCrossSection.py configFile.yml
'''

import sys
import argparse
import numpy as np
import yaml
from ROOT import TFile, TF1, TCanvas, TH1F, TLegend, TGraphAsymmErrors  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.AnalysisUtils import ComputeRatioDiffBins
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor, GetROOTMarker

SetGlobalStyle(padbottommargin=0.13)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCutSetFile:
    inCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)

inFile = TFile.Open(inCfg['infilename'])
hCrossSec = inFile.Get('hCrossSection')
hCrossSec.SetName('hCrossSecMeas')
hCrossSec.SetDirectory(0)
gCrossSec = inFile.Get('gCrossSectionSystTot')
gCrossSec.SetName('gCrossSectionSystTotMeas')
gCrossSecLumi = inFile.Get('gCrossSectionSystLumi')
systErr = inFile.Get('AliHFSystErr')
SetObjectStyle(hCrossSec, color=GetROOTColor('kRed+1'), markerstyle=GetROOTMarker('kFullSquare'), markersize=0.8)
SetObjectStyle(gCrossSec, color=GetROOTColor('kRed+1'), fillstyle=0)
inFile.Close()

ptMaxMeas = hCrossSec.GetXaxis().GetXmax()
ptLimsMeas = [hCrossSec.GetBinLowEdge(iPt) for iPt in range(1, hCrossSec.GetNbinsX()+1)]
ptLimsMeas.append(ptMaxMeas)
ptMinMeas = ptLimsMeas[0]

colors = {'cent': GetROOTColor('kBlack'), 'min': GetROOTColor('kRed+1'), 'max': GetROOTColor('kAzure+4')}
hFONLL, hDataOverFONLL, fDataOverFONLL = {}, {}, {}
inFileFONLL = TFile.Open(inCfg['extrap']['FONLLfilename'])
for var in ['cent', 'min', 'max']:
    hFONLL[var] = inFileFONLL.Get(inCfg['extrap']['hist'][var])
    hDataOverFONLL[var] = ComputeRatioDiffBins(hCrossSec, hFONLL[var])
    hDataOverFONLL[var].Scale(1.e6)
    SetObjectStyle(hDataOverFONLL[var], color=colors[var], fillstyle=0)
    fDataOverFONLL[var] = TF1(f'fDataOverFONLL{var}', 'pol0', ptMinMeas, ptMaxMeas)
    SetObjectStyle(fDataOverFONLL[var], color=colors[var])
    hDataOverFONLL[var].Fit(fDataOverFONLL[var], '0')

ptLimsForExtrap = inCfg['extrap']['ptlims']

if ptLimsForExtrap[0] != ptMaxMeas:
   print('ERROR: pT range of measurement and extrapolation do not match! Exit')
   sys.exit()

ptLimsAll = np.unique(ptLimsMeas + ptLimsForExtrap)
nPtBinsAll = len(ptLimsAll)-1
hFONLLReb = {}
for var in ['cent', 'min', 'max']:
    hFONLLReb[var] = hFONLL[var].Rebin(nPtBinsAll, f'hFONLLReb{var}', np.array(ptLimsAll, 'd'))
    hFONLLReb[var].Scale(hFONLL[var].GetBinWidth(1), 'width')

hCrossSecExtrap = TH1F('hCrossSection',
                       ';#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} #times BR (#mub GeV^{-1} #it{c})',
                       nPtBinsAll, np.array(ptLimsAll, 'd'))
gCrossSecExtrap = TGraphAsymmErrors(0) # put all systematics together
gCrossSecExtrap.SetName('gCrossSectionSystTot')
gCrossSecExtrapSyst = TGraphAsymmErrors(0)
gCrossSecExtrapSyst.SetName('gCrossSectionSystExtrap')
SetObjectStyle(hCrossSecExtrap, color=GetROOTColor('kBlack'), markersize=0.8)
SetObjectStyle(gCrossSecExtrap, color=GetROOTColor('kBlack'), fillstyle=0)
SetObjectStyle(gCrossSecExtrapSyst, color=GetROOTColor('kBlack'), fillstyle=0)

for iPt in range(1, hCrossSec.GetNbinsX()+1):
    ptCent = hCrossSec.GetBinCenter(iPt)
    crossSec = hCrossSec.GetBinContent(iPt)
    statUnc = hCrossSec.GetBinError(iPt)
    systLow = gCrossSec.GetErrorYlow(iPt-1)
    systHigh = gCrossSec.GetErrorYhigh(iPt-1)
    hCrossSecExtrap.SetBinContent(iPt, crossSec)
    hCrossSecExtrap.SetBinError(iPt, statUnc)
    gCrossSecExtrap.SetPoint(iPt-1, ptCent, crossSec)
    gCrossSecExtrapSyst.SetPoint(iPt-1, ptCent, crossSec)
    gCrossSecExtrap.SetPointError(iPt-1, 0.4, 0.4, systLow, systHigh)
    gCrossSecExtrapSyst.SetPointError(iPt-1, 0.4, 0.4, 0., 0.)
for iPt in range(hCrossSec.GetNbinsX()+1, hCrossSecExtrap.GetNbinsX()+1):
    ptCent = hCrossSecExtrap.GetBinCenter(iPt)
    extrapCross = hFONLLReb['cent'].GetBinContent(iPt) * fDataOverFONLL['cent'].GetParameter(0) * 1.e-6
    systLow = extrapCross - hFONLLReb['cent'].GetBinContent(iPt) * fDataOverFONLL['max'].GetParameter(0) * 1.e-6
    systHigh = hFONLLReb['cent'].GetBinContent(iPt) * fDataOverFONLL['min'].GetParameter(0) * 1.e-6 - extrapCross
    hCrossSecExtrap.SetBinContent(iPt, extrapCross)
    hCrossSecExtrap.SetBinError(iPt, 0.)
    gCrossSecExtrap.SetPoint(iPt-1, ptCent, extrapCross)
    gCrossSecExtrapSyst.SetPoint(iPt-1, ptCent, extrapCross)
    gCrossSecExtrap.SetPointError(iPt-1, 0.4, 0.4, systLow, systHigh)
    gCrossSecExtrapSyst.SetPointError(iPt-1, 0.4, 0.4, systLow, systHigh)

cFit = TCanvas('cFit', '', 500, 500)
hFrame = cFit.DrawFrame(ptMinMeas, 0., ptMaxMeas, hDataOverFONLL['min'].GetMaximum()*2,
                        ';#it{p}_{T} (GeV/#it{c});data/FONLL')
hFrame.GetYaxis().SetDecimals()
for var in ['cent', 'min', 'max']:
    hDataOverFONLL[var].Draw('same')
    fDataOverFONLL[var].Draw('same')

cCrossSec = TCanvas('cCrossSec', '', 500, 500)
cCrossSec.DrawFrame(ptMinMeas, hCrossSecExtrap.GetMinimum()/5, ptLimsForExtrap[-1], hCrossSecExtrap.GetMaximum()*5,
                    ';#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} #times BR (#mub GeV^{-1} #it{c})')
cCrossSec.SetLogy()
gCrossSec.Draw('2')
hCrossSec.Draw('same')
gCrossSecExtrap.Draw('2')
hCrossSecExtrap.Draw('same')

outFile = TFile(inCfg['outfilename'], 'recreate')
hCrossSecExtrap.Write()
gCrossSecExtrap.Write()
gCrossSecExtrapSyst.Write()
gCrossSecLumi.Write()
systErr.Write()
outFile.Close()

input('Press enter to exit')
