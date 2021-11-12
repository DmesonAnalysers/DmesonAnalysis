'''
Script for the computation of the data-driven fraction ratio
run: python ComputeDataDrivenFracRatio.py config.yml
'''

import os
import argparse
import numpy as np
import yaml
from ROOT import TCanvas, TFile, TLine, TLegend # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor, GetROOTMarker #pylint: disable=wrong-import-position,import-error,no-name-in-module

SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.16,
               padtopmargin=0.035, titleoffsety=1.4, maxdigits=3)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('config', metavar='text', default='config.yml',
                    help='yaml config file')
args = parser.parse_args()

with open(args.config, 'r') as ymlCfg: # pylint: disable=unspecified-encoding
    cfg = yaml.load(ymlCfg, yaml.FullLoader)

inFileNames, colors, correl, labels, titles = ([] for _ in range(5))
cfgDen = cfg['inputs']['denominator']
inFileNames.append(cfgDen['file'])
colors.append(cfgDen['color'])
labels.append(cfgDen['label'])
titles.append(cfgDen['title'])
correl.append(None)
cfgNum = cfg['inputs']['numerator']
for fileName, color, label, title, corr in zip(
        cfgNum['files'],
        cfgNum['colors'],
        cfgNum['labels'],
        cfgNum['titles'],
        cfgNum['corr']
    ):
    inFileNames.append(fileName)
    colors.append(color)
    labels.append(label)
    titles.append(title)
    correl.append(corr)

hFracPrompt, hFracNonPrompt, hRatioFracPrompt, hRatioFracNonPrompt = [], [], [], []
for iFile, fileName in enumerate(inFileNames):

    inFile = TFile.Open(fileName)

    hFracPrompt.append(inFile.Get('hPromptFracCorr'))
    hFracPrompt[iFile].SetDirectory(0)
    hFracPrompt[iFile].SetName(f'hFracPrompt_{labels[iFile]}')
    SetObjectStyle(hFracPrompt[iFile], color=GetROOTColor(colors[iFile]))

    hFracNonPrompt.append(inFile.Get('hFDFracCorr'))
    hFracNonPrompt[iFile].SetDirectory(0)
    hFracNonPrompt[iFile].SetName(f'hFracNonPrompt_{labels[iFile]}')
    SetObjectStyle(hFracNonPrompt[iFile], color=GetROOTColor(colors[iFile]))

    if iFile == 0:
        hRatioFracPrompt.append(None)
        hRatioFracNonPrompt.append(None)
    else:
        hRatioFracPrompt.append(
            hFracPrompt[iFile].Clone(f'hRatioFracPrompt_{labels[iFile]}')
        )
        hRatioFracNonPrompt.append(
            hFracNonPrompt[iFile].Clone(f'hRatioFracNonPrompt_{labels[iFile]}')
        )
        hRatioFracPrompt[iFile].Divide(hFracPrompt[0])
        hRatioFracNonPrompt[iFile].Divide(hFracNonPrompt[0])
        hRatioFracPrompt[iFile].SetDirectory(0)
        hRatioFracNonPrompt[iFile].SetDirectory(0)
        if correl[iFile]:
            for iPt in range(1, hFracPrompt[iFile].GetNbinsX()+1):
                uncMult = hFracPrompt[iFile].GetBinError(iPt)
                uncMB = hFracPrompt[0].GetBinError(iPt)
                valmult = hFracPrompt[iFile].GetBinContent(iPt)
                valMB = hFracPrompt[0].GetBinContent(iPt)
                rho = uncMult / uncMB
                if rho > 1.:
                    rho = 1. / rho
                uncRatio = np.sqrt(
                    uncMult**2 / valMB**2 + valmult**2 / valMB**4 * uncMB**2 -
                    2 * rho * valmult / valMB**3 * uncMult * uncMB
                )
                hRatioFracPrompt[iFile].SetBinError(iPt, uncRatio)

                uncMult = hFracNonPrompt[iFile].GetBinError(iPt)
                uncMB = hFracNonPrompt[0].GetBinError(iPt)
                valmult = hFracNonPrompt[iFile].GetBinContent(iPt)
                valMB = hFracNonPrompt[0].GetBinContent(iPt)
                rho = uncMult / uncMB
                if rho > 1.:
                    rho = 1. / rho
                uncRatio = np.sqrt(
                    uncMult**2 / valMB**2 + valmult**2 / valMB**4 * uncMB**2 -
                    2 * rho * valmult / valMB**3 * uncMult * uncMB
                )
                hRatioFracNonPrompt[iFile].SetBinError(iPt, uncRatio)

    inFile.Close()

nPtBins = hFracPrompt[0].GetNbinsX()
ptMin = hFracPrompt[0].GetBinLowEdge(1)
ptMax = hFracPrompt[0].GetXaxis().GetBinUpEdge(nPtBins)

leg = TLegend(0.2, 0.2, 0.4, 0.4)
leg.SetTextSize(0.045)
leg.SetBorderSize(0)
leg.SetFillStyle(0)

cFracs = TCanvas('cFrac', '', 800, 800)
cFracs.Divide(2, 2)
hFrame = cFracs.cd(1).DrawFrame(
    ptMin,
    0.7,
    ptMax,
    1.,
    ';#it{p}_{T} (GeV/#it{c});#it{f}_{prompt}'
)
hFrame.GetYaxis().SetDecimals()
for iHist, histo in enumerate(hFracPrompt):
    leg.AddEntry(histo, titles[iHist], 'lp')
    histo.DrawCopy('same')
leg.Draw()

hFrame = cFracs.cd(2).DrawFrame(
    ptMin,
    0.,
    ptMax,
    0.2,
    ';#it{p}_{T} (GeV/#it{c});#it{f}_{non-prompt}'
)
hFrame.GetYaxis().SetDecimals()
hFrame.GetYaxis().SetNdivisions(505)
for histo in hFracNonPrompt:
    histo.DrawCopy('same')

line = TLine(ptMin, 1., ptMax, 1.)
line.SetLineWidth(2)
line.SetLineColor(GetROOTColor('kGrey+2'))
line.SetLineStyle(9)

hFrame = cFracs.cd(3).DrawFrame(
    ptMin,
    0.9,
    ptMax,
    1.1,
    f';#it{{p}}_{{T}} (GeV/#it{{c}});#it{{f}}_{{prompt}} / #it{{f}}_{{prompt}}^{{{  titles[0]}}}'
)
hFrame.GetYaxis().SetDecimals()
line.Draw()
for histo in hRatioFracPrompt:
    if histo:
        histo.DrawCopy('same')

hFrame = cFracs.cd(4).DrawFrame(
    ptMin,
    0.,
    ptMax,
    2.,
    f';#it{{p}}_{{T}} (GeV/#it{{c}});#it{{f}}_{{non-prompt}} / #it{{f}}_{{non-prompt}}^{{{  titles[0]}}}'
)
hFrame.GetYaxis().SetDecimals()
hFrame.GetYaxis().SetNdivisions(505)
line.Draw()
for histo in hRatioFracNonPrompt:
    if histo:
        histo.DrawCopy('same')

cFracs.Modified()
cFracs.Update()

cFracs.SaveAs(os.path.join(cfg['output']['directory'], cfg['output']['filewoext']) + '.pdf')

input('Press enter to exit')
