'''
Script for the computation of the data-driven fraction ratio
run: python ComputeDataDrivenFracRatio.py config.yml
'''

import os
import argparse
import numpy as np
import yaml
from ROOT import TCanvas, TFile, TLine, TLegend, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor #pylint: disable=wrong-import-position,import-error,no-name-in-module

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

rawYieldSyst, cutSetSyst, multGenSyst, ptGenSyst = ([] for _ in range(4))
rawYieldSystRat, cutSetSystRat, multGenSystRat, ptGenSystRat = ([] for _ in range(4))
cfgSystDen = cfg['systematics']['fractions']['denominator']
rawYieldSyst.append(cfgSystDen['rawyield'])
cutSetSyst.append(cfgSystDen['cutsets'])
multGenSyst.append(cfgSystDen['multweights'])
ptGenSyst.append(cfgSystDen['ptweights'])
rawYieldSystRat.append(None)
cutSetSystRat.append(None)
multGenSystRat.append(None)
ptGenSystRat.append(None)
cfgSystNum = cfg['systematics']['fractions']['numerator']
cfgSystRat = cfg['systematics']['ratios']
for raw, cut, mult, pt, rawR, cutR, multR, ptR in zip(
        cfgSystNum['rawyield'],
        cfgSystNum['cutsets'],
        cfgSystNum['multweights'],
        cfgSystNum['ptweights'],
        cfgSystRat['rawyield'],
        cfgSystRat['cutsets'],
        cfgSystRat['multweights'],
        cfgSystRat['ptweights']
    ):
    rawYieldSyst.append(raw)
    cutSetSyst.append(cut)
    multGenSyst.append(mult)
    ptGenSyst.append(pt)
    rawYieldSystRat.append(rawR)
    cutSetSystRat.append(cutR)
    multGenSystRat.append(multR)
    ptGenSystRat.append(ptR)

hFracPrompt, hFracNonPrompt, hRatioFracPrompt, hRatioFracNonPrompt = ([] for _ in range(4))
gFracNonPrompt, gRatioFracNonPrompt = ([] for _ in range(2))
hFracNonPromptRawYieldSys, hFracNonPromptCutSetsSys, \
    hFracNonPromptMultWeightsSys, hFracNonPromptPtWeightsSys = ([] for _ in range(4))
hRatioFracNonPromptRawYieldSys, hRatioFracNonPromptCutSetsSys, \
    hRatioFracNonPromptMultWeightsSys, hRatioFracNonPromptPtWeightsSys = ([] for _ in range(4))

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

    hFracNonPromptRawYieldSys.append(
        hFracNonPrompt[iFile].Clone(f'hFracNonPromptRawYieldSys_{labels[iFile]}')
    )
    hFracNonPromptCutSetsSys.append(
        hFracNonPrompt[iFile].Clone(f'hFracNonPromptCutSetsSys_{labels[iFile]}')
    )
    hFracNonPromptMultWeightsSys.append(
        hFracNonPrompt[iFile].Clone(f'hFracNonPromptMultWeightsSys_{labels[iFile]}')
    )
    hFracNonPromptPtWeightsSys.append(
        hFracNonPrompt[iFile].Clone(f'hFracNonPromptPtWeightsSys_{labels[iFile]}')
    )
    hFracNonPromptRawYieldSys[iFile].SetDirectory(0)
    hFracNonPromptCutSetsSys[iFile].SetDirectory(0)
    hFracNonPromptMultWeightsSys[iFile].SetDirectory(0)
    hFracNonPromptPtWeightsSys[iFile].SetDirectory(0)

    gFracNonPrompt.append(TGraphAsymmErrors(1))
    gFracNonPrompt[iFile].SetName(f'gFracNonPrompt_{labels[iFile]}')
    SetObjectStyle(gFracNonPrompt[iFile], color=GetROOTColor(colors[iFile]), fillstyle=0)

    for iPt in range(0, hFracNonPrompt[iFile].GetNbinsX()):
        pt = hFracNonPrompt[iFile].GetBinCenter(iPt+1)
        frac = hFracNonPrompt[iFile].GetBinContent(iPt+1)
        totSyst = np.sqrt(
            rawYieldSyst[iFile][iPt]**2 + cutSetSyst[iFile][iPt]**2 + \
                multGenSyst[iFile][iPt]**2 + ptGenSyst[iFile][iPt]**2
        ) * frac
        hFracNonPromptRawYieldSys[iFile].SetBinContent(iPt+1, rawYieldSyst[iFile][iPt])
        hFracNonPromptCutSetsSys[iFile].SetBinContent(iPt+1, cutSetSyst[iFile][iPt])
        hFracNonPromptMultWeightsSys[iFile].SetBinContent(iPt+1, multGenSyst[iFile][iPt])
        hFracNonPromptPtWeightsSys[iFile].SetBinContent(iPt+1, ptGenSyst[iFile][iPt])
        gFracNonPrompt[iFile].SetPoint(iPt, pt, frac)
        gFracNonPrompt[iFile].SetPointError(iPt, 0.4, 0.4, totSyst, totSyst)

    if iFile == 0:
        hRatioFracPrompt.append(None)
        hRatioFracNonPrompt.append(None)
        gRatioFracNonPrompt.append(None)
        hRatioFracNonPromptRawYieldSys.append(None)
        hRatioFracNonPromptCutSetsSys.append(None)
        hRatioFracNonPromptMultWeightsSys.append(None)
        hRatioFracNonPromptPtWeightsSys.append(None)
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

        hRatioFracNonPromptRawYieldSys.append(
            hRatioFracNonPrompt[iFile].Clone(f'hRatioFracNonPromptRawYieldSys_{labels[iFile]}')
        )
        hRatioFracNonPromptCutSetsSys.append(
            hRatioFracNonPrompt[iFile].Clone(f'hRatioFracNonPromptCutSetsSys_{labels[iFile]}')
        )
        hRatioFracNonPromptMultWeightsSys.append(
            hRatioFracNonPrompt[iFile].Clone(f'hRatioFracNonPromptMultWeightsSys_{labels[iFile]}')
        )
        hRatioFracNonPromptPtWeightsSys.append(
            hRatioFracNonPrompt[iFile].Clone(f'hRatioFracNonPromptPtWeightsSys_{labels[iFile]}')
        )
        hRatioFracNonPromptRawYieldSys[iFile].SetDirectory(0)
        hRatioFracNonPromptCutSetsSys[iFile].SetDirectory(0)
        hRatioFracNonPromptMultWeightsSys[iFile].SetDirectory(0)
        hRatioFracNonPromptPtWeightsSys[iFile].SetDirectory(0)

        gRatioFracNonPrompt.append(TGraphAsymmErrors(1))
        gRatioFracNonPrompt[iFile].SetName(f'gRatioFracNonPrompt_{labels[iFile]}')
        SetObjectStyle(gRatioFracNonPrompt[iFile], color=GetROOTColor(colors[iFile]), fillstyle=0)

        for iPt in range(0, hRatioFracNonPrompt[iFile].GetNbinsX()):
            pt = hRatioFracNonPrompt[iFile].GetBinCenter(iPt+1)
            ratio = hRatioFracNonPrompt[iFile].GetBinContent(iPt+1)
            totSyst = np.sqrt(
                rawYieldSystRat[iFile][iPt]**2 + cutSetSystRat[iFile][iPt]**2 + \
                    multGenSystRat[iFile][iPt]**2 + ptGenSystRat[iFile][iPt]**2
            ) * ratio
            hRatioFracNonPromptRawYieldSys[iFile].SetBinContent(iPt+1, rawYieldSystRat[iFile][iPt])
            hRatioFracNonPromptCutSetsSys[iFile].SetBinContent(iPt+1, cutSetSystRat[iFile][iPt])
            hRatioFracNonPromptMultWeightsSys[iFile].SetBinContent(iPt+1, multGenSystRat[iFile][iPt])
            hRatioFracNonPromptPtWeightsSys[iFile].SetBinContent(iPt+1, ptGenSystRat[iFile][iPt])
            gRatioFracNonPrompt[iFile].SetPoint(iPt, pt, ratio)
            gRatioFracNonPrompt[iFile].SetPointError(iPt, 0.4, 0.4, totSyst, totSyst)

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
for histo, graph in zip(hFracNonPrompt, gFracNonPrompt):
    graph.Draw('2')
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
for histo, graph in zip(hRatioFracNonPrompt, gRatioFracNonPrompt):
    if histo:
        graph.Draw('2')
        histo.DrawCopy('same')

cFracs.Modified()
cFracs.Update()

cFracs.SaveAs(os.path.join(cfg['output']['directory'], cfg['output']['filewoext']) + '.pdf')

outFile = TFile(
    os.path.join(cfg['output']['directory'], cfg['output']['filewoext']) + '.root',
    'recreate'
)
cFracs.Write()
for iFile, _ in enumerate(hFracPrompt):
    hFracPrompt[iFile].Write()
    hFracNonPrompt[iFile].Write()
    gFracNonPrompt[iFile].Write()
    hFracNonPromptRawYieldSys[iFile].Write()
    hFracNonPromptCutSetsSys[iFile].Write()
    hFracNonPromptMultWeightsSys[iFile].Write()
    hFracNonPromptPtWeightsSys[iFile].Write()
    if iFile > 0:
        hRatioFracPrompt[iFile].Write()
        hRatioFracNonPrompt[iFile].Write()
        gRatioFracNonPrompt[iFile].Write()
        hRatioFracNonPromptRawYieldSys[iFile].Write()
        hRatioFracNonPromptCutSetsSys[iFile].Write()
        hRatioFracNonPromptMultWeightsSys[iFile].Write()
        hRatioFracNonPromptPtWeightsSys[iFile].Write()

outFile.Close()

input('Press enter to exit')
