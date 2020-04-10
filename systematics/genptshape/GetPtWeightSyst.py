'''
Script for the estimation of the MC pT shape systematic uncertainty
run: GetPtWeigthSyst.py cfgFileName.yml
'''

import sys
import argparse
import yaml
from ROOT import TCanvas, TFile, TLine, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kBlue, kOrange, kGreen, kRed, kAzure  # pylint: disable=import-error,no-name-in-module
from ROOT import kFullTriangleUp, kFullTriangleDown, kFullCross, kFullCircle, kFullSquare, kFullDiamond  # pylint: disable=import-error,no-name-in-module
sys.path.insert(0, '../..')
#pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgFileName', metavar='text', default='config_ptshape_syst.yml')
args = parser.parse_args()

SetGlobalStyle(padleftmargin=0.15, padbottommargin=0.14, titlesize=0.045, labelsize=0.04, titleoffsety=1.4)

modelNames = ['pythia', 'fonll', 'tamu', 'phsd', 'mc@shq', 'catania']
colors = {'pythia': kBlack, 'fonll': kBlue+1, 'tamu': kOrange+7, 'phsd': kGreen+2, 'catania': kRed+2, 'mc@shq':kAzure+4}
markers = {'pythia': kFullCircle, 'fonll': kFullSquare, 'tamu': kFullCross, 'phsd': kFullDiamond,
           'catania': kFullTriangleUp, 'mc@shq': kFullTriangleDown}
legNames = {'pythia': 'Pythia', 'fonll': 'FONLL', 'tamu': 'FONLL #times TAMU', 'phsd': 'FONLL #times PHSD',
            'catania': 'FONLL #times Catania', 'mc@shq': 'FONLL #times MC@sHQ'}
histoNames = {'pythia': 'hPtGen', 'fonll': 'hPtFONLL', 'tamu': 'hPtFONLLtimesTAMU', 'phsd': 'hPtFONLLtimesPHSD',
              'catania': 'hPtFONLLtimesCatania', 'mc@shq': 'hPtFONLLtimesGossiaux'}
histoWeightNames = {'fonll': 'hPtWeightsFONLL', 'tamu': 'hPtWeightsFONLLtimesTAMU', 'phsd': 'hPtWeightsFONLLtimesPHSD',
                    'catania': 'hPtWeightsFONLLtimesCatania', 'mc@shq': 'hPtWeightsFONLLtimesGossiaux'}

leg = TLegend(0.2, 0.2, 0.5, 0.4)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)

legEff = TLegend(0.6, 0.2, 0.9, 0.4)
legEff.SetFillStyle(0)
legEff.SetTextSize(0.04)

with open(args.cfgFileName, 'r') as ymlfitConfigFile:
    inputCfg = yaml.load(ymlfitConfigFile, yaml.FullLoader)

outDir = inputCfg['outputfile']['directory']
outSuffix = inputCfg['outputfile']['suffix']

inFilePtWeightsName = inputCfg['inputs']['ptweightsfile']
inFilePtWeights = TFile.Open(inFilePtWeightsName)
enableBweights = inputCfg['inputs']['enableBweights']

hPtDistrD, hPtWeightsD, hPtDistrB, hPtWeightsB, hEffPrompt, \
    hEffFD, hEffPromptRatio, hEffFDRatio = ({} for _ in range(8))

for model in modelNames:
    if model == 'pythia' or inputCfg['inputs']['shapes'][model]['enabled']:
        if model != 'pythia':
            hPtDistrD[model] = inFilePtWeights.Get(f'{histoNames[model]}Dcent')
            hPtWeightsD[model] = inFilePtWeights.Get(f'{histoWeightNames[model]}Dcent')
            SetObjectStyle(hPtWeightsD[model], linewidth=2, linecolor=colors[model],
                           markerstyle=markers[model], markercolor=colors[model])
        else:
            hPtDistrD[model] = inFilePtWeights.Get(f'{histoNames[model]}D')
        SetObjectStyle(hPtDistrD[model], linewidth=2, linecolor=colors[model],
                       markerstyle=markers[model], markercolor=colors[model])
        if enableBweights:
            if model != 'pythia':
                hPtDistrB[model] = inFilePtWeights.Get(f'{histoNames[model]}Bcent')
                hPtWeightsB[model] = inFilePtWeights.Get(f'{histoWeightNames[model]}Bcent')
                SetObjectStyle(hPtWeightsB[model], linewidth=2, linecolor=colors[model],
                               markerstyle=markers[model], markercolor=colors[model])
            else:
                hPtDistrB[model] = inFilePtWeights.Get(f'{histoNames[model]}B')
            SetObjectStyle(hPtDistrB[model], linewidth=2, linecolor=colors[model],
                           markerstyle=markers[model], markercolor=colors[model])
        inFileEff = TFile.Open(inputCfg['inputs']['shapes'][model]['efffile'])
        hEffPrompt[model] = inFileEff.Get('hEffPrompt')
        hEffFD[model] = inFileEff.Get('hEffFD')
        hEffPrompt[model].SetDirectory(0)
        hEffFD[model].SetDirectory(0)
        SetObjectStyle(hEffPrompt[model], linewidth=2, linecolor=colors[model],
                       markerstyle=markers[model], markercolor=colors[model])
        SetObjectStyle(hEffFD[model], linewidth=2, linecolor=colors[model],
                       markerstyle=markers[model], markercolor=colors[model])
        if model != 'pythia':
            hEffPromptRatio[model] = hEffPrompt[model].Clone(f'hEffPromptRatio{model}')
            hEffFDRatio[model] = hEffFD[model].Clone(f'hEffFDRatio{model}')
            hEffPromptRatio[model].SetDirectory(0)
            hEffFDRatio[model].SetDirectory(0)
            hEffPromptRatio[model].Divide(hEffPrompt[model], hEffPrompt['pythia'])
            hEffFDRatio[model].Divide(hEffFD[model], hEffFD['pythia'])
            for iPt in range(1, hEffPromptRatio[model].GetNbinsX()+1):
                hEffPromptRatio[model].SetBinError(iPt, 1.e-20)
                hEffFDRatio[model].SetBinError(iPt, 1.e-20)

        leg.AddEntry(hPtDistrD[model], legNames[model], 'l')
        legEff.AddEntry(hPtDistrD[model], legNames[model], 'lp')

ptMin = hEffPrompt['pythia'].GetBinLowEdge(1)
ptMax = hEffPrompt['pythia'].GetBinLowEdge(hEffPrompt['pythia'].GetNbinsX()) \
    + hEffPrompt['pythia'].GetBinWidth(hEffPrompt['pythia'].GetNbinsX())

lineatone = TLine(ptMin, 1., ptMax, 1.)
lineatone.SetLineWidth(2)
lineatone.SetLineStyle(9)
lineatone.SetLineColor(kBlack)

lineatoneWeight = TLine(0., 1., 50., 1.)
lineatoneWeight.SetLineWidth(2)
lineatoneWeight.SetLineStyle(9)
lineatoneWeight.SetLineColor(kBlack)

cShapesD = TCanvas('cShapesD', '', 1000, 500)
cShapesD.Divide(2, 1)
cShapesD.cd(1).DrawFrame(0., hPtDistrD['pythia'].GetMinimum()*0.1, ptMax, hPtDistrD['pythia'].GetMaximum()*10,
                         ';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T} (a. u.)')
cShapesD.cd(1).SetLogy()
for model in hPtDistrD:
    print(model)
    hPtDistrD[model].DrawCopy('chistsame')
leg.Draw()
cShapesD.cd(2).DrawFrame(0., 0.8, ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); #it{p}_{T} weights')
lineatoneWeight.Draw('same')
for model in hPtWeightsD:
    hPtWeightsD[model].DrawCopy('chistsame')

cShapesD.SaveAs(f'{outDir}/PtWeightsD{outSuffix}.pdf')

if enableBweights:
    cShapesB = TCanvas('cShapesB', '', 1000, 500)
    cShapesB.Divide(2, 1)
    cShapesB.cd(1).DrawFrame(0, hPtDistrB['pythia'].GetMinimum()*0.1, ptMax*2, hPtDistrD['pythia'].GetMaximum()*10,
                             ';#it{p}_{T}^{B} (GeV/#it{c}); d#it{N}/d#it{p}_{T}^{B} (a. u.)')
    cShapesB.cd(1).SetLogy()
    for model in hPtDistrB:
        hPtDistrB[model].DrawCopy('chistsame')
    leg.Draw()
    cShapesB.cd(2).DrawFrame(0, 0.8, ptMax*2, 1.2, ';#it{p}_{T}^{B} (GeV/#it{c}); #it{p}_{T}^{B} weights')
    lineatoneWeight.Draw('same')
    for model in hPtWeightsB:
        hPtWeightsB[model].DrawCopy('chistsame')

    cShapesB.SaveAs(f'{outDir}/PtWeightsB{outSuffix}.pdf')


cEffD = TCanvas('cEffD', '', 1000, 500)
cEffD.Divide(2, 1)
cEffD.cd(1).DrawFrame(ptMin, 1.e-4, ptMax, 1., ';#it{p}_{T} (GeV/#it{c}); prompt efficiency')
cEffD.cd(1).SetLogy()
for model in hEffPrompt:
    hEffPrompt[model].DrawCopy('same')
legEff.Draw()
cEffD.cd(2).DrawFrame(ptMin, 0.2, ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); prompt efficiency ratio')
lineatone.Draw('same')
for model in hEffPromptRatio:
    hEffPromptRatio[model].DrawCopy('same')
cEffD.SaveAs(f'{outDir}/SystPtWeightsPrompt{outSuffix}.pdf')

cEffB = TCanvas('cEffB', '', 1000, 500)
cEffB.Divide(2, 1)
cEffB.cd(1).DrawFrame(ptMin, 1.e-4, ptMax, 1., ';#it{p}_{T} (GeV/#it{c}); FD efficiency')
cEffB.cd(1).SetLogy()
for model in hEffFD:
    hEffFD[model].DrawCopy('same')
legEff.Draw()
cEffB.cd(2).DrawFrame(ptMin, 0.2, ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); FD efficiency ratio')
lineatone.Draw('same')
for model in hEffFDRatio:
    hEffFDRatio[model].DrawCopy('same')
cEffB.SaveAs(f'{outDir}/SystPtWeightsFD{outSuffix}.pdf')

input('Press enter to exit')
