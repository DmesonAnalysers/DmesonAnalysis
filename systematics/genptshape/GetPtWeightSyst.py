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

SetGlobalStyle(padleftmargin=0.14, titlesize=0.045, labelsize=0.04)

modelNames = ['pythia', 'fonll', 'tamu', 'phsd', 'mc@shq', 'catania']
colors = {'pythia': kBlack, 'fonll': kBlue+1, 'tamu': kOrange+7, 'phsd': kGreen+2, 'catania': kRed+2, 'mc@shq':kAzure+4}
markers = {'pythia': kFullCircle, 'fonll': kFullSquare, 'tamu': kFullCross, 'phsd': kFullDiamond,
           'catania': kFullTriangleUp, 'mc@shq': kFullTriangleDown}
legNames = {'pythia': 'Pythia', 'fonll': 'FONLL', 'tamu': 'FONLL #times TAMU', 'phsd': 'FONLL #times PHSD',
            'catania': 'FONLL #times Catania', 'mc@shq': 'FONLL #times MC@sHQ'}

leg = TLegend(0.3, 0.7, 0.7, 0.9)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)

legEff = TLegend(0.3, 0.2, 0.7, 0.4)
legEff.SetFillStyle(0)
legEff.SetTextSize(0.04)

with open(args.fitConfigFileName, 'r') as ymlfitConfigFile:
    inputCfg = yaml.load(ymlfitConfigFile, yaml.FullLoader)

outDir = inputCfg['output']['directory']
outSuffix = inputCfg['output']['suffix']

inFilePtWeightsName = inputCfg['inputs']['ptweightsfile']
inFilePtWeights = TFile.Open(inFilePtWeightsName)
enableBweights = inputCfg['inputs']['enableBweights']

hPtDistrD, hPtWeightsD, hPtDistrB, hPtWeightsB, hEffPrompt, \
    hEffFD, hEffPromptRatio, hEffFDRatio = ({} for _ in range(8))

for model in modelNames:
    if model == 'pythia' or inputCfg['inputs'][model]['enable']:
        hPtDistrD[model] = inFilePtWeights.Get(f'')
        SetObjectStyle(hPtDistrD[model], linewidth=2, linecolor=colors[model],
                       markerstyle=markers[model], markercolor=colors[model])
        if model != 'pythia':
            hPtWeightsD[model] = inFilePtWeights.Get(f'')
            SetObjectStyle(hPtWeightsD[model], linewidth=2, linecolor=colors[model],
                           markerstyle=markers[model], markercolor=colors[model])
        if enableBweights:
            hPtDistrB[model] = inFilePtWeights.Get(f'')
            SetObjectStyle(hPtDistrB[model], linewidth=2, linecolor=colors[model],
                           markerstyle=markers[model], markercolor=colors[model])
            if model != 'pythia':
                hPtWeightsB[model] = inFilePtWeights.Get(f'')
                SetObjectStyle(hPtWeightsB[model], linewidth=2, linecolor=colors[model],
                               markerstyle=markers[model], markercolor=colors[model])
        inFileEff = TFile.Open(inputCfg['inputs'][model]['efffile'])
        hEffPrompt[model] = inFileEff.Get('hEffPrompt')
        hEffFD[model] = inFileEff.Get('hEffFD')
        SetObjectStyle(hEffPrompt[model], linewidth=2, linecolor=colors[model],
                       markerstyle=markers[model], markercolor=colors[model])
        SetObjectStyle(hEffFD[model], linewidth=2, linecolor=colors[model],
                       markerstyle=markers[model], markercolor=colors[model])

        leg.AddEntry(hPtDistrD[model], legNames[model], 'l')
        legEff.AddEntry(hPtDistrD[model], legNames[model], 'lp')

ptMin = hPtDistrD['pythia'].GetBinLowEdge(1)
ptMax = hPtDistrD['pythia'].GetBinLowEdge(hPtDistrD['pythia'].GetNbinsX()) \
    + hPtDistrD['pythia'].GetBinWidth(hPtDistrD['pythia'].GetNbinsX())

lineatone = TLine(ptMin, 1., ptMax, 1.)
lineatone.SetLineWidth(2)
lineatone.SetLineStyle(9)
lineatone.SetLineColor(kBlack)

cShapesD = TCanvas('cShapes', '', 1000, 500)
cShapesD.Divide(2, 1)
cShapesD.cd(1).DrawFrame(ptMin, 1.e-8, ptMax, 1., ';#it{p}_{T} GeV/#it{c}; d#it{N}/d#it{p}_{T} (a. u.)')
cShapesD.cd(1).SetLogy()
for histo in hPtDistrD:
    histo.DrawCopy('chistsame')
leg.Draw()
cShapesD.cd(2).DrawFrame(ptMin, 1.e-3, ptMax, 1.e+1, ';#it{p}_{T} GeV/#it{c}; #it{p}_{T} weights')
cShapesD.cd(2).SetLogy()
lineatone.Draw('same')
for histo in hPtWeightsD:
    histo.DrawCopy('chistsame')

cShapesD.SaveAs(f'{outDir}/PtWeightsD{outSuffix}.pdf')

if enableBweights:
    cShapesB = TCanvas('cShapes', '', 1000, 500)
    cShapesB.Divide(2, 1)
    cShapesB.cd(1).DrawFrame(ptMin, 1.e-8, ptMax, 1., ';#it{p}_{T} GeV/#it{c}; d#it{N}/d#it{p}_{T}^{B} (a. u.)')
    cShapesB.cd(1).SetLogy()
    for histo in hPtDistrD:
        histo.DrawCopy('chistsame')
    leg.Draw()
    cShapesB.cd(2).DrawFrame(ptMin, 1.e-3, ptMax, 1.e+1, ';#it{p}_{T} GeV/#it{c}; #it{p}_{T}^{B} weights')
    cShapesB.cd(2).SetLogy()
    lineatone.Draw('same')
    for histo in hPtWeightsD:
        histo.DrawCopy('chistsame')

    cShapesB.SaveAs(f'{outDir}/PtWeightsB{outSuffix}.pdf')


cEffD = TCanvas('cEffD', '', 1000, 500)
cEffD.Divide(2, 1)
cEffD.cd(1).DrawFrame(ptMin, 1.e-4, ptMax, 1., ';#it{p}_{T} GeV/#it{c}; prompt efficiency')
cEffD.cd(1).SetLogy()
for histo in hEffPrompt:
    histo.DrawCopy('same')
legEff.Draw()
cEffD.cd(2).DrawFrame(ptMin, 0.8, ptMax, 1.2, ';#it{p}_{T} GeV/#it{c}; prompt efficiency ratio')
for histo in hEffPromptRatio:
    histo.DrawCopy('same')
cEffD.SaveAs(f'{outDir}/SystPtWeightsPrompt{outSuffix}.pdf')

cEffB = TCanvas('cEffB', '', 1000, 500)
cEffB.Divide(2, 1)
cEffB.cd(1).DrawFrame(ptMin, 1.e-4, ptMax, 1., ';#it{p}_{T} GeV/#it{c}; FD efficiency')
cEffB.cd(1).SetLogy()
for histo in hEffFD:
    histo.DrawCopy('same')
legEff.Draw()
cEffB.cd(2).DrawFrame(ptMin, 0.8, ptMax, 1.2, ';#it{p}_{T} GeV/#it{c}; FD efficiency ratio')
for histo in hEffFDRatio:
    histo.DrawCopy('same')
cEffB.SaveAs(f'{outDir}/SystPtWeightsFD{outSuffix}.pdf')

input('Press enter to exit')
