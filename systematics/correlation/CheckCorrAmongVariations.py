'''
Script for the evaluation of the correlation degree among variations of two variables
run: python CheckCorrAmongVariations.py cfgFileName.yml
'''

import sys
from os.path import join
import argparse
import yaml
from ROOT import TCanvas, TFile, TGraph, TLatex # pylint: disable=import-error,no-name-in-module
from ROOT import kAzure, kFullCircle # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error,no-name-in-module

# load inputs
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgFileName', metavar='text', default='config_comparison.yml')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

varName = {'var1': inputCfg['inputs']['var1']['name'], 'var2': inputCfg['inputs']['var2']['name']}
objNameVar = {'var1': inputCfg['inputs']['var1']['histoname'], 'var2': inputCfg['inputs']['var2']['histoname']}
inFileRefNameVar = {'var1': inputCfg['inputs']['var1']['reffilename'],
                    'var2': inputCfg['inputs']['var2']['reffilename']}
inDirNameVar = {'var1': inputCfg['inputs']['var1']['variations']['dirname'],
                'var2': inputCfg['inputs']['var2']['variations']['dirname']}
inFileNamesVar = {'var1': inputCfg['inputs']['var1']['variations']['filenames'],
                  'var2': inputCfg['inputs']['var2']['variations']['filenames']}

doRelVar = inputCfg['outputs']['variations']['relative']
doAbsVar = inputCfg['outputs']['variations']['relative']
outFileName = inputCfg['outputs']['filename']

if len(inFileNamesVar['var1']) != len(inFileNamesVar['var2']):
    print(f"ERROR: different number of files for {varName['var1']} and {varName['var2']}! Exit")
    sys.exit()

# set global style
SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14, titleoffsety=1.7)

hRefVar = {}
for var in inFileRefNameVar:
    inFileRefVar = TFile.Open(inFileRefNameVar[var])
    hRefVar[var] = inFileRefVar.Get(objNameVar[var])
    hRefVar[var].SetDirectory(0)

if hRefVar['var1'].GetNbinsX() != hRefVar['var2'].GetNbinsX():
    print(f"ERROR: different number of bins in {inFileRefNameVar['var1']} and {inFileRefNameVar['var2']}! Exit")
    sys.exit()

hVar = {}
for var in inFileRefNameVar:
    hVar[var] = []
    for iFile, inFileNameVar in enumerate(inFileNamesVar[var]):

        if inDirNameVar[var]:
            inFileNameVar = join(inDirNameVar[var], inFileNameVar)
        inFileVar = TFile.Open(inFileNameVar)
        hVar[var].append(inFileVar.Get(objNameVar[var]))
        hVar[var][-1].SetDirectory(0)

        if hVar[var][-1].GetNbinsX() != hRefVar[var].GetNbinsX():
            print(f"ERROR: different number of bins in {inFileNameVar} and {inFileRefNameVar}! Exit")
            sys.exit()

iPoint = 0
gRelVar, gAbsVar = TGraph(0), TGraph(0)
for iVar, (histoVar1, histoVar2) in enumerate(zip(hVar['var1'], hVar['var2'])):
    for iBin in range(1, hRefVar['var1'].GetNbinsX()+1):
        if doRelVar:
            gRelVar.SetPoint(iPoint, histoVar1.GetBinContent(iBin) / hRefVar['var1'].GetBinContent(iBin) - 1,
                             histoVar2.GetBinContent(iBin) / hRefVar['var2'].GetBinContent(iBin) - 1)
        if doAbsVar:
            gAbsVar.SetPoint(iPoint, histoVar1.GetBinContent(iBin) - hRefVar['var1'].GetBinContent(iBin),
                             histoVar2.GetBinContent(iBin) - hRefVar['var2'].GetBinContent(iBin))
        iPoint += 1

latRho = TLatex()
latRho.SetNDC()
latRho.SetTextSize(0.05)
latRho.SetTextFont(42)

if doRelVar:
    cRelVar = TCanvas('cRelVar', '', 800, 800)
    SetObjectStyle(gRelVar, color=kAzure+4, markerstyle=kFullCircle)
    gRelVar.GetYaxis().SetDecimals()
    gRelVar.SetTitle(f";relative variation of {varName['var1']};relative variation of {varName['var2']}")
    gRelVar.Draw('AP')
    latRho.DrawLatex(0.22, 0.89, f'#rho = {gRelVar.GetCorrelationFactor():0.3f}')
    cRelVar.SaveAs(f'{outFileName}_RelativeVariation.pdf')
if doAbsVar:
    cAbsVar = TCanvas('cAbsVar', '', 800, 800)
    SetObjectStyle(gAbsVar, color=kAzure+4, markerstyle=kFullCircle)
    gAbsVar.GetYaxis().SetDecimals()
    gAbsVar.SetTitle(f";absolute variation of {varName['var1']};absolute variation of {varName['var2']}")
    gAbsVar.Draw('AP')
    latRho.DrawLatex(0.22, 0.89, f'#rho = {gAbsVar.GetCorrelationFactor():0.3f}')
    gAbsVar.SaveAs(f'{outFileName}_AbsoluteVariation.pdf')

input('Press enter to exit')
