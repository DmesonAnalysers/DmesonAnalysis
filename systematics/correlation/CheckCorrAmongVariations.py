'''
Script for the evaluation of the correlation degree among variations of two variables
run: python CheckCorrAmongVariations.py cfgFileName.yml
'''

import sys
from os.path import join
import argparse
import numpy as np
from sklearn.cluster import dbscan
import yaml
from ROOT import TCanvas, TFile, TGraph, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kAzure, kRed, kWhite, kFullCircle # pylint: disable=import-error,no-name-in-module
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

doRelVar = inputCfg['outputs']['variations']['relative']['enable']
if doRelVar:
    removeOutliersRel = inputCfg['outputs']['variations']['relative']['outliers']['remove']
doAbsVar = inputCfg['outputs']['variations']['absolute']['enable']
if doAbsVar:
    removeOutliersAbs = inputCfg['outputs']['variations']['relative']['outliers']['remove']
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

relVar, absVar = [], []
for iVar, (histoVar1, histoVar2) in enumerate(zip(hVar['var1'], hVar['var2'])):
    for iBin in range(1, hRefVar['var1'].GetNbinsX()+1):
        if doRelVar:
            relVar.append([histoVar1.GetBinContent(iBin) / hRefVar['var1'].GetBinContent(iBin) - 1,
                           histoVar2.GetBinContent(iBin) / hRefVar['var2'].GetBinContent(iBin) - 1])
        if doAbsVar:
            absVar.append([histoVar1.GetBinContent(iBin) - hRefVar['var1'].GetBinContent(iBin),
                           histoVar2.GetBinContent(iBin) - hRefVar['var2'].GetBinContent(iBin)])
relVar = np.array(relVar)
absVar = np.array(absVar)

if doRelVar and removeOutliersRel:
    radius = inputCfg['outputs']['variations']['relative']['outliers']['DBSCAN']['radius']
    minsample = inputCfg['outputs']['variations']['relative']['outliers']['DBSCAN']['minsample']
    _, clLabelsRel = dbscan(relVar, eps=radius, min_samples=minsample, metric='euclidean')
    gRelVarOutliers = TGraph(0)
if doAbsVar and removeOutliersAbs:
    radius = inputCfg['outputs']['variations']['absolute']['outliers']['DBSCAN']['radius']
    minsample = inputCfg['outputs']['variations']['absolute']['outliers']['DBSCAN']['minsample']
    _, clLabelsAbs = dbscan(relVar, eps=radius, min_samples=minsample, metric='euclidean')
    gAbsVarOutliers = TGraph(0)

iPointOutlier, iPointCore = 0, 0
gRelVar, gAbsVar = TGraph(0), TGraph(0)
if doRelVar:
    for iPoint, point in enumerate(relVar):
        if not removeOutliersRel or (removeOutliersRel and clLabelsRel[iPoint] != -1):
            gRelVar.SetPoint(iPointCore, point[0], point[1])
            iPointCore += 1
        else:
            gRelVarOutliers.SetPoint(iPointOutlier, point[0], point[1])
            iPointOutlier += 1

iPointOutlier, iPointCore = 0, 0
if doAbsVar:
    for iPoint, point in enumerate(absVar):
        if not removeOutliersRel or (removeOutliersRel and clLabelsAbs[iPoint] != -1):
            gAbsVar.SetPoint(iPointCore, point[0], point[1])
            iPointCore += 1
        else:
            gAbsVarOutliers.SetPoint(iPointOutlier, point[0], point[1])
            iPointOutlier += 1

if doRelVar:
    legRel = TLegend(0.18, 0.815, 0.45, 0.965)
    legRel.SetTextSize(0.045)
    legRel.SetBorderSize(1)
    legRel.SetHeader(f'#rho = {gRelVar.GetCorrelationFactor():0.3f}')
    if removeOutliersRel:
        legRel.AddEntry(gRelVar, 'considered', 'p')
        legRel.AddEntry(gRelVarOutliers, 'outliers', 'p')

if doAbsVar:
    legAbs = TLegend(0.18, 0.815, 0.45, 0.965)
    legAbs.SetTextSize(0.045)
    legAbs.SetBorderSize(1)
    legAbs.SetFillColor(kWhite)
    legAbs.SetHeader(f'#rho = {gAbsVar.GetCorrelationFactor():0.3f}')
    if removeOutliersRel:
        legAbs.AddEntry(gRelVar, 'considered', 'p')
        legAbs.AddEntry(gRelVarOutliers, 'outliers', 'p')

if doRelVar:
    cRelVar = TCanvas('cRelVar', '', 800, 800)
    minVals = np.min(relVar, axis=0)
    maxVals = np.max(relVar, axis=0)
    hFrameRel = cRelVar.DrawFrame(minVals[0]-0.2*abs(minVals[0]), minVals[1]-0.2*abs(minVals[1]),
                                  maxVals[0]+0.2*abs(maxVals[0]), maxVals[1]+0.2*abs(maxVals[1]),
                                  f";relative variation of {varName['var1']};relative variation of {varName['var2']}")
    hFrameRel.GetYaxis().SetDecimals()
    hFrameRel.GetYaxis().SetNdivisions(505)
    hFrameRel.GetXaxis().SetNdivisions(505)
    SetObjectStyle(gRelVar, color=kAzure+4, markerstyle=kFullCircle)
    gRelVar.Draw('P')
    if removeOutliersRel:
        SetObjectStyle(gRelVarOutliers, color=kRed+1, markerstyle=kFullCircle)
        gRelVarOutliers.Draw('P')
        legRel.Draw()
    cRelVar.SaveAs(f'{outFileName}_RelativeVariation.pdf')
if doAbsVar:
    cAbsVar = TCanvas('cAbsVar', '', 800, 800)
    minVals = np.min(absVar, axis=0)
    maxVals = np.max(absVar, axis=0)
    hFrameAbs = cAbsVar.DrawFrame(minVals[0]-0.2*abs(minVals[0]), minVals[1]-0.2*abs(minVals[1]),
                                  maxVals[0]+0.2*abs(maxVals[0]), maxVals[1]+0.2*abs(maxVals[1]),
                                  f";absolute variation of {varName['var1']};absolute variation of {varName['var2']}")
    hFrameAbs.GetYaxis().SetDecimals()
    hFrameAbs.GetYaxis().SetNdivisions(505)
    hFrameAbs.GetXaxis().SetNdivisions(505)
    SetObjectStyle(gAbsVar, color=kAzure+4, markerstyle=kFullCircle)
    gAbsVar.Draw('P')
    if removeOutliersAbs:
        SetObjectStyle(gAbsVarOutliers, color=kRed+1, markerstyle=kFullCircle)
        gAbsVarOutliers.Draw('P')
        legAbs.Draw()
    cAbsVar.SaveAs(f'{outFileName}_AbsoluteVariation.pdf')

input('Press enter to exit')
