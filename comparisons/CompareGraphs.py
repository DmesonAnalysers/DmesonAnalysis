'''
Script for the comparison of ROOT TH1s or TGraphs
run: python CompareGraphs.py cfgFileName.yml
'''

import sys
from os.path import join
import argparse
import numpy as np
import yaml
from ROOT import TCanvas, TFile, TLegend, TLine # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
#pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor, GetROOTMarker
from utils.AnalysisUtils import ComputeRatioDiffBins

# load inputs
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgFileName', metavar='text', default='config_comparison.yml')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inDirName = inputCfg['inputs']['dirname']
inFileNames = inputCfg['inputs']['filenames']
objNames = inputCfg['inputs']['objectnames']

outFileName = inputCfg['output']['filename']
outExtensions = inputCfg['output']['extensions']

objTypes = inputCfg['options']['ROOTobject']
scales = inputCfg['options']['scale']
colors = inputCfg['options']['colors']
markers = inputCfg['options']['markers']
drawOptions = inputCfg['options']['drawopt']

doRatio = inputCfg['options']['ratio']['enable']
drawRatioUnc = inputCfg['options']['ratio']['uncertainties']['enable']
ratioUncCorr = inputCfg['options']['ratio']['uncertainties']['corr']
displayRMS = inputCfg['options']['ratio']['displayRMS']

doCompareUnc = inputCfg['options']['errcomp']['enable']
compareRelUnc = inputCfg['options']['errcomp']['relative']

wCanv = inputCfg['options']['canvas']['width']
hCanv = inputCfg['options']['canvas']['heigth']
xLimits = inputCfg['options']['canvas']['xlimits']
yLimits = inputCfg['options']['canvas']['ylimits']
yLimitsRatio = inputCfg['options']['canvas']['ylimitsratio']
yLimitsUnc = inputCfg['options']['canvas']['ylimitserr']
xTitle = inputCfg['options']['canvas']['xaxistitle']
yTitle = inputCfg['options']['canvas']['yaxistitle']
logX = inputCfg['options']['canvas']['logx']
logY = inputCfg['options']['canvas']['logy']
ratioLogX = inputCfg['options']['canvas']['ratio']['logx']
ratioLogY = inputCfg['options']['canvas']['ratio']['logy']
uncCompLogX = inputCfg['options']['canvas']['errcomp']['logx']
uncCompLogY = inputCfg['options']['canvas']['errcomp']['logy']

xLegLimits = inputCfg['options']['legend']['xlimits']
yLegLimits = inputCfg['options']['legend']['ylimits']
legNames = inputCfg['options']['legend']['titles']
legOpt = inputCfg['options']['legend']['options']
legTextSize = inputCfg['options']['legend']['textsize']

# set global style
SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14, titleoffsety=1.5)

leg = TLegend(xLegLimits[0], yLegLimits[0], xLegLimits[1], yLegLimits[1])
leg.SetFillStyle(0)
leg.SetTextSize(legTextSize)

hToCompare, hRatioToCompare, hUncToCompare = [], [], []
for iFile, (inFileName, objName, objType, scale, color, marker) in \
    enumerate(zip(inFileNames, objNames, objTypes, scales, colors, markers)):
    if inDirName:
        inFileName = join(inDirName, inFileName)
    inFile = TFile.Open(inFileName)
    hToCompare.append(inFile.Get(objName))
    if 'TH' in objType:
        hToCompare[iFile].SetName(f'h{iFile}')
        hToCompare[iFile].SetStats(0)
    else:
        hToCompare[iFile].SetName(f'g{iFile}')
    SetObjectStyle(hToCompare[iFile], color=GetROOTColor(color), markerstyle=GetROOTMarker(marker), fillstyle=0)
    if 'TH' in objType:
        hToCompare[iFile].SetDirectory(0)
        hToCompare[iFile].SetStats(0)
        hToCompare[iFile].Scale(scale)
    #TODO: else: scale graph --> add util function in AnalysisUtils
    if doRatio:
        if 'TH' in objType:
            if drawRatioUnc:
                if ratioUncCorr:
                    hRatioToCompare.append(ComputeRatioDiffBins(hToCompare[iFile], hToCompare[0], 'B'))
                else:
                    hRatioToCompare.append(ComputeRatioDiffBins(hToCompare[iFile], hToCompare[0]))
            else:
                hRatioToCompare.append(ComputeRatioDiffBins(hToCompare[iFile], hToCompare[0]))
                print(hRatioToCompare[0].GetBinContent(1))
                for iBin in range(1, hRatioToCompare[iFile].GetNbinsX()+1):
                    hRatioToCompare[iFile].SetBinError(iBin, 1.e-20)
            hRatioToCompare[iFile].SetDirectory(0)
        #TODO: add util function in AnalysisUtils to manage ratios between graphs or graph and histo
        SetObjectStyle(hRatioToCompare[iFile], color=GetROOTColor(color), markerstyle=GetROOTMarker(marker),
                       fillstyle=0)
    if doCompareUnc:
        if 'TH' in objType:
            hUncToCompare.append(hToCompare[iFile].Clone(f'hUncToCompare{iFile}'))
            for iBin in range(1, hUncToCompare[iFile].GetNbinsX()+1):
                unc = hUncToCompare[iFile].GetBinError(iBin)
                cent = hUncToCompare[iFile].GetBinContent(iBin)
                if compareRelUnc:
                    hUncToCompare[iFile].SetBinContent(iBin, unc/cent)
                else:
                    hUncToCompare[iFile].SetBinContent(iBin, unc)
                hUncToCompare[iFile].SetBinError(iBin, 1.e-20)
            hUncToCompare[iFile].SetDirectory(0)
            SetObjectStyle(hUncToCompare[iFile], color=GetROOTColor(color), markerstyle=GetROOTMarker(marker),
                           fillstyle=0)
        else:
            #TODO: add uncertainty comparison for TGraphs
            print('WARNING: uncertainty comparison for TGraphs not implemented. Switching off')
            doCompareUnc = False

    leg.AddEntry(hToCompare[iFile], legNames[iFile], legOpt[iFile])

ratios, RMS, shift = [], [], []
if doRatio and displayRMS:
    for iBin in range(hRatioToCompare[1].GetNbinsX()):
        ratios.append([])
        for iFile, _ in enumerate(inFileNames):
            if iFile == 0:
                continue
            ratios[iBin].append(hRatioToCompare[iFile].GetBinContent(iBin+1))
        aRatios = np.array(ratios[iBin])
        RMS.append(np.std(aRatios))
        shift.append(np.mean(aRatios))
print('\033[92mRMS values:', np.around(RMS, decimals=3), '\033[0m')
print('\033[92mshift values:', np.around(shift, decimals=3) - 1., '\033[0m')

cOut = TCanvas('cOutput', '', wCanv, hCanv)

if doRatio or doCompareUnc:
    if doRatio and doCompareUnc:
        cOut.Divide(3, 1)
        ratioPad = 2
        uncPad = 3
    else:
        cOut.Divide(2, 1)
        if doRatio:
            ratioPad = 2
        else:
            uncPad = 2

    hFrame = cOut.cd(1).DrawFrame(xLimits[0], yLimits[0], xLimits[1], yLimits[1], f';{xTitle};{yTitle}')
    if logX:
        cOut.cd(1).SetLogx()
    if logY:
        cOut.cd(1).SetLogy()
else:
    hFrame = cOut.cd().DrawFrame(xLimits[0], yLimits[0], xLimits[1], yLimits[1], f';{xTitle};{yTitle}')
    if logX:
        cOut.cd().SetLogx()
    if logY:
        cOut.cd().SetLogy()
hFrame.GetYaxis().SetDecimals()

for histo, objType, drawOpt in zip(hToCompare, objTypes, drawOptions):
    if 'TH' in objType:
        histo.DrawCopy(f'{drawOpt}same')
    else:
        histo.Draw(drawOpt)
leg.Draw()

if doRatio:
    hFrameRatio = cOut.cd(ratioPad).DrawFrame(xLimits[0], yLimitsRatio[0], xLimits[1], yLimitsRatio[1],
                                              f';{xTitle};Ratio')
    hFrameRatio.GetYaxis().SetDecimals()
    if ratioLogX:
        cOut.cd(ratioPad).SetLogx()
    if ratioLogY:
        cOut.cd(ratioPad).SetLogy()
    lineAtOne = TLine(xLimits[0], 1., xLimits[1], 1.)
    lineAtOne.SetLineColor(GetROOTColor('kBlack'))
    lineAtOne.SetLineWidth(2)
    lineAtOne.SetLineStyle(9)
    lineAtOne.Draw()
    for iHisto, (histo, objType, drawOpt) in enumerate(zip(hRatioToCompare, objTypes, drawOptions)):
        if iHisto > 0:
            if 'TH' in objType:
                histo.DrawCopy(f'{drawOpt}same')
            else:
                histo.Draw(drawOpt)

if doCompareUnc:
    if compareRelUnc:
        uncTitle = 'Relative uncertainties'
    else:
        uncTitle = 'Absolute uncertainties'
    hFrameUnc = cOut.cd(uncPad).DrawFrame(xLimits[0], yLimitsUnc[0], xLimits[1], yLimitsUnc[1],
                                          f';{xTitle};{uncTitle}')
    hFrameUnc.GetYaxis().SetDecimals()
    if uncCompLogX:
        cOut.cd(uncPad).SetLogx()
    if uncCompLogY:
        cOut.cd(uncPad).SetLogy()
    for iHisto, (histo, drawOpt) in enumerate(zip(hUncToCompare, drawOptions)):
        histo.DrawCopy(f'{drawOpt}same')

for ext in outExtensions:
    if 'root' in ext:
        outFile = TFile(f'{outFileName}.root', 'recreate')
        cOut.Write()
        for histo in hToCompare:
            histo.Write()
        if doRatio:
            for histo in hRatioToCompare:
                histo.Write()
        if doCompareUnc:
            for histo in hUncToCompare:
                histo.Write()
        outFile.Close()
    else:
        cOut.SaveAs(f'{outFileName}.{ext}')

input()
