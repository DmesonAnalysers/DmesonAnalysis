'''
Script for the comparison of ROOT TH1s or TGraphs
run: python CompareGraphs.py cfgFileName.yml
'''

import sys
from os.path import join
import argparse
import numpy as np
import yaml
from ROOT import TCanvas, TFile, TLegend, TLine, gROOT # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor, GetROOTMarker #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ComputeRatioDiffBins, ScaleGraph, ComputeRatioGraph #pylint: disable=wrong-import-position,import-error
from utils.DfUtils import GetObjectFromFile #pylint: disable=wrong-import-position,import-error

# load inputs
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgFileName', metavar='text', default='config_comparison.yml')
parser.add_argument('-b', default=False, action='store_true')

args = parser.parse_args()
gROOT.SetBatch(args.b)

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inDirName = inputCfg['inputs']['dirname']
inFileNames = inputCfg['inputs']['filenames']
objNames = inputCfg['inputs']['objectnames']
nObj = len(objNames)

outFileName = inputCfg['output']['filename']
outExtensions = inputCfg['output']['extensions']

objTypes = inputCfg['options']['ROOTobject']
if not isinstance(objTypes, list):
    objTypes = [objTypes] * nObj

scales = inputCfg['options']['scale']
if not isinstance(scales, list):
    scales = [scales] * nObj

lambdaParams = inputCfg['options']['lambdaParams']
if not isinstance(lambdaParams, list):
    lambdaParams = [lambdaParams] * nObj
normalizes = inputCfg['options']['normalize']
if not isinstance(normalizes, list):
    normalizes = [normalizes] * nObj
colors = inputCfg['options']['colors']
markers = inputCfg['options']['markers']
markersize = inputCfg['options']['markersize']
linewidth = inputCfg['options']['linewidth']
fillstyles = inputCfg['options']['fillstyle']
if not isinstance(fillstyles, list):
    fillstyles = [fillstyles] * nObj
fillalphas = inputCfg['options']['fillalpha']
if not isinstance(fillalphas, list):
    fillalphas = [fillalphas] * nObj
drawOptions = inputCfg['options']['drawopt']
if not isinstance(drawOptions, list):
    drawOptions = [drawOptions] * nObj

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

avoidLeg = inputCfg['options']['legend']['avoid']
xLegLimits = inputCfg['options']['legend']['xlimits']
yLegLimits = inputCfg['options']['legend']['ylimits']
legHeader = inputCfg['options']['legend']['header']
legNames = inputCfg['options']['legend']['titles']
legOpt = inputCfg['options']['legend']['options']
legTextSize = inputCfg['options']['legend']['textsize']
ncolumns = inputCfg['options']['legend']['ncolumns']

# set global style
SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14, titleoffsety=1.5)

leg = TLegend(xLegLimits[0], yLegLimits[0], xLegLimits[1], yLegLimits[1])
leg.SetFillStyle(0)
leg.SetTextSize(legTextSize)
leg.SetNColumns(ncolumns)
if legHeader is not None:
    leg.SetHeader(legHeader, 'C')

hToCompare, hRatioToCompare, hUncToCompare = [], [], []
for iFile, (inFileName, objName, objType, scale, lambdaParam, normalize, color, marker, fillstyle, fillalpha) in \
    enumerate(
        zip(inFileNames, objNames, objTypes, scales, lambdaParams, normalizes, colors, markers, fillstyles, fillalphas)
    ):
    if inDirName:
        inFileName = join(inDirName, inFileName)
    inFile = TFile.Open(inFileName)
    if inFile == None:
        print(f"ERROR: cannot open {inFileName}. Check your config. Exit!")
        sys.exit()
    objToCompare = GetObjectFromFile(inFile, objName)
    if objToCompare == None:
        print(f"ERROR: couldn't load the histogram \'{objName}\' in \'{inFileName}\'. Check your config. Exit! ")
        sys.exit()
    hToCompare.append(objToCompare)
    if 'TH' in objType:
        hToCompare[iFile].SetName(f'h{iFile}')
        hToCompare[iFile].SetStats(0)
    else:
        hToCompare[iFile].SetName(f'g{iFile}')
    SetObjectStyle(hToCompare[iFile],
                   color=GetROOTColor(color),
                   markerstyle=GetROOTMarker(marker),
                   markersize=markersize,
                   linewidth=linewidth,
                   fillstyle=fillstyle,
                   fillalpha=fillalpha)
    if 'TH' in objType:
        hToCompare[iFile].SetDirectory(0)
        hToCompare[iFile].SetStats(0)
        if normalize:
            if scale != 1.:
                print('WARNING: you are both scaling and normalizing the histogram, check if it makes sense!')
            hToCompare[iFile].Scale(1. / hToCompare[iFile].Integral())
        hToCompare[iFile].Scale(scale)

        # apply scaling using lambda parameters
        for iBin in range(hToCompare[iFile].GetNbinsX()):
            bc = hToCompare[iFile].GetBinContent(iBin+1)
            hToCompare[iFile].SetBinContent(iBin+1, (bc -1) * lambdaParam +1)
            hToCompare[iFile].SetBinError(iBin+1, hToCompare[iFile].GetBinError(iBin+1) * lambdaParam)
    else:
        ScaleGraph(hToCompare[iFile], scale)
        
        # apply scaling using lambda parameters
        for iBin in range(hToCompare[iFile].GetN()):
            bc = hToCompare[iFile].GetPointY(iBin)
            hToCompare[iFile].SetPointY(iBin, (bc -1) * lambdaParam +1)
            hToCompare[iFile].SetPointError(iBin,
                hToCompare[iFile].GetErrorXlow(iBin),
                hToCompare[iFile].GetErrorXhigh(iBin),
                hToCompare[iFile].GetErrorYlow(iBin) * lambdaParam,
                hToCompare[iFile].GetErrorYhigh(iBin) * lambdaParam)
    if doRatio:
        if 'TH' in objType:
            if drawRatioUnc:
                if ratioUncCorr:
                    hRatioToCompare.append(ComputeRatioDiffBins(hToCompare[iFile], hToCompare[0], 'B'))
                else:
                    hRatioToCompare.append(ComputeRatioDiffBins(hToCompare[iFile], hToCompare[0]))
            else:
                hRatioToCompare.append(ComputeRatioDiffBins(hToCompare[iFile], hToCompare[0]))
                for iBin in range(1, hRatioToCompare[iFile].GetNbinsX()+1):
                    hRatioToCompare[iFile].SetBinError(iBin, 1.e-20)
            hRatioToCompare[iFile].SetDirectory(0)
        else:
            if drawRatioUnc:
                if ratioUncCorr:
                    print('WARNING: correlated uncertainty in ratio for TGraphs not implemented. Switching off')
                    ratioUncCorr = False
                     #TODO: extend ComputeRatioGraph to account for correlated uncertainties
                else:
                    hRatioToCompare.append(ComputeRatioGraph(hToCompare[iFile], hToCompare[0]))
            else:
                hRatioToCompare.append(ComputeRatioGraph(hToCompare[iFile], hToCompare[0]))
                for iBin in range(hRatioToCompare[iFile].GetN()):
                    hRatioToCompare[iFile].SetPointEYlow(iBin, 1.e-20)
                    hRatioToCompare[iFile].SetPointEYhigh(iBin, 1.e-20)
        #TODO: add case to manage ratio between graph and histo (utility function already available in AnalysisUtils)
        hRatioToCompare[iFile].SetName(f'hRatio{iFile}')
        SetObjectStyle(hRatioToCompare[iFile],
                       color=GetROOTColor(color),
                       markerstyle=GetROOTMarker(marker),
                       markersize=markersize,
                       linewidth=linewidth,
                       fillstyle=fillstyle,
                       fillalpha=fillalpha)
    if doCompareUnc:
        if 'TH' in objType:
            hUncToCompare.append(hToCompare[iFile].Clone(f'hUncToCompare{iFile}'))
            for iBin in range(1, hUncToCompare[iFile].GetNbinsX()+1):
                unc = hUncToCompare[iFile].GetBinError(iBin)
                cent = hUncToCompare[iFile].GetBinContent(iBin)
                if compareRelUnc:
                    unctocomp = unc/cent if cent != 0 else 0
                    hUncToCompare[iFile].SetBinContent(iBin, unctocomp)
                else:
                    hUncToCompare[iFile].SetBinContent(iBin, unc)
                hUncToCompare[iFile].SetBinError(iBin, 1.e-20)
            hUncToCompare[iFile].SetDirectory(0)
            SetObjectStyle(hUncToCompare[iFile],
                           color=GetROOTColor(color),
                           markerstyle=GetROOTMarker(marker),
                           markersize=markersize,
                           linewidth=linewidth,
                           fillstyle=fillstyle,
                           fillalpha=fillalpha)
        else:
            #TODO: add uncertainty comparison for TGraphs
            print('WARNING: uncertainty comparison for TGraphs not implemented. Switching off')
            doCompareUnc = False

    leg.AddEntry(hToCompare[iFile], legNames[iFile], legOpt[iFile])

ratios, RMS, shift = [], [], []
if doRatio and displayRMS:
    if 'TH' in objType:
        nPoints = hRatioToCompare[1].GetNbinsX()
    else:
        nPoints = hRatioToCompare[1].GetN()
    for iBin in range(nPoints):
        ratios.append([])
        for iFile, _ in enumerate(inFileNames):
            if iFile == 0:
                continue
            if 'TH' in objType:
                ratios[iBin].append(hRatioToCompare[iFile].GetBinContent(iBin+1))
            else:
                ratios[iBin].append(hRatioToCompare[iFile].GetPointY(iBin))

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
if  not avoidLeg:
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
        if inputCfg.get('DrawGrid'):
            cOut.SetGridx()
            cOut.SetGridy()
        cOut.SaveAs(f'{outFileName}.{ext}')

input("Press enter to exit")
