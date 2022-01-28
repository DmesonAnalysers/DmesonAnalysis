'''
Script to do template fits to MC and estimate the fraction of primary/secondary/material pions and kaons
'''

import argparse
import sys
from typing import Collection
import ctypes
import numpy as np
import yaml
import os
from ROOT import TCanvas, TFractionFitter, TFile, TObjArray, TGraphErrors
from ROOT import THStack, kRed, gROOT, TLatex, TLegend, TH1D
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor, GetROOTMarker, DivideCanvas

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', help='yaml config file name')
parser.add_argument('-b', action='store_true', help='Run un batch mode')
args = parser.parse_args()

gROOT.SetBatch(args.b)

with open(args.cfgFileName, 'r') as ymlConfig:
    cfg = yaml.load(ymlConfig, yaml.FullLoader)

channel = cfg['input']['channel']
dataFile = TFile(cfg['input']['data'], 'read')
mcFile = TFile(cfg['input']['mc'], 'read')
sourceList = cfg['sources']
sourceLabelList = cfg['sourcelabels']
colors = [GetROOTColor(color) for color in cfg['colors']]
tolerance = cfg['tolerance']
rebinPt = cfg['rebinpt']
dcaFitRange = cfg['dcafitrange']
dcaFracRange = cfg['dcafracrange']
ptRange = cfg['ptrange']
useLogY = cfg['output']['logy']
oFileName = cfg['output']['filename']

# fractions-vs-pt plot settings
fracNLegCols = cfg['output']['fractionsVsPt']['legcols']
fracLegPos = cfg['output']['fractionsVsPt']['legposition']
fracYRange = cfg['output']['fractionsVsPt']['yrange']

# data-and-prediction plot settings
dataPredNLegCols = cfg['output']['dataprediction']['legcols']
dataPredLegPos = cfg['output']['dataprediction']['legposition']
dataPredYRange = cfg['output']['dataprediction']['yrange']

# data-and-MC-templates plot settings
dataTemplatesNLegCols = cfg['output']['mctemplates']['legcols']
dataTemplLegPos = cfg['output']['mctemplates']['legposition']
dataTemplatesYRange = cfg['output']['mctemplates']['yrange']
dataTemplatesYRebin = cfg['output']['mctemplates']['rebin']

if channel == 'DK':
    buddy = 'kaon'
elif channel == 'Dpi':
    buddy = 'pion'
else:
    print(f'\033[91mError:\033[0m the {channel} channel is not supported. Exit!')
    sys.exit()

oFile = TFile(oFileName + f'_{buddy}' + '.root', 'recreate')

dirNameDataList = [f'HM_CharmFemto_D{buddy}_TrackCuts0/HM_CharmFemto_D{buddy}_TrackCuts0',
                   f'HM_CharmFemto_D{buddy}_AntiTrackCuts0/HM_CharmFemto_D{buddy}_AntiTrackCuts0']
dirNameMCList = [f'HM_CharmFemto_D{buddy}_TrackCutsMC0/HM_CharmFemto_D{buddy}_TrackCutsMC0',
                 f'HM_CharmFemto_D{buddy}_AntiTrackCutsMC0/HM_CharmFemto_D{buddy}_AntiTrackCutsMC0']


for dirNameData, dirNameMC in zip(dirNameDataList, dirNameMCList):
    # Update buddy name according to particle/antiparticle
    if '_AntiTrackCuts' in dirNameData:
        buddy = f'anti{buddy}'

    hDcaPtList_data = dataFile.Get(dirNameData)
    hDcaPtList_mc = mcFile.Get(dirNameMC)

    hDcaPt_data = hDcaPtList_data.FindObject('DCAXYPtBinningTot')
    hDcaPt_data.GetXaxis().SetRangeUser(ptRange[0], ptRange[1])
    hDcaPt_data.RebinX(rebinPt)
    hDcaPt_data.GetYaxis().SetRangeUser(-dcaFitRange, dcaFitRange)

    nPtBins = hDcaPt_data.GetNbinsX()

    ptMins = [hDcaPt_data.GetXaxis().GetBinLowEdge(iPtBin + 1) for iPtBin in range(0, nPtBins)]
    ptMaxs = [hDcaPt_data.GetXaxis().GetBinLowEdge(iPtBin + 2) for iPtBin in range(0, nPtBins)]

    ptMins = list(filter(lambda pt: ptRange[0] < pt and pt < ptRange[1], ptMins))
    ptMaxs = ptMaxs[:len(ptMins)]
    nPtBins = len(ptMins)
    hMcToProject = []

    # Make plots
    for iSource, source in enumerate(sourceList):
        hDcaPt_mc = hDcaPtList_mc.FindObject('DCAPtBinning').FindObject(f'DCAPtBinning{source}')
        hDcaPt_mc.RebinX(rebinPt)
        hDcaPt_mc.GetYaxis().SetRangeUser(-dcaFitRange, dcaFitRange)

        if hDcaPt_data.GetNbinsX() != hDcaPt_mc.GetNbinsX():
            print(f'\033[31mError\033[0m: number of pT bins are not the same for data ({hDcaPt_data.GetNbinsX()}) \
                and MC ({hDcaPt_mc.GetNbinsX()}). Exit!')
            sys.exit()

        for iPtBin in range(nPtBins):
            if hDcaPt_data.GetXaxis().GetBinLowEdge(iPtBin + 1) - hDcaPt_mc.GetXaxis().GetBinLowEdge(iPtBin + 1) > 1e-6:
                print('\033[31mError\033[0m: low edge of the pT bins are not the same for data and MC. Exit!')
                sys.exit()
        hMcToProject.append(hDcaPt_mc)

    fractionsInFracRange, fractionsUncInFracRange = [], []
    dataToDraw = []
    cocktail = []
    mcToDraw = [[] for _ in range(nPtBins)]
    fractionsInMC = []

    for iPtBin, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        hDataToFit = hDcaPt_data.ProjectionY(f'hDca_{iPtBin+1}_data', iPtBin + 1, iPtBin + 1)
        dataToDraw.append(hDataToFit.Clone(f'hDataToDraw_{iPtBin}_data'))

        nSources = len(sourceList)
        # Create list of MC templates
        hMcToFit = TObjArray(nSources)
        for iSource, source in enumerate(sourceList):
            hMC = hMcToProject[iSource].ProjectionY(f'hDca_{source}_{iPtBin+1}_mc', iPtBin + 1, iPtBin + 1)
            hMcToFit.Add(hMC)
            mcToDraw[iPtBin].append(hMC.Clone(f'hMCToDraw_{source}_{iPtBin+1}'))
        cocktail.append(TFractionFitter(hDataToFit, hMcToFit, 'Q'))

        # Compute fractions in MC, used as a reference to constrain the fit parameters
        totalMCIntegral = sum([hMcToFit.At(iSource).Integral() for iSource in range(nSources)])
        fractionsInMC.append([hMcToFit.At(iSource).Integral() / totalMCIntegral for iSource in range(nSources)])

        for iSource in range(nSources):
            fracMin = fractionsInMC[iPtBin][iSource] * (1 - tolerance)
            fracMax = min(1, fractionsInMC[iPtBin][iSource] / (1 - tolerance))
            cocktail[iPtBin].Constrain(iSource, fracMin, fracMax)
        status = cocktail[iPtBin].Fit()

        fractionsInData = [ctypes.c_double() for iSource in range(nSources)]
        fractionsInDataUnc = [ctypes.c_double() for iSource in range(nSources)]

        for iSource in range(nSources):
            cocktail[iPtBin].GetResult(iSource, fractionsInData[iSource], fractionsInDataUnc[iSource])
        totalDataIntegral = sum([float(fractionsInData[iSource].value) for iSource in range(nSources)])

        print('      From TFractionFitter          |      in MC')
        for source, fData, fDataUnc, fMC in zip(sourceList, fractionsInData, fractionsInDataUnc, fractionsInMC[iPtBin]):
            print(f'{source.ljust(10)}: {fData.value:.6f} +/- {fDataUnc.value:.6f},   {fMC:.6f}')
        print('sum frac form TFractionFitter: ', totalDataIntegral)

        uncList = []
        contributions = []
        for iSource in range(nSources):
            predSource = cocktail[iPtBin].GetMCPrediction(iSource)
            minBin = int(predSource.GetXaxis().FindBin(-dcaFitRange))
            maxBin = int(predSource.GetXaxis().FindBin(+dcaFitRange))

            uncCounts = ctypes.c_double()
            contributions.append(predSource.IntegralAndError(minBin, maxBin, uncCounts))
            uncList.append(uncCounts.value)

        fractionsInFracRange.append([contrib / sum(contributions) for contrib in contributions])
        fractionsUncInFracRange.append([unc / sum(contributions) for unc in uncList])

    # Make drawings
    oFile.mkdir(f'HM_CharmFemto_D{buddy}')
    oFile.cd(f'HM_CharmFemto_D{buddy}')

    SetGlobalStyle(
        textsize=0.04,
        labelsize=0.035,
        titlesize=0.035,
        titleoffsety=1.8,
        padleftmargin=0.12,
        padtopmargin=0.06,
        padbottommargin=0.12,
        padrightmargin=0.1)

    cDataTemplatesStack = TCanvas(f'cDataTemplatesStack_{buddy}', f'MC DCA templates - {buddy}', 1920, 1080)
    DivideCanvas(cDataTemplatesStack, nPtBins)
    cDataTemplatesStack.cd(iPtBin + 1)

    cPred = TCanvas(f'cPred_{buddy}', f'Data and prediction - {buddy}', 1920, 1080)
    DivideCanvas(cPred, nPtBins)

    hDataTemplatesStack = []

    for iPtBin, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        if useLogY:
            cPred.cd(iPtBin + 1).SetLogy()

        # Draw Data and prediction
        cPred.cd(iPtBin + 1)
        cPred.SetBottomMargin(0)
        cPred.SetTopMargin(0)
        hPred = cocktail[iPtBin].GetPlot()
        SetObjectStyle(hPred, linecolor=kRed, linewidth=1)
        hPred.SetLineColor(GetROOTColor('kRed'))
        SetObjectStyle(dataToDraw[iPtBin])

        if isinstance(dataPredYRange, list):
            dataToDraw[iPtBin].GetYaxis().SetRangeUser(dataPredYRange[0], dataPredYRange[1])
        elif dataPredYRange == 'auto':
            pass
        else:
            print(f'\033[91mError:\033[0m the range option {fracYRange} is not supported. Exit!')
            sys.exit()

        dataToDraw[iPtBin].SetTitle(';DCA_{xy} (cm);Entries (a.u.)')
        dataToDraw[iPtBin].DrawCopy('pe')
        dataToDraw[iPtBin].Write(f'hData_pT_{ptMin*1000:.0f}_{ptMax*1000:.0f}')
        hPred.SetTitle(';DCA_{xy} (cm);Entries (a.u.)')
        hPred.DrawCopy('same hist')
        hPred.Write(f'hPrediction_{source}_pT_{ptMin*1000:.0f}_{ptMax*1000:.0f}')

        if iPtBin == 0:
            legPred = TLegend(dataPredLegPos[0], dataPredLegPos[1], dataPredLegPos[2], dataPredLegPos[3])
            legPred.AddEntry(dataToDraw[iPtBin], 'Data', 'lp')
            legPred.AddEntry(hPred, 'Prediction', 'l')
            legPred.Draw('same')
        latPred = TLatex()
        latPred.SetTextSize(0.04)
        chi2NDF = cocktail[iPtBin].GetChisquare() / cocktail[iPtBin].GetNDF()
        latPred.DrawLatexNDC(0.65, 0.85, f'#chi^{{2}}/NDF = {chi2NDF:.2f}')
        latPred.DrawLatexNDC(0.65, 0.8, f'{ptMin:.2f} < #it{{p}}_{{T}} < {ptMax:.2f} GeV/c')

        if useLogY:
            cDataTemplatesStack.cd(iPtBin + 1).SetLogy()
        hDataTemplatesStack.append(THStack(f'hDataTemplatesStack_{buddy}',
                                           'MC DCA templates;DCA_{xy} (cm); Entries (a.u.)'))

        # stack with MC templates and data
        legStack = TLegend(dataTemplLegPos[0], dataTemplLegPos[1], dataTemplLegPos[2], dataTemplLegPos[3])
        hDataToDraw = dataToDraw[iPtBin].Clone('dataToDraw_stack')  # Data
        SetObjectStyle(hDataToDraw, markerstyle=20, SetPadLeftMargin=0.3)
        hDataToDraw.Rebin(dataTemplatesYRebin)
        hDataToDraw.Scale(1. / hDataToDraw.Integral())
        if isinstance(dataTemplatesYRange, list):
            hDataToDraw.GetYaxis().SetRangeUser(dataTemplatesYRange[0], dataTemplatesYRange[1])
        elif dataTemplatesYRange == 'auto':
            pass
        else:
            print(f'\033[91mError:\033[0m the range option {fracYRange} is not supported. Exit!')
            sys.exit()
        legStack.AddEntry(hDataToDraw, 'Data')
        hDataTemplatesStack[-1].Add(hDataToDraw)

        for iSource, (source, sourceLabel, color) in enumerate(zip(sourceList, sourceLabelList, colors, )):  # MC
            hMCTemplateToDraw = mcToDraw[iPtBin][iSource]
            hMCTemplateToDraw.Write(f'hMCTemplate_{source}_pT_{ptMin*1000:.0f}_{ptMax*1000:.0f}')
            SetObjectStyle(hMCTemplateToDraw, markercolor=color, linecolor=color, markerstyle=21 + iSource)
            hMCTemplateToDraw.Rebin(dataTemplatesYRebin)
            hMCTemplateToDraw.Scale(1. / hMCTemplateToDraw.Integral() if hMCTemplateToDraw.Integral() > 0 else 1)
            legStack.AddEntry(hMCTemplateToDraw, sourceLabel)
            hDataTemplatesStack[-1].Add(hMCTemplateToDraw)

        cDataTemplatesStack.cd(iPtBin + 1)
        hDataTemplatesStack[-1].Draw('nostack hist')

        legStack.Draw('same l')
        latStack = TLatex()
        latStack.SetTextSize(0.04)

        latStack.DrawLatexNDC(0.65, 0.85, f'{ptMin:.2f} < #it{{p}}_{{T}} < {ptMax:.2f} GeV/c')
        for iSource, (source, frac, fracUnc) in enumerate(zip(sourceList, fractionsInData, fractionsInDataUnc)):
            latStack.DrawLatexNDC(0.65, 0.8 - 0.05 * iSource,
                                  f'f_{{{source}}} = {frac.value:.4f} #pm {fracUnc.value:.4f}')

    cDataTemplatesStack.Modified()
    cDataTemplatesStack.Update()

    cPred.Modified()
    cPred.Update()

    cPred.SaveAs(f'{oFileName}_{buddy}_pred.pdf')
    cDataTemplatesStack.SaveAs(f'{oFileName}_{buddy}_shapes.pdf')

    # Fractions vs pT
    cFracVsPt = TCanvas(f'cFracVsPt_{buddy}', f'Fractions vs pT - {buddy}', 600, 600)
    gFraction = [TGraphErrors(1) for _ in range(nSources)]
    for iPtBin, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        SetObjectStyle(gFraction[iSource])
        x = (ptMax + ptMin) / 2
        xUnc = (ptMax - ptMin) / 2
        for iSource in range(nSources):
            y = fractionsInFracRange[iPtBin][iSource]
            yUnc = fractionsUncInFracRange[iPtBin][iSource]
            gFraction[iSource].SetPoint(iPtBin, x, y)
            gFraction[iSource].SetPointError(iPtBin, xUnc, yUnc)
    gFraction[0].SetTitle(f';#it{{p}}_{{T}} (GeV/c);fraction')
    if isinstance(fracYRange, list):
        gFraction[0].GetYaxis().SetRangeUser(fracYRange[0], fracYRange[1])
    elif fracYRange == 'auto':
        pass
    else:
        print(f'\033[91mError:\033[0m the range option {fracYRange} is not supported. Exit!')
        sys.exit()

    if useLogY:
        cFracVsPt.SetLogy()
    legFrac = TLegend(fracLegPos[0], fracLegPos[1], fracLegPos[2], fracLegPos[3])
    legFrac.SetNColumns(fracNLegCols)
    if buddy == 'kaon':
        legTitle = 'K^{+}'
    elif buddy == 'antikaon':
        legTitle = 'K^{-}'
    elif buddy == 'pion':
        legTitle = '#pi^{+}'
    elif buddy == 'antipion':
        legTitle = '#pi^{-}'
    legFrac.SetHeader(legTitle, 'C')
    legFrac.SetTextSize(0.04)
    for iSource, (source, sourceLabel, color)in enumerate(zip(sourceList, sourceLabelList, colors)):
        gFraction[iSource].SetName(f'gFraction_{sourceList[iSource]}')
        SetObjectStyle(gFraction[iSource], markercolor=color, linecolor=color, markerstyle=21 + iSource)
        legFrac.AddEntry(gFraction[iSource], sourceLabel, 'lp')
        if iSource == 0:
            gFraction[iSource].Draw('ape')
        gFraction[iSource].Draw('pe same')
        gFraction[iSource].Write()
    legFrac.Draw('same')
    cFracVsPt.SaveAs(f'{oFileName}_{buddy}_fractions.pdf')
    cFracVsPt.Write()

    # Compute fraction averages
    hFracAritmAverage = TH1D('hFracAritmAverage', 'Fractions - aritmetic average;;fraction', nSources, 0, nSources)
    for iSource, (sourceLabel, gFrac) in enumerate(zip(sourceLabelList, gFraction)):
        hFracAritmAverage.GetXaxis().SetBinLabel(iSource + 1, sourceLabel)
        fracAverage = sum([gFrac.GetPointY(iPoint) for iPoint in range(gFrac.GetN())])/gFrac.GetN()
        
        hFracAritmAverage.SetBinContent(iSource + 1, fracAverage)
    hFracAritmAverage.Write()
