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
from ROOT import THStack, kRed, gROOT, TLatex, TLegend
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor, GetROOTMarker, DivideCanvas

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', help='yaml config file name')
parser.add_argument('--onlyMC', action='store_true', help='Compute the fractions only in the MC sample')
parser.add_argument('-b', action='store_true', help='Run un batch mode')
args = parser.parse_args()

gROOT.SetBatch(args.b)

with open(args.cfgFileName, 'r') as ymlConfig:
    cfg = yaml.load(ymlConfig, yaml.FullLoader)

channel = cfg['input']['channel']
dataFile = TFile(cfg['input']['data'], 'read')
mcFile = TFile(cfg['input']['mc'], 'read')
sourceList = cfg['sources']
colors = cfg['colors']
colors = [GetROOTColor(color) for color in colors]
tolerance = cfg['tolerance']
rebin = cfg['rebin']
dcaFitRange = cfg['dcafitrange']
dcaFracRange = cfg['dcafracrange']
ptRange = cfg['ptrange']
fractoinRange = cfg['fractionrange']
useLogY = cfg['output']['logy']
sourceLegend = cfg['sourcelegend']
oFileName = cfg['output']['filename']

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
    hDcaPt_data.RebinX(rebin)
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
        hDcaPt_mc.RebinX(rebin)
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

        # Create list of MC templates
        hMcToFit = TObjArray(len(sourceList))
        for iSource, source in enumerate(sourceList):
            hMC = hMcToProject[iSource].ProjectionY(f'hDca_{source}_{iPtBin+1}_mc', iPtBin + 1, iPtBin + 1)
            hMcToFit.Add(hMC)
            mcToDraw[iPtBin].append(hMC.Clone(f'hMCToDraw_{source}_{iPtBin+1}'))
        cocktail.append(TFractionFitter(hDataToFit, hMcToFit, 'Q'))

        # Compute fractions in MC, used as a reference to constrain the fit parameters
        totalMCIntegral = sum([hMcToFit.At(iSource).Integral() for iSource in range(len(sourceList))])
        fractionsInMC.append([hMcToFit.At(iSource).Integral() / totalMCIntegral for iSource in range(len(sourceList))])
        if args.onlyMC:
            continue
        for iSource in range(len(sourceList)):
            cocktail[iPtBin].Constrain(iSource, fractionsInMC[iPtBin][iSource] * (1 - tolerance), min(1, fractionsInMC[iPtBin][iSource] / (1 - tolerance)))
        status = cocktail[iPtBin].Fit()

        fractionsInData = [ctypes.c_double() for iSource in range(len(sourceList))]
        fractionsInDataUnc = [ctypes.c_double() for iSource in range(len(sourceList))]

        for iSource in range(len(sourceList)):
            cocktail[iPtBin].GetResult(iSource, fractionsInData[iSource], fractionsInDataUnc[iSource])
        totalDataIntegral = sum([float(fractionsInData[iSource].value) for iSource in range(len(sourceList))])

        print('      From TFractionFitter          |      in MC')
        for source, fData, fDataUnc, fMC in zip(sourceList, fractionsInData, fractionsInDataUnc, fractionsInMC[iPtBin]):
            print(f'{source.ljust(10)}: {fData.value:.6f} +/- {fDataUnc.value:.6f},   {fMC:.6f}')
        print('sum frac form TFractionFitter: ', totalDataIntegral)

        uncList = []
        contributions = []
        for iSource in range(len(sourceList)):
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
    # gStyle.SetPalette(52)
    SetGlobalStyle(
        textsize=0.03,
        labelsize=0.035,
        titlesize=0.035,
        titleoffsety=1.8,
        padleftmargin=0.12,
        padtopmargin=0.06,
        padbottommargin=0.12,
        padrightmargin=0.1)

    # todo: fix this part
    # cDataTemplatesStack = TCanvas('cDataTemplatesStack', 'MC DCA templates', 1920, 1080)
    # DivideCanvas(cDataTemplatesStack, nPtBins)

    cPred = TCanvas('cPred', 'Data and prediction', 1920, 1080)
    DivideCanvas(cPred, nPtBins)

    hDataTemplatesStack = []

    for iPtBin, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        if useLogY:
            cPred.cd(iPtBin + 1).SetLogy()
            # cDataTemplatesStack.cd(iPtBin + 1).SetLogy() # todo: fix this part
        if not args.onlyMC:
            # Draw Data and prediction
            cPred.cd(iPtBin + 1)
            cPred.SetBottomMargin(0)
            cPred.SetTopMargin(0)
            hPred = cocktail[iPtBin].GetPlot()
            SetObjectStyle(hPred, linecolor=kRed, linewidth=1)
            hPred.SetTitle(';DCA_{xy} (cm); Entries (a.u.)')
            hPred.SetLineColor(GetROOTColor('kRed'))
            SetObjectStyle(dataToDraw[iPtBin])
            dataToDraw[iPtBin].DrawCopy('pe')
            dataToDraw[iPtBin].Write(f'hData_pT_{ptMin*1000:.0f}_{ptMax*1000:.0f}')
            if iPtBin == 0:
                legPred = TLegend(0.15, 0.8, 0.4, 0.9)
                legPred.AddEntry(dataToDraw[iPtBin], 'Data', 'lp')
                legPred.AddEntry(hPred, 'Prediction', 'l')
                legPred.Draw('same')
            hPred.DrawCopy('same hist')
            hPred.Write(f'hPrediction_{source}_pT_{ptMin*1000:.0f}_{ptMax*1000:.0f}')
            latPred = TLatex()
            latPred.SetTextSize(0.03)
            latPred.DrawLatexNDC(0.65, 0.85, f'#chi^{{2}}/NDF = {cocktail[iPtBin].GetChisquare()/cocktail[iPtBin].GetNDF():.2f}')
            latPred.DrawLatexNDC(0.65, 0.8, f'{ptMin:.2f} < #it{{p}}_{{T}} < {ptMax:.2f} GeV/c')

        # stack templates # todo: fix this part
        # cDataTemplatesStack.cd(iPtBin + 1)
        hDataTemplatesStack.append(THStack(f'hDataTemplatesStack_{buddy}', 'MC DCA templates;DCA_{xy} (cm); Entries (a.u.)'))

        # legStack = TLegend(0.15, 0.55, 0.4, 0.92)
        # if not args.onlyMC:
        #     hDataToDraw = dataToDraw[iPtBin].Clone('dataToDraw_stack')  # Data
        #     SetObjectStyle(hDataToDraw, markerstyle=20, SetPadLeftMargin=0.3)
        #     hDataToDraw.Scale(1. / hDataToDraw.Integral())
        #     legStack.AddEntry(hDataToDraw, 'Data')
        #     hDataTemplatesStack[-1].Add(hDataToDraw)
        for iSource, source in enumerate(sourceList):  # MC
            hMCTemplateToDraw = mcToDraw[iPtBin][iSource]
            hMCTemplateToDraw.Write(f'hMCTemplate_{source}_pT_{ptMin*1000:.0f}_{ptMax*1000:.0f}')
        #     SetObjectStyle(hMCTemplateToDraw, markercolor=colors[iSource], linecolor=colors[iSource], markerstyle=21 + iSource)
        #     hMCTemplateToDraw.Scale(1. / hMCTemplateToDraw.Integral() if hMCTemplateToDraw.Integral() > 0 else 1)
        #     legStack.AddEntry(hMCTemplateToDraw, sourceList[iSource])
        #     hDataTemplatesStack[-1].Add(hMCTemplateToDraw)

        # hDataTemplatesStack[-1].DrawClone('nostack plm')
        # hDataTemplatesStack[-1].Write()
        # legStack.Draw('same plm')
        # latStack = TLatex()
        # latStack.SetTextSize(0.03)
        # latStack.DrawLatexNDC(0.65, 0.85, f'{ptMin:.2f} < #it{{p}}_{{T}} < {ptMax:.2f} GeV/c')
        # for iSource, source in enumerate(sourceList):
        #     latStack.DrawLatexNDC(
        #         0.65, 0.8 - 0.05 * iSource, f'f_{{{source}}} = {fractionsInData[iSource].value if not args.onlyMC else fractionsInMC[iPtBin][iSource]:.4f} #pm {fractionsInDataUnc[iSource].value if not args.onlyMC else 0:.4f}')

    # cDataTemplatesStack.Modified() # todo: fix this partali
    # cDataTemplatesStack.Update() # todo: fix this part

    cPred.SaveAs(f'{oFileName}_{buddy}_pred.pdf')
    # cDataTemplatesStack.SaveAs(f'{oFileName}_{buddy}_shapes.pdf')  # todo: fix this part

    # Fractions vs pT

    cFracVsPt = TCanvas(f'cFracVsPt_{buddy}', f'Fractions vs pT - {buddy}', 600, 600)
    gFraction = [TGraphErrors(1) for _ in range(len(sourceList))]
    for iPtBin, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
        SetObjectStyle(gFraction[iSource])
        x = (ptMax + ptMin) / 2
        xUnc = (ptMax - ptMin) / 2
        for iSource in range(len(sourceList)):
            if args.onlyMC:
                y = fractionsInMC[iPtBin][iSource]
                yUnc = 0
            else:
                y = fractionsInFracRange[iPtBin][iSource]
                yUnc = fractionsUncInFracRange[iPtBin][iSource]
            gFraction[iSource].SetPoint(iPtBin, x, y)
            gFraction[iSource].SetPointError(iPtBin, xUnc, yUnc)

    gFraction[0].GetYaxis().SetRangeUser(fractoinRange[0], fractoinRange[1])
    gFraction[0].SetTitle(f'Contribution to the sample with {-dcaFitRange} #leq DCA_{{xy}} #leq {-dcaFitRange} (cm);#it{{p}}_{{T}} (GeV/c);fraction')
    if useLogY:
        cFracVsPt.SetLogy()
    legFrac = TLegend(0.35, 0.2)
    legFrac.SetNColumns(2)
    if buddy == 'kaon':
        legTitle = 'K^{+}'
    elif buddy == 'antikaon':
        legTitle = 'K^{-}'
    elif buddy == 'pion':
        legTitle = '#pi^{+}'
    elif buddy == 'antipion':
        legTitle = '#pi^{-}'
    legFrac.SetHeader(legTitle, 'C')
    for iSource in range(len(sourceList)):
        gFraction[iSource].SetName(f'gFraction_{buddy}_{sourceList[iSource]}')
        SetObjectStyle(gFraction[iSource], markercolor=colors[iSource], linecolor=colors[iSource], markerstyle=21 + iSource)
        legFrac.AddEntry(gFraction[iSource], sourceLegend[iSource], 'lp')
        if iSource == 0:
            gFraction[iSource].Draw('ape')
        gFraction[iSource].Draw('pe same')
        gFraction[iSource].Write()
    legFrac.Draw('same')
    cFracVsPt.SaveAs(f'{oFileName}_{buddy}_fractions.pdf')
    cFracVsPt.Write()
oFile.Close()
