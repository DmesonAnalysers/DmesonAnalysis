'''
python script for the computation of the prompt / non-prompt fraction with the cut-variation method
run: python ComputeCutVarPromptFrac.py cfgFileName.yml outFileName.root
'''

import argparse
import os
from itertools import product
import numpy as np
import yaml
from ROOT import TFile, TH1F, TH2F, TCanvas, TLegend, TGraphAsymmErrors, TLatex  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure, kGreen, kRainBow, kFullCircle, kFullSquare, kOpenSquare, kOpenCircle  # pylint: disable=import-error,no-name-in-module
from utils.AnalysisUtils import GetPromptFDYieldsAnalyticMinimisation, GetPromptFDFractionFc, GetFractionNb
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('outFileName', metavar='text', default='outFile.root',
                    help='output root file name')
args = parser.parse_args()

outFileNameEffPDF = args.outFileName.replace('.root', '_Eff.pdf')
outFileNameDistrPDF = args.outFileName.replace('.root', '_Distr.pdf')
outFileNameFracPDF = args.outFileName.replace('.root', '_Frac.pdf')
outFileNameCorrMatrixPDF = args.outFileName.replace('.root', '_CorrMatrix.pdf')

with open(args.cfgFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)

if len(cutSetCfg['rawyields']['inputfiles']) != len(cutSetCfg['efficiencies']['inputfiles']):
    print('ERROR: number or raw yield files and efficiency files not consistent! Please check your config file. Exit')
    exit()

nSets = len(cutSetCfg['rawyields']['inputfiles'])
inputDirRaw = cutSetCfg['rawyields']['inputdir']
inputDirEff = cutSetCfg['efficiencies']['inputdir']

areCorrSets = cutSetCfg['minimisation']['correlated']

hRawYields, hEffPrompt, hEffFD, hEv = [], [], [], []
for iCutSet, (inFileNameRawYield, inFileNameEff) in enumerate(
        zip(cutSetCfg['rawyields']['inputfiles'], cutSetCfg['efficiencies']['inputfiles'])):

    inFileNameRawYield = os.path.join(inputDirRaw, inFileNameRawYield)
    inFileRawYield = TFile.Open(inFileNameRawYield)
    hRawYields.append(inFileRawYield.Get(cutSetCfg['rawyields']['histoname']))
    hRawYields[-1].SetDirectory(0)
    hEv.append(inFileRawYield.Get('hEvForNorm'))
    hEv[-1].SetDirectory(0)

    inFileNameEff = os.path.join(inputDirEff, inFileNameEff)
    inFileEff = TFile.Open(inFileNameEff)
    hEffPrompt.append(inFileEff.Get(cutSetCfg['efficiencies']['histonames']['prompt']))
    hEffFD.append(inFileEff.Get(cutSetCfg['efficiencies']['histonames']['feeddown']))
    hEffPrompt[-1].SetDirectory(0)
    hEffFD[-1].SetDirectory(0)

# load inputs for theory-driven methods
compareToFc = cutSetCfg['theorydriven']['enableFc']
compareToNb = cutSetCfg['theorydriven']['enableNb']
hCrossSecPrompt, hCrossSecFD = [], []
if compareToFc or compareToNb:
    inFileCrossSec = TFile.Open(cutSetCfg['theorydriven']['predictions']['crosssec']['inputfile'])

    hCrossSecPrompt.append(inFileCrossSec.Get(
        f"{cutSetCfg['theorydriven']['predictions']['crosssec']['histonames']['prompt']}_central"))
    hCrossSecPrompt.append(inFileCrossSec.Get(
        f"{cutSetCfg['theorydriven']['predictions']['crosssec']['histonames']['prompt']}_min"))
    hCrossSecPrompt.append(inFileCrossSec.Get(
        f"{cutSetCfg['theorydriven']['predictions']['crosssec']['histonames']['prompt']}_max"))

    hCrossSecFD.append(inFileCrossSec.Get(
        f"{cutSetCfg['theorydriven']['predictions']['crosssec']['histonames']['feeddown']}_central_corr"))
    hCrossSecFD.append(inFileCrossSec.Get(
        f"{cutSetCfg['theorydriven']['predictions']['crosssec']['histonames']['feeddown']}_min_corr"))
    hCrossSecFD.append(inFileCrossSec.Get(
        f"{cutSetCfg['theorydriven']['predictions']['crosssec']['histonames']['feeddown']}_max_corr"))

    if compareToNb:
        BR = cutSetCfg['theorydriven']['BR']
        sigmaMB = cutSetCfg['theorydriven']['sigmaMB']

SetGlobalStyle(padleftmargin=0.15, padtopmargin=0.08, titleoffsetx=1., titleoffsety=1.4, opttitle=1, palette=kRainBow)

legDistr = TLegend(0.45, 0.69, 0.75, 0.89)
legDistr.SetFillStyle(0)
legDistr.SetBorderSize(0)
legDistr.SetTextSize(0.045)

legEff = TLegend(0.2, 0.2, 0.4, 0.4)
legEff.SetFillStyle(0)
legEff.SetBorderSize(0)
legEff.SetTextSize(0.045)

legFrac = TLegend(0.2, 0.79, 0.4, 0.89)
legFrac.SetFillStyle(0)
legFrac.SetBorderSize(0)
legFrac.SetTextSize(0.045)

latInfo = TLatex()
latInfo.SetNDC()
latInfo.SetTextSize(0.045)
latInfo.SetTextFont(42)
latInfo.SetTextColor(1)

hCorrYieldPrompt = hRawYields[0].Clone('hCorrYieldPrompt')
hCorrYieldPrompt.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{N}_{prompt}')
SetObjectStyle(hCorrYieldPrompt, color=kRed+1, fillcolor=kRed+1, markerstyle=kFullCircle)

hCorrYieldFD = hRawYields[0].Clone('hCorrYieldFD')
hCorrYieldFD.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{N}_{non-prompt}')
SetObjectStyle(hCorrYieldFD, color=kAzure+4, fillcolor=kAzure+4, markerstyle=kFullSquare)

hCovCorrYields = [[hRawYields[0].Clone('hCovPromptPrompt'), hRawYields[0].Clone('hCovPromptFD')],
                  [hRawYields[0].Clone('hCovFDPrompt'), hRawYields[0].Clone('hCovFDFD')]]
for iRow, row in enumerate(hCovCorrYields):
    for iCol, hCov in enumerate(row):
        SetObjectStyle(hCov, linecolor=kBlack)
        if iRow == 0:
            rowName = '#it{N}_{prompt}'
        else:
            rowName = '#it{N}_{non-prompt}'
        if iCol == 0:
            colName = '#it{N}_{prompt}'
        else:
            colName = '#it{N}_{non-prompt}'

        hCov.SetTitle(f';#it{{p}}_{{T}} (GeV/#it{{c}}); #sigma({rowName}, {colName})')

hRawYieldsVsCut, hRawYieldsVsCutReSum, hRawYieldPromptVsCut, hRawYieldFDVsCut, cDistr = ([] for _ in range(5))
hEffPromptVsCut, hEffFDVsCut, cEff = ([] for _ in range(3))
hPromptFracVsCut, hFDFracVsCut, gPromptFracFcVsCut, gFDFracFcVsCut, gPromptFracNbVsCut, gFDFracNbVsCut, cFrac = (
    [] for _ in range(7))
hCorrMatrixCutSets, cCorrMatrix = ([] for _ in range(2))

for iPt in range(hRawYields[0].GetNbinsX()):
    ptMin = hRawYields[0].GetBinLowEdge(iPt+1)
    ptMax = ptMin + hRawYields[0].GetBinWidth(iPt+1)

    listRawYield = [hRaw.GetBinContent(iPt+1) for hRaw in hRawYields]
    listRawYieldUnc = [hRaw.GetBinError(iPt+1) for hRaw in hRawYields]

    listEffPrompt = [hEffP.GetBinContent(iPt+1) for hEffP in hEffPrompt]
    listEffPromptUnc = [hEffP.GetBinError(iPt+1) for hEffP in hEffPrompt]

    listEffFD = [hEffF.GetBinContent(iPt+1) for hEffF in hEffFD]
    listEffFDUnc = [hEffF.GetBinError(iPt+1) for hEffF in hEffFD]

    corrYields, covMatrixCorrYields, chiSquare, matrices = GetPromptFDYieldsAnalyticMinimisation(\
        listEffPrompt, listEffFD, listRawYield, listEffPromptUnc, listEffFDUnc, listRawYieldUnc, areCorrSets)

    hCorrYieldPrompt.SetBinContent(iPt+1, corrYields.item(0))
    hCorrYieldPrompt.SetBinError(iPt+1, np.sqrt(covMatrixCorrYields.item(0, 0)))
    hCorrYieldFD.SetBinContent(iPt+1, corrYields.item(1))
    hCorrYieldFD.SetBinError(iPt+1, np.sqrt(covMatrixCorrYields.item(1, 1)))
    for covElem in product(range(2), range(2)):
        hCovCorrYields[covElem[0]][covElem[1]].SetBinContent(iPt+1, covMatrixCorrYields.item(covElem))
        hCovCorrYields[covElem[0]][covElem[1]].SetBinError(iPt+1, 0.)

    hRawYieldsVsCut.append(TH1F(f'hRawYieldsVsCutPt_pT{ptMin:.0f}_{ptMax:.0f}',
                                f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;raw yield',
                                nSets, 0.5, nSets+0.5))
    hRawYieldsVsCutReSum.append(TH1F(f'hRawYieldsVsCutReSum_pT{ptMin:.0f}_{ptMax:.0f}',
                                     f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;raw yield',
                                     nSets, 0.5, nSets+0.5))
    hRawYieldPromptVsCut.append(TH1F(f'hRawYieldPromptVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                                     f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;raw yield',
                                     nSets, 0.5, nSets+0.5))
    hRawYieldFDVsCut.append(TH1F(f'hRawYieldFDVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                                 f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;raw yield',
                                 nSets, 0.5, nSets+0.5))
    hEffPromptVsCut.append(TH1F(f'hEffPromptVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                                f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;efficiency',
                                nSets, 0.5, nSets+0.5))
    hEffFDVsCut.append(TH1F(f'hEffFDVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                            f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;efficiency',
                            nSets, 0.5, nSets+0.5))
    hPromptFracVsCut.append(TH1F(f'hPromptFracVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                                 f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;#it{{f}}_{{prompt}}',
                                 nSets, 0.5, nSets+0.5))
    hFDFracVsCut.append(TH1F(f'hFDFracVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                             f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;#it{{f}}_{{FD}}',
                             nSets, 0.5, nSets+0.5))

    SetObjectStyle(hRawYieldsVsCut[iPt], linecolor=kBlack, markercolor=kBlack, markerstyle=kFullSquare)
    SetObjectStyle(hRawYieldPromptVsCut[iPt], color=kRed+1, fillcolor=kRed+1, markerstyle=kFullCircle, fillalpha=0.3)
    SetObjectStyle(hRawYieldFDVsCut[iPt], color=kAzure+4, fillcolor=kAzure+4, markerstyle=kFullSquare, fillalpha=0.3)
    SetObjectStyle(hRawYieldsVsCutReSum[iPt], linecolor=kGreen+2)
    SetObjectStyle(hEffPromptVsCut[iPt], color=kRed+1, markerstyle=kFullCircle)
    SetObjectStyle(hEffFDVsCut[iPt], color=kAzure+4, markerstyle=kFullSquare)
    SetObjectStyle(hPromptFracVsCut[iPt], color=kRed+1, markerstyle=kFullCircle)
    SetObjectStyle(hFDFracVsCut[iPt], color=kAzure+4, markerstyle=kFullSquare)

    hCorrMatrixCutSets.append(TH2F(f'hCorrMatrixCutSets_pT{ptMin:.0f}_{ptMax:.0f}',
                                   f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;cut set',
                                   nSets, 0.5, nSets+0.5, nSets, 0.5, nSets+0.5))

    for iCutSetRow in range(nSets):
        for iCutSetCol in range(nSets):
            hCorrMatrixCutSets[iPt].SetBinContent(
                iCutSetRow+1, iCutSetCol+1, matrices['corrMatrix'].item(iCutSetRow, iCutSetCol))

    # cross sections from theory if comparison enabled
    if compareToFc or compareToNb:
        crossSecPrompt = [h.Integral(h.GetXaxis().FindBin(
            ptMin*1.0001), h.GetXaxis().FindBin(ptMax*0.9999)) / (ptMax-ptMin) for h in hCrossSecPrompt]
        crossSecFD = [h.Integral(h.GetXaxis().FindBin(
            ptMin*1.0001), h.GetXaxis().FindBin(ptMax*0.9999)) / (ptMax-ptMin) for h in hCrossSecFD]

        if compareToFc:
            gPromptFracFcVsCut.append(TGraphAsymmErrors(nSets))
            gPromptFracFcVsCut[iPt].SetNameTitle(
                f'gPromptFracFcVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;#it{{f}}_{{prompt}}')

            gFDFracFcVsCut.append(TGraphAsymmErrors(nSets))
            gFDFracFcVsCut[iPt].SetNameTitle(
                f'gFDFracFcVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;#it{{f}}_{{FD}}')

            SetObjectStyle(gPromptFracFcVsCut[iPt], color=kRed+3, markerstyle=kOpenCircle)
            SetObjectStyle(gFDFracFcVsCut[iPt], color=kAzure+3, markerstyle=kOpenSquare)

        if compareToNb:
            gPromptFracNbVsCut.append(TGraphAsymmErrors(nSets))
            gPromptFracNbVsCut[iPt].SetNameTitle(
                f'gPromptFracNbVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;#it{{f}}_{{prompt}}')

            gFDFracNbVsCut.append(TGraphAsymmErrors(nSets))
            gFDFracNbVsCut[iPt].SetNameTitle(
                f'gFDFracNbVsCut_pT{ptMin:.0f}_{ptMax:.0f}',
                f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;#it{{f}}_{{FD}}')

            SetObjectStyle(gPromptFracNbVsCut[iPt], color=kRed-7, markerstyle=kOpenCircle)
            SetObjectStyle(gFDFracNbVsCut[iPt], color=kAzure+5, markerstyle=kOpenSquare)

    for iCutSet, (rawY, effP, effF, rawYunc, effPunc, effFunc) in enumerate(zip(
            listRawYield, listEffPrompt, listEffFD, listRawYieldUnc, listEffPromptUnc, listEffFDUnc)):

        # efficiency
        hEffPromptVsCut[iPt].SetBinContent(iCutSet+1, effP)
        hEffPromptVsCut[iPt].SetBinError(iCutSet+1, effPunc)

        hEffFDVsCut[iPt].SetBinContent(iCutSet+1, effF)
        hEffFDVsCut[iPt].SetBinError(iCutSet+1, effFunc)

        # raw yields (including prompt and non-prompt raw yields)
        hRawYieldsVsCut[iPt].SetBinContent(iCutSet+1, rawY)
        hRawYieldsVsCut[iPt].SetBinError(iCutSet+1, rawYunc)

        hRawYieldPromptVsCut[iPt].SetBinContent(iCutSet+1, corrYields.item(0)*effP)
        hRawYieldPromptVsCut[iPt].SetBinError(iCutSet+1, np.sqrt(covMatrixCorrYields.item(0, 0))*effP)

        hRawYieldFDVsCut[iPt].SetBinContent(iCutSet+1, corrYields.item(1)*effF)
        hRawYieldFDVsCut[iPt].SetBinError(iCutSet+1, np.sqrt(covMatrixCorrYields.item(1, 1))*effF)

        hRawYieldsVsCutReSum[iPt].SetBinContent(iCutSet+1, \
            hRawYieldPromptVsCut[iPt].GetBinContent(iCutSet+1)+hRawYieldFDVsCut[iPt].GetBinContent(iCutSet+1))

        # prompt fraction
        fPrompt = effP * corrYields.item(0) / (effP * corrYields.item(0) + effF * corrYields.item(1))
        defPdeNP = (effP * (effP * corrYields.item(0) + effF * corrYields.item(1)) - effP**2
                    * corrYields.item(0)) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
        defPdeNF = - effF * effP * corrYields.item(0) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
        fPromptUnc = np.sqrt(defPdeNP**2 * covMatrixCorrYields.item(0, 0) + \
            defPdeNF**2 * covMatrixCorrYields.item(1, 1) + \
            2 * defPdeNP * defPdeNF * covMatrixCorrYields.item(1, 0))

        # feed-down fraction
        fFD = effF * corrYields.item(1) / (effP * corrYields.item(0) + effF * corrYields.item(1))
        defFdeNF = (effF * (effF * corrYields.item(1) + effP * corrYields.item(0)) - effF**2
                    * corrYields.item(1)) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
        defFdeNP = - effF * effP * corrYields.item(1) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
        fFDUnc = np.sqrt(defFdeNF**2 * covMatrixCorrYields.item(1, 1) + \
            defFdeNP**2 * covMatrixCorrYields.item(0, 0) + \
            2 * defFdeNF * defFdeNP * covMatrixCorrYields.item(1, 0))

        hPromptFracVsCut[iPt].SetBinContent(iCutSet+1, fPrompt)
        hPromptFracVsCut[iPt].SetBinError(iCutSet+1, fPromptUnc)
        hFDFracVsCut[iPt].SetBinContent(iCutSet+1, fFD)
        hFDFracVsCut[iPt].SetBinError(iCutSet+1, fFDUnc)

        # theory-driven, if enabled
        if compareToFc:
            fPromptFc, fFDFc = GetPromptFDFractionFc(effP, effF, crossSecPrompt, crossSecFD)
            gPromptFracFcVsCut[iPt].SetPoint(iCutSet, iCutSet+1, fPromptFc[0])
            gPromptFracFcVsCut[iPt].SetPointError(iCutSet, 0.5, 0.5, fPromptFc[0]-fPromptFc[1],
                                                  fPromptFc[2]-fPromptFc[0])
            gFDFracFcVsCut[iPt].SetPoint(iCutSet, iCutSet+1, fFDFc[0])
            gFDFracFcVsCut[iPt].SetPointError(iCutSet, 0.5, 0.5, fFDFc[0]-fFDFc[1], fFDFc[2]-fFDFc[0])

        if compareToNb:
            fPromptNb = GetFractionNb(rawY, effP, effF, crossSecFD, ptMax-ptMin, 1., BR,
                                      hEv[iCutSet].GetBinContent(1), sigmaMB)
            fFDNb = [1-fPromptNb[0], 1-fPromptNb[2], 1-fPromptNb[1]] #inverse Nb method not reliable

            gPromptFracNbVsCut[iPt].SetPoint(iCutSet, iCutSet+1, fPromptNb[0])
            gPromptFracNbVsCut[iPt].SetPointError(iCutSet, 0.5, 0.5, fPromptNb[0]-fPromptNb[1],
                                                  fPromptNb[2]-fPromptNb[0])
            gFDFracNbVsCut[iPt].SetPoint(iCutSet, iCutSet+1, fFDNb[0])
            gFDFracNbVsCut[iPt].SetPointError(iCutSet, 0.5, 0.5, fFDNb[0]-fFDNb[1], fFDNb[2]-fFDNb[0])

    if iPt == 0:
        legDistr.AddEntry(hRawYieldsVsCut[iPt], 'Measured raw yield', 'lpe')
        legDistr.AddEntry(hRawYieldPromptVsCut[iPt], 'Prompt', 'f')
        legDistr.AddEntry(hRawYieldFDVsCut[iPt], 'Non-prompt', 'f')
        legDistr.AddEntry(hRawYieldsVsCutReSum[iPt], 'Prompt + non-prompt', 'l')

        legEff.AddEntry(hEffPromptVsCut[iPt], 'Prompt', 'lpe')
        legEff.AddEntry(hEffFDVsCut[iPt], 'Non-prompt', 'lpe')

        legFrac.AddEntry(hPromptFracVsCut[iPt], 'Prompt', 'lpe')
        legFrac.AddEntry(hFDFracVsCut[iPt], 'Non-prompt', 'lpe')

        deltaY = 0.
        if compareToFc:
            legFrac.AddEntry(gPromptFracFcVsCut[iPt], 'Prompt #it{f}_{c}', 'lpe')
            legFrac.AddEntry(gFDFracFcVsCut[iPt], 'Non-prompt #it{f}_{c}', 'lpe')
            deltaY += 0.1
            legFrac.SetY1(0.83-deltaY)
        if compareToNb:
            legFrac.AddEntry(gPromptFracNbVsCut[iPt], 'Prompt #it{N}_{b}', 'lpe')
            legFrac.AddEntry(gFDFracNbVsCut[iPt], 'Non-prompt #it{N}_{b}', 'lpe')
            deltaY += 0.1
            legFrac.SetY1(0.83-deltaY)

    cEff.append(TCanvas(f'cEff_pT{ptMin:.0f}_{ptMax:.0f}', '', 800, 800))
    cEff[iPt].DrawFrame(0.5, 1.e-5, nSets+0.5, 1.,
                        f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;efficiency')
    cEff[iPt].SetLogy()
    hEffPromptVsCut[iPt].DrawCopy('same')
    hEffFDVsCut[iPt].DrawCopy('same')
    legEff.Draw()

    cDistr.append(TCanvas(f'cDistr_pT{ptMin:.0f}_{ptMax:.0f}', '', 800, 800))
    cDistr[iPt].DrawFrame(0.5, 0., nSets+0.5, hRawYieldsVsCut[iPt].GetMaximum()*1.2,
                          f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;raw yield')
    hRawYieldsVsCut[iPt].Draw('same')
    hRawYieldPromptVsCut[iPt].DrawCopy('histsame')
    hRawYieldFDVsCut[iPt].DrawCopy('histsame')
    hRawYieldsVsCutReSum[iPt].Draw('same')
    legDistr.Draw()
    latInfo.DrawLatex(0.47, 0.65, f'#chi^{{2}} / ndf = {chiSquare:.3f}')

    cFrac.append(TCanvas(f'cFrac_pT{ptMin:.0f}_{ptMax:.0f}', '', 800, 800))
    cFrac[iPt].DrawFrame(0.5, 0., nSets+0.5, 1.8,
                         f'{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}};cut set;fraction')
    hPromptFracVsCut[iPt].DrawCopy('Esame')
    hFDFracVsCut[iPt].DrawCopy('Esame')
    if compareToFc:
        gPromptFracFcVsCut[iPt].Draw('PZ')
        gFDFracFcVsCut[iPt].Draw('PZ')
    if compareToNb:
        gPromptFracNbVsCut[iPt].Draw('PZ')
        gFDFracNbVsCut[iPt].Draw('PZ')
    legFrac.Draw()

    cCorrMatrix.append(TCanvas(f'cCorrMatrix_pT{ptMin:.0f}_{ptMax:.0f}', '', 800, 800))
    cCorrMatrix[-1].cd().SetRightMargin(0.14)
    hCorrMatrixCutSets[iPt].Draw('colz')

nPtBins = hCorrYieldPrompt.GetNbinsX()
cCorrYield = TCanvas('cCorrYield', '', 800, 800)
cCorrYield.DrawFrame(hCorrYieldPrompt.GetBinLowEdge(1), 1.,
                     hCorrYieldPrompt.GetBinLowEdge(nPtBins)+hCorrYieldPrompt.GetBinWidth(nPtBins),
                     hCorrYieldPrompt.GetMaximum()*1.2, ';#it{p}_{T} (GeV/#it{c});corrected yield')
cCorrYield.SetLogy()
hCorrYieldPrompt.Draw('same')
hCorrYieldFD.Draw('same')
legEff.Draw()

outFile = TFile(args.outFileName, 'recreate')
cCorrYield.Write()
hCorrYieldPrompt.Write()
hCorrYieldFD.Write()
for covElem in product(range(2), range(2)):
    hCovCorrYields[covElem[0]][covElem[1]].Write()
for iPt in range(hRawYields[0].GetNbinsX()):
    cDistr[iPt].Write()
    cEff[iPt].Write()
    cFrac[iPt].Write()
    hRawYieldsVsCut[iPt].Write()
    hRawYieldPromptVsCut[iPt].Write()
    hRawYieldFDVsCut[iPt].Write()
    hRawYieldsVsCutReSum[iPt].Write()
    hEffPromptVsCut[iPt].Write()
    hEffFDVsCut[iPt].Write()
    hPromptFracVsCut[iPt].Write()
    hFDFracVsCut[iPt].Write()
    hCorrMatrixCutSets[iPt].Write()
outFile.Close()

for iPt in range(hRawYields[0].GetNbinsX()):
    if iPt == 0:
        cEff[iPt].SaveAs(f'{outFileNameEffPDF}[')
        cDistr[iPt].SaveAs(f'{outFileNameDistrPDF}[')
        cFrac[iPt].SaveAs(f'{outFileNameFracPDF}[')
        cCorrMatrix[iPt].SaveAs(f'{outFileNameCorrMatrixPDF}[')
    cEff[iPt].SaveAs(outFileNameEffPDF)
    cDistr[iPt].SaveAs(outFileNameDistrPDF)
    cFrac[iPt].SaveAs(outFileNameFracPDF)
    cCorrMatrix[iPt].SaveAs(outFileNameCorrMatrixPDF)
    if iPt == hRawYields[0].GetNbinsX()-1:
        cEff[iPt].SaveAs(f'{outFileNameEffPDF}]')
        cDistr[iPt].SaveAs(f'{outFileNameDistrPDF}]')
        cFrac[iPt].SaveAs(f'{outFileNameFracPDF}]')
        cCorrMatrix[iPt].SaveAs(f'{outFileNameCorrMatrixPDF}]')

input('Press enter to exit')
