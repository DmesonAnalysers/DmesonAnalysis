'''
python script for the computation of the prompt / non-prompt fraction with the cut-variation method
run: python ComputeCutVarPromptFrac.py cfgFileName.yml outFileName.root
'''

import argparse
import numpy as np
import yaml
from ROOT import TFile, TH1F, TCanvas, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure, kGreen, kFullCircle, kFullSquare, kOpenSquare  # pylint: disable=import-error,no-name-in-module
from utils.AnalysisUtils import GetPromptFDYieldsAnalyticMinimisation
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

with open(args.cfgFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)

if len(cutSetCfg['rawyields']['inputfiles']) != len(cutSetCfg['efficiencies']['inputfiles']):
    print('ERROR: number or raw yield files and efficiency files not consistent! Please check your config file. Exit')
    exit()

SetGlobalStyle(padleftmargin=0.15, titleoffsetx=1., titleoffsety=1.4)

nSets = len(cutSetCfg['rawyields']['inputfiles'])

hRawYields, hEffPrompt, hEffFD = [], [], []
for iCutSet, (inFileNameRawYield, inFileNameEff) in enumerate(
        zip(cutSetCfg['rawyields']['inputfiles'], cutSetCfg['efficiencies']['inputfiles'])):

    inFileRawYield = TFile.Open(inFileNameRawYield)
    hRawYields.append(inFileRawYield.Get(cutSetCfg['rawyields']['histoname']))
    hRawYields[-1].SetDirectory(0)

    inFileEff = TFile.Open(inFileNameEff)
    hEffPrompt.append(inFileEff.Get(cutSetCfg['efficiencies']['histonames']['prompt']))
    hEffFD.append(inFileEff.Get(cutSetCfg['efficiencies']['histonames']['feeddown']))
    hEffPrompt[-1].SetDirectory(0)
    hEffFD[-1].SetDirectory(0)

legDistr = TLegend(0.4, 0.7, 0.7, 0.9)
legDistr.SetFillStyle(0)
legDistr.SetBorderSize(0)
legDistr.SetTextSize(0.05)

legEff = TLegend(0.2, 0.2, 0.4, 0.4)
legEff.SetFillStyle(0)
legEff.SetBorderSize(0)
legEff.SetTextSize(0.05)

hCorrYieldPrompt = hRawYields[0].Clone('hCorrYieldPrompt')
hCorrYieldPrompt.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{N}_{prompt}')
SetObjectStyle(hCorrYieldPrompt, color=kRed+1, fillcolor=kRed + 1, markerstyle=kFullCircle, fillalpha=0.3)

hCorrYieldFD = hRawYields[0].Clone('hCorrYieldFD')
hCorrYieldFD.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{N}_{non-prompt}')
SetObjectStyle(hCorrYieldFD, color=kAzure+4, fillcolor=kAzure + 4, markerstyle=kOpenSquare, fillalpha=0.3)

hRawYieldsVsCut, hRawYieldsVsCutReSum, hRawYieldPromptVsCut, hRawYieldFDVsCut, cDistr = ([] for _ in range(5))
hEffPromptVsCut, hEffFDVsCut, cEff = ([] for _ in range(3))
hPromptFracVsCut, hFDFracVsCut, cFrac = ([] for _ in range(3))
for iPt in range(hRawYields[0].GetNbinsX()):
    ptMin = hRawYields[0].GetBinLowEdge(iPt+1)
    ptMax = ptMin + hRawYields[0].GetBinWidth(iPt+1)

    listRawYield = [hRaw.GetBinContent(iPt+1) for hRaw in hRawYields]
    listRawYieldUnc = [hRaw.GetBinError(iPt+1) for hRaw in hRawYields]

    listEffPrompt = [hEffP.GetBinContent(iPt+1) for hEffP in hEffPrompt]
    listEffPromptUnc = [hEffP.GetBinError(iPt+1) for hEffP in hEffPrompt]

    listEffFD = [hEffF.GetBinContent(iPt+1) for hEffF in hEffFD]
    listEffFDUnc = [hEffF.GetBinError(iPt+1) for hEffF in hEffFD]

    corrYields, covMatrixCorrYields = GetPromptFDYieldsAnalyticMinimisation(\
        listEffPrompt, listEffFD, listRawYield, listEffPromptUnc, listEffFDUnc, listRawYieldUnc)

    hCorrYieldPrompt.SetBinContent(iPt+1, corrYields.item(0))
    hCorrYieldPrompt.SetBinError(iPt+1, np.sqrt(covMatrixCorrYields.item(0, 0)))
    hCorrYieldFD.SetBinContent(iPt+1, corrYields.item(1))
    hCorrYieldFD.SetBinError(iPt+1, np.sqrt(covMatrixCorrYields.item(1, 1)))

    hRawYieldsVsCut.append(TH1F('hRawYieldsVsCutPt_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), ';cut set;raw yield',
                                nSets, 0.5, nSets+0.5))
    hRawYieldsVsCutReSum.append(TH1F('hRawYieldsVsCutReSum_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax),
                                     ';cut set;raw yield', nSets, 0.5, nSets+0.5))
    hRawYieldPromptVsCut.append(TH1F('hRawYieldPromptVsCut_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax),
                                     ';cut set;raw yield', nSets, 0.5, nSets+0.5))
    hRawYieldFDVsCut.append(TH1F('hRawYieldFDVsCut_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), ';cut set;raw yield',
                                 nSets, 0.5, nSets+0.5))
    hEffPromptVsCut.append(TH1F('hEffPromptVsCut_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), ';cut set;efficiency',
                                nSets, 0.5, nSets+0.5))
    hEffFDVsCut.append(TH1F('hEffFDVsCut_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), ';cut set;efficiency',
                            nSets, 0.5, nSets+0.5))
    hPromptFracVsCut.append(TH1F('hPromptFracVsCut_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), ';cut set;#it{f}_{prompt}',
                            nSets, 0.5, nSets+0.5))
    hFDFracVsCut.append(TH1F('hFDFracVsCut_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), ';cut set;#it{f}_{FD}',
                        nSets, 0.5, nSets+0.5))

    SetObjectStyle(hRawYieldsVsCut[iPt], linecolor=kBlack, markercolor=kBlack, markerstyle=kFullSquare)
    SetObjectStyle(hRawYieldPromptVsCut[iPt], color=kRed+1, fillcolor=kRed+1, markerstyle=kFullCircle, fillalpha=0.3)
    SetObjectStyle(hRawYieldFDVsCut[iPt], color=kAzure+4, fillcolor=kAzure+4, markerstyle=kOpenSquare, fillalpha=0.3)
    SetObjectStyle(hRawYieldsVsCutReSum[iPt], linecolor=kGreen+2)
    SetObjectStyle(hEffPromptVsCut[iPt], color=kRed+1, markerstyle=kFullCircle, fillalpha=0.3)
    SetObjectStyle(hEffFDVsCut[iPt], color=kAzure+4, markerstyle=kOpenSquare, fillalpha=0.3)
    SetObjectStyle(hPromptFracVsCut[iPt], color=kRed+1, markerstyle=kFullCircle, fillalpha=0.3)
    SetObjectStyle(hFDFracVsCut[iPt], color=kAzure+4, markerstyle=kOpenSquare, fillalpha=0.3)

    for iCutSet, (rawY, effP, effF, rawYunc, effPunc, effFunc) in enumerate(zip(
            listRawYield, listEffPrompt, listEffFD, listRawYieldUnc, listEffPromptUnc, listEffFDUnc)):

        #efficiency
        hEffPromptVsCut[iPt].SetBinContent(iCutSet+1, effP)
        hEffPromptVsCut[iPt].SetBinError(iCutSet+1, effPunc)

        hEffFDVsCut[iPt].SetBinContent(iCutSet+1, effF)
        hEffFDVsCut[iPt].SetBinError(iCutSet+1, effFunc)

        #raw yields (including prompt and non-prompt raw yields)
        hRawYieldsVsCut[iPt].SetBinContent(iCutSet+1, rawY)
        hRawYieldsVsCut[iPt].SetBinError(iCutSet+1, rawYunc)

        hRawYieldPromptVsCut[iPt].SetBinContent(iCutSet+1, corrYields.item(0)*effP)
        hRawYieldPromptVsCut[iPt].SetBinError(iCutSet+1, np.sqrt(covMatrixCorrYields.item(0, 0))*effP)

        hRawYieldFDVsCut[iPt].SetBinContent(iCutSet+1, corrYields.item(1)*effF)
        hRawYieldFDVsCut[iPt].SetBinError(iCutSet+1, np.sqrt(covMatrixCorrYields.item(1, 1))*effF)

        hRawYieldsVsCutReSum[iPt].SetBinContent(iCutSet+1, \
            hRawYieldPromptVsCut[iPt].GetBinContent(iCutSet+1)+hRawYieldFDVsCut[iPt].GetBinContent(iCutSet+1))

        #prompt fraction
        fPrompt = effP * corrYields.item(0) / (effP * corrYields.item(0) + effF * corrYields.item(1))
        defPdeNP = (effP * (effP * corrYields.item(0) + effF * corrYields.item(1)) - effP**2
                    * corrYields.item(0)) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
        defPdeNF = - effF * effP * corrYields.item(0) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
        fPromptUnc = np.sqrt(defPdeNP**2 * covMatrixCorrYields.item(0, 0) + \
            defPdeNF**2 * covMatrixCorrYields.item(1, 1) + \
            2 * defPdeNP * defPdeNF * covMatrixCorrYields.item(1, 0))


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

    if iPt == 0:
        legDistr.AddEntry(hRawYieldsVsCut[iPt], 'Measured raw yield', 'lpe')
        legDistr.AddEntry(hRawYieldPromptVsCut[iPt], 'Prompt', 'f')
        legDistr.AddEntry(hRawYieldFDVsCut[iPt], 'Non-prompt', 'f')
        legDistr.AddEntry(hRawYieldsVsCutReSum[iPt], 'Prompt + non-prompt', 'l')

        legEff.AddEntry(hEffPromptVsCut[iPt], 'Prompt', 'lpe')
        legEff.AddEntry(hEffFDVsCut[iPt], 'Non-prompt', 'lpe')

    cEff.append(TCanvas('cEff_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), '', 800, 800))
    cEff[iPt].DrawFrame(0.5, 1.e-5, nSets+0.5, 1., ';cut set;efficiency')
    cEff[iPt].SetLogy()
    hEffPromptVsCut[iPt].DrawCopy('same')
    hEffFDVsCut[iPt].DrawCopy('same')
    legEff.Draw()

    cDistr.append(TCanvas('cDistr_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), '', 800, 800))
    cDistr[iPt].DrawFrame(0.5, 0., nSets+0.5, hRawYieldsVsCut[iPt].GetMaximum()*1.2, ';cut set;raw yield')
    hRawYieldsVsCut[iPt].Draw('same')
    hRawYieldPromptVsCut[iPt].DrawCopy('histsame')
    hRawYieldFDVsCut[iPt].DrawCopy('histsame')
    hRawYieldsVsCutReSum[iPt].Draw('same')
    legDistr.Draw()

    cFrac.append(TCanvas('cFrac_pT{0:.0f}_{1:.0f}'.format(ptMin, ptMax), '', 800, 800))
    cFrac[iPt].DrawFrame(0.5, 0., nSets+0.5, 1., ';cut set;fraction')
    hPromptFracVsCut[iPt].DrawCopy('Esame')
    hFDFracVsCut[iPt].DrawCopy('Esame')
    legEff.Draw()

nPtBins = hCorrYieldPrompt.GetNbinsX()
cCorrYield = TCanvas('cCorrYield', '', 800, 800)
cCorrYield.DrawFrame(hCorrYieldPrompt.GetBinLowEdge(1), 1.,
                     hCorrYieldPrompt.GetBinLowEdge(nPtBins)+hCorrYieldPrompt.GetBinWidth(nPtBins),
                     hCorrYieldPrompt.GetMaximum()*1.2, ';#it{p}_{T} (GeV/#it{c});corrected yield')
cCorrYield.SetLogy()
hCorrYieldPrompt.Draw('same')
hCorrYieldFD.Draw('same')

outFile = TFile(args.outFileName, 'recreate')
cCorrYield.Write()
hCorrYieldPrompt.Write()
hCorrYieldFD.Write()
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
outFile.Close()

for iPt in range(hRawYields[0].GetNbinsX()):
    if iPt == 0:
        cEff[iPt].SaveAs('{}['.format(outFileNameEffPDF))
        cDistr[iPt].SaveAs('{}['.format(outFileNameDistrPDF))
        cFrac[iPt].SaveAs('{}['.format(outFileNameFracPDF))
    cEff[iPt].SaveAs(outFileNameEffPDF)
    cDistr[iPt].SaveAs(outFileNameDistrPDF)
    cFrac[iPt].SaveAs(outFileNameFracPDF)
    if iPt == hRawYields[0].GetNbinsX()-1:
        cEff[iPt].SaveAs('{}]'.format(outFileNameEffPDF))
        cDistr[iPt].SaveAs('{}]'.format(outFileNameDistrPDF))
        cFrac[iPt].SaveAs('{}]'.format(outFileNameFracPDF))

input('Press enter to exit')
