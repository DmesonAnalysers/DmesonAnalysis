'''
python script for the evaluation of the D-meson S/B as a function of k*
'''

import sys
import os
import argparse
import ctypes
import numpy as np
import yaml
from ROOT import TFile, TCanvas, TGraphAsymmErrors, TLatex, TF1, TH1F, TGaxis, TLegend, TDirectoryFile, gPad, gRandom # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kAzure, kBlack, kOpenCircle # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle # pylint: disable=import-error,wrong-import-position

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='yaml config file name')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlFile:
    cfg = yaml.load(ymlFile, yaml.FullLoader)

rawYiedlFileName = cfg['rawyields']['filename']
SoverBHistoName = cfg['rawyields']['histoname']
inFileName = cfg['input']['filename']
suffix = cfg['input']['suffix']
prefix = cfg['input']['prefix']
outDirName = cfg['output']['directory']
outSuffix = cfg['output']['suffix']
addTDirectory = cfg['output']['addTDirectory']

# setting global style
SetGlobalStyle(padbottommargin=0.13, padtopmargin=0.075, titleoffsety=1.2, maxdigits=2)

# load S/B histo and compute purity
inFile = TFile.Open(rawYiedlFileName)
hSoverB = inFile.Get(SoverBHistoName)
SetObjectStyle(hSoverB)
hSignal = inFile.Get(SoverBHistoName.replace('SoverB', 'Signal'))
hBkg = inFile.Get(SoverBHistoName.replace('SoverB', 'Bkg'))
hSoverB.SetDirectory(0)
hSignal.SetDirectory(0)
hBkg.SetDirectory(0)
hSoverB.SetTitle(';#it{p}_{T}(D^{#pm}) (GeV/#it{c});S/B')
hPurity = hSoverB.Clone('hPurity')
hPurity.GetYaxis().SetTitle('S/(S+B)')
hPurity.SetDirectory(0)
for iPt in range(hSoverB.GetNbinsX()):
    bkg = hBkg.GetBinContent(iPt+1)
    bkgErr = hBkg.GetBinError(iPt+1)
    sgn = hSignal.GetBinContent(iPt+1)
    sgnErr = hSignal.GetBinError(iPt+1)
    sgnPlusBkg = sgn+bkg
    purity = sgn/sgnPlusBkg
    purityErr = np.sqrt((bkg/sgnPlusBkg**2)**2 * sgnErr**2 + (sgn/sgnPlusBkg**2)**2 * bkgErr**2)
    hPurity.SetBinContent(iPt+1, purity)
    hPurity.SetBinError(iPt+1, purityErr)
ptMax = hSoverB.GetXaxis().GetBinUpEdge(hSoverB.GetNbinsX())
inFile.Close()

# load corr histos
corrNames = ['Particle0_Particle2', 'Particle1_Particle3', 'Particle0_Particle3', 'Particle1_Particle2']
corrTitles = {'Particle0_Particle2': 'p - D^{+}',
              'Particle0_Particle3': 'p - D^{#font[122]{-}}',
              'Particle1_Particle2': '#bar{p} - D^{+}',
              'Particle1_Particle3': '#bar{p} - D^{#font[122]{-}}'}

inFile = TFile.Open(inFileName)
listName = f'{prefix}_CharmFemto_ResultQA{suffix}/{prefix}_CharmFemto_ResultQA{suffix}'
inList = inFile.Get(listName)

kStarDelta = 0.2
kStarMins = [iK * kStarDelta for iK in range(15)]
kStarMaxs = [(iK+1) * kStarDelta for iK in range(15)]

hSEPairVsPtVsKstar, hSEPairVsPt = {}, {}
for corrName in corrNames:
    listCorr = inList.FindObject(f'QA_{corrName}')
    hSEPairVsPtVsKstar[corrName], hSEPairVsPt[corrName] = [], []
    for kStarMin, kStarMax in zip(kStarMins, kStarMaxs):
        hSEPairVsPtVsKstar[corrName].append(listCorr.FindObject(f'KstarPtSEPartTwo_{corrName}'))
        hSEPairVsPtVsKstar[corrName][-1].SetDirectory(0)
        kStarMaxBin = hSEPairVsPtVsKstar[corrName][-1].GetXaxis().FindBin(kStarMax)
        # assuming constant pT binning
        ptRebin = round(hSoverB.GetBinWidth(1) / hSEPairVsPtVsKstar[corrName][-1].GetYaxis().GetBinWidth(1))
        hSEPairVsPt[corrName].append(
            hSEPairVsPtVsKstar[corrName][-1].ProjectionY(f'hSEPairVsPt_{corrName}', 1, kStarMaxBin))
        hSEPairVsPt[corrName][-1].SetDirectory(0)
        hSEPairVsPt[corrName][-1].Rebin(ptRebin)
        SetObjectStyle(hSEPairVsPt[corrName][-1], color=kRed+1, markerstyle=kOpenCircle, fillstyle=0)
        if 'Particle2' in corrName:
            hSEPairVsPt[corrName][-1].SetTitle(';#it{p}_{T}(D^{+}) (GeV/#it{c});SE pairs')
        else:
            hSEPairVsPt[corrName][-1].SetTitle(';#it{p}_{T}(D^{#font[122]{-}}) (GeV/#it{c});SE pairs')
        hSEPairVsPt[corrName][-1].GetXaxis().SetTitleSize(0.05)
        hSEPairVsPt[corrName][-1].GetXaxis().SetLabelSize(0.05)
        hSEPairVsPt[corrName][-1].GetYaxis().SetTitleSize(0.05)
        hSEPairVsPt[corrName][-1].GetYaxis().SetLabelSize(0.05)
        hSEPairVsPt[corrName][-1].GetYaxis().SetDecimals()
        ptMaxPair = hSEPairVsPt[corrName][-1].GetXaxis().GetBinUpEdge(hSEPairVsPt[corrName][-1].GetNbinsX())
        if ptMaxPair < ptMax:
            ptMax = ptMaxPair
inFile.Close()

for corrName in ['Particle0_Particle2_plus_Particle1_Particle3', 'Particle0_Particle3_plus_Particle1_Particle2']:
    corrNamePart = corrName.split('_plus_')
    hSEPairVsPt[corrName] = []
    for iK, (kStarMin, kStarMax) in enumerate(zip(kStarMins, kStarMaxs)):
        hSEPairVsPt[corrName].append(hSEPairVsPt[corrNamePart[0]][iK].Clone(f'hSEPairVsPt_{corrName}'))
        hSEPairVsPt[corrName][iK].Add(hSEPairVsPt[corrNamePart[1]][iK])
        hSEPairVsPt[corrName][iK].SetDirectory(0)
        hSEPairVsPt[corrName][iK].SetTitle(';#it{p}_{T}(D^{#pm}) (GeV/#it{c});SE pairs')
        if iK == 0:
            corrTitles[corrName] = corrTitles[corrNamePart[0]] + ' #oplus ' + corrTitles[corrNamePart[1]]

# compute < S/B > and < purity >
gAverageSoverB, gSoverBAvPt, gAveragePurity, gPurityAvPt = ({} for _ in range(4))
for corrName in ['Particle0_Particle2_plus_Particle1_Particle3', 'Particle0_Particle3_plus_Particle1_Particle2']:
    if corrName == 'Particle0_Particle2_plus_Particle1_Particle3':
        suffix = 'DpPr'
    else:
        suffix = 'DmPr'
    gAverageSoverB[corrName] = TGraphAsymmErrors(0)
    gAverageSoverB[corrName].SetNameTitle(f'gAverageSoverB_{suffix}', ';#it{k}* (GeV/#it{c}); #LT S/B #GT')
    gSoverBAvPt[corrName] = TGraphAsymmErrors(0)
    gSoverBAvPt[corrName].SetNameTitle(f'gSoverBAvPt_{suffix}',
                                       ';#it{k}* (GeV/#it{c}); S/B (#LT #it{p}_{T}(D^{#pm}) #GT)')
    gAveragePurity[corrName] = TGraphAsymmErrors(0)
    gAveragePurity[corrName].SetNameTitle(f'gAveragePurity_{suffix}', ';#it{k}* (GeV/#it{c}); #LT S/(S+B) #GT')
    gPurityAvPt[corrName] = TGraphAsymmErrors(0)
    gPurityAvPt[corrName].SetNameTitle(f'gPurityAvPt_{suffix}',
                                       ';#it{k}* (GeV/#it{c}); S/(S+B) (#LT #it{p}_{T}(D^{#pm}) #GT)')
    SetObjectStyle(gSoverBAvPt[corrName])
    SetObjectStyle(gAverageSoverB[corrName], color=kRed+1, fillstyle=0, markerstyle=kOpenCircle)
    SetObjectStyle(gPurityAvPt[corrName])
    SetObjectStyle(gAveragePurity[corrName], color=kRed+1, fillstyle=0, markerstyle=kOpenCircle)

    hDistrAvSoverB = TH1F(f'hDistrAvSoverB{corrName}', '', 100, 0., hSoverB.GetMaximum())
    hDistrAvPurity = TH1F(f'hDistrAvPurity{corrName}', '', 100, 0., hPurity.GetMaximum())
    for iK, (kStarMin, kStarMax) in enumerate(zip(kStarMins, kStarMaxs)):
        weighAvSoverB, weighAvPurity, totCounts = 0., 0., 0.
        hDistrAvSoverB.Reset()
        hDistrAvPurity.Reset()
        for iGen in range(100): # estimate uncertainty due to uncertainty on weights
            weighAvSoverBSmeared, weighAvPuritySmeared, totCountsSmeared = 0., 0., 0.
            for iPt in range(hSoverB.GetXaxis().FindBin(ptMax)):
                ptCent = hSoverB.GetBinCenter(iPt+1)
                ptBinPair = hSEPairVsPt[corrName][iK].GetXaxis().FindBin(ptCent)
                nPairs = hSEPairVsPt[corrName][iK].GetBinContent(ptBinPair)
                nPairsSmeared = gRandom.PoissonD(hSEPairVsPt[corrName][iK].GetBinContent(ptBinPair))
                weighAvSoverBSmeared += hSoverB.GetBinContent(iPt+1) * nPairsSmeared
                weighAvPuritySmeared += hPurity.GetBinContent(iPt+1) * nPairsSmeared
                totCountsSmeared += nPairsSmeared
                if iGen == 0:
                    weighAvSoverB += hSoverB.GetBinContent(iPt+1) * nPairs
                    weighAvPurity += hPurity.GetBinContent(iPt+1) * nPairs
                    totCounts += nPairs
            hDistrAvSoverB.Fill(weighAvSoverBSmeared/totCountsSmeared)
            hDistrAvPurity.Fill(weighAvPuritySmeared/totCountsSmeared)

        weighAvSoverB /= totCounts
        weighAvSoverBUnc = hDistrAvSoverB.GetRMS()
        weighAvPurity /= totCounts
        weighAvPurityUnc = hDistrAvPurity.GetRMS()

        avPt = hSEPairVsPt[corrName][iK].GetMean()
        avPtUnc = hSEPairVsPt[corrName][iK].GetMeanError()
        # fit data points around mean value
        fSoverB = TF1(f'fSoverB_kstar{iK}', 'pol2', 0., 10.)
        fSoverB.SetLineColor(kAzure+4)
        hSoverB.Fit(f'fSoverB_kstar{iK}', 'Q')
        fPurity = TF1(f'fPurity_kstar{iK}', 'pol4', 0., 10.)
        fPurity.SetLineColor(kAzure+4)
        hPurity.Fit(f'fPurity_kstar{iK}', 'Q')

        kStarCent = (kStarMax + kStarMin) / 2
        gAverageSoverB[corrName].SetPoint(iK, kStarCent, weighAvSoverB)
        gSoverBAvPt[corrName].SetPoint(iK, kStarCent, fSoverB.Eval(avPt))
        gAverageSoverB[corrName].SetPointError(iK, kStarDelta/2, kStarDelta/2, weighAvSoverBUnc, weighAvSoverBUnc)
        gSoverBAvPt[corrName].SetPointError(iK, kStarDelta/2, kStarDelta/2,
                                            fSoverB.Eval(avPt)-fSoverB.Eval(avPt-avPtUnc),
                                            fSoverB.Eval(avPt)-fSoverB.Eval(avPt+avPtUnc))
        gAveragePurity[corrName].SetPoint(iK, kStarCent, weighAvPurity)
        gPurityAvPt[corrName].SetPoint(iK, kStarCent, fPurity.Eval(avPt))
        gAveragePurity[corrName].SetPointError(iK, kStarDelta/2, kStarDelta/2, weighAvPurityUnc, weighAvPurityUnc)
        gPurityAvPt[corrName].SetPointError(iK, kStarDelta/2, kStarDelta/2,
                                            fPurity.Eval(avPt)-fPurity.Eval(avPt-avPtUnc),
                                            fPurity.Eval(avPt)-fPurity.Eval(avPt+avPtUnc))

# plots
legSoverB = TLegend(0.18, 0.7, 0.5, 0.85)
legSoverB.SetTextSize(0.045)
legSoverB.SetBorderSize(0)
legSoverB.SetFillStyle(0)
legSoverB.AddEntry(gAverageSoverB[corrName], 'weighted average', 'lp')
legSoverB.AddEntry(gSoverBAvPt[corrName], 'S/B (#LT #it{p}_{T}(D^{#pm}) #GT) - pol2 parm', 'lp')

legPurity = TLegend(0.18, 0.2, 0.5, 0.45)
legPurity.SetTextSize(0.045)
legPurity.SetBorderSize(0)
legPurity.SetFillStyle(0)
legPurity.AddEntry(gAverageSoverB[corrName], 'weighted average', 'lp')
legPurity.AddEntry(gSoverBAvPt[corrName], 'S/(S+B) (#LT #it{p}_{T}(D^{#pm}) #GT) - pol4 parm', 'lp')

cSoverBvsKstar = TCanvas('cSoverBvsKstar', '', 1000, 500)
cSoverBvsKstar.Divide(2, 1)
for iPad, corrName in enumerate(gAverageSoverB):
    cSoverBvsKstar.cd(iPad+1).DrawFrame(0., 0., 3., hSoverB.GetMaximum()*1.2, ';#it{k}* (GeV/#it{c});#LT S/B #GT')
    gAverageSoverB[corrName].Draw('p')
    gSoverBAvPt[corrName].Draw('p')
    legSoverB.Draw()
cSoverBvsKstar.Modified()
cSoverBvsKstar.Update()

cPurityvsKstar = TCanvas('cPurityvsKstar', '', 1000, 500)
cPurityvsKstar.Divide(2, 1)
for iPad, corrName in enumerate(gAveragePurity):
    cPurityvsKstar.cd(iPad+1).DrawFrame(0., 0., 3., hPurity.GetMaximum()*1.2, ';#it{k}* (GeV/#it{c});#LT S/(S+B) #GT')
    gAveragePurity[corrName].Draw('p')
    gPurityAvPt[corrName].Draw('p')
    legPurity.Draw()
cPurityvsKstar.Modified()
cPurityvsKstar.Update()

info = TLatex()
info.SetTextColor(kBlack)
info.SetTextSize(0.045)
info.SetTextFont(42)
info.SetNDC()

padOrder = {'Particle0_Particle2': 1,
            'Particle1_Particle3': 2,
            'Particle0_Particle2_plus_Particle1_Particle3': 3,
            'Particle0_Particle3': 4,
            'Particle1_Particle2': 5,
            'Particle0_Particle3_plus_Particle1_Particle2': 6}

cSEDistr = TCanvas('cSEDistr', '', 1200, 900)
cSEDistr.Divide(3, 2)
for corrName in hSEPairVsPt:
    cSEDistr.cd(padOrder[corrName])
    hSEPairVsPt[corrName][0].GetYaxis().SetRangeUser(0., hSEPairVsPt[corrName][0].GetMaximum()*2)
    hSEPairVsPt[corrName][0].DrawCopy('e')
    info.DrawLatex(0.2, 0.86, corrTitles[corrName])
    DName = ''
    if 'Particle2' in corrName and 'Particle3' in corrName:
        DName = 'D^{#pm}'
    elif 'Particle2' in corrName:
        DName = 'D^{+}'
    elif 'Particle3' in corrName:
        DName = 'D^{#font[122]{-}}'

    info.DrawLatex(0.2, 0.80, f'#LT #it{{p}}_{{T}} ({DName}) #GT = '
                   f'{hSEPairVsPt[corrName][0].GetMean():0.2f} #pm {hSEPairVsPt[corrName][0].GetMeanError():0.2f}')
    info.DrawLatex(0.2, 0.74, '#it{k}* < 200 MeV/#it{c}')
cSEDistr.Modified()
cSEDistr.Update()

cSoverBvsDistr = TCanvas('cSoverBvsDistr', '', 1000, 500)
cSoverBvsDistr.Divide(2, 1)
axisPairSoverB = {}
for iPad, corrName in enumerate(['Particle0_Particle2_plus_Particle1_Particle3',
                                 'Particle0_Particle3_plus_Particle1_Particle2']):

    maxSoverB = hSoverB.GetMaximum()
    maxSEPair = hSEPairVsPt[corrName][0].GetMaximum()

    hFrameSoverB = cSoverBvsDistr.cd(iPad+1).DrawFrame(0., 0., ptMax, maxSoverB,
                                                       ';#it{p}_{T}(D^{#pm}) (GeV/#it{c});S/B')
    hFrameSoverB.GetYaxis().SetDecimals()
    hSoverB.Draw('same')
    hPair = hSEPairVsPt[corrName][0].Clone()
    hPair.Scale(maxSoverB / maxSEPair)
    hPair.DrawCopy('esame')
    kStar, avSoverB, SoverBavPt = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    gAverageSoverB[corrName].GetPoint(0, kStar, avSoverB)
    gSoverBAvPt[corrName].GetPoint(0, kStar, SoverBavPt)
    info.DrawLatex(0.2, 0.86, corrTitles[corrName])
    info.DrawLatex(0.2, 0.80, '#it{k}* < 200 MeV/#it{c}')
    info.DrawLatex(0.2, 0.74, f'#LT S/B #GT = '
                   f'{avSoverB.value:0.2f} #pm {gAverageSoverB[corrName].GetErrorYlow(0):0.2f}')
    info.DrawLatex(0.2, 0.68, 'S/B (#LT #it{p}_{T}(D^{#pm}) #GT) = '
                   f'{SoverBavPt.value:0.2f} #pm {gSoverBAvPt[corrName].GetErrorYlow(0):0.2f}')
    cSoverBvsDistr.cd(iPad+1).SetRightMargin(0.14)
    cSoverBvsDistr.cd(iPad+1).SetTicky(0)
    cSoverBvsDistr.cd(iPad+1).Modified()
    cSoverBvsDistr.cd(iPad+1).Update()
    axisPairSoverB[corrName] = TGaxis(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),
                                0., maxSEPair, 510, "+L")
    axisPairSoverB[corrName].SetLineColor(kRed+1)
    axisPairSoverB[corrName].SetLabelColor(kRed+1)
    axisPairSoverB[corrName].SetLabelFont(42)
    axisPairSoverB[corrName].SetLabelSize(0.045)
    axisPairSoverB[corrName].SetTitle('SE pairs')
    axisPairSoverB[corrName].SetTitleOffset(1.4)
    axisPairSoverB[corrName].SetLabelOffset(0.012)
    axisPairSoverB[corrName].SetTitleColor(kRed+1)
    axisPairSoverB[corrName].SetTitleFont(42)
    axisPairSoverB[corrName].SetTitleSize(0.05)
    axisPairSoverB[corrName].SetMaxDigits(3)
    axisPairSoverB[corrName].Draw()
cSoverBvsDistr.Modified()
cSoverBvsDistr.Update()

cPurityvsDistr = TCanvas('cPurityvsDistr', '', 1000, 500)
cPurityvsDistr.Divide(2, 1)
axisPairPurity = {}
for iPad, corrName in enumerate(['Particle0_Particle2_plus_Particle1_Particle3',
                                 'Particle0_Particle3_plus_Particle1_Particle2']):

    maxPurity = 2.
    maxSEPair = hSEPairVsPt[corrName][0].GetMaximum()

    hFramePurity = cPurityvsDistr.cd(iPad+1).DrawFrame(0., 0., ptMax, maxPurity,
                                                       ';#it{p}_{T}(D^{#pm}) (GeV/#it{c});S/(S+B)')
    hFramePurity.GetYaxis().SetDecimals()
    hPurity.Draw('same')
    hPair = hSEPairVsPt[corrName][0].Clone()
    hPair.Scale(maxPurity / maxSEPair)
    hPair.DrawCopy('esame')
    kStar, avPurity, PurityavPt = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
    gAveragePurity[corrName].GetPoint(0, kStar, avPurity)
    gPurityAvPt[corrName].GetPoint(0, kStar, PurityavPt)
    info.DrawLatex(0.2, 0.86, corrTitles[corrName])
    info.DrawLatex(0.2, 0.80, '#it{k}* < 200 MeV/#it{c}')
    info.DrawLatex(0.2, 0.74, f'#LT S/(S+B) #GT = '
                   f'{avPurity.value:0.2f} #pm {gAveragePurity[corrName].GetErrorYlow(0):0.2f}')
    info.DrawLatex(0.2, 0.68, 'S/(S+B) (#LT #it{p}_{T}(D^{#pm}) #GT) = '
                   f'{PurityavPt.value:0.2f} #pm {gPurityAvPt[corrName].GetErrorYlow(0):0.2f}')
    cPurityvsDistr.cd(iPad+1).SetRightMargin(0.14)
    cPurityvsDistr.cd(iPad+1).SetTicky(0)
    cPurityvsDistr.cd(iPad+1).Modified()
    cPurityvsDistr.cd(iPad+1).Update()
    axisPairPurity[corrName] = TGaxis(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(),
                                0., maxSEPair, 510, "+L")
    axisPairPurity[corrName].SetLineColor(kRed+1)
    axisPairPurity[corrName].SetLabelColor(kRed+1)
    axisPairPurity[corrName].SetLabelFont(42)
    axisPairPurity[corrName].SetLabelSize(0.045)
    axisPairPurity[corrName].SetTitle('SE pairs')
    axisPairPurity[corrName].SetTitleOffset(1.4)
    axisPairPurity[corrName].SetLabelOffset(0.012)
    axisPairPurity[corrName].SetTitleColor(kRed+1)
    axisPairPurity[corrName].SetTitleFont(42)
    axisPairPurity[corrName].SetTitleSize(0.05)
    axisPairPurity[corrName].SetMaxDigits(3)
    axisPairPurity[corrName].Draw()
cPurityvsDistr.Modified()
cPurityvsDistr.Update()

outFileName = 'Purity_vs_kstar' + outSuffix + '.root'
outFileName = os.path.join(outDirName, outFileName)
outFile = TFile(outFileName, 'recreate')
if addTDirectory:
    outDir = TDirectoryFile(outSuffix, outSuffix)
    outDir.Write()
    outDir.cd()
cSEDistr.Write()
cSoverBvsDistr.Write()
cSoverBvsKstar.Write()
cPurityvsDistr.Write()
cPurityvsKstar.Write()
for corrName in gAverageSoverB:
    gAverageSoverB[corrName].Write()
    gSoverBAvPt[corrName].Write()
    gAveragePurity[corrName].Write()
    gPurityAvPt[corrName].Write()
if addTDirectory:
    outDir.Close()
outFile.Close()
print('Output file saved: ', outFileName)

outFileNamePDF = outFileName.replace('.root', '.pdf')
cSEDistr.SaveAs(outFileNamePDF.replace('Purity_vs_kstar', 'PairDistr'))
cSoverBvsDistr.SaveAs(outFileNamePDF.replace('Purity_vs_kstar', 'SoverB_vs_PairDistr_kstar200'))
cPurityvsDistr.SaveAs(outFileNamePDF.replace('Purity_vs_kstar', 'Purity_vs_PairDistr_kstar200'))
cSoverBvsKstar.SaveAs(outFileNamePDF.replace('Purity', 'SoverB'))
cPurityvsKstar.SaveAs(outFileNamePDF)

input('Press enter to exit')
