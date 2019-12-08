'''
python script for the expected significance of D+ and Ds+ mesons using THnSparses
and manage sparses for v2 analyses with ML
run: python GetExpectedSignificance.py cfgFileName.yml cutSetFile.yml outFileName
with option --suffix a suffix is added to the output file name, otherwise the input file is overwritten
'''

import sys
import array
import math
import six
import yaml
from ROOT import TFile, TCanvas, TH1F, TSpline3, TF1, TGraphAsymmErrors, TLine, TLegend  # pylint: disable=import-error, no-name-in-module
from ROOT import kWhite, kBlack, kOrange, kRed, kGreen, kBlue, kAzure, kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp  # pylint: disable=import-error, no-name-in-module, line-too-long
from utils.ReadModel import ReadFONLL, ReadTAMU, ReadPHSD, ReadGossiaux, ReadCatania
from utils.ReadHepData import ReadHepDataROOT
from utils.StyleFormatter import SetGlobalStyle

#pylint: disable=redefined-outer-name
def GetExpectedBackgroundFromSB(hSB, mean, sigma, Nexp, Nanal):
    fMassBkg = TF1('fMassBkg', 'expo', 1.8, 2.12)
    hSB.Fit('fMassBkg', 'QR0')
    B = fMassBkg.Integral(mean-3*sigma, mean+3*sigma) / \
        hSB.GetBinWidth(1) * Nexp / Nanal
    return B


def GetExpectedBackgroundFromMC(hBkgMC, mean, sigma, Nexp, Nanal):
    binmin = hBkgMC.GetXaxis().FindBin(mean-3*sigma)
    binmax = hBkgMC.GetXaxis().FindBin(mean+3*sigma)
    B = hBkgMC.Integral(binmin, binmax) / hBkgMC.GetBinWidth(1) * Nexp / Nanal
    return B


def GetExpectedSignal(DeltaPt, sigma, Raa, Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp):
    corrYield = sigma * fractoD * Raa * Taa / DeltaPt
    rawYield = 2 * DeltaPt * corrYield * effPrompt * Acc * BR * Nexp / fprompt
    return rawYield, corrYield


def GetSideBandHisto(hMassData, mean, sigma, name):
    hMassDataReb = hMassData.Clone('hSB')
    hMassDataReb.Rebin(5)
    hSB = hMassDataReb.Clone(name)
    massbinmin = hSB.GetXaxis().FindBin(mean-4*sigma)
    massbinmax = hSB.GetXaxis().FindBin(mean+4*sigma)
    for iMassBin in range(1, hSB.GetNbinsX()):
        if iMassBin >= massbinmin and iMassBin <= massbinmax:
            hSB.SetBinContent(iMassBin, 0.)
            hSB.SetBinError(iMassBin, 0.)
        else:
            hSB.SetBinContent(iMassBin, hMassDataReb.GetBinContent(iMassBin))
            hSB.SetBinError(iMassBin, hMassDataReb.GetBinError(iMassBin))
    return hSB


def ComputeExpectedFprompt(ptmin, ptmax, effPrompt, hPredPrompt, hPredPromptMin, hPredPromptMax, effFD, hPredFD, hPredFDMin, hPredFDMax, RatioRaaFDPrompt):
    ptbinmin = hPredPrompt.GetXaxis().FindBin(ptmin*1.001)
    ptbinmax = hPredFD.GetXaxis().FindBin(ptmax*0.999)
    predPrompt, predPromptMin, predPromptMax, predFD, predFDMin, predFDMax = (
        0 for iPred in range(6))
    for iPt in range(ptbinmin, ptbinmax+1):
        predPrompt += hPredPrompt.GetBinContent(iPt) * \
            hPredPrompt.GetBinWidth(iPt)
        predPromptMin += hPredPromptMin.GetBinContent(
            iPt) * hPredPromptMin.GetBinWidth(iPt)
        predPromptMax += hPredPromptMax.GetBinContent(
            iPt) * hPredPromptMax.GetBinWidth(iPt)
        predFD += hPredFD.GetBinContent(iPt) * hPredFD.GetBinWidth(iPt)
        predFDMin += hPredFDMin.GetBinContent(iPt) * \
            hPredFDMin.GetBinWidth(iPt)
        predFDMax += hPredFDMax.GetBinContent(iPt) * \
            hPredFDMax.GetBinWidth(iPt)
    predPrompt /= (ptmax-ptmin)
    predPromptMin /= (ptmax-ptmin)
    predPromptMax /= (ptmax-ptmin)
    predFD /= (ptmax-ptmin)
    predFDMin /= (ptmax-ptmin)
    predFDMax /= (ptmax-ptmin)
    PromptSignal = predPrompt * effPrompt
    PromptSignalMin = predPromptMin * effPrompt
    PromptSignalMax = predPromptMax * effPrompt
    FDSignal = predFD * effFD
    FDSignalMin = predFDMin * effFD
    FDSignalMax = predFDMax * effFD
    fprompt = 1. / (1 + RatioRaaFDPrompt * FDSignal / PromptSignal)
    fpromptmin = 1. / (1 + RatioRaaFDPrompt * FDSignalMax / PromptSignalMin)
    fpromptmax = 1. / (1 + RatioRaaFDPrompt * FDSignalMin / PromptSignalMax)
    return fprompt, fpromptmin, fpromptmax


def SetGraphStyle(graph, color, colorfill, marker, linewidth=2):
    if colorfill == kWhite:
        graph.SetFillStyle(0)
    else:
        graph.SetFillColorAlpha(colorfill, 0.25)
    graph.SetMarkerColor(color)
    graph.SetLineColor(color)
    graph.SetLineWidth(linewidth)
    graph.SetMarkerStyle(marker)


def SetSplineStyle(spline, color):
    spline.SetLineColor(color)
    spline.SetLineWidth(2)


colors = {'TAMU': kOrange+10, 'PHSD': kGreen +
          3, 'Gossiaux': kAzure+4, 'Catania': kRed+2}
colorsFill = {'TAMU': kOrange+7, 'PHSD': kGreen +
              2, 'Gossiaux': kAzure+2, 'Catania': kRed}
markers = {'TAMU': kFullSquare, 'PHSD': kFullCircle,
           'Gossiaux': kFullDiamond, 'Catania': kFullTriangleUp}

####################################################################################################
# main function

inputCfg = sys.argv[1]
cutSetFile = sys.argv[2]
outFileName = sys.argv[3]

with open(inputCfg, 'r') as ymlinputCfg:
    if six.PY2:
        inputCfg = yaml.load(ymlinputCfg)
    else:
        inputCfg = yaml.safe_load(ymlinputCfg)

if inputCfg['bkgConfiguration']['getbkgfromMC'] is False:
    infileDataLowPt = TFile(inputCfg['fileDataLowPt']['filename'])
    indirDataLowPt = infileDataLowPt.Get(inputCfg['fileDataLowPt']['dirname'])
    inlistDataLowPt = indirDataLowPt.Get(inputCfg['fileDataLowPt']['listname'])
    sMassPtCutVarsLowPt = inlistDataLowPt.FindObject(
        inputCfg['fileDataLowPt']['sparsenameAll'])
    hEvLowPt = inlistDataLowPt.FindObject(
        inputCfg['fileDataLowPt']['histoevname'])

    infileDataHighPt = TFile(inputCfg['fileDataHighPt']['filename'])
    indirDataHighPt = infileDataHighPt.Get(
        inputCfg['fileDataHighPt']['dirname'])
    inlistDataHighPt = indirDataHighPt.Get(
        inputCfg['fileDataHighPt']['listname'])
    sMassPtCutVarsHighPt = inlistDataHighPt.FindObject(
        inputCfg['fileDataHighPt']['sparsenameAll'])
    hEvHighPt = inlistDataHighPt.FindObject(
        inputCfg['fileDataHighPt']['histoevname'])

infileMCLowPt = TFile(inputCfg['fileMCLowPt']['filename'])
indirMCLowPt = infileMCLowPt.Get(inputCfg['fileMCLowPt']['dirname'])
inlistMCLowPt = indirMCLowPt.Get(inputCfg['fileMCLowPt']['listname'])
sMassPtCutVarsPromptLowPt = inlistMCLowPt.FindObject(
    inputCfg['fileMCLowPt']['sparsenamePrompt'])
sMassPtCutVarsFDLowPt = inlistMCLowPt.FindObject(
    inputCfg['fileMCLowPt']['sparsenameFD'])
sGenPromptLowPt = inlistMCLowPt.FindObject(
    inputCfg['fileMCLowPt']['sparsenameGenPrompt'])
sGenFDLowPt = inlistMCLowPt.FindObject(
    inputCfg['fileMCLowPt']['sparsenameGenFD'])
if inputCfg['bkgConfiguration']['getbkgfromMC']:
    sMassPtCutVarsBkgLowPt = inlistMCLowPt.FindObject(
        inputCfg['fileMCLowPt']['sparsenameBkg'])
    hEvLowPt = inlistMCLowPt.FindObject(inputCfg['fileMCLowPt']['histoevname'])

infileMCHighPt = TFile(inputCfg['fileMCHighPt']['filename'])
indirMCHighPt = infileMCHighPt.Get(inputCfg['fileMCHighPt']['dirname'])
inlistMCHighPt = indirMCHighPt.Get(inputCfg['fileMCHighPt']['listname'])
sMassPtCutVarsPromptHighPt = inlistMCHighPt.FindObject(
    inputCfg['fileMCHighPt']['sparsenamePrompt'])
sMassPtCutVarsFDHighPt = inlistMCHighPt.FindObject(
    inputCfg['fileMCHighPt']['sparsenameFD'])
sGenPromptHighPt = inlistMCHighPt.FindObject(
    inputCfg['fileMCHighPt']['sparsenameGenPrompt'])
sGenFDHighPt = inlistMCLowPt.FindObject(
    inputCfg['fileMCLowPt']['sparsenameGenFD'])
if inputCfg['bkgConfiguration']['getbkgfromMC']:
    sMassPtCutVarsBkgHighPt = inlistMCLowPt.FindObject(
        inputCfg['fileMCHighPt']['sparsenameBkg'])
    hEvHighPt = inlistMCLowPt.FindObject(
        inputCfg['fileMCHighPt']['histoevname'])

if inputCfg['bkgConfiguration']['applyCorrFactor']:
    infileCorrFactor = TFile.Open(
        inputCfg['bkgConfiguration']['bkgCorrFactorfilename'])
    hBkgCorrFactor = infileCorrFactor.Get(
        inputCfg['bkgConfiguration']['bkgCorrFactorhistoname'])

PtThreshold = inputCfg['PtThreshold']

Nexp = inputCfg['nExpectedEvents']
Taa = inputCfg['Taa']
BR = inputCfg['BR']
fractoD = inputCfg['fractoD']

if inputCfg['PredForFprompt']['estimateFprompt']:
    infilePred = TFile.Open(inputCfg['PredForFprompt']['filename'])
    hPredPrompt = infilePred.Get(inputCfg['PredForFprompt']['histonamePrompt'])
    hPredPromptMin = infilePred.Get(
        inputCfg['PredForFprompt']['histonamePrompt'].replace('central', 'min'))
    hPredPromptMax = infilePred.Get(
        inputCfg['PredForFprompt']['histonamePrompt'].replace('central', 'max'))
    hPredFD = infilePred.Get(inputCfg['PredForFprompt']['histonameFD'])
    hPredFDMin = infilePred.Get(
        inputCfg['PredForFprompt']['histonameFD'].replace('central', 'min'))
    hPredFDMax = infilePred.Get(
        inputCfg['PredForFprompt']['histonameFD'].replace('central', 'max'))
    RatioRaaFDPrompt = inputCfg['PredForFprompt']['RatioRaaFDPrompt']

infileAcc = TFile(inputCfg['filenameAcc'])
hPtGenAcc = infileAcc.Get('hPtGenAcc')
hPtGenLimAcc = infileAcc.Get('hPtGenLimAcc')

FONLL = ReadFONLL(inputCfg['filenameFONLL'])
TAMU = ReadTAMU(inputCfg['filenameRaaPredTAMU'])
PHSD = ReadPHSD(inputCfg['filenameRaaPredPHSD'])
Gossiaux = ReadGossiaux(inputCfg['filenameRaaPredGossiaux'])
Catania = ReadCatania(inputCfg['filenameRaaPredCatania'])

outfile = TFile('%s.root' % outFileName, 'recreate')

sTAMU, sGossiaux, sPHSD, sCatania = ({} for iDic in range(4))
for iRaaTAMU in TAMU:
    if iRaaTAMU != 'PtCent':
        sTAMU[iRaaTAMU] = TSpline3('sTAMU%s' % iRaaTAMU, array.array(
            'd', TAMU['PtCent']), array.array('d', TAMU[iRaaTAMU]), len(TAMU['PtCent']))
for iRaaPHSD in PHSD:
    if iRaaPHSD != 'PtCent':
        sPHSD[iRaaPHSD] = TSpline3('sPHSD%s' % iRaaPHSD, array.array(
            'd', PHSD['PtCent']), array.array('d', PHSD[iRaaPHSD]), len(PHSD['PtCent']))
for iRaaGossiaux in Gossiaux:
    if iRaaGossiaux != 'PtCent':
        sGossiaux[iRaaGossiaux] = TSpline3('sGossiaux%s' % iRaaGossiaux, array.array(
            'd', Gossiaux['PtCent']), array.array('d', Gossiaux[iRaaGossiaux]), len(Gossiaux['PtCent']))
for iRaaCatania in Catania:
    if iRaaCatania != 'PtCent':
        sCatania[iRaaCatania] = TSpline3('sCatania%s' % iRaaCatania, array.array(
            'd', Catania['PtCent']), array.array('d', Catania[iRaaCatania]), len(Catania['PtCent']))

with open(cutSetFile, 'r') as ymlCutSetFile:
    if six.PY2:
        cutSetFile = yaml.load(ymlCutSetFile)
    else:
        cutSetFile = yaml.safe_load(ymlCutSetFile)
cutVars = cutSetFile['cutvars']

cutVars['Pt']['limits'] = list(cutVars['Pt']['min'])
cutVars['Pt']['limits'].append(
    cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1])

names = ['TAMU', 'PHSD', 'Gossiaux', 'Catania']
gSignifTAMU, gSignifPHSD, gSignifGossiaux, gSignifCatania = (
    TGraphAsymmErrors(0) for iGraph in range(4))
hSignifTAMU, hSignifPHSD, hSignifGossiaux, hSignifCatania = (TH1F('hSignificance%s' % names[iGraph], ';#it{p}_{T} (GeV/#it{c});significance', len(
    cutVars['Pt']['limits'])-1, array.array('f', cutVars['Pt']['limits'])) for iGraph in range(4))
gCorrYieldTAMU, gCorrYieldPHSD, gCorrYieldGossiaux, gCorrYieldCatania = (
    TGraphAsymmErrors(0) for iGraph in range(4))
gRawYieldTAMU, gRawYieldPHSD, gRawYieldGossiaux, gRawYieldCatania = (
    TGraphAsymmErrors(0) for iGraph in range(4))
hRawYieldTAMU, hRawYieldPHSD, hRawYieldGossiaux, hRawYieldCatania = (TH1F('hRawYield%s' % names[iGraph], ';#it{p}_{T} (GeV/#it{c});raw yield', len(
    cutVars['Pt']['limits'])-1, array.array('f', cutVars['Pt']['limits'])) for iGraph in range(4))

gSignifTAMU.SetTitle(';#it{p}_{T} (GeV/#it{c});significance')
gSignifPHSD.SetTitle(';#it{p}_{T} (GeV/#it{c});significance')
gSignifGossiaux.SetTitle(';#it{p}_{T} (GeV/#it{c});significance')
gSignifCatania.SetTitle(';#it{p}_{T} (GeV/#it{c});significance')

gCorrYieldTAMU.SetTitle(';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}')
gCorrYieldPHSD.SetTitle(';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}')
gCorrYieldGossiaux.SetTitle(';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}')
gCorrYieldCatania.SetTitle(';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}')

gRawYieldTAMU.SetTitle(';#it{p}_{T} (GeV/#it{c}); raw yield')
gRawYieldPHSD.SetTitle(';#it{p}_{T} (GeV/#it{c}); raw yield')
gRawYieldGossiaux.SetTitle(';#it{p}_{T} (GeV/#it{c}); raw yield')
gRawYieldCatania.SetTitle(';#it{p}_{T} (GeV/#it{c}); raw yield')

hBackground = TH1F('hBackground', ';#it{p}_{T} (GeV/#it{c});B(3#sigma)', len(
    cutVars['Pt']['limits'])-1, array.array('f', cutVars['Pt']['limits']))
SetGraphStyle(hBackground, kBlack, kWhite, kFullSquare)

gFprompt = TGraphAsymmErrors(0)
gFprompt.SetTitle(';#it{p}_{T} (GeV/#it{c});#it{f}_{prompt} (#it{f}_{c})')
SetGraphStyle(gFprompt, kBlack, kWhite, kFullSquare)

hAccEffPrompt = TH1F('hAccEffPrompt', ';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon) #times 2#it{y}_{fid}', len(
    cutVars['Pt']['limits'])-1, array.array('f', cutVars['Pt']['limits']))
hAccEffFD = TH1F('hAccEffFD', ';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon) #times 2#it{y}_{fid}', len(
    cutVars['Pt']['limits'])-1, array.array('f', cutVars['Pt']['limits']))
SetGraphStyle(hAccEffPrompt, kRed+1, kWhite, kFullSquare)
SetGraphStyle(hAccEffFD, kBlue+1, kWhite, kFullCircle)

for iPt, ptmin in enumerate(cutVars['Pt']['min']):

    ptmax = cutVars['Pt']['max'][iPt]
    PtCent = (ptmin+ptmax) / 2
    PtUnc = (ptmax-ptmin) / 2

    # check if low or high pt
    if ptmax <= PtThreshold:
        nEvBkg = hEvLowPt.GetBinContent(5)
        sMassPtCutVarsPrompt = sMassPtCutVarsPromptLowPt.Clone(
            'sMassPtCutVarsPrompt')
        sMassPtCutVarsFD = sMassPtCutVarsFDLowPt.Clone('sMassPtCutVarsFD')
        sGenPrompt = sGenPromptLowPt.Clone('sGenPrompt')
        sGenFD = sGenFDLowPt.Clone('sGenFD')
        if inputCfg['bkgConfiguration']['getbkgfromMC']:
            sMassPtCutVarsBkg = sMassPtCutVarsBkgLowPt.Clone(
                'sMassPtCutVarsBkg')
        else:
            sMassPtCutVars = sMassPtCutVarsLowPt.Clone('sMassPtCutVars')
    else:
        nEvBkg = hEvHighPt.GetBinContent(5)
        sMassPtCutVarsPrompt = sMassPtCutVarsPromptHighPt.Clone(
            'sMassPtCutVarsPrompt')
        sMassPtCutVarsFD = sMassPtCutVarsFDHighPt.Clone('sMassPtCutVarsFD')
        sGenPrompt = sGenPromptHighPt.Clone('sGenPrompt')
        sGenFD = sGenFDHighPt.Clone('sGenFD')
        if inputCfg['bkgConfiguration']['getbkgfromMC']:
            sMassPtCutVarsBkg = sMassPtCutVarsBkgHighPt.Clone(
                'sMassPtCutVarsBkg')
        else:
            sMassPtCutVars = sMassPtCutVarsHighPt.Clone('sMassPtCutVars')

    # gen for efficiency
    binGenPtMin = sGenPrompt.GetAxis(0).FindBin(
        cutVars['Pt']['min'][iPt]*1.0001)
    binGenPtMax = sGenPrompt.GetAxis(0).FindBin(
        cutVars['Pt']['max'][iPt]*0.9999)
    hPtGenPrompt = sGenPrompt.GetAxis(0).SetRange(binGenPtMin, binGenPtMax)
    hPtGenPrompt = sGenPrompt.Projection(0)
    hPtGenPrompt.SetName('hPtPrompt_pT_%0.f_%0.f' % (
        cutVars['Pt']['min'][iPt], cutVars['Pt']['max'][iPt]))
    nGenPrompt = hPtGenPrompt.Integral()
    hPtGenFD = sGenFD.GetAxis(0).SetRange(binGenPtMin, binGenPtMax)
    hPtGenFD = sGenFD.Projection(0)
    hPtGenFD.SetName('hPtFD_pT_%0.f_%0.f' %
                     (cutVars['Pt']['min'][iPt], cutVars['Pt']['max'][iPt]))
    nGenFD = hPtGenFD.Integral()

    # acceptance from TOY MC
    binAccPtMin = hPtGenAcc.GetXaxis().FindBin(
        cutVars['Pt']['min'][iPt]*1.0001)
    binAccPtMax = hPtGenAcc.GetXaxis().FindBin(
        cutVars['Pt']['max'][iPt]*0.9999)
    Acc = hPtGenAcc.Integral(binAccPtMin, binAccPtMax) / \
        hPtGenLimAcc.Integral(binAccPtMin, binAccPtMax)

    # FONLL
    binFONLL = FONLL['PtMin'].index(cutVars['Pt']['min'][iPt])
    sigmaFONLL = []
    sigmaFONLL.append(FONLL['Cent'][binFONLL])
    sigmaFONLL.append(FONLL['Min'][binFONLL])
    sigmaFONLL.append(FONLL['Max'][binFONLL])

    # RAA
    RaaTAMU, RaaPHSD, RaaGossiaux, RaaCatania = ({} for iDic in range(4))
    for iRaaTAMU in sTAMU:
        RaaTAMU[iRaaTAMU] = sTAMU[iRaaTAMU].Eval(
            (cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2)
    for iRaaPHSD in sPHSD:
        RaaPHSD[iRaaPHSD] = sPHSD[iRaaPHSD].Eval(
            (cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2)
    for iRaaGossiaux in sGossiaux:
        RaaGossiaux[iRaaGossiaux] = sGossiaux[iRaaGossiaux].Eval(
            (cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2)
    for iRaaCatania in sCatania:
        RaaCatania[iRaaCatania] = sCatania[iRaaCatania].Eval(
            (cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2)

    for iVar in cutVars:
        if iVar == 'InvMass':
            continue
        binMin = sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).FindBin(
            cutVars[iVar]['min'][iPt]*1.0001)
        binMax = sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).FindBin(
            cutVars[iVar]['max'][iPt]*0.9999)
        sMassPtCutVarsPrompt.GetAxis(
            cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
        sMassPtCutVarsFD.GetAxis(
            cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
        if inputCfg['bkgConfiguration']['getbkgfromMC']:
            sMassPtCutVarsBkg.GetAxis(
                cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
        else:
            sMassPtCutVars.GetAxis(
                cutVars[iVar]['axisnum']).SetRange(binMin, binMax)

    if inputCfg['bkgConfiguration']['getbkgfromMC']:
        hMassBkg = sMassPtCutVarsBkg.Projection(0)
        if inputCfg['bkgConfiguration']['applyCorrFactor']:
            binCorrFactor = hBkgCorrFactor.GetXaxis().FindBin(PtCent)
            corrfactor = hBkgCorrFactor.GetBinContent(binCorrFactor)
            hMassBkg.Scale(corrfactor)
        hMassBkg.SetName('hMassBkgpT_%0.f_%0.f' %
                         (cutVars['Pt']['min'][iPt], cutVars['Pt']['max'][iPt]))
    else:
        hMassData = sMassPtCutVars.Projection(0)
        if inputCfg['bkgConfiguration']['applyCorrFactor']:
            binCorrFactor = hBkgCorrFactor.GetXaxis().FindBin(PtCent)
            corrfactor = hBkgCorrFactor.GetBinContent(binCorrFactor)
            hMassData.Scale(corrfactor)
        hMassData.SetName('hMassData_pT_%0.f_%0.f' % (
            cutVars['Pt']['min'][iPt], cutVars['Pt']['max'][iPt]))

    hMassPrompt = sMassPtCutVarsPrompt.Projection(0)
    hMassPrompt.SetName('hMassPrompt_pT_%0.f_%0.f' %
                        (cutVars['Pt']['min'][iPt], cutVars['Pt']['max'][iPt]))
    hMassFD = sMassPtCutVarsFD.Projection(0)
    hMassFD.SetName('hMassFD_pT_%0.f_%0.f' %
                    (cutVars['Pt']['min'][iPt], cutVars['Pt']['max'][iPt]))
    hMassSgn = hMassPrompt.Clone('hMassSgn_pT_%0.f_%0.f' % (
        cutVars['Pt']['min'][iPt], cutVars['Pt']['max'][iPt]))
    hMassSgn.Add(hMassFD)

    nRecoPrompt = hMassPrompt.Integral()
    nRecoFD = hMassFD.Integral()
    effPrompt = nRecoPrompt / nGenPrompt
    effFD = nRecoFD / nGenFD
    effPromptUnc = math.sqrt(effPrompt*(1-effPrompt)/nGenPrompt)
    effFDUnc = math.sqrt(effFD*(1-effFD)/nGenFD)

    hAccEffPrompt.SetBinContent(iPt+1, effPrompt*Acc)
    hAccEffPrompt.SetBinError(iPt+1, effPromptUnc*Acc)
    hAccEffFD.SetBinContent(iPt+1, effFD*Acc)
    hAccEffFD.SetBinError(iPt+1, effFDUnc*Acc)

    fMassSgn = TF1('fMassSgn', 'gaus', 1.7, 2.15)
    hMassSgn.Fit('fMassSgn', 'Q0')
    mean = fMassSgn.GetParameter(1)
    sigma = fMassSgn.GetParameter(2)

    if inputCfg['bkgConfiguration']['getbkgfromMC']:
        B = GetExpectedBackgroundFromMC(hMassBkg, mean, sigma, Nexp, nEvBkg)
    else:
        hMassSB = GetSideBandHisto(hMassData, mean, sigma, 'hMassSB_pT_%0.f_%0.f' % (
            cutVars['Pt']['min'][iPt], cutVars['Pt']['max'][iPt]))
        B = GetExpectedBackgroundFromSB(hMassSB, mean, sigma, Nexp, nEvBkg)

    if inputCfg['PredForFprompt']['estimateFprompt']:
        fprompt, fpromptmin, fpromptmax = ComputeExpectedFprompt(
            ptmin, ptmax, effPrompt, hPredPrompt, hPredPromptMin, hPredPromptMax, effFD, hPredFD, hPredFDMin, hPredFDMax, RatioRaaFDPrompt)
    else:
        fprompt, fpromptmin, fpromptmax = (
            inputCfg['fprompt'] for iF in range(3))

    gFprompt.SetPoint(iPt, PtCent, fprompt)
    gFprompt.SetPointError(iPt, PtUnc, PtUnc, fprompt -
                           fpromptmin, fpromptmax-fprompt)

    SignifTAMU, SignifPHSD, SignifGossiaux, SignifCatania = (
        {} for iDic in range(4))
    SignifUncTAMU, SignifUncPHSD, SignifUncGossiaux, SignifUncCatania = (
        {} for iDic in range(4))
    CorrYieldTAMU, CorrYieldPHSD, CorrYieldGossiaux, CorrYieldCatania = (
        {} for iDic in range(4))
    RawYieldTAMU, RawYieldPHSD, RawYieldGossiaux, RawYieldCatania = (
        {} for iDic in range(4))
    for iFONLL in FONLL:
        if iFONLL == 'Cent' or iFONLL == 'Min' or iFONLL == 'Max':

            SignifTAMU[iFONLL], SignifUncTAMU[iFONLL], CorrYieldTAMU[iFONLL], RawYieldTAMU[iFONLL] = (
                dict(RaaTAMU) for iDic in range(4))
            SignifPHSD[iFONLL], SignifUncPHSD[iFONLL], CorrYieldPHSD[iFONLL], RawYieldPHSD[iFONLL] = (
                dict(RaaPHSD) for iDic in range(4))
            SignifGossiaux[iFONLL], SignifUncGossiaux[iFONLL], CorrYieldGossiaux[iFONLL], RawYieldGossiaux[iFONLL] = (
                dict(RaaGossiaux) for iDic in range(4))
            SignifCatania[iFONLL], SignifUncCatania[iFONLL], CorrYieldCatania[iFONLL], RawYieldCatania[iFONLL] = (
                dict(RaaCatania) for iDic in range(4))

            for iRaaTAMU in RaaTAMU:
                rawyield, corryield = GetExpectedSignal(cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt],
                                                        FONLL[iFONLL][binFONLL], RaaTAMU[iRaaTAMU], Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp)
                SignifTAMU[iFONLL][iRaaTAMU] = rawyield/math.sqrt(rawyield+B)
                SignifUncTAMU[iFONLL][iRaaTAMU] = rawyield / \
                    math.sqrt(rawyield+B)
                CorrYieldTAMU[iFONLL][iRaaTAMU] = corryield
                RawYieldTAMU[iFONLL][iRaaTAMU] = rawyield
            for iRaaPHSD in RaaPHSD:
                rawyield, corryield = GetExpectedSignal(cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt],
                                                        FONLL[iFONLL][binFONLL], RaaPHSD[iRaaPHSD], Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp)
                SignifPHSD[iFONLL][iRaaPHSD] = rawyield/math.sqrt(rawyield+B)
                SignifUncPHSD[iFONLL][iRaaPHSD] = rawyield / \
                    math.sqrt(rawyield+B)
                CorrYieldPHSD[iFONLL][iRaaPHSD] = corryield
                RawYieldPHSD[iFONLL][iRaaPHSD] = rawyield
            for iRaaGossiaux in RaaGossiaux:
                rawyield, corryield = GetExpectedSignal(cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt], FONLL[iFONLL]
                                                        [binFONLL], RaaGossiaux[iRaaGossiaux], Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp)
                SignifGossiaux[iFONLL][iRaaGossiaux] = rawyield / \
                    math.sqrt(rawyield+B)
                SignifUncGossiaux[iFONLL][iRaaGossiaux] = rawyield / \
                    math.sqrt(rawyield+B)
                CorrYieldGossiaux[iFONLL][iRaaGossiaux] = corryield
                RawYieldGossiaux[iFONLL][iRaaGossiaux] = rawyield
            for iRaaCatania in RaaCatania:
                rawyield, corryield = GetExpectedSignal(cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt],
                                                        FONLL[iFONLL][binFONLL], RaaCatania[iRaaCatania], Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp)
                SignifCatania[iFONLL][iRaaCatania] = rawyield / \
                    math.sqrt(rawyield+B)
                SignifUncCatania[iFONLL][iRaaCatania] = rawyield / \
                    math.sqrt(rawyield+B)
                CorrYieldCatania[iFONLL][iRaaCatania] = corryield
                RawYieldCatania[iFONLL][iRaaCatania] = rawyield

    gSignifTAMU.SetPoint(iPt, PtCent, SignifTAMU['Max']['Min'])
    gSignifTAMU.SetPointError(iPt, PtUnc, PtUnc, SignifTAMU['Max']['Min']-SignifTAMU['Min'][min(
        SignifTAMU['Min'], key=SignifTAMU['Min'].get)], SignifTAMU['Max'][max(SignifTAMU['Max'], key=SignifTAMU['Max'].get)]-SignifTAMU['Max']['Min'])
    gSignifPHSD.SetPoint(iPt, PtCent, SignifPHSD['Max']['Cent'])
    gSignifPHSD.SetPointError(iPt, PtUnc, PtUnc, SignifPHSD['Max']['Cent']-SignifPHSD['Min'][min(
        SignifPHSD['Min'], key=SignifPHSD['Min'].get)], SignifPHSD['Max'][max(SignifPHSD['Max'], key=SignifPHSD['Max'].get)]-SignifPHSD['Max']['Cent'])
    gSignifGossiaux.SetPoint(iPt, PtCent, SignifGossiaux['Max']['ColRad'])
    gSignifGossiaux.SetPointError(iPt, PtUnc, PtUnc, SignifGossiaux['Max']['ColRad']-SignifGossiaux['Min'][min(
        SignifGossiaux['Min'], key=SignifGossiaux['Min'].get)], SignifGossiaux['Max'][max(SignifGossiaux['Max'], key=SignifGossiaux['Max'].get)]-SignifGossiaux['Max']['ColRad'])
    gSignifCatania.SetPoint(iPt, PtCent, SignifCatania['Max']['Cent'])
    gSignifCatania.SetPointError(iPt, PtUnc, PtUnc, SignifCatania['Max']['Cent']-SignifCatania['Min'][min(
        SignifCatania['Min'], key=SignifCatania['Min'].get)], SignifCatania['Max'][max(SignifCatania['Max'], key=SignifCatania['Max'].get)]-SignifCatania['Max']['Cent'])

    gCorrYieldTAMU.SetPoint(iPt, PtCent, CorrYieldTAMU['Max']['Min'])
    gCorrYieldTAMU.SetPointError(iPt, PtUnc, PtUnc, CorrYieldTAMU['Max']['Min']-CorrYieldTAMU['Min'][min(
        CorrYieldTAMU['Min'], key=CorrYieldTAMU['Min'].get)], CorrYieldTAMU['Max'][max(CorrYieldTAMU['Max'], key=CorrYieldTAMU['Max'].get)]-CorrYieldTAMU['Max']['Min'])
    gCorrYieldPHSD.SetPoint(iPt, PtCent, CorrYieldPHSD['Max']['Cent'])
    gCorrYieldPHSD.SetPointError(iPt, PtUnc, PtUnc, CorrYieldPHSD['Max']['Cent']-CorrYieldPHSD['Min'][min(
        CorrYieldPHSD['Min'], key=CorrYieldPHSD['Min'].get)], CorrYieldPHSD['Max'][max(CorrYieldPHSD['Max'], key=CorrYieldPHSD['Max'].get)]-CorrYieldPHSD['Max']['Cent'])
    gCorrYieldGossiaux.SetPoint(
        iPt, PtCent, CorrYieldGossiaux['Max']['ColRad'])
    gCorrYieldGossiaux.SetPointError(iPt, PtUnc, PtUnc, CorrYieldGossiaux['Max']['ColRad']-CorrYieldGossiaux['Min'][min(
        CorrYieldGossiaux['Min'], key=CorrYieldGossiaux['Min'].get)], CorrYieldGossiaux['Max'][max(CorrYieldGossiaux['Max'], key=CorrYieldGossiaux['Max'].get)]-CorrYieldGossiaux['Max']['ColRad'])
    gCorrYieldCatania.SetPoint(iPt, PtCent, CorrYieldCatania['Max']['Cent'])
    gCorrYieldCatania.SetPointError(iPt, PtUnc, PtUnc, CorrYieldCatania['Max']['Cent']-CorrYieldCatania['Min'][min(
        CorrYieldCatania['Min'], key=CorrYieldCatania['Min'].get)], CorrYieldCatania['Max'][max(CorrYieldCatania['Max'], key=CorrYieldCatania['Max'].get)]-CorrYieldCatania['Max']['Cent'])

    gRawYieldTAMU.SetPoint(iPt, PtCent, RawYieldTAMU['Max']['Min'])
    gRawYieldTAMU.SetPointError(iPt, PtUnc, PtUnc, RawYieldTAMU['Max']['Min']-RawYieldTAMU['Min'][min(
        RawYieldTAMU['Min'], key=RawYieldTAMU['Min'].get)], RawYieldTAMU['Max'][max(RawYieldTAMU['Max'], key=RawYieldTAMU['Max'].get)]-RawYieldTAMU['Max']['Min'])
    gRawYieldPHSD.SetPoint(iPt, PtCent, RawYieldPHSD['Max']['Cent'])
    gRawYieldPHSD.SetPointError(iPt, PtUnc, PtUnc, RawYieldPHSD['Max']['Cent']-RawYieldPHSD['Min'][min(
        RawYieldPHSD['Min'], key=RawYieldPHSD['Min'].get)], RawYieldPHSD['Max'][max(RawYieldPHSD['Max'], key=RawYieldPHSD['Max'].get)]-RawYieldPHSD['Max']['Cent'])
    gRawYieldGossiaux.SetPoint(iPt, PtCent, RawYieldGossiaux['Max']['ColRad'])
    gRawYieldGossiaux.SetPointError(iPt, PtUnc, PtUnc, RawYieldGossiaux['Max']['ColRad']-RawYieldGossiaux['Min'][min(
        RawYieldGossiaux['Min'], key=RawYieldGossiaux['Min'].get)], RawYieldGossiaux['Max'][max(RawYieldGossiaux['Max'], key=RawYieldGossiaux['Max'].get)]-RawYieldGossiaux['Max']['ColRad'])
    gRawYieldCatania.SetPoint(iPt, PtCent, RawYieldCatania['Max']['Cent'])
    gRawYieldCatania.SetPointError(iPt, PtUnc, PtUnc, RawYieldCatania['Max']['Cent']-RawYieldCatania['Min'][min(
        RawYieldCatania['Min'], key=RawYieldCatania['Min'].get)], RawYieldCatania['Max'][max(RawYieldCatania['Max'], key=RawYieldCatania['Max'].get)]-RawYieldCatania['Max']['Cent'])

    sigTAMU = RawYieldTAMU['Max']['Min']
    sigPHSD = RawYieldPHSD['Max']['Cent']
    sigGossiaux = RawYieldGossiaux['Max']['ColRad']
    sigCatania = RawYieldCatania['Max']['Cent']

    siguncTAMU = math.sqrt(sigTAMU)
    siguncPHSD = math.sqrt(sigPHSD)
    siguncGossiaux = math.sqrt(sigGossiaux)
    siguncCatania = math.sqrt(sigCatania)

    hRawYieldTAMU.SetBinContent(iPt+1, sigTAMU)
    hRawYieldTAMU.SetBinError(iPt+1, siguncTAMU)
    hRawYieldPHSD.SetBinContent(iPt+1, sigPHSD)
    hRawYieldPHSD.SetBinError(iPt+1, siguncPHSD)
    hRawYieldGossiaux.SetBinContent(iPt+1, sigGossiaux)
    hRawYieldGossiaux.SetBinError(iPt+1, siguncGossiaux)
    hRawYieldCatania.SetBinContent(iPt+1, sigCatania)
    hRawYieldCatania.SetBinError(iPt+1, siguncCatania)

    Bunc = math.sqrt(B)
    hBackground.SetBinContent(iPt+1, B)
    hBackground.SetBinError(iPt+1, Bunc)

    hSignifTAMU.SetBinContent(iPt+1, SignifTAMU['Max']['Min'])
    hSignifTAMU.SetBinError(iPt+1, math.sqrt(SignifTAMU['Max']['Min'])*math.sqrt(
        (siguncTAMU**2+Bunc**2)/(4.*(sigTAMU+B)**2)+(B/(sigTAMU+B))*(siguncTAMU/sigTAMU)**2))
    hSignifPHSD.SetBinContent(iPt+1, SignifPHSD['Max']['Cent'])
    hSignifPHSD.SetBinError(iPt+1, math.sqrt(SignifPHSD['Max']['Cent'])*math.sqrt(
        (siguncPHSD**2+Bunc**2)/(4.*(sigPHSD+B)**2)+(B/(sigPHSD+B))*(siguncPHSD/sigPHSD)**2))
    hSignifGossiaux.SetBinContent(iPt+1, SignifGossiaux['Max']['ColRad'])
    hSignifGossiaux.SetBinError(iPt+1, math.sqrt(SignifGossiaux['Max']['ColRad'])*math.sqrt(
        (siguncGossiaux**2+Bunc**2)/(4.*(sigGossiaux+B)**2)+(B/(sigGossiaux+B))*(siguncGossiaux/sigGossiaux)**2))
    hSignifCatania.SetBinContent(iPt+1, SignifCatania['Max']['Cent'])
    hSignifCatania.SetBinError(iPt+1, math.sqrt(SignifCatania['Max']['Cent'])*math.sqrt(
        (siguncCatania**2+Bunc**2)/(4.*(sigCatania+B)**2)+(B/(sigCatania+B))*(siguncCatania/sigCatania)**2))

    hMassData.Write()
    hMassPrompt.Write()
    hMassFD.Write()
    hMassSB.Write()

# load published result
if inputCfg['PublishedResult']['plotPublished']:
    hDRaaStat, gDRaaSyst = ReadHepDataROOT(
        inputCfg['PublishedResult']['filenameHepData'], inputCfg['PublishedResult']['tableRaa'])
    hDCorrYieldStat, gDCorrYieldSyst = ReadHepDataROOT(
        inputCfg['PublishedResult']['filenameHepData'], inputCfg['PublishedResult']['tableCorrYield'])
    SetGraphStyle(gDRaaSyst, kBlack, kWhite, kFullSquare, 1)
    SetGraphStyle(gDCorrYieldSyst, kBlack, kWhite, kFullSquare, 1)
    SetGraphStyle(hDRaaStat, kBlack, kWhite, kFullSquare, 1)
    SetGraphStyle(hDCorrYieldStat, kBlack, kWhite, kFullSquare, 1)

# plots
SetGlobalStyle(padleftmargin=0.14, titlesize=0.045, labelsize=0.04)

legTAMU, legPHSD, legGossiaux, legCatania = ({'CorrYield': TLegend(0.45, 0.7, 0.65, 0.85), 'Raa': TLegend(
    0.45, 0.7, 0.65, 0.85), 'Signif': TLegend(0.45, 0.75, 0.65, 0.85)} for iDic in range(4))

for iLeg in legTAMU:
    legTAMU[iLeg].SetFillStyle(0)
    legTAMU[iLeg].SetTextSize(0.045)
    legPHSD[iLeg].SetFillStyle(0)
    legPHSD[iLeg].SetTextSize(0.045)
    legGossiaux[iLeg].SetFillStyle(0)
    legGossiaux[iLeg].SetTextSize(0.045)
    legCatania[iLeg].SetFillStyle(0)
    legCatania[iLeg].SetTextSize(0.045)
    if iLeg != 'Signif' and inputCfg['PublishedResult']['plotPublished']:
        legTAMU[iLeg].AddEntry(
            hDCorrYieldStat, 'ALICE - JHEP 10 (2018) 174', 'pl')
        legPHSD[iLeg].AddEntry(
            hDCorrYieldStat, 'ALICE - JHEP 10 (2018) 174', 'pl')
        legGossiaux[iLeg].AddEntry(
            hDCorrYieldStat, 'ALICE - JHEP 10 (2018) 174', 'pl')
        legCatania[iLeg].AddEntry(
            hDCorrYieldStat, 'ALICE - JHEP 10 (2018) 174', 'pl')

legTAMU['CorrYield'].AddEntry(gCorrYieldTAMU, 'FONLL x TAMU', 'fp')
legTAMU['Raa'].AddEntry(sTAMU['Min'], 'TAMU', 'l')
legTAMU['Signif'].AddEntry(gCorrYieldTAMU, 'FONLL x TAMU', 'fp')
legPHSD['CorrYield'].AddEntry(gCorrYieldPHSD, 'FONLL x PHSD', 'fp')
legPHSD['Raa'].AddEntry(sPHSD['Cent'], 'PHSD', 'l')
legPHSD['Signif'].AddEntry(gCorrYieldPHSD, 'FONLL x PHSD', 'fp')
legGossiaux['CorrYield'].AddEntry(
    gCorrYieldGossiaux, 'FONLL x MC@sHQ+EPOS2', 'fp')
legGossiaux['Raa'].AddEntry(sGossiaux['Col'], 'MC@sHQ+EPOS2', 'l')
legGossiaux['Signif'].AddEntry(
    gCorrYieldGossiaux, 'FONLL x MC@sHQ+EPOS2', 'fp')
legCatania['CorrYield'].AddEntry(
    gCorrYieldCatania, 'FONLL x Catania (frag+coal)', 'fp')
legCatania['Raa'].AddEntry(sCatania['Cent'], 'Catania (frag+coal)', 'l')
legCatania['Signif'].AddEntry(
    gCorrYieldCatania, 'FONLL x Catania (frag+coal)', 'fp')

lineAtOne = TLine(cutVars['Pt']['min'][0]-1, 1., cutVars['Pt']
                  ['max'][len(cutVars['Pt']['max'])-1]+1, 1.)
lineAtOne.SetLineColor(kBlack)
lineAtOne.SetLineStyle(9)

outfile.cd()

cRaa = TCanvas('cRaa', '', 1000, 1000)
cRaa.Divide(2, 2)
cRaa.cd(1).DrawFrame(cutVars['Pt']['min'][0]-1, 0, cutVars['Pt']['max']
                     [len(cutVars['Pt']['max'])-1]+1, 2, ';#it{p}_{T} (GeV/#it{c});#it{R}_{AA}')
lineAtOne.Draw('same')
for iRaaTAMU in sTAMU:
    SetSplineStyle(sTAMU[iRaaTAMU], colorsFill['TAMU'])
    sTAMU[iRaaTAMU].Draw('same')
    sTAMU[iRaaTAMU].Write('sTAMU%s' % iRaaTAMU)
if inputCfg['PublishedResult']['plotPublished']:
    hDRaaStat.Draw('same')
    gDRaaSyst.Draw('2')
legTAMU['Raa'].Draw('same')
cRaa.cd(2).DrawFrame(cutVars['Pt']['min'][0]-1, 0, cutVars['Pt']['max']
                     [len(cutVars['Pt']['max'])-1]+1, 2, ';#it{p}_{T} (GeV/#it{c});#it{R}_{AA}')
lineAtOne.Draw('same')
for iRaaPHSD in sPHSD:
    SetSplineStyle(sPHSD[iRaaPHSD], colorsFill['PHSD'])
    sPHSD[iRaaPHSD].Draw('same')
    sPHSD[iRaaPHSD].Write('sPHSD%s' % iRaaPHSD)
if inputCfg['PublishedResult']['plotPublished']:
    hDRaaStat.Draw('same')
    gDRaaSyst.Draw('2')
legPHSD['Raa'].Draw('same')
cRaa.cd(3).DrawFrame(cutVars['Pt']['min'][0]-1, 0, cutVars['Pt']['max']
                     [len(cutVars['Pt']['max'])-1]+1, 2, ';#it{p}_{T} (GeV/#it{c});#it{R}_{AA}')
lineAtOne.Draw('same')
for iRaaGossiaux in sGossiaux:
    SetSplineStyle(sGossiaux[iRaaGossiaux], colorsFill['Gossiaux'])
    sGossiaux[iRaaGossiaux].Draw('same')
    sGossiaux[iRaaGossiaux].Write('sGossiaux%s' % iRaaGossiaux)
if inputCfg['PublishedResult']['plotPublished']:
    hDRaaStat.Draw('same')
    gDRaaSyst.Draw('2')
legGossiaux['Raa'].Draw('same')
cRaa.cd(4).DrawFrame(cutVars['Pt']['min'][0]-1, 0, cutVars['Pt']['max']
                     [len(cutVars['Pt']['max'])-1]+1, 2, ';#it{p}_{T} (GeV/#it{c});#it{R}_{AA}')
lineAtOne.Draw('same')
for iRaaCatania in sCatania:
    SetSplineStyle(sCatania[iRaaCatania], colorsFill['Catania'])
    sCatania[iRaaCatania].Draw('same')
    sCatania[iRaaCatania].Write('sCatania%s' % iRaaCatania)
if inputCfg['PublishedResult']['plotPublished']:
    hDRaaStat.Draw('same')
    gDRaaSyst.Draw('2')
legCatania['Raa'].Draw('same')
cRaa.Write()
cRaa.SaveAs('%s_Raa.pdf' % outFileName)

cCorrYield = TCanvas('cCorrYield', '', 1000, 1000)
cCorrYield.Divide(2, 2)
cCorrYield.cd(1).DrawFrame(cutVars['Pt']['min'][0]-1, 3.e-6, cutVars['Pt']['max'][len(
    cutVars['Pt']['max'])-1]+1, 2., ';#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}')
cCorrYield.cd(1).SetLogy()
SetGraphStyle(gCorrYieldTAMU, colors['TAMU'], kWhite, markers['TAMU'])
gCorrYieldTAMU.Draw('P2')
if inputCfg['PublishedResult']['plotPublished']:
    hDCorrYieldStat.Draw('same')
    gDCorrYieldSyst.Draw('2')
legTAMU['CorrYield'].Draw('same')
gCorrYieldTAMU.Write('gCorrYieldTAMU')
cCorrYield.cd(2).SetLogy()
cCorrYield.cd(2).DrawFrame(cutVars['Pt']['min'][0]-1, 3.e-6, cutVars['Pt']['max'][len(
    cutVars['Pt']['max'])-1]+1, 2., ';#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}')
SetGraphStyle(gCorrYieldPHSD, colors['PHSD'], kWhite, markers['PHSD'])
gCorrYieldPHSD.Draw('P2')
if inputCfg['PublishedResult']['plotPublished']:
    hDCorrYieldStat.Draw('same')
    gDCorrYieldSyst.Draw('2')
legPHSD['CorrYield'].Draw('same')
gCorrYieldPHSD.Write('gCorrYieldPHSD')
cCorrYield.cd(3).SetLogy()
cCorrYield.cd(3).DrawFrame(cutVars['Pt']['min'][0]-1, 3.e-6, cutVars['Pt']['max'][len(
    cutVars['Pt']['max'])-1]+1, 2., ';#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}')
SetGraphStyle(gCorrYieldGossiaux,
              colors['Gossiaux'], kWhite, markers['Gossiaux'])
gCorrYieldGossiaux.Draw('P2')
if inputCfg['PublishedResult']['plotPublished']:
    hDCorrYieldStat.Draw('same')
    gDCorrYieldSyst.Draw('2')
legGossiaux['CorrYield'].Draw('same')
gCorrYieldGossiaux.Write('gCorrYieldGossiaux')
cCorrYield.cd(4).SetLogy()
cCorrYield.cd(4).DrawFrame(cutVars['Pt']['min'][0]-1, 3.e-6, cutVars['Pt']['max'][len(
    cutVars['Pt']['max'])-1]+1, 2., ';#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}')
SetGraphStyle(gCorrYieldCatania, colors['Catania'], kWhite, markers['Catania'])
gCorrYieldCatania.Draw('P2')
if inputCfg['PublishedResult']['plotPublished']:
    hDCorrYieldStat.Draw('same')
    gDCorrYieldSyst.Draw('2')
gCorrYieldCatania.Write('gCorrYieldCatania')
legCatania['CorrYield'].Draw('same')
cCorrYield.Write()
cCorrYield.SaveAs('%s_CorrYield.pdf' % outFileName)

cRawYield = TCanvas('cRawYield', '', 1000, 1000)
cRawYield.Divide(2, 2)
cRawYield.cd(1).DrawFrame(cutVars['Pt']['min'][0]-1, hRawYieldTAMU.GetMinimum()*0.5, cutVars['Pt']['max']
                          [len(cutVars['Pt']['max'])-1]+1, hRawYieldTAMU.GetMaximum()*1.5, ';#it{p}_{T} (GeV/#it{c});raw yield')
cRawYield.cd(1).SetLogy()
SetGraphStyle(gRawYieldTAMU, colors['TAMU'], kWhite, markers['TAMU'])
SetGraphStyle(hRawYieldTAMU, colors['TAMU'], kWhite, markers['TAMU'])
gRawYieldTAMU.Draw('2')
hRawYieldTAMU.Draw('same')
legTAMU['Signif'].Draw('same')
gRawYieldTAMU.Write('gRawYieldTAMU')
hRawYieldTAMU.Write('hRawYieldTAMU')
cRawYield.cd(2).SetLogy()
cRawYield.cd(2).DrawFrame(cutVars['Pt']['min'][0]-1, hRawYieldPHSD.GetMinimum()*0.5, cutVars['Pt']['max']
                          [len(cutVars['Pt']['max'])-1]+1, hRawYieldPHSD.GetMaximum()*1.5, ';#it{p}_{T} (GeV/#it{c});raw yield')
SetGraphStyle(gRawYieldPHSD, colors['PHSD'], kWhite, markers['PHSD'])
SetGraphStyle(hRawYieldPHSD, colors['PHSD'], kWhite, markers['PHSD'])
gRawYieldPHSD.Draw('2')
hRawYieldPHSD.Draw('same')
legPHSD['Signif'].Draw('same')
gRawYieldPHSD.Write('gRawYieldPHSD')
hRawYieldPHSD.Write('hRawYieldPHSD')
cRawYield.cd(3).SetLogy()
cRawYield.cd(3).DrawFrame(cutVars['Pt']['min'][0]-1, hRawYieldGossiaux.GetMinimum()*0.5, cutVars['Pt']['max']
                          [len(cutVars['Pt']['max'])-1]+1, hRawYieldGossiaux.GetMaximum()*1.5, ';#it{p}_{T} (GeV/#it{c});raw yield')
SetGraphStyle(gRawYieldGossiaux,
              colors['Gossiaux'], kWhite, markers['Gossiaux'])
SetGraphStyle(hRawYieldGossiaux,
              colors['Gossiaux'], kWhite, markers['Gossiaux'])
gRawYieldGossiaux.Draw('2')
hRawYieldGossiaux.Draw('same')
legGossiaux['Signif'].Draw('same')
gRawYieldGossiaux.Write('gRawYieldGossiaux')
hRawYieldGossiaux.Write('hRawYieldGossiaux')
cRawYield.cd(4).SetLogy()
cRawYield.cd(4).DrawFrame(cutVars['Pt']['min'][0]-1, hRawYieldCatania.GetMinimum()*0.5, cutVars['Pt']['max']
                          [len(cutVars['Pt']['max'])-1]+1, hRawYieldCatania.GetMaximum()*1.5, ';#it{p}_{T} (GeV/#it{c});raw yield')
SetGraphStyle(gRawYieldCatania, colors['Catania'], kWhite, markers['Catania'])
SetGraphStyle(hRawYieldCatania, colors['Catania'], kWhite, markers['Catania'])
gRawYieldCatania.Draw('2')
hRawYieldCatania.Draw('same')
gRawYieldCatania.Write('gRawYieldCatania')
hRawYieldCatania.Write('hRawYieldCatania')
legCatania['Signif'].Draw('same')
cRawYield.Write()
cRawYield.SaveAs('%s_RawYield.pdf' % outFileName)

cSignif = TCanvas('cSignif', '', 1000, 1000)
cSignif.Divide(2, 2)
cSignif.cd(1).DrawFrame(cutVars['Pt']['min'][0]-1, 0, cutVars['Pt']['max'][len(
    cutVars['Pt']['max'])-1]+1, hSignifTAMU.GetMaximum()*1.5, ';#it{p}_{T} (GeV/#it{c});significance')
SetGraphStyle(gSignifTAMU, colors['TAMU'], kWhite, markers['TAMU'])
SetGraphStyle(hSignifTAMU, colors['TAMU'], kWhite, markers['TAMU'])
gSignifTAMU.Draw('2')
hSignifTAMU.Draw('same')
legTAMU['Signif'].Draw('same')
gSignifTAMU.Write('gSignifTAMU')
hSignifTAMU.Write('hSignifTAMU')
cSignif.cd(2).DrawFrame(cutVars['Pt']['min'][0]-1, 0, cutVars['Pt']['max'][len(
    cutVars['Pt']['max'])-1]+1, hSignifPHSD.GetMaximum()*1.5, ';#it{p}_{T} (GeV/#it{c});significance')
SetGraphStyle(gSignifPHSD, colors['PHSD'], kWhite, markers['PHSD'])
SetGraphStyle(hSignifPHSD, colors['PHSD'], kWhite, markers['PHSD'])
gSignifPHSD.Draw('P2')
hSignifPHSD.Draw('same')
legPHSD['Signif'].Draw('same')
gSignifPHSD.Write('gSignifPHSD')
cSignif.cd(3).DrawFrame(cutVars['Pt']['min'][0]-1, 0, cutVars['Pt']['max'][len(cutVars['Pt']
                                                                               ['max'])-1]+1, gSignifGossiaux.GetMaximum()*1.5, ';#it{p}_{T} (GeV/#it{c});significance')
SetGraphStyle(gSignifGossiaux, colors['Gossiaux'], kWhite, markers['Gossiaux'])
SetGraphStyle(hSignifGossiaux, colors['Gossiaux'], kWhite, markers['Gossiaux'])
gSignifGossiaux.Draw('P2')
hSignifGossiaux.Draw('same')
legGossiaux['Signif'].Draw('same')
gSignifGossiaux.Write('gSignifGossiaux')
hSignifGossiaux.Write('hSignifGossiaux')
cSignif.cd(4).DrawFrame(cutVars['Pt']['min'][0]-1, 0, cutVars['Pt']['max'][len(cutVars['Pt']
                                                                               ['max'])-1]+1, hSignifCatania.GetMaximum()*1.5, ';#it{p}_{T} (GeV/#it{c});significance')
SetGraphStyle(gSignifCatania, colors['Catania'], kWhite, markers['Catania'])
SetGraphStyle(hSignifCatania, colors['Catania'], kWhite, markers['Catania'])
gSignifCatania.Draw('P2')
hSignifCatania.Draw('same')
legCatania['Signif'].Draw('same')
gSignifCatania.Write('gSignifCatania')
hSignifCatania.Write('hSignifCatania')
cSignif.Write()
cSignif.SaveAs('%s_Significance.pdf' % outFileName)

cBackground = TCanvas('cBackground', '', 500, 500)
cBackground.SetLogy()
hBackground.Draw()
hBackground.Write('hBackground')
cBackground.SaveAs('%s_Background.pdf' % outFileName)

cFprompt = TCanvas('cFprompt', '', 500, 500)
cFprompt.DrawFrame(cutVars['Pt']['min'][0], 0., cutVars['Pt']['max'][len(
    cutVars['Pt']['max'])-1], 1., ';#it{p}_{T} (GeV/#it{c});#it{f}_{prompt} (#it{f}_{c})')
gFprompt.Draw('p')
gFprompt.Write('gFprompt')
cFprompt.SaveAs('%s_Fprompt.pdf' % outFileName)

cEff = TCanvas('cEff', '', 500, 500)
cEff.DrawFrame(cutVars['Pt']['min'][0], 1.e-4, cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1],
               1., ';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon) #times 2#it{y}_{fid}')
cEff.SetLogy()
hAccEffPrompt.Draw('same')
hAccEffFD.Draw('same')
hAccEffPrompt.Write()
hAccEffFD.Write()
cEff.Write()
cEff.SaveAs('%s_Efficiency.pdf' % outFileName)

if six.PY2:
    raw_input('Press enter to exit')
elif six.PY3:
    input('Press enter to exit')
