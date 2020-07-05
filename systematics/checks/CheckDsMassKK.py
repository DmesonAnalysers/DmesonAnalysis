'''
python script to perform checks of massKK distributions in data and wd
run: python CheckDsMassKK.py cfgFileNameMC.yml cfgFileNameData.yml
'''

import sys
import argparse
import ctypes
import numpy as np
import yaml
from ROOT import AliHFInvMassFitter, AliVertexingHFUtils  # pylint: disable=import-error,no-name-in-module
from ROOT import TFile, TCanvas, TH1F, TF1, TVirtualFitter, TDatabasePDG, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kGreen, kAzure, kRed, kRainBow, kFullCircle, kFullSquare, kFullDiamond # pylint: disable=import-error,no-name-in-module
from ROOT import kFullCross, kOpenCircle, kOpenSquare, kOpenDiamond, kOpenCross  # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
# pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle, DivideCanvas
from utils.AnalysisUtils import VoigtFunc, ExpoPowLaw

# set global style
SetGlobalStyle(padrightmargin=0.1, padleftmargin=0.15, padbottommargin=0.14, padtopmargin=0.075,
               titleoffsety=1.5, palette=kRainBow, maxdigits=2)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileNameMC', metavar='text', default='cfgFileNameMC.yml',
                    help='config file name with MC root input files')
parser.add_argument('cfgFileNameData', metavar='text', default='cfgFileNameData.yml',
                    help='config file name with data root input files')
parser.add_argument('outputPath', metavar='text', default='outputPath',
                    help='path for outputs')
parser.add_argument('--subSB', action='store_true', default=False,
                    help='flag to subtract bkg from SB instead than fit')
parser.add_argument('--nSigma', type=float, default=3.,
                    help='number of sigmas for signal region')
args = parser.parse_args()

# load histos from input files
hMassKKvsMassKKpiMC, hMassKKvsMassKKpiData = [], []
# MC
with open(args.cfgFileNameMC, 'r') as ymlCfgFileMC:
    inputCfgMC = yaml.load(ymlCfgFileMC, yaml.FullLoader)
inFileNamesMC = inputCfgMC['filename']
inDirNameMC = inputCfgMC['dirname']
inListNameMC = inputCfgMC['listname']
if not isinstance(inFileNamesMC, list):
    inFileNamesMC = [inFileNamesMC]
for iFile, inFileName in enumerate(inFileNamesMC):
    inFileMC = TFile.Open(inFileName)
    if iFile == 0:
        inListMC = inFileMC.Get(f'{inDirNameMC}/{inListNameMC}')
    else:
        inListMC.Add(inFileMC.Get(f'{inDirNameMC}/{inListNameMC}'))
iPt = 0
while inListMC.FindObject(f'hMassKKVsKKpiPt{iPt}'):
    hMassKKvsMassKKpiMC.append(inListMC.FindObject(f'hMassKKVsKKpiPt{iPt}'))
    hMassKKvsMassKKpiMC[-1].SetNameTitle(f'hMassKKVsKKpiPt{iPt}_MC',
                                         ';#it{M}(KK#pi} (GeV/#it{c}^{2});#it{M}(KK) (GeV/#it{c}^{2})')
    iPt += 1
nPtBinsMC = len(hMassKKvsMassKKpiMC)
# data
with open(args.cfgFileNameData, 'r') as ymlCfgFileData:
    inputCfgData = yaml.load(ymlCfgFileData, yaml.FullLoader)
inFileNamesData = inputCfgData['filename']
inDirNameData = inputCfgData['dirname']
inListNameData = inputCfgData['listname']
inCutObjName = inputCfgData['listname'].replace('coutputDs', 'coutputDsCuts')

if not isinstance(inFileNamesData, list):
    inFileNamesData = [inFileNamesData]
for iFile, inFileName in enumerate(inFileNamesData):
    inFileData = TFile.Open(inFileName)
    if iFile == 0:
        inListData = inFileData.Get(f'{inDirNameData}/{inListNameData}')
    else:
        inListData.Merge(inFileData.Get(f'{inDirNameData}/{inListNameData}'))
    if iFile == 0:
        cutObj = inFileData.Get(f'{inDirNameData}/{inCutObjName}').FindObject('AnalysisCuts')
iPt = 0
while inListData.FindObject(f'hMassKKVsKKpiPt{iPt}'):
    hMassKKvsMassKKpiData.append(inListData.FindObject(f'hMassKKVsKKpiPt{iPt}'))
    hMassKKvsMassKKpiData[-1].SetNameTitle(f'hMassKKVsKKpiPt{iPt}_data',
                                           ';#it{M}(KK#pi} (GeV/#it{c}^{2});#it{M}(KK) (GeV/#it{c}^{2})')
    iPt += 1
nPtBinsData = len(hMassKKvsMassKKpiData)

if nPtBinsMC != nPtBinsData:
    print('ERROR: number of pT bins different in data and MC! Exit')
    sys.exit()

# get bin limits from cut object
ptLimsArray = cutObj.GetPtBinLimits()
ptLims = []
for iPt in range(nPtBinsMC+1):
    ptLims.append(ptLimsArray[iPt])

# perform fit to massKKpi to infer mean, sigma, S and B
hMassForFitMC, hMassForFitData, massFitterMC, massFitterData, mean, sigma, signalMC, bkgMC, \
    signalData, bkgData = ([] for _ in range(10))
for iPt in range(nPtBinsMC):
    minBinMassKK = hMassKKvsMassKKpiMC[iPt].GetYaxis().FindBin(0.98)
    maxBinMassKK = hMassKKvsMassKKpiMC[iPt].GetYaxis().FindBin(1.08)
    hMassForFitMC.append(TH1F())
    hMassForFitData.append(TH1F())

    # to cast TH1D to TH1F
    AliVertexingHFUtils.RebinHisto(hMassKKvsMassKKpiMC[iPt].ProjectionX(
        f'hMassKKMC{iPt}', minBinMassKK, maxBinMassKK), 1).Copy(hMassForFitMC[iPt])
    AliVertexingHFUtils.RebinHisto(hMassKKvsMassKKpiData[iPt].ProjectionX(
        f'hMassKKData{iPt}', minBinMassKK, maxBinMassKK), 1).Copy(hMassForFitData[iPt])

    massFitterMC.append(AliHFInvMassFitter(hMassForFitMC[iPt], 1.9, 2.15,
                                           AliHFInvMassFitter.kExpo, AliHFInvMassFitter.kGaus))
    massFitterMC[iPt].SetUseLikelihoodFit()
    massFitterMC[iPt].SetInitialGaussianMean(1.970)
    massFitterMC[iPt].SetInitialGaussianSigma(0.008)
    massFitterMC[iPt].MassFitter(False)

    mean.append(massFitterMC[iPt].GetMean())
    sigma.append(massFitterMC[iPt].GetSigma())
    Serr = ctypes.c_double()
    massLowLimit = mean[iPt]-args.nSigma*sigma[iPt]
    massHighLimit = mean[iPt]+args.nSigma*sigma[iPt]
    binLow = hMassKKvsMassKKpiMC[iPt].GetXaxis().FindBin(massLowLimit*1.0001)
    binHigh = hMassKKvsMassKKpiMC[iPt].GetXaxis().FindBin(massHighLimit*0.9999)
    S = massFitterMC[iPt].GetRawYieldBinCounting(Serr, massLowLimit, massHighLimit)
    B = hMassForFitMC[iPt].Integral(binLow, binHigh) - S
    signalMC.append(S)
    bkgMC.append(B)

    massFitterData.append(AliHFInvMassFitter(hMassForFitData[iPt], 1.9, 2.15,
                                             AliHFInvMassFitter.kExpo, AliHFInvMassFitter.kGaus))
    massFitterData[iPt].SetUseLikelihoodFit()
    massFitterData[iPt].SetFixGaussianMean(mean[iPt])
    massFitterData[iPt].SetFixGaussianSigma(sigma[iPt])
    massFitterData[iPt].MassFitter(False)

    Sdata = massFitterData[iPt].GetRawYieldBinCounting(Serr, massLowLimit, massHighLimit)
    Bdata = hMassForFitData[iPt].Integral(binLow, binHigh) - Sdata
    signalData.append(Sdata)
    bkgData.append(Bdata)

hMassKKSigRegionMC, hMassKKSigRegionData, hSideBandsMC, hSideBandsData, hMassKKwoPhiMC, hMassKKwoPhiData, \
    fMassKKwoPhiMC, fMassKKwoPhiData, hMassKKBkgSubMC, hMassKKBkgSubData, fMassKKMC, fMassKKData, \
        hConfIntMassKKMC, hConfIntMassKKData = ([] for _ in range(14))

hMeanMC = TH1F('hMeanMC', ';#it{p}_{T} (GeV/#it{c}); mean (GeV/#it{c}^{2})', nPtBinsMC, np.array(ptLims, 'f'))
hMeanData = TH1F('hMeanData', ';#it{p}_{T} (GeV/#it{c}); mean (GeV/#it{c}^{2})', nPtBinsMC, np.array(ptLims, 'f'))
hFWHMMC = TH1F('hFWHMMC', ';#it{p}_{T} (GeV/#it{c}); sigma (GeV/#it{c}^{2})', nPtBinsMC, np.array(ptLims, 'f'))
hFWHMData = TH1F('hFWHMData', ';#it{p}_{T} (GeV/#it{c}); sigma (GeV/#it{c}^{2})', nPtBinsMC, np.array(ptLims, 'f'))
SetObjectStyle(hMeanMC, color=kRed+1, markerstyle=kFullSquare)
SetObjectStyle(hMeanData, color=kAzure+4, markerstyle=kFullCircle)
SetObjectStyle(hFWHMMC, color=kRed+1, markerstyle=kFullSquare)
SetObjectStyle(hFWHMData, color=kAzure+4, markerstyle=kFullCircle)
hEffMC = {5: TH1F('hEffMC05', ';#it{p}_{T} (GeV/#it{c}); efficiency', nPtBinsMC, np.array(ptLims, 'f')),
          10: TH1F('hEffMC10', ';#it{p}_{T} (GeV/#it{c}); efficiency', nPtBinsMC, np.array(ptLims, 'f')),
          15: TH1F('hEffMC15', ';#it{p}_{T} (GeV/#it{c}); efficiency', nPtBinsMC, np.array(ptLims, 'f')),
          20: TH1F('hEffMC20', ';#it{p}_{T} (GeV/#it{c}); efficiency', nPtBinsMC, np.array(ptLims, 'f'))}
hEffData = {5: TH1F('hEffData05', ';#it{p}_{T} (GeV/#it{c}); efficiency', nPtBinsMC, np.array(ptLims, 'f')),
            10: TH1F('hEffData10', ';#it{p}_{T} (GeV/#it{c}); efficiency', nPtBinsMC, np.array(ptLims, 'f')),
            15: TH1F('hEffData15', ';#it{p}_{T} (GeV/#it{c}); efficiency', nPtBinsMC, np.array(ptLims, 'f')),
            20: TH1F('hEffData20', ';#it{p}_{T} (GeV/#it{c}); efficiency', nPtBinsMC, np.array(ptLims, 'f'))}
markersData = {5: kFullSquare, 10: kFullCircle, 15: kFullDiamond, 20: kFullCross}
markersMC = {5: kOpenSquare, 10: kOpenCircle, 15: kOpenDiamond, 20: kOpenCross}
for cut in hEffMC:
    SetObjectStyle(hEffMC[cut], color=kRed+1, fillstyle=0, markerstyle=markersMC[cut])
    SetObjectStyle(hEffData[cut], color=kAzure+4, fillstyle=0, markerstyle=markersData[cut])

massPhi = TDatabasePDG.Instance().GetParticle(333).Mass()
massK = TDatabasePDG.Instance().GetParticle(321).Mass()
widthPhi = TDatabasePDG.Instance().GetParticle(333).Width()

for iPt in range(nPtBinsMC):
    # get massKK distributions in signal region
    massLowLimit = mean[iPt]-args.nSigma*sigma[iPt]
    massHighLimit = mean[iPt]+args.nSigma*sigma[iPt]
    binLow = hMassKKvsMassKKpiMC[iPt].GetXaxis().FindBin(massLowLimit*1.0001)
    binHigh = hMassKKvsMassKKpiMC[iPt].GetXaxis().FindBin(massHighLimit*0.9999)
    hMassKKvsMassKKpiMC[iPt].GetYaxis().SetRangeUser(0.98, 1.08)
    hMassKKSigRegionMC.append(hMassKKvsMassKKpiMC[iPt].ProjectionY(f'hMassKKSigRegionMCPt{iPt}', binLow, binHigh))
    SetObjectStyle(hMassKKSigRegionMC[iPt], markersize=0.5)

    massLowLimit = mean[iPt]-args.nSigma*sigma[iPt]
    massHighLimit = mean[iPt]+args.nSigma*sigma[iPt]
    hMassKKvsMassKKpiData[iPt].GetYaxis().SetRangeUser(0.98, 1.08)
    binLow = hMassKKvsMassKKpiData[iPt].GetXaxis().FindBin(massLowLimit*1.0001)
    binHigh = hMassKKvsMassKKpiData[iPt].GetXaxis().FindBin(massHighLimit*0.9999)
    hMassKKSigRegionData.append(hMassKKvsMassKKpiData[iPt].ProjectionY(
        f'hMassKKSigRegionDataPt{iPt}', binLow, binHigh))
    SetObjectStyle(hMassKKSigRegionData[iPt], markersize=0.5)

    # get massKK distributions from massKKpi sidebands
    massLowLimit = 1.8
    massHighLimit = 2.0
    binLow = hMassKKvsMassKKpiMC[iPt].GetXaxis().FindBin(massLowLimit*0.9999)
    binHigh = hMassKKvsMassKKpiMC[iPt].GetXaxis().FindBin(massHighLimit*1.0001)
    hMassKKvsMassKKpiMC[iPt].GetYaxis().SetRangeUser(0.98, 1.08)
    hSBMC = hMassKKvsMassKKpiMC[iPt].ProjectionY(f'hMassKKSBMCPt{iPt}', 0, binLow) # left SB
    hSBMC.Add(hMassKKvsMassKKpiMC[iPt].ProjectionY(f'hMassKKSBMCPt{iPt}_R', binHigh, 1000000)) # right SB
    hSideBandsMC.append(hSBMC)
    hSideBandsMC[iPt].SetDirectory(0)
    hSideBandsMC[iPt].Scale(bkgMC[iPt]/hSideBandsMC[iPt].Integral())
    SetObjectStyle(hSideBandsMC[iPt], color=kGreen+2, linewidth=0, fillalpha=0.3)

    binLow = hMassKKvsMassKKpiData[iPt].GetXaxis().FindBin(massLowLimit*0.9999)
    binHigh = hMassKKvsMassKKpiData[iPt].GetXaxis().FindBin(massHighLimit*1.0001)
    hMassKKvsMassKKpiData[iPt].GetYaxis().SetRangeUser(0.98, 1.08)
    hSBData = hMassKKvsMassKKpiData[iPt].ProjectionY(f'hMassKKSBDataPt{iPt}', 0, binLow) # left SB
    hSBData.Add(hMassKKvsMassKKpiData[iPt].ProjectionY(f'hMassKKSBDataPt{iPt}_R', binHigh, 1000000)) # right SB
    hSideBandsData.append(hSBData)
    hSideBandsData[iPt].SetDirectory(0)
    hSideBandsData[iPt].Scale(bkgData[iPt]/hSideBandsData[iPt].Integral())
    SetObjectStyle(hSideBandsData[iPt], color=kGreen+2, linewidth=0, fillalpha=0.3)

    # fit massKK ditribution w/o phi peak
    hMassKKwoPhiMC.append(hMassKKSigRegionMC[iPt].Clone(f'hMassKKwoPhiMCPt{iPt}'))
    for iBin in range(1, hMassKKwoPhiMC[iPt].GetNbinsX()+1):
        binCentre = hMassKKwoPhiMC[iPt].GetBinCenter(iBin)
        if massPhi-0.020 < binCentre < massPhi+0.020:
            hMassKKwoPhiMC[iPt].SetBinContent(iBin, 0)
    fMassKKwoPhiMC.append(TF1(f'fMassKKwoPhiMCPt{iPt}', ExpoPowLaw, 2*massK, 1.08, 3))
    fMassKKwoPhiMC[iPt].SetParameters(hSideBandsMC[iPt].Integral(), 2*massK, 0.01)
    fMassKKwoPhiMC[iPt].FixParameter(1, 2*massK)
    hMassKKwoPhiMC[iPt].Fit(fMassKKwoPhiMC[iPt], 'ER0SQ')

    hMassKKwoPhiData.append(hMassKKSigRegionData[iPt].Clone(f'hMassKKwoPhiDataPt{iPt}'))
    for iBin in range(1, hMassKKwoPhiData[iPt].GetNbinsX()+1):
        binCentre = hMassKKwoPhiData[iPt].GetBinCenter(iBin)
        if massPhi-0.020 < binCentre < massPhi+0.020:
            hMassKKwoPhiData[iPt].SetBinContent(iBin, 0)
    fMassKKwoPhiData.append(TF1(f'fMassKKwoPhiDataPt{iPt}', ExpoPowLaw, 2*massK, 1.08, 3))
    fMassKKwoPhiData[iPt].SetParameters(hSideBandsData[iPt].Integral(), 2*massK, 0.01)
    fMassKKwoPhiData[iPt].FixParameter(1, 2*massK)
    hMassKKwoPhiData[iPt].Fit(fMassKKwoPhiData[iPt], 'ER0SQ')

    # subtract bkg from inclusive histo
    hMassKKBkgSubData.append(hMassKKSigRegionData[iPt].Clone(f'hMassKKBkgSubDataPt{iPt}'))
    hMassKKBkgSubData[iPt].Sumw2()
    hMassKKBkgSubMC.append(hMassKKSigRegionMC[iPt].Clone(f'hMassKKBkgSubMCPt{iPt}'))
    hMassKKBkgSubMC[iPt].Sumw2()
    if args.subSB:
        hMassKKBkgSubMC[iPt].Add(hSideBandsMC[iPt], -1)
        hMassKKBkgSubData[iPt].Add(hSideBandsData[iPt], -1)
    else:
        hMassKKBkgSubMC[iPt].Add(fMassKKwoPhiMC[iPt], -1)
        hMassKKBkgSubData[iPt].Add(fMassKKwoPhiData[iPt], -1)

    hMassKKBkgSubMC[iPt].Scale(1./hMassKKBkgSubMC[iPt].Integral())
    SetObjectStyle(hMassKKBkgSubMC[iPt], color=kRed+1, markersize=0.5, markerstyle=kFullSquare)
    hMassKKBkgSubData[iPt].Scale(1./hMassKKBkgSubData[iPt].Integral())
    SetObjectStyle(hMassKKBkgSubData[iPt], color=kAzure+4, markersize=0.5, markerstyle=kFullCircle)

    # fit distribution with voigt function
    fMassKKMC.append(TF1(f'fMassKKMCPt{iPt}', VoigtFunc, 0.98, 1.08, 4))
    fMassKKMC[iPt].SetLineColor(kRed+1)
    fMassKKMC[iPt].SetParameters(hMassKKBkgSubMC[iPt].GetBinWidth(1), massPhi, 0.0015, widthPhi)
    fMassKKMC[iPt].FixParameter(0, hMassKKBkgSubMC[iPt].GetBinWidth(1))
    fMassKKMC[iPt].SetParLimits(1, massPhi*0.5, massPhi*2)
    fMassKKMC[iPt].SetParLimits(2, 0., 0.1)
    fMassKKMC[iPt].FixParameter(3, widthPhi)
    resMC = hMassKKBkgSubMC[iPt].Fit(fMassKKMC[iPt], 'ER0SQ')
    hConfIntMassKKMC.append(TH1F(f'hConfIntMassKKMC_{iPt}', '', 500, 0.98, 1.08))
    TVirtualFitter.GetFitter().GetConfidenceIntervals(hConfIntMassKKMC[iPt], 0.68)
    SetObjectStyle(hConfIntMassKKMC[iPt], color=kRed+1, alpha=0.3, markerstyle=0)

    fMassKKData.append(TF1(f'fMassKKDataPt{iPt}', VoigtFunc, 0.98, 1.08, 4))
    fMassKKData[iPt].SetLineColor(kAzure+4)
    fMassKKData[iPt].SetParameters(hMassKKBkgSubData[iPt].GetBinWidth(1), massPhi, 0.0015, widthPhi)
    fMassKKData[iPt].FixParameter(0, hMassKKBkgSubData[iPt].GetBinWidth(1))
    fMassKKData[iPt].SetParLimits(1, massPhi*0.5, massPhi*2)
    fMassKKData[iPt].SetParLimits(2, 0., 0.1)
    fMassKKData[iPt].FixParameter(3, widthPhi)
    resData = hMassKKBkgSubData[iPt].Fit(fMassKKData[iPt], 'ER0SQ')
    hConfIntMassKKData.append(TH1F(f'hConfIntMassKKData_{iPt}', '', 500, 0.98, 1.08))
    TVirtualFitter.GetFitter().GetConfidenceIntervals(hConfIntMassKKData[iPt], 0.68)
    SetObjectStyle(hConfIntMassKKData[iPt], color=kAzure+4, alpha=0.3, markerstyle=0)

    # fill histos
    hMeanMC.SetBinContent(iPt+1, fMassKKMC[iPt].GetParameter(1))
    hMeanMC.SetBinError(iPt+1, fMassKKMC[iPt].GetParError(1))
    hFWHMMC.SetBinContent(iPt+1, widthPhi+np.sqrt(widthPhi**2/4 + fMassKKMC[iPt].GetParameter(2)**2))
    hFWHMMC.SetBinError(iPt+1, fMassKKMC[iPt].GetParError(2))
    if (resData.Status() == 0 or resData.Status() == 4000) and resData.Chi2()/resData.Ndf() < 2:
        hMeanData.SetBinContent(iPt+1, fMassKKData[iPt].GetParameter(1))
        hMeanData.SetBinError(iPt+1, fMassKKData[iPt].GetParError(1))
        hFWHMData.SetBinContent(iPt+1, widthPhi+np.sqrt(widthPhi**2/4 + fMassKKData[iPt].GetParameter(2)**2))
        hFWHMData.SetBinError(iPt+1, fMassKKData[iPt].GetParError(2))
    for cut in hEffMC:
        totInt = hMassKKBkgSubMC[iPt].GetBinWidth(1)
        effMC = fMassKKMC[iPt].Integral(massPhi-cut/1000, massPhi+cut/1000, 1.e-9) / totInt
        effMCUnc = fMassKKMC[iPt].IntegralError(massPhi-cut/1000, massPhi+cut/1000) / totInt
        if effMCUnc == 0:
            effMCUnc = 1.e-20
        hEffMC[cut].SetBinContent(iPt+1, effMC)
        hEffMC[cut].SetBinError(iPt+1, effMCUnc)
        if (resData.Status() == 0 or resData.Status() == 4000) and resData.Chi2()/resData.Ndf() < 2:
            effData = fMassKKData[iPt].Integral(massPhi-cut/1000, massPhi+cut/1000, 1.e-9) / totInt
            effDataUnc = fMassKKData[iPt].IntegralError(massPhi-cut/1000, massPhi+cut/1000) / totInt
            hEffData[cut].SetBinContent(iPt+1, effData)
            hEffData[cut].SetBinError(iPt+1, effDataUnc)

legMassKK = TLegend(0.55, 0.75, 0.99, 0.95)
legMassKK.SetBorderSize(1)
legMassKK.SetTextSize(0.05)
legMassKK.AddEntry(hMassKKSigRegionMC[0], 'D_{s}^{+} signal region', 'lp')
legMassKK.AddEntry(hSideBandsMC[0], 'D_{s}^{+} sidebands', 'f')
legMassKK.AddEntry(fMassKKwoPhiMC[0], '#it{M}(KK) bkg fit', 'l')

legMassKKsub = TLegend(0.7, 0.825, 0.9, 0.925)
legMassKKsub.SetBorderSize(1)
legMassKKsub.SetTextSize(0.05)
legMassKKsub.AddEntry(hMassKKBkgSubMC[0], 'MC', 'lp')
legMassKKsub.AddEntry(hMassKKBkgSubData[0], 'Data', 'lp')

legMassKKsubFunc = TLegend(0.7, 0.825, 0.9, 0.925)
legMassKKsubFunc.SetBorderSize(1)
legMassKKsubFunc.SetTextSize(0.05)
legMassKKsubFunc.AddEntry(hConfIntMassKKMC[0], 'MC', 'fl')
legMassKKsubFunc.AddEntry(hConfIntMassKKData[0], 'Data', 'fl')

legPars = TLegend(0.6, 0.7, 0.85, 0.85)
legPars.SetBorderSize(1)
legPars.SetTextSize(0.05)
legPars.AddEntry(hMeanMC, 'MC', 'lp')
legPars.AddEntry(hMeanData, 'Data', 'lp')

legEff = TLegend(0.15, 0.14, 0.75, 0.44)
legEff.SetBorderSize(1)
legEff.SetTextSize(0.04)
legEff.SetNColumns(2)
legEff.SetHeader('MC       Data')
for cut in hEffMC:
    legEff.AddEntry(hEffMC[cut], '    ', 'lp')
    legEff.AddEntry(hEffData[cut], f'|#Delta#it{{M}}(KK)| < {cut:.0f} MeV/#it{{c}}^{{2}}', 'lp')

cMassKKvsMassKKpiMC = TCanvas('cMassKKvsMassKKpiMC', 'MC', 1920, 1080)
DivideCanvas(cMassKKvsMassKKpiMC, iPt)
for iPt in range(nPtBinsMC):
    cMassKKvsMassKKpiMC.cd(iPt+1).SetLogz()
    hMassKKvsMassKKpiMC[iPt].Draw('colz')
cMassKKvsMassKKpiMC.Modified()
cMassKKvsMassKKpiMC.Update()

cMassKKvsMassKKpiData = TCanvas('cMassKKvsMassKKpiData', 'Data', 1920, 1080)
DivideCanvas(cMassKKvsMassKKpiData, iPt)
for iPt in range(nPtBinsData):
    cMassKKvsMassKKpiData.cd(iPt+1).SetLogz()
    hMassKKvsMassKKpiData[iPt].Draw('colz')
cMassKKvsMassKKpiData.Modified()
cMassKKvsMassKKpiData.Update()

cMassKKMC = TCanvas('cMassKKMC', 'MC', 1920, 1080)
DivideCanvas(cMassKKMC, iPt)
for iPt in range(nPtBinsMC):
    cMassKKMC.cd(iPt+1)
    hMassKKSigRegionMC[iPt].DrawCopy('E')
    hSideBandsMC[iPt].DrawCopy('hist same')
    fMassKKwoPhiMC[iPt].Draw('same')
    legMassKK.Draw()
cMassKKMC.Modified()
cMassKKMC.Update()

cMassKKData = TCanvas('cMassKKData', 'Data', 1920, 1080)
DivideCanvas(cMassKKData, iPt)
for iPt in range(nPtBinsData):
    cMassKKData.cd(iPt+1)
    hMassKKSigRegionData[iPt].DrawCopy('E')
    hSideBandsData[iPt].DrawCopy('hist same')
    fMassKKwoPhiData[iPt].Draw('same')
    legMassKK.Draw()
cMassKKData.Modified()
cMassKKData.Update()

cMassKKBkgSub = TCanvas('cMassKKBkgSub', '', 1920, 1080)
DivideCanvas(cMassKKBkgSub, iPt)
for iPt in range(nPtBinsMC):
    cMassKKBkgSub.cd(iPt+1)
    hMassKKBkgSubMC[iPt].DrawCopy('E')
    hMassKKBkgSubData[iPt].DrawCopy('Esame')
    legMassKKsub.Draw()
cMassKKBkgSub.Modified()
cMassKKBkgSub.Update()

cMassKKFunc = TCanvas('cMassKKFunc', '', 1920, 1080)
DivideCanvas(cMassKKFunc, iPt)
for iPt in range(nPtBinsMC):
    hFrame = cMassKKFunc.cd(iPt+1).DrawFrame(1., 0., 1.04, 0.15, ';#it{M}(KK) (GeV/#it{c}^{2}); a.u.')
    hFrame.GetXaxis().SetNdivisions(505)
    hConfIntMassKKMC[iPt].DrawCopy('e3same')
    fMassKKMC[iPt].Draw('same')
    hConfIntMassKKData[iPt].DrawCopy('e3same')
    fMassKKData[iPt].Draw('same')
    legMassKKsubFunc.Draw()
cMassKKFunc.Modified()
cMassKKFunc.Update()

cMeanSigma = TCanvas('cMeanSigma', '', 1000, 500)
cMeanSigma.Divide(2, 1)
hFrameMean = cMeanSigma.cd(1).DrawFrame(ptLims[0], 1.01, ptLims[-1], 1.03,
                                        ';#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})')
hFrameMean.GetYaxis().SetNdivisions(505)
hMeanMC.DrawCopy('same')
hMeanData.DrawCopy('same')
legPars.Draw()
cMeanSigma.cd(2).DrawFrame(ptLims[0], 0, ptLims[-1], 0.015, ';#it{p}_{T} (GeV/#it{c});FWHM (GeV/#it{c}^{2})')
hFWHMMC.DrawCopy('same')
hFWHMData.DrawCopy('same')
legPars.Draw()
cMeanSigma.Modified()
cMeanSigma.Update()

cEfficiency = TCanvas('cEfficiency', '', 800, 800)
cEfficiency.DrawFrame(ptLims[0], 0., ptLims[-1], 1., ';#it{p}_{T} (GeV/#it{c});efficiency')
for cut in hEffMC:
    hEffMC[cut].DrawCopy('epsame')
    hEffData[cut].DrawCopy('epsame')
legEff.Draw()
cEfficiency.Modified()
cEfficiency.Update()

input('Press Enter to exit')
