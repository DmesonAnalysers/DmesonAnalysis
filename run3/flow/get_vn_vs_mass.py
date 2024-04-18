'''
Script for fitting D+, D0 and Ds+ invariant-mass spectra
run: python get_vn_vs_mass.py fitConfigFileName.yml centClass inputFileName.root outFileName.root
            [--refFileName][--isMC][--batch]
'''

import sys
import argparse
import ctypes
import numpy as np
import yaml
from ROOT import TLatex, TFile, TCanvas, TH1D, TH1F, TDatabasePDG, AliHFInvMassFitter, AliVertexingHFUtils, AliHFVnVsMassFitter, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, gPad, kBlack, kRed, kFullCircle, kFullSquare # pylint: disable=import-error,no-name-in-module
from flow_analysis_utils import get_centrality_bins
sys.path.append('../../')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, DivideCanvas
from utils.FitUtils import SingleGaus, DoubleGaus, DoublePeakSingleGaus, DoublePeakDoubleGaus

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('fitConfigFileName', metavar='text', default='config_Ds_Fit.yml')
parser.add_argument('centClass', metavar='text', default='')
parser.add_argument('inFileName', metavar='text', default='')
parser.add_argument("--outputdir", "-o", metavar="text",
                    default=".", help="output directory")
parser.add_argument("--suffix", "-s", metavar="text",
                    default="", help="suffix for output files")
parser.add_argument('--refFileName', metavar='text', default='')
parser.add_argument('--isOutOfPlane', action='store_true', default=False)
parser.add_argument('--isInPlane', action='store_true', default=False)
parser.add_argument('--doEP', action='store_true', default=False)
parser.add_argument('--isMC', action='store_true', default=False)
parser.add_argument('--batch', help='suppress video output', action='store_true')
args = parser.parse_args()

cent, _ = get_centrality_bins(args.centClass)

with open(args.fitConfigFileName, 'r', encoding='utf8') as ymlfitConfigFile:
    fitConfig = yaml.load(ymlfitConfigFile, yaml.FullLoader)

gROOT.SetBatch(args.batch)
SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, padtopmargin=0.12, opttitle=1)

ptMins = fitConfig['ptmins']
ptMaxs = fitConfig['ptmaxs']
fixSigma = fitConfig['FixSigma']
fixMean = fitConfig['FixMean']
harmonic = fitConfig['harmonic']
if 'EnableRef' not in fitConfig:
    enableRef = False
else:
    enableRef = fitConfig['EnableRef']
if not isinstance(fixSigma, list):
    fixSigma = [fixSigma for _ in ptMins]
if not isinstance(fixMean, list):
    fixMean = [fixMean for _ in ptMins]
ptLims = list(ptMins)
nPtBins = len(ptMins)
ptLims.append(ptMaxs[-1])

particleName = fitConfig['Dmeson']
inclSecPeak = fitConfig['InclSecPeak']

SgnFuncStr = fitConfig['SgnFunc']
if not isinstance(SgnFuncStr, list):
    SgnFuncStr = [SgnFuncStr] * nPtBins

BkgFuncStr = fitConfig['BkgFunc']
if not isinstance(BkgFuncStr, list):
    BkgFuncStr = [BkgFuncStr] * nPtBins

BkgFuncVnStr = fitConfig['BkgFuncVn']
if not isinstance(BkgFuncVnStr, list):
    BkgFuncVn = [BkgFuncVnStr] * nPtBins

SgnFunc, BkgFunc, BkgFuncVn, degPol = [], [], [], []
for iPt, (bkgStr, sgnStr, bkgVnStr) in enumerate(zip(BkgFuncStr, SgnFuncStr, BkgFuncVnStr)):
    degPol.append(-1)
    if bkgStr == 'kExpo':
        BkgFunc.append(AliHFInvMassFitter.kExpo)
    elif bkgStr == 'kLin':
        BkgFunc.append(AliHFInvMassFitter.kLin)
    elif bkgStr == 'kPol2':
        BkgFunc.append(AliHFInvMassFitter.kPol2)
    elif bkgStr == 'kPol3':
        BkgFunc.append(6)
        degPol[-1] = 3
        if len(ptMins) > 1 and inclSecPeak[iPt] == 1:
            print('ERROR: Pol3 and Pol4 fits work only with one bin if you have the secondary peak! Exit!')
            sys.exit()
    elif bkgStr == 'kPol4':
        BkgFunc.append(6)
        degPol[-1] = 4
        if len(ptMins) > 1 and inclSecPeak[iPt] == 1:
            print('ERROR: Pol3 and Pol4 fits work only with one bin if you have the secondary peak! Exit!')
            sys.exit()
    elif bkgStr == 'kPow':
        BkgFunc.append(AliHFInvMassFitter.kPow)
    elif bkgStr == 'kPowEx':
        BkgFunc.append(AliHFInvMassFitter.kPowEx)
    else:
        print('ERROR: only kExpo, kLin, kPol2, kPol3, kPol4, kPow, and kPowEx background functions supported! Exit')
        sys.exit()

    if bkgVnStr == 'kExpo':
        BkgFuncVn.append(AliHFInvMassFitter.kExpo)
    elif bkgVnStr == 'kLin':
        BkgFuncVn.append(AliHFInvMassFitter.kLin)
    elif bkgVnStr == 'kPol2':
        BkgFuncVn.append(AliHFInvMassFitter.kPol2)
    else:
        print('ERROR: only kExpo, kLin, and kPol2 background functions supported for vn! Exit')
        sys.exit()

    if sgnStr == 'kGaus':
        SgnFunc.append(AliHFInvMassFitter.kGaus)
    elif sgnStr == 'k2Gaus':
        SgnFunc.append(AliHFInvMassFitter.k2Gaus)
    elif sgnStr == 'k2GausSigmaRatioPar':
        SgnFunc.append(AliHFInvMassFitter.k2GausSigmaRatioPar)
    else:
        print('ERROR: only kGaus, k2Gaus and k2GausSigmaRatioPar signal functions supported! Exit!')
        sys.exit()

if particleName == 'Dplus':
    massAxisTit = '#it{M}(K#pi#pi) (GeV/#it{c}^{2})'
elif particleName == 'Ds':
    massAxisTit = '#it{M}(KK#pi) (GeV/#it{c}^{2})'
elif particleName == 'LctopKpi':
    massAxisTit = '#it{M}(pK#pi) (GeV/#it{c}^{2})'
elif particleName == 'LctopK0s':
    massAxisTit = '#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})'
elif particleName == 'Dstar':
    massAxisTit = '#it{M}(K#pi#pi) - #it{M}(K#pi) (GeV/#it{c}^{2})'
elif particleName == 'D0':
    massAxisTit = '#it{M}(K#pi) (GeV/#it{c}^{2})'
else:
    print(f'ERROR: the particle "{particleName}" is not supported! Choose between Dplus, Ds, Dstar, and Lc. Exit!')
    sys.exit()

# load inv-mass histos
infile = TFile.Open(args.inFileName)
if not infile or not infile.IsOpen():
    print(f'ERROR: file "{args.inFileName}" cannot be opened! Exit!')
    sys.exit()
if enableRef:
    infileref = TFile.Open(args.refFileName)
    if not (infileref and infileref.IsOpen()):
        print(f'ERROR: file "{args.refFileName}" cannot be opened! Exit!')
        sys.exit()

hRel, hSig, hMassForRel, hMassForSig  = [], [], [], []
hMass, hMassForFit, hVn, hVnForFit = [], [], [], []
fTotFuncMass, fTotFuncVn = [], []
inclSecPeak = [inclSecPeak] * len(ptMins) if not isinstance(inclSecPeak, list) else inclSecPeak
for iPt, (ptMin, ptMax, secPeak) in enumerate(zip(ptMins, ptMaxs, inclSecPeak)):
    if not args.isMC:
        print(f'loading: cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_mass_cent{cent}_pt{ptMin}_{ptMax}')
        hMass.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_mass_cent{cent}_pt{ptMin}_{ptMax}'))
        hMass[iPt].SetDirectory(0)
        if args.doEP:
            print(f'loading: cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_vn_ep_pt{ptMin}_{ptMax}')
            hVn.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_vn_ep_pt{ptMin}_{ptMax}'))
        else:
            print(f'loading: cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_vn_sp_pt{ptMin}_{ptMax}')
            hVn.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_vn_sp_pt{ptMin}_{ptMax}'))
        hVn[iPt].SetDirectory(0)
        if enableRef:
            hRel.append(infileref.Get(f'hVarReflMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
            hSig.append(infileref.Get(f'hFDMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
            hSig[iPt].Add(infileref.Get(f'hPromptMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
            hRel[iPt].SetDirectory(0)
            hSig[iPt].SetDirectory(0)
            hRel[iPt].Sumw2()
            hSig[iPt].Sumw2()
    else:
        hMass.append(infile.Get(f'hPromptMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
        hMass[iPt].Add(infile.Get(f'hFDMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
        if secPeak:
            hMass[iPt].Add(infile.Get(f'hPromptSecPeakMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
            hMass[iPt].Add(infile.Get(f'hFDSecPeakMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
        hMass[iPt].SetDirectory(0)
    hMass[iPt].Sumw2()
    SetObjectStyle(hMass[iPt], color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hVn[iPt], color=kBlack, markerstyle=kFullCircle)
infile.Close()

if fitConfig['FixSigmaRatio']:
    # load sigma of first gaussian
    infileSigma = TFile.Open(fitConfig['SigmaRatioFile'])
    if not infileSigma:
        print(f'ERROR: file "{infileSigma}" cannot be opened! Exit!')
        sys.exit()
    hSigmaToFix = infileSigma.Get('hRawYieldsSigma')
    hSigmaToFix.SetDirectory(0)
    if hSigmaToFix.GetNbinsX() != nPtBins:
        print('WARNING: Different number of bins for this analysis and histo for fix sigma')
    infileSigma.Close()
    # load sigma of second gaussian
    infileSigma2 = TFile.Open(fitConfig['SigmaRatioFile'])
    if not infileSigma2:
        print(f'ERROR: file "{infileSigma2}" cannot be opened! Exit!')
        sys.exit()
    hSigmaToFix2 = infileSigma2.Get('hRawYieldsSigma2')
    hSigmaToFix2.SetDirectory(0)
    if hSigmaToFix2.GetNbinsX() != nPtBins:
        print('WARNING: Different number of bins for this analysis and histo for fix sigma')
    infileSigma2.Close()

hMeanToFix = None
if sum(fixMean) > 0:
    infileMean = TFile.Open(fitConfig['MeanFile'])
    if not infileMean:
        print(f'ERROR: file "{infileMean}" cannot be opened! Exit!')
        sys.exit()
    hMeanToFix = infileMean.Get('hRawYieldsMean')
    hMeanToFix.SetDirectory(0)
    if hMeanToFix.GetNbinsX() != nPtBins:
        print('WARNING: Different number of bins for this analysis and histo for fix mean')
    infileMean.Close()

hSigmaFirstPeakMC = None
hSigmaToFixSecPeak = None
infileSigmaSecPeak = None
if fitConfig['SigmaFileSecPeak']:
    infileSigmaSecPeak = TFile.Open(fitConfig['SigmaFileSecPeak'])
if fitConfig['FixSigmaToFirstPeak'] and not infileSigmaSecPeak:
    print(f'ERROR: file "{fitConfig["SigmaFileSecPeak"]}" cannot be opened! Exit!')
    sys.exit()
if infileSigmaSecPeak:
    hSigmaFirstPeakMC = infileSigmaSecPeak.Get("hRawYieldsSigma")
    hSigmaToFixSecPeak = infileSigmaSecPeak.Get("hRawYieldsSigmaSecondPeak")
    hSigmaFirstPeakMC.SetDirectory(0)
    hSigmaToFixSecPeak.SetDirectory(0)
    if hSigmaFirstPeakMC.GetNbinsX() != nPtBins or hSigmaToFixSecPeak.GetNbinsX() != nPtBins:
        print('WARNING: Different number of bins for this analysis and histoss for fix mean')
    infileSigmaSecPeak.Close()

ptBinsArr = np.asarray(ptLims, 'd')
ptTit = '#it{p}_{T} (GeV/#it{c})'

hSigmaSimFit = TH1D('hSigmaSimFit', f';{ptTit};#sigma', nPtBins, ptBinsArr)
hMeanSimFit = TH1D('hMeanSimFit', f';{ptTit};mean', nPtBins, ptBinsArr)
hMeanSecPeakFitMass = TH1D('hMeanSecondPeakFitMass', f';{ptTit};mean second peak mass fit', nPtBins, ptBinsArr)
hMeanSecPeakFitVn = TH1D('hMeanSecondPeakFitVn', f';{ptTit};mean second peak vn fit', nPtBins, ptBinsArr)
hSigmaSecPeakFitMass = TH1D('hSigmaSecondPeakFitMass', f';{ptTit};width second peak mass fit', nPtBins, ptBinsArr)
hSigmaSecPeakFitVn = TH1D('hSigmaSecondPeakFitVn', f';{ptTit};width second peak vn fit', nPtBins, ptBinsArr)
hRawYieldsSimFit = TH1D('hRawYieldsSimFit', f';{ptTit};raw yield', nPtBins, ptBinsArr)
hRawYieldsSecPeakSimFit = TH1D('hRawYieldsSecondPeakSimFit', f';{ptTit};raw yield second peak', nPtBins, ptBinsArr)
hRawYieldsSignificanceSimFit = TH1D('hRawYieldsSignificanceSimFit', f';{ptTit};significance', nPtBins, ptBinsArr)
hRawYieldsSoverBSimFit = TH1D('hRawYieldsSoverBSimFit', f';{ptTit};S/B', nPtBins, ptBinsArr)
hRedChi2SimFit = TH1D('hRedChi2SimFit', f';{ptTit};#chi^{{2}}/#it{{ndf}}', nPtBins, ptBinsArr)
hProbSimFit = TH1D('hProbSimFit', f';{ptTit};prob', nPtBins, ptBinsArr)
hRedChi2SBVnPrefit = TH1D('hRedChi2SBVnPrefit', f';{ptTit};#chi^{{2}}/#it{{ndf}}', nPtBins, ptBinsArr)
hProbSBVnPrefit = TH1D('hProbSBVnPrefit', f';{ptTit};prob', nPtBins, ptBinsArr)
gvnSimFit = TGraphAsymmErrors(1)
gvnSimFit.SetName('gvnSimFit')
gvnSimFitSecPeak = TGraphAsymmErrors(1)
gvnSimFitSecPeak.SetName('gvnSimFitSecPeak')

SetObjectStyle(hSigmaSimFit, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hMeanSimFit, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hMeanSecPeakFitMass, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hSigmaSecPeakFitMass, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hMeanSecPeakFitVn, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hSigmaSecPeakFitVn, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSimFit, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSecPeakSimFit, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSignificanceSimFit, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSoverBSimFit, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRedChi2SimFit, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hProbSimFit, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRedChi2SBVnPrefit, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hProbSBVnPrefit, color=kRed, markerstyle=kFullSquare)

# additional S, B, S/B, and significance histos for different Nsigma values (filled only in case of data)
nSigma4SandB = [1.5, 1.75, 2., 2.25, 2.5, 2.75]
hRawYieldsSignalDiffSigma, hRawYieldsBkgDiffSigma, hRawYieldsSoverBDiffSigma, \
    hRawYieldsSignifDiffSigma = ([] for _ in range(4))

for iS, sigma in enumerate(nSigma4SandB):
    hRawYieldsSignalDiffSigma.append(TH1D(f'hRawYieldsSignal_{sigma:0.2f}sigma',
                                          f';{ptTit};Signal ({sigma:0.2f}#sigma)', nPtBins, ptBinsArr))
    hRawYieldsBkgDiffSigma.append(TH1D(f'hRawYieldsBkg_{sigma:0.2f}sigma',
                                       f';{ptTit};Background ({sigma:0.2f}#sigma)', nPtBins, ptBinsArr))
    hRawYieldsSoverBDiffSigma.append(TH1D(f'hRawYieldsSoverB_{sigma:0.2f}sigma',
                                          f';{ptTit};S/B ({sigma:0.2f}#sigma)', nPtBins, ptBinsArr))
    hRawYieldsSignifDiffSigma.append(TH1D(f'hRawYieldsSignif_{sigma:0.2f}sigma',
                                          f';{ptTit};significance ({sigma:0.2f}#sigma)', nPtBins, ptBinsArr))
    SetObjectStyle(hRawYieldsSignalDiffSigma[iS], color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRawYieldsBkgDiffSigma[iS], color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRawYieldsSoverBDiffSigma[iS], color=kBlack, markerstyle=kFullCircle)
    SetObjectStyle(hRawYieldsSignifDiffSigma[iS], color=kBlack, markerstyle=kFullCircle)

# fit histos
massDplus = TDatabasePDG.Instance().GetParticle(411).Mass()
massDs = TDatabasePDG.Instance().GetParticle(431).Mass()
massLc = TDatabasePDG.Instance().GetParticle(4122).Mass()
massDstar = TDatabasePDG.Instance().GetParticle(413).Mass() - TDatabasePDG.Instance().GetParticle(421).Mass()
massD0 = TDatabasePDG.Instance().GetParticle(421).Mass()

if particleName == 'Dplus':
    massForFit=massDplus
elif particleName == 'Ds':
    massForFit = massDs
elif particleName == 'Dstar':
    massForFit = massDstar
elif particleName == 'D0':
    massForFit = massD0
else:
    massForFit = massLc

canvSizes = [1920, 1080]
if nPtBins == 1:
    canvSizes = [500, 500]

cSimFit = []
for i in range(nPtBins):
    ptLow = ptMins[i]
    ptHigh = ptMaxs[i]
    cSimFit.append(TCanvas(f'cSimFit_Pt{ptLow}_{ptHigh}', f'cSimFit_Pt{ptLow}_{ptHigh}', 400, 900))

massFitter, vnFitter = [], []

rebins = fitConfig['Rebin']
if not isinstance(rebins, list):
    rebins = [rebins] * len(ptMins)

massMins = fitConfig['MassMin']
if not isinstance(massMins, list):
    massMins = [massMins] * len(ptMins)

massMaxs = fitConfig['MassMax']
if not isinstance(massMaxs, list):
    massMaxs = [massMaxs] * len(ptMins)

for iPt, (hM, hV, ptMin, ptMax, reb, sgnEnum, bkgEnum, bkgVnEnum, secPeak, massMin, massMax) in enumerate(
        zip(hMass, hVn, ptMins, ptMaxs, rebins, SgnFunc, BkgFunc, BkgFuncVn, inclSecPeak, massMins, massMaxs)):
    iCanv = iPt
    hMassForFit.append(TH1F())
    hVnForFit.append(TH1F())
    AliVertexingHFUtils.RebinHisto(hM, reb).Copy(hMassForFit[iPt]) #to cast TH1D to TH1F
    hMassForFit[iPt].SetDirectory(0)
    xbins = np.asarray(hV.GetXaxis().GetXbins())
    hDummy = TH1F('hDummy', '', len(xbins)-1, xbins)
    for iBin in range(1, hV.GetNbinsX()+1):
        hDummy.SetBinContent(iBin, hV.GetBinContent(iBin))
        hDummy.SetBinError(iBin, hV.GetBinError(iBin))
    hVnForFit[iPt] = hDummy
    hVnForFit[iPt].SetDirectory(0)
    hVnForFit[iPt].GetXaxis().SetTitle(massAxisTit)
    hVnForFit[iPt].GetYaxis().SetTitle(f'#it{{v}}{harmonic}')
    binWidth = hMassForFit[iPt].GetBinWidth(1)
    hMassForFit[iPt].SetTitle((f'{ptMin:0.1f} < #it{{p}}_{{T}} < {ptMax:0.1f} GeV/#it{{c}};{massAxisTit};'
                               f'Counts per {binWidth*1000:.0f} MeV/#it{{c}}^{{2}}'))
    hMassForFit[iPt].SetName(f'MassForFit{iPt}')
    if not args.isMC and enableRef:
        hMassForRel.append(TH1F())
        hMassForSig.append(TH1F())
        AliVertexingHFUtils.RebinHisto(hRel[iPt], reb).Copy(hMassForRel[iPt])
        AliVertexingHFUtils.RebinHisto(hSig[iPt], reb).Copy(hMassForSig[iPt])
    if nPtBins < 15:
        markerSize = 1.
    else:
        markerSize = 0.5
    SetObjectStyle(hMassForFit[iPt], color=kBlack, markerstyle=kFullCircle, markersize=markerSize)
    SetObjectStyle(hVnForFit[iPt], color=kBlack, markerstyle=kFullCircle, markersize=markerSize)

    print(f'Fitting {ptMin} - {ptMax} GeV/c')
    vnFitter.append(AliHFVnVsMassFitter(hMassForFit[iPt], hVnForFit[iPt], massMin, 2.15, bkgEnum, sgnEnum, bkgVnEnum))
    vnFitter[iPt].SetHarmonic(harmonic)

    # set the parameters for the fit
    # Mean
    vnFitter[iPt].SetInitialGaussianMean(massForFit, 1)
    if fixMean[iPt]:
        vnFitter[iPt].FixMeanFromMassFit(hMeanToFix.GetBinContent(iPt+1))
    # Sigma
    vnFitter[iPt].SetInitialGaussianSigma(fitConfig['Sigma'][iPt], 1)
    if fixSigma[iPt]:
        vnFitter[iPt].FixSigmaFromMassFit()
    # Second peak (Ds specific)
    if secPeak and particleName == 'Ds':
        vnFitter[iPt].IncludeSecondGausPeak(massDplus, False, fitConfig['SigmaSecPeak'][iPt], False, 1)
        if fixSigma[iPt]:
            vnFitter[iPt].FixSigma2GausFromMassFit()
    vnFitter[iPt].FixFrac2GausFromMassFit()
    # TODO: Add the reflections
    vnFitter[iPt].SimultaneusFit(False)

    # collect the results
    vn = vnFitter[iPt].GetVn()
    vnUnc = vnFitter[iPt].GetVnUncertainty()
    mean = vnFitter[iPt].GetMean()
    meanUnc = vnFitter[iPt].GetMeanUncertainty()
    sigma = vnFitter[iPt].GetSigma()
    sigmaUnc = vnFitter[iPt].GetSigmaUncertainty()
    ry = vnFitter[iPt].GetRawYield()
    ryUnc = vnFitter[iPt].GetRawYieldUncertainty()
    chi2 = vnFitter[iPt].GetReducedChiSquare()
    prob = vnFitter[iPt].GetFitProbability()
    fTotFuncMass.append(vnFitter[iPt].GetMassTotFitFunc())
    # _________________________
    # fTotFuncVn parameters:
    # 0: BkgInt
    # 1: BkgSlope
    # 2: SgnInt
    # 3: Mean
    # 4: Sigma
    # 5: SecPeakInt
    # 6: SecPeakMean
    # 7: SecPeakSigma
    # 8: ConstVnBkg
    # 9: SlopeVnBkg
    # 10: v2Sgn
    # 11: v2SecPeak
    # _________________________
    fTotFuncVn.append(vnFitter[iPt].GetVnVsMassTotFitFunc())

    hSigmaSimFit.SetBinContent(iPt+1, sigma)
    hSigmaSimFit.SetBinError(iPt+1, sigmaUnc)
    hMeanSimFit.SetBinContent(iPt+1, mean)
    hMeanSimFit.SetBinError(iPt+1, meanUnc)
    hRedChi2SimFit.SetBinContent(iPt+1, chi2)
    hRedChi2SimFit.SetBinError(iPt+1, 1.e-20)
    hProbSimFit.SetBinContent(iPt+1, prob)
    hProbSimFit.SetBinError(iPt+1, 1.e-20)
    hRawYieldsSimFit.SetBinContent(iPt+1, ry)
    hRawYieldsSimFit.SetBinError(iPt+1, ryUnc)
    gvnSimFit.SetPoint(iPt, (ptMin+ptMax)/2, vn)
    gvnSimFit.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, vnUnc, vnUnc)

    if secPeak:
        secPeakMeanMass = fTotFuncMass[iPt].GetParameter(7)
        secPeakMeanMassUnc = fTotFuncMass[iPt].GetParError(7)
        secPeakSigmaMass = fTotFuncMass[iPt].GetParameter(8)
        secPeakSigmaMassUnc = fTotFuncMass[iPt].GetParError(8)
        secPeakMeanVn = fTotFuncVn[iPt].GetParameter(7)
        secPeakMeanVnUnc = fTotFuncVn[iPt].GetParError(7)
        secPeakSigmaVn = fTotFuncVn[iPt].GetParameter(8)
        secPeakSigmaVnUnc = fTotFuncVn[iPt].GetParError(8)
        vnSecPeak = fTotFuncVn[iPt].GetParameter(11)
        vnSecPeakUnc = fTotFuncVn[iPt].GetParError(11)

        hMeanSecPeakFitMass.SetBinContent(iPt+1, secPeakMeanMass)
        hMeanSecPeakFitMass.SetBinError(iPt+1, secPeakMeanMassUnc)
        hSigmaSecPeakFitMass.SetBinContent(iPt+1, secPeakSigmaMass)
        hSigmaSecPeakFitMass.SetBinError(iPt+1, secPeakSigmaMassUnc)
        hMeanSecPeakFitVn.SetBinContent(iPt+1, secPeakMeanVn)
        hMeanSecPeakFitVn.SetBinError(iPt+1, secPeakMeanVnUnc)
        hSigmaSecPeakFitVn.SetBinContent(iPt+1, secPeakSigmaVn)
        hSigmaSecPeakFitVn.SetBinError(iPt+1, secPeakSigmaVnUnc)
        gvnSimFitSecPeak.SetPoint(iPt, (ptMin+ptMax)/2, vnSecPeak)
        gvnSimFitSecPeak.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, vnSecPeakUnc, vnSecPeakUnc)

    if vn != 0:
        cSimFit[iPt].cd()
        vnFitter[iPt].DrawHere(gPad)
        latex = TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.04)
        latex.DrawLatex(0.5, 0.15, f'{ptMin:.1f} < #it{{p}}_{{T}} < {ptMax:.1f} GeV/#it{{c}}')
    cSimFit[iCanv].Modified()
    cSimFit[iCanv].Update()

#save output histos
print(f'Saving output to {args.outputdir}')
for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
    if iPt == 0:
        suffix = '('
    elif iPt == nPtBins-1:
        suffix = ')'
    else:
        suffix = ''
    cSimFit[iPt].SaveAs(f'{args.outputdir}/SimFit{args.suffix}_{particleName}.pdf{suffix}')
outfile_name = f'{args.outputdir}/raw_yields{args.suffix}.root'
outFile = TFile(outfile_name, 'recreate')
for canv in cSimFit:
    canv.Write()
for hist in hMass:
    hist.Write()
for hist in hVn:
    hist.Write()
for ipt, (ptmin, ptmax) in enumerate(zip(ptMins, ptMaxs)):
    fTotFuncMass[ipt].Write(f'fTotFuncMass_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
for ipt, (ptmin, ptmax) in enumerate(zip(ptMins, ptMaxs)):
    fTotFuncVn[ipt].Write(f'fTotFuncVn_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
hSigmaSimFit.Write()
hMeanSimFit.Write()
hMeanSecPeakFitMass.Write()
hMeanSecPeakFitVn.Write()
hSigmaSecPeakFitMass.Write()
hSigmaSecPeakFitVn.Write()
hRawYieldsSimFit.Write()
hRawYieldsSecPeakSimFit.Write()
hRawYieldsSignificanceSimFit.Write()
hRawYieldsSoverBSimFit.Write()
hRedChi2SimFit.Write()
hProbSimFit.Write()
hRedChi2SBVnPrefit.Write()
hProbSBVnPrefit.Write()
gvnSimFit.Write()
gvnSimFitSecPeak.Write()

outFile.Close()

if not args.batch:
    input('Press enter to exit')