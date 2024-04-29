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
from ROOT import gROOT, gPad, kBlack, kRed, kAzure, kGray, kFullCircle, kFullSquare, kOpenCircle # pylint: disable=import-error,no-name-in-module
from flow_analysis_utils import get_centrality_bins, get_vnfitter_results, get_ep_vn # pylint: disable=import-error,no-name-in-module
sys.path.append('../../')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, DivideCanvas
from utils.FitUtils import SingleGaus, DoubleGaus, DoublePeakSingleGaus, DoublePeakDoubleGaus

def get_vn_vs_mass(fitConfigFileName, centClass, inFileName,
                   outputdir, suffix, refFileName, vn_method, batch):


    with open(fitConfigFileName, 'r', encoding='utf8') as ymlfitConfigFile:
        fitConfig = yaml.load(ymlfitConfigFile, yaml.FullLoader)
    
    gROOT.SetBatch(batch)
    SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, padtopmargin=0.12, opttitle=1)
    cent, _ = get_centrality_bins(centClass)

    # read global configuration
    ptMins = fitConfig['ptmins']
    ptMaxs = fitConfig['ptmaxs']
    ptLims = list(ptMins)
    nPtBins = len(ptMins)
    ptLims.append(ptMaxs[-1])
    ptBinsArr = np.asarray(ptLims, 'd')
    ptTit = '#it{p}_{T} (GeV/#it{c})'
    fixSigma = fitConfig['FixSigma']
    fixMean = fitConfig['FixMean']
    harmonic = fitConfig['harmonic']
    particleName = fitConfig['Dmeson']
    inclSecPeak = fitConfig['InclSecPeak']
    rebins = fitConfig['Rebin']
    if not isinstance(rebins, list):
        rebins = [rebins] * len(ptMins)
    massMins = fitConfig['MassMin']
    if not isinstance(massMins, list):
        massMins = [massMins] * len(ptMins)
    massMaxs = fitConfig['MassMax']
    if not isinstance(massMaxs, list):
        massMaxs = [massMaxs] * len(ptMins)
    
    # read fit configuration
    if not isinstance(fixSigma, list):
        fixSigma = [fixSigma for _ in ptMins]
    if not isinstance(fixMean, list):
        fixMean = [fixMean for _ in ptMins]
    SgnFuncStr = fitConfig['SgnFunc']
    if not isinstance(SgnFuncStr, list):
        SgnFuncStr = [SgnFuncStr] * nPtBins
    BkgFuncStr = fitConfig['BkgFunc']
    if not isinstance(BkgFuncStr, list):
        BkgFuncStr = [BkgFuncStr] * nPtBins
    BkgFuncVnStr = fitConfig['BkgFuncVn']
    if not isinstance(BkgFuncVnStr, list):
        BkgFuncVn = [BkgFuncVnStr] * nPtBins

    # sanity check of fit configuration
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

    # set particle configuration
    massDplus = TDatabasePDG.Instance().GetParticle(411).Mass()
    massDs = TDatabasePDG.Instance().GetParticle(431).Mass()
    massLc = TDatabasePDG.Instance().GetParticle(4122).Mass()
    massDstar = TDatabasePDG.Instance().GetParticle(413).Mass() - TDatabasePDG.Instance().GetParticle(421).Mass()
    massD0 = TDatabasePDG.Instance().GetParticle(421).Mass()
    if particleName == 'Dplus':
        massAxisTit = '#it{M}(K#pi#pi) (GeV/#it{c}^{2})'
        massForFit=massDplus
    elif particleName == 'Ds':
        massAxisTit = '#it{M}(KK#pi) (GeV/#it{c}^{2})'
        massForFit = massDs
    elif particleName == 'LctopKpi':
        massAxisTit = '#it{M}(pK#pi) (GeV/#it{c}^{2})'
        massForFit = massLc
    elif particleName == 'LctopK0s':
        massAxisTit = '#it{M}(pK^{0}_{s}) (GeV/#it{c}^{2})'
        massForFit = massLc
    elif particleName == 'Dstar':
        massAxisTit = '#it{M}(K#pi#pi) - #it{M}(K#pi) (GeV/#it{c}^{2})'
        massForFit = massDstar
    elif particleName == 'D0':
        massAxisTit = '#it{M}(K#pi) (GeV/#it{c}^{2})'
        massForFit = massD0
    else:
        print(f'ERROR: the particle "{particleName}" is not supported! Choose between Dplus, Ds, Dstar, and Lc. Exit!')
        sys.exit()

    # load histos
    infile = TFile.Open(inFileName)
    if not infile or not infile.IsOpen():
        print(f'ERROR: file "{inFileName}" cannot be opened! Exit!')
        sys.exit()
    hRel, hSig, hMassForRel, hMassForSig  = [], [], [], []
    hMass, hMassForFit, hVn, hVnForFit = [], [], [], []
    hMassIns, hMassOuts, hMassInsForFit, hMassOutsForFit = [], [], [], []
    fTotFuncMass, fTotFuncVn = [], []
    inclSecPeak = [inclSecPeak] * len(ptMins) if not isinstance(inclSecPeak, list) else inclSecPeak
    for iPt, (ptMin, ptMax, secPeak) in enumerate(zip(ptMins, ptMaxs, inclSecPeak)):
        if not vn_method == 'sp' and not vn_method == 'ep':
            print(f'loading: cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_mass_cent{cent}_pt{ptMin}_{ptMax}')
            hMassIns.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_mass_inplane_cent{cent}_pt{ptMin}_{ptMax}'))
            hMassOuts.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_mass_outplane_cent{cent}_pt{ptMin}_{ptMax}'))
            hMassIns[iPt].SetDirectory(0)
            hMassOuts[iPt].SetDirectory(0)
        else:
            print(f'loading: cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_vn_{vn_method}_pt{ptMin}_{ptMax}')
            hMass.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_mass_cent{cent}_pt{ptMin}_{ptMax}'))
            hVn.append(infile.Get(f'cent_bins{cent}/pt_bins{ptMin}_{ptMax}/hist_vn_{vn_method}_pt{ptMin}_{ptMax}'))
            hVn[iPt].SetDirectory(0)
            hMass[iPt].SetDirectory(0)
            SetObjectStyle(hMass[iPt], color=kBlack, markerstyle=kFullCircle)
            SetObjectStyle(hVn[iPt], color=kBlack, markerstyle=kFullCircle)
    infile.Close()

    hSigmaToFix = None
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

    # create histos for fit results
    if vn_method == 'sp' or vn_method == 'ep':
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
        SetObjectStyle(gvnSimFit, color=kBlack, markerstyle=kFullCircle)
        SetObjectStyle(gvnSimFitSecPeak, color=kBlack, markerstyle=kFullCircle)
    else:
        hRawYieldsIn = TH1D('hRawYieldsIn', f';{ptTit};raw yield in-plane', nPtBins, ptBinsArr)
        hRawYieldsOut = TH1D('hRawYieldsOut', f';{ptTit};raw yield out-of-plane', nPtBins, ptBinsArr)
        hSigmaIn = TH1D('hSigmaIn', f';{ptTit};#sigma in-plane', nPtBins, ptBinsArr)
        hSigmaOut = TH1D('hSigmaOut', f';{ptTit};#sigma out-of-plane', nPtBins, ptBinsArr)
        hMeanIn = TH1D('hMeanIn', f';{ptTit};mean in-plane', nPtBins, ptBinsArr)
        hMeanOut = TH1D('hMeanOut', f';{ptTit};mean out-of-plane', nPtBins, ptBinsArr)
        hRedChi2In = TH1D('hRedChi2In', f';{ptTit};#chi^{{2}}/#it{{ndf}} in-plane', nPtBins, ptBinsArr)
        hRedChi2Out = TH1D('hRedChi2Out', f';{ptTit};#chi^{{2}}/#it{{ndf}} out-of-plane', nPtBins, ptBinsArr)
        hProbIn = TH1D('hProbIn', f';{ptTit};prob in-plane', nPtBins, ptBinsArr)
        hProbOut = TH1D('hProbOut', f';{ptTit};prob out-of-plane', nPtBins, ptBinsArr)
        hRawYieldsSignificanceIn = TH1D('hRawYieldsSignificanceIn', f';{ptTit};significance in-plane',
                                        nPtBins, ptBinsArr)
        hRawYieldsSignificanceOut = TH1D('hRawYieldsSignificanceOut', f';{ptTit};significance out-of-plane',
                                         nPtBins, ptBinsArr)
        hRawYieldsSoverBIn = TH1D('hRawYieldsSoverBIn', f';{ptTit};S/B in-plane', nPtBins, ptBinsArr)
        hRawYieldsSoverBOut = TH1D('hRawYieldsSoverBOut', f';{ptTit};S/B out-of-plane', nPtBins, ptBinsArr)
        hSigmaSecPeakFitIn = TH1D('hSigmaSecPeakFitIn', f';{ptTit};width second peak in-plane', nPtBins, ptBinsArr)
        hSigmaSecPeakFitOut = TH1D('hSigmaSecPeakFitOut', f';{ptTit};width second peak out-of-plane',
                                   nPtBins, ptBinsArr)
        gvnSimFit = TGraphAsymmErrors(1)
        gvnSimFit.SetName('gvnSimFit')
        gvnSimFitSecPeak = TGraphAsymmErrors(1)
        gvnSimFitSecPeak.SetName('gvnSimFitSecPeak')

        SetObjectStyle(hRawYieldsIn, color=kRed, markerstyle=kFullCircle)
        SetObjectStyle(hRawYieldsOut, color=kAzure, markerstyle=kOpenCircle)
        SetObjectStyle(hSigmaIn, color=kRed, markerstyle=kFullCircle)
        SetObjectStyle(hSigmaOut, color=kAzure, markerstyle=kOpenCircle)
        SetObjectStyle(hMeanIn, color=kRed, markerstyle=kFullCircle)
        SetObjectStyle(hMeanOut, color=kAzure, markerstyle=kOpenCircle)
        SetObjectStyle(hRedChi2In, color=kRed, markerstyle=kFullCircle)
        SetObjectStyle(hRedChi2Out, color=kAzure, markerstyle=kOpenCircle)
        SetObjectStyle(hProbIn, color=kRed, markerstyle=kFullCircle)
        SetObjectStyle(hProbOut, color=kAzure, markerstyle=kOpenCircle)
        SetObjectStyle(hRawYieldsSignificanceIn, color=kRed, markerstyle=kFullCircle)
        SetObjectStyle(hRawYieldsSignificanceOut, color=kAzure, markerstyle=kOpenCircle)
        SetObjectStyle(hRawYieldsSoverBIn, color=kRed, markerstyle=kFullCircle)
        SetObjectStyle(hRawYieldsSoverBOut, color=kAzure, markerstyle=kOpenCircle)
        SetObjectStyle(hSigmaSecPeakFitIn, color=kRed, markerstyle=kFullCircle)
        SetObjectStyle(hSigmaSecPeakFitOut, color=kAzure, markerstyle=kOpenCircle)
        SetObjectStyle(gvnSimFit, color=kBlack, markerstyle=kFullCircle)
        SetObjectStyle(gvnSimFitSecPeak, color=kBlack, markerstyle=kOpenCircle)

    # create canvases
    canvSizes = [1920, 1080]
    nMaxCanvases = 10 # do not put more than 20 bins per canvas to make them visible
    nCanvases = 1
    if nPtBins == 1:
        canvSizes = [500, 500]
    latex = TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.04)

    if vn_method == 'sp' or vn_method == 'ep':
        cSimFit = []
        for i in range(nPtBins):
            ptLow = ptMins[i]
            ptHigh = ptMaxs[i]
            cSimFit.append(TCanvas(f'cSimFit_Pt{ptLow}_{ptHigh}', f'cSimFit_Pt{ptLow}_{ptHigh}', 400, 900))
    else:
        cMass, cResiduals = [], []
        cMass = TCanvas('cMass', 'cMass', canvSizes[0], canvSizes[1])
        nPads = nPtBins if nCanvases == 1 else nMaxCanvases
        DivideCanvas(cMass, nPads)

    #_____________________________________________________
    # Vn estimation with Scalar Product / Event Plane
    if vn_method == 'sp' or vn_method == 'ep':
        massFitter, vnFitter = [], []
    else:
        massFitterIns, massFitterOuts = [], []
    if vn_method == 'sp' or vn_method == 'ep':
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
            SetObjectStyle(hMassForFit[iPt], color=kBlack, markerstyle=kFullCircle, markersize=1)
            SetObjectStyle(hVnForFit[iPt], color=kBlack, markerstyle=kOpenCircle, markersize=1)

            print(f'Fitting {ptMin} - {ptMax} GeV/c')
            vnFitter.append(AliHFVnVsMassFitter(hMassForFit[iPt], hVnForFit[iPt], massMin, 2.15, bkgEnum, sgnEnum, bkgVnEnum))
            vnFitter[iPt].SetHarmonic(harmonic)

            #_____________________________________________________
            # set the parameters for the fit
            # Mean
            vnFitter[iPt].SetInitialGaussianMean(massForFit, 1)
            if fixMean[iPt]:
                vnFitter[iPt].FixMeanFromMassFit()
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
            # TODO: Add reflections for D0
            # collect fit results
            vnFitter[iPt].SimultaneusFit(False)
            vnResults = get_vnfitter_results(vnFitter[iPt], secPeak)
            fTotFuncMass.append(vnResults['fTotFuncMass'])
            fTotFuncVn.append(vnResults['fTotFuncVn'])

            hSigmaSimFit.SetBinContent(iPt+1, vnResults['sigma'])
            hSigmaSimFit.SetBinError(iPt+1, vnResults['sigmaUnc'])
            hMeanSimFit.SetBinContent(iPt+1, vnResults['mean'])
            hMeanSimFit.SetBinError(iPt+1, vnResults['meanUnc'])
            hRedChi2SimFit.SetBinContent(iPt+1, vnResults['chi2'])
            hRedChi2SimFit.SetBinError(iPt+1, 1.e-20)
            hProbSimFit.SetBinContent(iPt+1, vnResults['prob'])
            hProbSimFit.SetBinError(iPt+1, 1.e-20)
            hRawYieldsSimFit.SetBinContent(iPt+1, vnResults['ry'])
            hRawYieldsSimFit.SetBinError(iPt+1, vnResults['ryUnc'])
            gvnSimFit.SetPoint(iPt, (ptMin+ptMax)/2, vnResults['vn'])
            gvnSimFit.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, vnResults['vnUnc'], vnResults['vnUnc'])

            if secPeak:
                hMeanSecPeakFitMass.SetBinContent(iPt+1, vnResults['secPeakMeanMass'])
                hMeanSecPeakFitMass.SetBinError(iPt+1, vnResults['secPeakMeanMassUnc'])
                hSigmaSecPeakFitMass.SetBinContent(iPt+1, vnResults['secPeakSigmaMass'])
                hSigmaSecPeakFitMass.SetBinError(iPt+1, vnResults['secPeakSigmaMassUnc'])
                hMeanSecPeakFitVn.SetBinContent(iPt+1, vnResults['secPeakMeanVn'])
                hMeanSecPeakFitVn.SetBinError(iPt+1, vnResults['secPeakMeanVnUnc'])
                hSigmaSecPeakFitVn.SetBinContent(iPt+1, vnResults['secPeakSigmaVn'])
                hSigmaSecPeakFitVn.SetBinError(iPt+1, vnResults['secPeakSigmaVnUnc'])
                gvnSimFitSecPeak.SetPoint(iPt, (ptMin+ptMax)/2, vnResults['vnSecPeak'])
                gvnSimFitSecPeak.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2,
                                               vnResults['vnSecPeakUnc'],
                                               vnResults['vnSecPeakUnc'])

            if vnResults['vn'] != 0:
                cSimFit[iPt].cd()
                vnFitter[iPt].DrawHere(gPad)
                latex.DrawLatex(0.5, 0.15, f'{ptMin:.1f} < #it{{p}}_{{T}} < {ptMax:.1f} GeV/#it{{c}}')
                cSimFit[iCanv].Modified()
                cSimFit[iCanv].Update()

    #_____________________________________________________
    # Mass fit
    else:
        for iPt, (hMassIn, hMassOut,\
                  ptMin, ptMax,\
                  reb, sgnEnum, bkgEnum, bkgVnEnum,\
                  secPeak, massMin, massMax) in enumerate(zip(hMassIns, hMassOuts,
                                                              ptMins, ptMaxs,
                                                              rebins, SgnFunc, BkgFunc,
                                                              BkgFuncVn, inclSecPeak,
                                                              massMins, massMaxs)):
            iCanv = iPt
            hMassInsForFit.append(TH1F())
            hMassOutsForFit.append(TH1F())
            AliVertexingHFUtils.RebinHisto(hMassIn, reb).Copy(hMassInsForFit[iPt])
            AliVertexingHFUtils.RebinHisto(hMassOut, reb).Copy(hMassOutsForFit[iPt])
            hMassInsForFit[iPt].SetDirectory(0)
            hMassOutsForFit[iPt].SetDirectory(0)
            binWidth = hMassInsForFit[iPt].GetBinWidth(1)
            hMassInsForFit[iPt].SetTitle((f'{ptMin:0.1f} < #it{{p}}_{{T}} < {ptMax:0.1f} GeV/#it{{c}};{massAxisTit};'
                                            f'Counts per {binWidth*1000:.0f} MeV/#it{{c}}^{{2}}'))
            hMassInsForFit[iPt].SetName(f'MassInForFit{iPt}')
            hMassOutsForFit[iPt].SetTitle((f'{ptMin:0.1f} < #it{{p}}_{{T}} < {ptMax:0.1f} GeV/#it{{c}};{massAxisTit};'
                                            f'Counts per {binWidth*1000:.0f} MeV/#it{{c}}^{{2}}'))
            hMassOutsForFit[iPt].SetName(f'MassOutForFit{iPt}')
            SetObjectStyle(hMassInsForFit[iPt], color=kRed-3, markerstyle=kFullCircle, markersize=0.8)
            SetObjectStyle(hMassOutsForFit[iPt], color=kAzure-3, markerstyle=kOpenCircle, markersize=0.8)

            print(f'Fitting {ptMin} - {ptMax} GeV/c')
            massFitterIns.append(AliHFInvMassFitter(hMassInsForFit[iPt],  massMin, massMax, bkgEnum, sgnEnum))
            massFitterOuts.append(AliHFInvMassFitter(hMassOutsForFit[iPt],  massMin, massMax, bkgEnum, sgnEnum))
            if degPol[iPt] > 0:
                massFitterIns[iPt].SetPolDegreeForBackgroundFit(degPol[iPt])
                massFitterOuts[iPt].SetPolDegreeForBackgroundFit(degPol[iPt])
            massFitterIns[iPt].SetUseLikelihoodFit()
            massFitterOuts[iPt].SetUseLikelihoodFit()
            if fitConfig['BoundMean']:
                massFitterIns[iPt].SetBoundGaussianMean(massForFit, massMin, massMax)
                massFitterOuts[iPt].SetBoundGaussianMean(massForFit, massMin, massMax)
            else:
                massFitterIns[iPt].SetInitialGaussianMean(massForFit)
                massFitterOuts[iPt].SetInitialGaussianMean(massForFit)
            if fitConfig['FixSigmaRatio']:
                massFitterIns[iPt].SetFixRatio2GausSigma(
                    hSigmaToFix.GetBinContent(iPt+1)/hSigmaToFix2.GetBinContent(iPt+1))
                massFitterOuts[iPt].SetFixRatio2GausSigma(
                    hSigmaToFix.GetBinContent(iPt+1)/hSigmaToFix2.GetBinContent(iPt+1))

            if fixSigma[iPt]:
                if isinstance(fitConfig['SigmaMultFactor'], (float, int)):
                    massFitterIns[iPt].SetFixGaussianSigma(
                        hSigmaToFix.GetBinContent(iPt+1)*fitConfig['SigmaMultFactor'])
                    massFitterOuts[iPt].SetFixGaussianSigma(
                        hSigmaToFix.GetBinContent(iPt+1)*fitConfig['SigmaMultFactor'])
                else:
                    if fitConfig['SigmaMultFactor'] == 'MinusUnc':
                        massFitterIns[iPt].SetFixGaussianSigma(
                            hSigmaToFix.GetBinContent(iPt+1)-hSigmaToFix.GetBinError(iPt+1))
                        massFitterOuts[iPt].SetFixGaussianSigma(
                            hSigmaToFix.GetBinContent(iPt+1)-hSigmaToFix.GetBinError(iPt+1))
                    elif fitConfig['SigmaMultFactor'] == 'PlusUnc':
                        massFitterIns[iPt].SetFixGaussianSigma(
                            hSigmaToFix.GetBinContent(iPt+1)+hSigmaToFix.GetBinError(iPt+1))
                        massFitterOuts[iPt].SetFixGaussianSigma(
                            hSigmaToFix.GetBinContent(iPt+1)+hSigmaToFix.GetBinError(iPt+1))
                    else:
                        print('WARNING: impossible to fix sigma! Wrong mult factor set in config file!')
            else:
                if hSigmaToFix:
                    massFitterIns[iPt].SetInitialGaussianSigma(
                        hSigmaToFix.GetBinContent(iPt+1)*fitConfig['SigmaMultFactor'])
                    massFitterOuts[iPt].SetInitialGaussianSigma(
                        hSigmaToFix.GetBinContent(iPt+1)*fitConfig['SigmaMultFactor'])
                else:
                    if particleName == 'Dstar':
                        massFitterIns[iPt].SetInitialGaussianSigma(0.001)
                        massFitterOuts[iPt].SetInitialGaussianSigma(0.001)
                    else:
                        massFitterIns[iPt].SetInitialGaussianSigma(0.008)
                        massFitterOuts[iPt].SetInitialGaussianSigma(0.008)

            if secPeak and particleName == 'Ds':
                massFitterIns[iPt].IncludeSecondGausPeak(massDplus, False, fitConfig['SigmaSecPeak'][iPt], True)
                massFitterOuts[iPt].IncludeSecondGausPeak(massDplus, False, fitConfig['SigmaSecPeak'][iPt], True)
            # TODO: Add reflections for D0
            massFitterIns[iPt].MassFitter(False)
            massFitterOuts[iPt].MassFitter(False)

            # collect fit results
            rawyield_in = massFitterIns[iPt].GetRawYield()
            rawyielderr_in = massFitterIns[iPt].GetRawYieldError()
            sigma_in = massFitterIns[iPt].GetSigma()
            sigmaerr_in = massFitterIns[iPt].GetSigmaUncertainty()
            mean_in = massFitterIns[iPt].GetMean()
            meanerr_in = massFitterIns[iPt].GetMeanUncertainty()
            redchi2_in = massFitterIns[iPt].GetReducedChiSquare()
            signif_in, signiferr_in = ctypes.c_double(), ctypes.c_double()
            sgn_in, sgnerr_in = ctypes.c_double(), ctypes.c_double()
            bkg_in, bkgerr_in = ctypes.c_double(), ctypes.c_double()
            massFitterIns[iPt].Significance(3, signif_in, signiferr_in)
            massFitterIns[iPt].Signal(3, sgn_in, sgnerr_in)
            massFitterIns[iPt].Background(3, bkg_in, bkgerr_in)

            rawyield_out = massFitterOuts[iPt].GetRawYield()
            rawyielderr_out = massFitterOuts[iPt].GetRawYieldError()
            sigma_out = massFitterOuts[iPt].GetSigma()
            sigmaerr_out = massFitterOuts[iPt].GetSigmaUncertainty()
            mean_out = massFitterOuts[iPt].GetMean()
            meanerr_out = massFitterOuts[iPt].GetMeanUncertainty()
            redchi2_out = massFitterOuts[iPt].GetReducedChiSquare()
            signif_out, signiferr_out = ctypes.c_double(), ctypes.c_double()
            sgn_out, sgnerr_out = ctypes.c_double(), ctypes.c_double()
            bkg_out, bkgerr_out = ctypes.c_double(), ctypes.c_double()
            massFitterOuts[iPt].Significance(3, signif_out, signiferr_out)
            massFitterOuts[iPt].Signal(3, sgn_out, sgnerr_out)
            massFitterOuts[iPt].Background(3, bkg_out, bkgerr_out)

            cMass.cd(iPt+1)
            hMassInsForFit[iPt].GetYaxis().SetRangeUser(0, hMassInsForFit[iPt].GetMaximum()*2)
            hMassInsForFit[iPt].GetXaxis().SetRangeUser(massMin, massMax)
            hMassInsForFit[iPt].Draw('PE')
            hMassOutsForFit[iPt].Draw('PE same')
            fTotFuncMassIn = massFitterIns[iPt].GetMassFunc()
            fTotFuncMassOut = massFitterOuts[iPt].GetMassFunc()
            fBkgFuncMassIn = massFitterIns[iPt].GetBackgroundRecalcFunc()
            fBkgFuncMassOut = massFitterOuts[iPt].GetBackgroundRecalcFunc()
            SetObjectStyle(fTotFuncMassIn, color=kRed, linestyle=2, linewidth=3)
            SetObjectStyle(fTotFuncMassOut, color=kAzure, linestyle=9, linewidth=3)
            SetObjectStyle(fBkgFuncMassIn, color=kGray+1, linestyle=9, linewidth=1)
            SetObjectStyle(fBkgFuncMassOut, color=kGray+1, linestyle=9, linewidth=1)
            fBkgFuncMassIn.Draw('same')
            fBkgFuncMassOut.Draw('same')
            fTotFuncMassIn.Draw('same')
            fTotFuncMassOut.Draw('same')
            latex.SetTextColor(kRed)
            latex.DrawLatex(0.15, 0.80, f'#mu = {mean_in:.3f} #pm {meanerr_in:.3f} GeV/c^{2}')
            latex.DrawLatex(0.15, 0.75, f'#sigma = {sigma_in:.3f} #pm {sigmaerr_in:.3f} GeV/c^{2}')
            latex.DrawLatex(0.15, 0.70, f'S = {rawyield_in:.0f} #pm {rawyielderr_in:.0f}')
            latex.DrawLatex(0.15, 0.65, f'S/B = {sgn_in.value/bkg_in.value:.2f}')
            latex.DrawLatex(0.15, 0.60, f'#chi^{{2}}/ndf = {redchi2_in:.2f}')
            latex.DrawLatex(0.15, 0.55, f'Signif. = {signif_in.value:.2f} #pm {signiferr_in.value:.2f}')
            latex.SetTextColor(kAzure)
            latex.DrawLatex(0.6, 0.80, f'#mu = {mean_out:.3f} #pm {meanerr_out:.3f} GeV/c^{2}')
            latex.DrawLatex(0.6, 0.75, f'#sigma = {sigma_out:.3f} #pm {sigmaerr_out:.3f} GeV/c^{2}')
            latex.DrawLatex(0.6, 0.70, f'S = {rawyield_out:.0f} #pm {rawyielderr_out:.0f}')
            latex.DrawLatex(0.6, 0.65, f'S/B = {sgn_out.value/bkg_out.value:.2f}')
            latex.DrawLatex(0.6, 0.60, f'#chi^{{2}}/ndf = {redchi2_out:.2f}')
            latex.DrawLatex(0.6, 0.55, f'Signif. = {signif_out.value:.2f} #pm {signiferr_out.value:.2f}')
            cMass.Modified()
            cMass.Update()

            hRawYieldsIn.SetBinContent(iPt+1, rawyield_in)
            hRawYieldsIn.SetBinError(iPt+1, rawyielderr_in)
            hSigmaIn.SetBinContent(iPt+1, sigma_in)
            hSigmaIn.SetBinError(iPt+1, sigmaerr_in)
            hMeanIn.SetBinContent(iPt+1, mean_in)
            hMeanIn.SetBinError(iPt+1, meanerr_in)
            hRedChi2In.SetBinContent(iPt+1, redchi2_in)
            hRedChi2In.SetBinError(iPt+1, 1.e-20)
            hRawYieldsSignificanceIn.SetBinContent(iPt+1, signif_in.value)
            hRawYieldsSignificanceIn.SetBinError(iPt+1, signiferr_in.value)
            hRawYieldsSoverBIn.SetBinContent(iPt+1, sgn_in.value/bkg_in.value)
            hRawYieldsSoverBIn.SetBinError(iPt+1, 1.e-20)
            hRawYieldsOut.SetBinContent(iPt+1, rawyield_out)
            hRawYieldsOut.SetBinError(iPt+1, rawyielderr_out)
            hSigmaOut.SetBinContent(iPt+1, sigma_out)
            hSigmaOut.SetBinError(iPt+1, sigmaerr_out)
            hMeanOut.SetBinContent(iPt+1, mean_out)
            hMeanOut.SetBinError(iPt+1, meanerr_out)
            hRedChi2Out.SetBinContent(iPt+1, redchi2_out)
            hRedChi2Out.SetBinError(iPt+1, 1.e-20)
            hRawYieldsSignificanceOut.SetBinContent(iPt+1, signif_out.value)
            hRawYieldsSignificanceOut.SetBinError(iPt+1, signiferr_out.value)
            hRawYieldsSoverBOut.SetBinContent(iPt+1, sgn_out.value/bkg_out.value)
            hRawYieldsSoverBOut.SetBinError(iPt+1, 1.e-20)

            # vn
            vn, vnUnc = get_ep_vn(harmonic,
                                  rawyield_in, rawyielderr_in,
                                  rawyield_out, rawyielderr_out)
            gvnSimFit.SetPoint(iPt, (ptMin+ptMax)/2, vn)
            gvnSimFit.SetPointError(iPt, (ptMax-ptMin)/2, (ptMax-ptMin)/2, vnUnc, vnUnc)

    #save output histos
    print(f'Saving output to {outputdir}')
    if vn_method == 'sp' or vn_method == 'ep':
        for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):
            if iPt == 0:
                suffix_pdf = '('
            elif iPt == nPtBins-1:
                suffix_pdf = ')'
            else:
                suffix_pdf = ''
            cSimFit[iPt].SaveAs(f'{outputdir}/SimFit{suffix}_{particleName}.pdf{suffix_pdf}')
    outfile_name = f'{outputdir}/raw_yields{suffix}.root'
    outFile = TFile(outfile_name, 'recreate')
    if vn_method == 'sp' or vn_method == 'ep':
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
    else:
        cMass.SaveAs(f'{outputdir}/MassFit{suffix}_{particleName}.pdf')
        cMass.Write()
        for hist in hMassInsForFit:
            hist.Write()
        for hist in hMassOutsForFit:
            hist.Write()
        for ipt, (ptmin, ptmax) in enumerate(zip(ptMins, ptMaxs)):
            print(f'Writing {ptmin} - {ptmax} GeV/c')
            fTotFuncMassIn = massFitterIns[ipt].GetMassFunc()
            fTotFuncMassOut = massFitterOuts[ipt].GetMassFunc()
            fTotFuncMassIn.Write(f'fTotFuncMassIn_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
            fTotFuncMassOut.Write(f'fTotFuncMassOut_pt{ptmin*10:.0f}_{ptmax*10:.0f}')
        hRawYieldsIn.Write()
        hRawYieldsOut.Write()
        hSigmaIn.Write()
        hSigmaOut.Write()
        hMeanIn.Write()
        hMeanOut.Write()
        hRedChi2In.Write()
        hRedChi2Out.Write()
        hRawYieldsSignificanceIn.Write()
        hRawYieldsSignificanceOut.Write()
        hRawYieldsSoverBIn.Write()
        hRawYieldsSoverBOut.Write()
        hSigmaSecPeakFitIn.Write()
        hSigmaSecPeakFitOut.Write()

    gvnSimFit.Write()
    gvnSimFitSecPeak.Write()

    outFile.Close()

    if not batch:
        input('Press enter to exit')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('fitConfigFileName', metavar='text', default='config_Ds_Fit.yml')
    parser.add_argument('centClass', metavar='text', default='')
    parser.add_argument('inFileName', metavar='text', default='')
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument('--refFileName', metavar='text', default='')
    parser.add_argument('--vn_method', '-vn', metavar='text', default='sp')
    parser.add_argument('--batch', help='suppress video output', action='store_true')
    args = parser.parse_args()

    get_vn_vs_mass(
        args.fitConfigFileName,
        args.centClass,
        args.inFileName,
        args.outputdir,
        args.suffix,
        args.refFileName,
        args.vn_method,
        args.batch
    )