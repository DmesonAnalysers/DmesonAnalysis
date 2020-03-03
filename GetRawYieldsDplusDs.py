'''
Script for fitting D+ and Ds+ invariant-mass spectra
run: python GetRawYieldsDsDplus.py fitConfigFileName.yml centClass inputFileName.root outFileName.root
'''

import argparse
import numpy as np
import yaml
from ROOT import TFile, TCanvas, TH1D, TH1F, TF1, TDatabasePDG, AliHFInvMassFitter, AliVertexingHFUtils  # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, gPad, Double, kBlack, kRed, kFullCircle, kFullSquare  # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, DivideCanvas
from utils.AnalysisUtils import SingleGaus, DoubleGaus, DoublePeakSingleGaus, DoublePeakDoubleGaus

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('fitConfigFileName', metavar='text', default='config_Ds_Fit.yml')
parser.add_argument('centClass', metavar='text', default='')
parser.add_argument('inFileName', metavar='text', default='')
parser.add_argument('outFileName', metavar='text', default='')
parser.add_argument('--isMC', action='store_true', default=False)
parser.add_argument('--batch', help='suppress video output', action='store_true')
args = parser.parse_args()

cent = ''
if args.centClass == 'k010':
    cent = 'Cent010'
elif args.centClass == 'k3050':
    cent = 'Cent3050'
elif args.centClass == 'k6080':
    cent = 'Cent6080'
elif args.centClass == 'kpp5TeVPrompt':
    cent = 'pp5TeVPrompt'
elif args.centClass == 'kpp5TeVFD':
    cent = 'pp5TeVFD'

with open(args.fitConfigFileName, 'r') as ymlfitConfigFile:
    fitConfig = yaml.load(ymlfitConfigFile, yaml.FullLoader)

gROOT.SetBatch(args.batch)
SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, padtopmargin=0.12, opttitle=1)

ptMins = fitConfig[cent]['PtMin']
ptMaxs = fitConfig[cent]['PtMax']
ptLims = list(ptMins)
nPtBins = len(ptMins)
ptLims.append(ptMaxs[-1])

mesonName = fitConfig[cent]['Meson']
inclSecPeak = fitConfig[cent]['InclSecPeak']

SgnFunc, BkgFunc = [], []
for iPt, (bkg, sgn) in enumerate(zip(fitConfig[cent]['BkgFunc'], fitConfig[cent]['SgnFunc'])):
    if bkg == 'kExpo':
        BkgFunc.append(AliHFInvMassFitter.kExpo)
    elif bkg == 'kLin':
        SgnFunc.append(AliHFInvMassFitter.kLin)
    elif bkg == 'kPol2':
        SgnFunc.append(AliHFInvMassFitter.kPol2)
    else:
        print('ERROR: only kExpo, kLin, and kPol2 background functions supported! Exit')
        exit()

    if sgn == 'kGaus':
        SgnFunc.append(AliHFInvMassFitter.kGaus)
    elif sgn == 'k2Gaus':
        SgnFunc.append(AliHFInvMassFitter.k2Gaus)
    else:
        print('ERROR: only kGaus and k2Gaus signal functions supported! Exit')
        exit()

if mesonName == 'Dplus':
    massAxisTit = '#it{M}(K#pi#pi) (GeV/#it{c}^{2})'
elif mesonName == 'Ds':
    massAxisTit = '#it{M}(KK#pi) (GeV/#it{c}^{2})'

# load inv-mass histos
infile = TFile.Open(args.inFileName)
if not infile or not infile.IsOpen():
    exit()

hMass, hMassForFit = [], []
for iPt, (ptMin, ptMax, secPeak) in enumerate(zip(ptMins, ptMaxs, inclSecPeak)):
    if not args.isMC:
        hMass.append(infile.Get('hMass_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10)))
        hMass[iPt].SetDirectory(0)
    else:
        hMass.append(infile.Get('hPromptMass_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10)))
        hMass[iPt].Add(infile.Get('hFDMass_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10)))
        if secPeak:
            hMass[iPt].Add(infile.Get('hPromptSecPeakMass_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10)))
            hMass[iPt].Add(infile.Get('hFDSecPeakMass_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10)))
            hMass[iPt].SetDirectory(0)
    hMass[iPt].Sumw2()
    SetObjectStyle(hMass[iPt], color=kBlack, markerstyle=kFullCircle)

hEv = infile.Get('hEvForNorm')
hEv.SetDirectory(0)
hEv.Sumw2()
SetObjectStyle(hEv, color=kBlack, markerstyle=kFullCircle)

infile.Close()

if fitConfig[cent]['FixSigma']:
    infileSigma = TFile.Open(fitConfig[cent]['SigmaFile'])
    if not infileSigma:
        exit()
    hSigmaToFix = infileSigma.Get("hRawYieldsSigma")
    if hSigmaToFix.GetNbinsX() != nPtBins:
        print("WARNING: Different number of bins for this analysis and histo for fix sigma")

if fitConfig[cent]['FixMean']:
    infileMean = TFile.Open(fitConfig[cent]['MeanFile'])
    if not infileMean or not infileMean.IsOpen():
        exit()
    hMeanToFix = infileMean.Get("hRawYieldsMean")
    if hMeanToFix.GetNbinsX() != nPtBins:
        print("WARNING: Different number of bins for this analysis and histo for fix mean")

hRawYields = TH1D('hRawYields', ';#it{p}_{T} (GeV/#it{c});raw yield', nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSigma = TH1D('hRawYieldsSigma', ';#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSigma2 = TH1D('hRawYieldsSigma2', ';#it{p}_{T} (GeV/#it{c});width (GeV/#it{c}^{2})', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsMean = TH1D('hRawYieldsMean', ';#it{p}_{T} (GeV/#it{c});mean (GeV/#it{c}^{2})', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsFracGaus2 = TH1D('hRawYieldsFracGaus2', ';#it{p}_{T} (GeV/#it{c});second-gaussian fraction', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSignificance = TH1D('hRawYieldsSignificance', ';#it{p}_{T} (GeV/#it{c});significance (3#sigma)', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSoverB = TH1D('hRawYieldsSoverB', ';#it{p}_{T} (GeV/#it{c});S/B (3#sigma)', nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSignal = TH1D('hRawYieldsSignal', ';#it{p}_{T} (GeV/#it{c});Signal (3#sigma)', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsBkg = TH1D('hRawYieldsBkg', ';#it{p}_{T} (GeV/#it{c});Background (3#sigma)', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsChiSquare = TH1D('hRawYieldsChiSquare', ';#it{p}_{T} (GeV/#it{c});#chi^{2}/#it{ndf}', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSecPeak = TH1D('hRawYieldsSecPeak', ';#it{p}_{T} (GeV/#it{c});raw yield second peak', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsMeanSecPeak = TH1D('hRawYieldsMeanSecPeak', ';#it{p}_{T} (GeV/#it{c});mean second peak (GeV/#it{c}^{2})', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSigmaSecPeak = TH1D('hRawYieldsSigmaSecPeak', \
    ';#it{p}_{T} (GeV/#it{c});width second peak (GeV/#it{c}^{2})', nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSignificanceSecPeak = TH1D('hRawYieldsSignificanceSecPeak', \
    ';#it{p}_{T} (GeV/#it{c});signficance second peak (3#sigma)', nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSigmaRatioSecondFirstPeak = TH1D('hRawYieldsSigmaRatioSecondFirstPeak', \
    ';#it{p}_{T} (GeV/#it{c});width second peak / width first peak', nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSoverBSecPeak = TH1D('hRawYieldsSoverBSecPeak', ';#it{p}_{T} (GeV/#it{c});S/B second peak (3#sigma)', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSignalSecPeak = TH1D('hRawYieldsSignalSecPeak', ';#it{p}_{T} (GeV/#it{c});Signal second peak (3#sigma)', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsBkgSecPeak = TH1D('hRawYieldsBkgSecPeak', ';#it{p}_{T} (GeV/#it{c});Background second peak (3#sigma)', \
    nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsTrue = TH1D('hRawYieldsTrue', ';#it{p}_{T} (GeV/#it{c});true signal', nPtBins, np.asarray(ptLims, 'd'))
hRawYieldsSecPeakTrue = TH1D('hRawYieldsSecPeakTrue', ';#it{p}_{T} (GeV/#it{c});true signal second peak', \
    nPtBins, np.asarray(ptLims, 'd'))
hRelDiffRawYieldsFitTrue = TH1D('hRelDiffRawYieldsFitTrue', \
    ';#it{p}_{T} (GeV/#it{c}); (Y_{fit} - Y_{true}) / Y_{true}', nPtBins, np.asarray(ptLims, 'd'))
hRelDiffRawYieldsSecPeakFitTrue = TH1D('hRelDiffRawYieldsSecPeakFitTrue', \
    ';#it{p}_{T} (GeV/#it{c});(Y_{fit} - Y_{true}) / Y_{true} second peak', nPtBins, np.asarray(ptLims, 'd'))

SetObjectStyle(hRawYields, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSigma, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSigma2, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsMean, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsFracGaus2, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSignificance, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSoverB, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSignal, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsBkg, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsChiSquare, color=kBlack, markerstyle=kFullCircle)
SetObjectStyle(hRawYieldsSecPeak, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsMeanSecPeak, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsSigmaSecPeak, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsSignificanceSecPeak, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsSigmaRatioSecondFirstPeak, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsSoverBSecPeak, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsSignalSecPeak, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsBkgSecPeak, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsTrue, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRawYieldsSecPeakTrue, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRelDiffRawYieldsFitTrue, color=kRed, markerstyle=kFullSquare)
SetObjectStyle(hRelDiffRawYieldsSecPeakFitTrue, color=kRed, markerstyle=kFullSquare)

# fit histos
massDplus = TDatabasePDG.Instance().GetParticle(411).Mass()
massDs = TDatabasePDG.Instance().GetParticle(431).Mass()
massForFit = massDplus if fitConfig[cent]['Meson'] == 'Dplus' else massDs

cMass = TCanvas("cMass", "cMass", 1920, 1080)
DivideCanvas(cMass, nPtBins)
cResiduals = TCanvas("cResiduals", "cResiduals", 1920, 1080)
DivideCanvas(cResiduals, nPtBins)

massFitter = []
for iPt, (hM, ptMin, ptMax, reb, sgn, bkg, secPeak, massMin, massMax) in enumerate(zip(
        hMass, ptMins, ptMaxs, fitConfig[cent]['Rebin'], SgnFunc, BkgFunc, \
            inclSecPeak, fitConfig[cent]['MassMin'], fitConfig[cent]['MassMax'])):

    hMassForFit.append(TH1F())
    AliVertexingHFUtils.RebinHisto(hM, reb).Copy(hMassForFit[iPt]) #to cast TH1D to TH1F
    hMassForFit[iPt].SetDirectory(0)
    hMassForFit[iPt].SetTitle(
        f"{ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};{massAxisTit};"
        f"Counts per {hMassForFit[iPt].GetBinWidth(1)*1000:.0f} MeV/#it{{c}}^{{2}}")
    hMassForFit[iPt].SetName(f"MassForFit{iPt}")
    SetObjectStyle(hMassForFit[iPt], color=kBlack, markerstyle=kFullCircle)

    # MC
    if args.isMC:
        parRawYield, parMean, parSigma1 = 0, 1, 2 # always the same
        parSigma2, parFrac2Gaus, parRawYieldSecPeak, parMeanSecPeak, parSigmaSecPeak = (-1 for _ in range(5))

        if sgn == AliHFInvMassFitter.kGaus:
            if not (secPeak and mesonName == 'Ds'):
                massFunc = TF1("massFunc{0}".format(iPt), SingleGaus, massMin, massMax, 3)
                massFunc.SetParameters(hMassForFit[iPt].Integral()*hMassForFit[iPt].GetBinWidth(1), massForFit, 0.010)
            else:
                massFunc = TF1("massFunc{0}".format(iPt), DoublePeakSingleGaus, massMin, massMax, 6)
                massFunc.SetParameters(hMassForFit[iPt].Integral()*hMassForFit[iPt].GetBinWidth(1), massForFit, 0.010, \
                    hMassForFit[iPt].Integral()*hMassForFit[iPt].GetBinWidth(1), massDplus, 0.010)
                parRawYieldSecPeak = 3
                parMeanSecPeak = 4
                parSigmaSecPeak = 5

        elif sgn == AliHFInvMassFitter.k2Gaus:
            parSigma2 = 3
            parFrac2Gaus = 4
            if not (secPeak and mesonName == 'Ds'):
                massFunc = TF1("massFunc{0}".format(iPt), DoubleGaus, massMin, massMax, 5)
                massFunc.SetParameters(hMassForFit[iPt].Integral()*hMassForFit[iPt].GetBinWidth(1), massForFit, \
                    0.010, 0.030, 0.9)
            else:
                massFunc = TF1("massFunc{0}".format(iPt), DoublePeakDoubleGaus, massMin, massMax, 8)
                massFunc.SetParameters(hMassForFit[iPt].Integral()*hMassForFit[iPt].GetBinWidth(1), massForFit, \
                    0.010, 0.030, 0.9, hMassForFit[iPt].Integral()*hMassForFit[iPt].GetBinWidth(1), massDplus, 0.010)
                parRawYieldSecPeak = 5
                parMeanSecPeak = 6
                parSigmaSecPeak = 7

        if nPtBins > 1:
            cMass.cd(iPt+1)
        else:
            cMass.cd()
        hMassForFit[iPt].Fit(massFunc, 'E')  # fit with chi2

        rawyield = massFunc.GetParameter(parRawYield)
        rawyielderr = massFunc.GetParError(parRawYield)
        sigma = massFunc.GetParameter(parSigma1)
        sigmaerr = massFunc.GetParError(parSigma1)
        mean = massFunc.GetParameter(parMean)
        meanerr = massFunc.GetParError(parMean)
        redchi2 = massFunc.GetChisquare() / massFunc.GetNDF()

        hRawYields.SetBinContent(iPt+1, rawyield)
        hRawYields.SetBinError(iPt+1, rawyielderr)
        hRawYieldsSigma.SetBinContent(iPt+1, sigma)
        hRawYieldsSigma.SetBinError(iPt+1, sigmaerr)
        hRawYieldsMean.SetBinContent(iPt+1, mean)
        hRawYieldsMean.SetBinError(iPt+1, meanerr)
        hRawYieldsChiSquare.SetBinContent(iPt+1, redchi2)
        hRawYieldsChiSquare.SetBinError(iPt+1, 0.)

        hRawYieldsTrue.SetBinContent(iPt+1, hMassForFit[iPt].Integral())
        hRawYieldsTrue.SetBinError(iPt+1, np.sqrt(hMassForFit[iPt].Integral()))
        hRelDiffRawYieldsFitTrue.SetBinContent(iPt+1, rawyield-hMassForFit[iPt].Integral())
        hRelDiffRawYieldsFitTrue.SetBinError(iPt+1, np.sqrt(rawyielderr*rawyielderr+hMassForFit[iPt].Integral()))

        if secPeak and mesonName == 'Ds':
            rawyieldSecPeak = massFunc.GetParameter(parRawYieldSecPeak)
            rawyieldSecPeakerr = massFunc.GetParError(parRawYieldSecPeak)
            sigmaSecPeak = massFunc.GetParameter(parSigmaSecPeak)
            sigmaSecPeakerr = massFunc.GetParError(parSigmaSecPeak)
            meanSecPeak = massFunc.GetParameter(parMeanSecPeak)
            meanSecPeakerr = massFunc.GetParError(parMeanSecPeak)
            hRawYieldsSecPeak.SetBinContent(iPt+1, rawyieldSecPeak)
            hRawYieldsSecPeak.SetBinError(iPt+1, rawyieldSecPeakerr)
            hRawYieldsMeanSecPeak.SetBinContent(iPt+1, meanSecPeak)
            hRawYieldsMeanSecPeak.SetBinError(iPt+1, meanSecPeakerr)
            hRawYieldsSigmaSecPeak.SetBinContent(iPt+1, sigmaSecPeak)
            hRawYieldsSigmaSecPeak.SetBinError(iPt+1, sigmaSecPeakerr)
            hRawYieldsSigmaRatioSecondFirstPeak.SetBinContent(iPt+1, sigmaSecPeak/sigma)
            hRawYieldsSigmaRatioSecondFirstPeak.SetBinError(iPt+1, \
                np.sqrt(sigmaerr**2/sigma**2+sigmaSecPeakerr**2/sigmaSecPeak**2)*sigmaSecPeak/sigma)

            hRawYieldsSecPeakTrue.SetBinContent(iPt+1, rawyield)
            hRelDiffRawYieldsSecPeakFitTrue.SetBinContent(iPt+1, rawyield)

        if sgn == AliHFInvMassFitter.k2Gaus:
            sigma2 = massFunc.GetParameter(parSigma2)
            sigma2err = massFunc.GetParError(parSigma2)
            frac2gaus = massFunc.GetParameter(parFrac2Gaus)
            frac2gauserr = massFunc.GetParError(parFrac2Gaus)
            hRawYieldsSigma2.SetBinContent(iPt+1, sigma2)
            hRawYieldsSigma2.SetBinError(iPt+1, sigma2err)
            hRawYieldsFracGaus2.SetBinContent(iPt+1, frac2gaus)
            hRawYieldsFracGaus2.SetBinError(iPt+1, frac2gauserr)

    else:  # data
        massFitter.append(AliHFInvMassFitter(hMassForFit[iPt], massMin, massMax, bkg, sgn))

        if fitConfig[cent]['UseLikelihood']:
            massFitter[iPt].SetUseLikelihoodFit()
        if fitConfig[cent]['FixMean']:
            massFitter[iPt].SetFixGaussianMean(hMeanToFix.GetBinContent(iPt+1))
        else:
            massFitter[iPt].SetInitialGaussianMean(massForFit)

        if fitConfig[cent]['FixSigma']:
            massFitter[iPt].SetFixGaussianSigma(
                hSigmaToFix.GetBinContent(iPt+1)*fitConfig[cent]['SigmaMultFactor'])

        if secPeak and mesonName == 'Ds':
            # TODO: add possibility to fix D+ peak to sigmaMC(D+)/sigmaMC(Ds+)*sigmaData(Ds+)
            massFitter[iPt].IncludeSecondGausPeak(massDplus, False, 0.008, True)
        massFitter[iPt].MassFitter(False)

        rawyield = massFitter[iPt].GetRawYield()
        rawyielderr = massFitter[iPt].GetRawYieldError()
        sigma = massFitter[iPt].GetSigma()
        sigmaerr = massFitter[iPt].GetSigmaUncertainty()
        mean = massFitter[iPt].GetMean()
        meanerr = massFitter[iPt].GetMeanUncertainty()
        redchi2 = massFitter[iPt].GetReducedChiSquare()
        signif, signiferr = Double(), Double()
        sgn, sgnerr = Double(), Double()
        bkg, bkgerr = Double(), Double()
        massFitter[iPt].Significance(3, signif, signiferr)
        massFitter[iPt].Signal(3, sgn, sgnerr)
        massFitter[iPt].Background(3, bkg, bkgerr)

        hRawYields.SetBinContent(iPt+1, rawyield)
        hRawYields.SetBinError(iPt+1, rawyielderr)
        hRawYieldsSigma.SetBinContent(iPt+1, sigma)
        hRawYieldsSigma.SetBinError(iPt+1, sigmaerr)
        hRawYieldsMean.SetBinContent(iPt+1, mean)
        hRawYieldsMean.SetBinError(iPt+1, meanerr)
        hRawYieldsSignificance.SetBinContent(iPt+1, signif)
        hRawYieldsSignificance.SetBinError(iPt+1, signiferr)
        hRawYieldsSoverB.SetBinContent(iPt+1, sgn/bkg)
        hRawYieldsSoverB.SetBinError(iPt+1, sgn/bkg*np.sqrt(sgnerr**2/sgn**2+bkgerr**2/bkg**2))
        hRawYieldsSignal.SetBinContent(iPt+1, sgn)
        hRawYieldsSignal.SetBinError(iPt+1, sgnerr)
        hRawYieldsBkg.SetBinContent(iPt+1, bkg)
        hRawYieldsBkg.SetBinError(iPt+1, bkgerr)
        hRawYieldsChiSquare.SetBinContent(iPt+1, redchi2)
        hRawYieldsChiSquare.SetBinError(iPt+1, 1.e-20)

        fTotFunc = massFitter[iPt].GetMassFunc()
        fBkgFunc = massFitter[iPt].GetBackgroundRecalcFunc()

        parFrac2Gaus, parsecondsigma = -1, -1
        if sgn == AliHFInvMassFitter.k2Gaus:
            if not (inclSecPeak and mesonName == 'Ds'):
                parFrac2Gaus = fTotFunc.GetNpar()-2
                parsecondsigma = fTotFunc.GetNpar()-1
            else:
                parFrac2Gaus = fTotFunc.GetNpar()-5
                parsecondsigma = fTotFunc.GetNpar()-4

            sigma2 = fTotFunc.GetParameter(parsecondsigma)
            sigma2err = fTotFunc.GetParError(parsecondsigma)
            frac2gaus = fTotFunc.GetParameter(parFrac2Gaus)
            frac2gauserr = fTotFunc.GetParError(parFrac2Gaus)
            hRawYieldsSigma2.SetBinContent(iPt+1, sigma2)
            hRawYieldsSigma2.SetBinError(iPt+1, sigma2err)
            hRawYieldsFracGaus2.SetBinContent(iPt+1, frac2gaus)
            hRawYieldsFracGaus2.SetBinError(iPt+1, frac2gauserr)

        if inclSecPeak and mesonName == 'Ds':
            paryieldSecPeak = fTotFunc.GetNpar()-3
            parMeanSecPeak = fTotFunc.GetNpar()-2
            parSigmaSecPeak = fTotFunc.GetNpar()-1

            rawyieldSecPeak = fTotFunc.GetParameter(paryieldSecPeak) / hMassForFit[iPt].GetBinWidth(1)
            rawyieldSecPeakerr = fTotFunc.GetParError(paryieldSecPeak) / hMassForFit[iPt].GetBinWidth(1)
            meanSecPeak = fTotFunc.GetParameter(parMeanSecPeak)
            meanSecPeakerr = fTotFunc.GetParError(parMeanSecPeak)
            sigmaSecPeak = fTotFunc.GetParameter(parSigmaSecPeak)
            sigmaSecPeakerr = fTotFunc.GetParError(parSigmaSecPeak)

            bkgSecPeak = fBkgFunc.Integral(\
                meanSecPeak-3*sigmaSecPeak, meanSecPeak+3*sigmaSecPeak) / hMassForFit[iPt].GetBinWidth(1)
            bkgSecPeakerr = np.sqrt(bkgSecPeak)
            signalSecPeak = fTotFunc.Integral(\
                meanSecPeak-3*sigmaSecPeak, meanSecPeak+3*sigmaSecPeak) / hMassForFit[iPt].GetBinWidth(1)-bkgSecPeak
            signalSecPeakerr = np.sqrt(signalSecPeak+bkgSecPeak)
            signifSecPeak, signifSecPeakerr = -1., 1.
            AliVertexingHFUtils.ComputeSignificance(signalSecPeak, signalSecPeakerr, bkgSecPeak, \
                bkgSecPeakerr, signifSecPeak, signifSecPeakerr)

            hRawYieldsSecPeak.SetBinContent(iPt+1, rawyieldSecPeak)
            hRawYieldsSecPeak.SetBinError(iPt+1, rawyieldSecPeakerr)
            hRawYieldsMeanSecPeak.SetBinContent(iPt+1, meanSecPeak)
            hRawYieldsMeanSecPeak.SetBinError(iPt+1, meanSecPeakerr)
            hRawYieldsSigmaSecPeak.SetBinContent(iPt+1, sigmaSecPeak)
            hRawYieldsSigmaSecPeak.SetBinError(iPt+1, sigmaSecPeakerr)
            hRawYieldsSignificanceSecPeak.SetBinContent(iPt+1, signifSecPeak)
            hRawYieldsSignificanceSecPeak.SetBinError(iPt+1, signifSecPeakerr)
            hRawYieldsSigmaRatioSecondFirstPeak.SetBinContent(iPt+1, sigmaSecPeak/sigma)
            hRawYieldsSigmaRatioSecondFirstPeak.SetBinError(iPt+1, \
                np.sqrt(sigmaerr**2/sigma**2+sigmaSecPeakerr**2/sigmaSecPeak**2)*sigmaSecPeak/sigma)
            hRawYieldsSoverBSecPeak.SetBinContent(iPt+1, signalSecPeak/bkgSecPeak)
            hRawYieldsSoverBSecPeak.SetBinError(iPt+1, \
                signalSecPeak/bkgSecPeak*np.sqrt(signalSecPeakerr**2/signalSecPeak**2+bkgSecPeakerr**2/bkgSecPeak**2))
            hRawYieldsSignalSecPeak.SetBinContent(iPt+1, signalSecPeak)
            hRawYieldsSignalSecPeak.SetBinError(iPt+1, signalSecPeakerr)
            hRawYieldsBkgSecPeak.SetBinContent(iPt+1, bkgSecPeak)
            hRawYieldsBkgSecPeak.SetBinError(iPt+1, bkgSecPeakerr)

        if nPtBins > 1:
            cMass.cd(iPt+1)
        else:
            cMass.cd()

        hMassForFit[iPt].GetYaxis().SetRangeUser(hMassForFit[iPt].GetMinimum()*0.95, hMassForFit[iPt].GetMaximum()*1.2)
        massFitter[iPt].DrawHere(gPad)

        # residuals
        if nPtBins > 1:
            cResiduals.cd(iPt+1)
        else:
            cResiduals.cd()
        massFitter[iPt].DrawHistoMinusFit(gPad)

    cMass.Modified()
    cMass.Update()

    cResiduals.Modified()
    cResiduals.Update()

#save output histos
outFile = TFile(args.outFileName, "recreate")
cMass.Write()
if not args.isMC:
    cResiduals.Write()
for iPt in range(nPtBins):
    hMass[iPt].Write()
hRawYields.Write()
hRawYieldsSigma.Write()
hRawYieldsMean.Write()
hRawYieldsSignificance.Write()
hRawYieldsSoverB.Write()
hRawYieldsSignal.Write()
hRawYieldsBkg.Write()
hRawYieldsChiSquare.Write()
hRawYieldsSigma2.Write()
hRawYieldsFracGaus2.Write()
hRawYieldsSecPeak.Write()
hRawYieldsMeanSecPeak.Write()
hRawYieldsSigmaSecPeak.Write()
hRawYieldsSignificanceSecPeak.Write()
hRawYieldsSigmaRatioSecondFirstPeak.Write()
hRawYieldsSoverBSecPeak.Write()
hRawYieldsSignalSecPeak.Write()
hRawYieldsBkgSecPeak.Write()
hRawYieldsTrue.Write()
hRawYieldsSecPeakTrue.Write()
hRelDiffRawYieldsFitTrue.Write()
hRelDiffRawYieldsSecPeakFitTrue.Write()
hEv.Write()
outFile.Close()

if not args.batch:
    outFileNamePDF = args.outFileName.replace('.root', '.pdf')
    cMass.SaveAs(outFileNamePDF)
    outFileNamePDF = args.outFileName.replace(".root", "_Residuals.pdf")
    cResiduals.SaveAs(outFileNamePDF)
    input('Press enter to exit')
