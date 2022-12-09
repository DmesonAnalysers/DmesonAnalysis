'''
Script for fitting D+, D0 and Ds+ invariant-mass spectra
run: python GetRawYieldsDsDplus.py fitConfigFileName.yml centClass inputFileName.root outFileName.root
            [--refFileName][--isMC][--batch]
'''

import sys
import argparse
import ctypes
import numpy as np
import yaml
from ROOT import TFile, TCanvas, TH1D, TH1F, TF1, TDatabasePDG, TDirectoryFile, AliHFInvMassFitter, AliVertexingHFUtils # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, gPad, kBlack, kRed, kFullCircle, kFullSquare # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, DivideCanvas
from utils.FitUtils import SingleGaus, DoubleGaus, DoublePeakSingleGaus, DoublePeakDoubleGaus

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('fitConfigFileName', metavar='text', default='config_Ds_Fit.yml')
parser.add_argument('centClass', metavar='text', default='')
parser.add_argument('inFileName', metavar='text', default='')
parser.add_argument('outFileName', metavar='text', default='')
parser.add_argument('--refFileName', metavar='text', default='')
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
elif args.centClass == 'kpp13TeVFD':
    cent = 'pp13TeVFD'
elif args.centClass == 'kpp13TeVPrompt':
    cent = 'pp13TeVPrompt'
else:
    print(f"ERROR: cent class \'{args.centClass}\' is not supported! Exit")
    sys.exit()


with open(args.fitConfigFileName, 'r', encoding='utf8') as ymlfitConfigFile:
    fitConfig = yaml.load(ymlfitConfigFile, yaml.FullLoader)

gROOT.SetBatch(args.batch)
SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, padtopmargin=0.12, opttitle=1)

ptMins = fitConfig[cent]['PtMin']
ptMaxs = fitConfig[cent]['PtMax']
fixSigma = fitConfig[cent]['FixSigma']
fixMean = fitConfig[cent]['FixMean']
if 'EnableRef' not in fitConfig[cent]:
    enableRef = False
else:
    enableRef = fitConfig[cent]['EnableRef']
if not isinstance(fixSigma, list):
    fixSigma = [fixSigma for _ in ptMins]
if not isinstance(fixMean, list):
    fixMean = [fixMean for _ in ptMins]
ptLims = list(ptMins)
nPtBins = len(ptMins)
ptLims.append(ptMaxs[-1])

particleName = fitConfig[cent]['Particle']
inclSecPeak = fitConfig[cent]['InclSecPeak']

SgnFuncStr = fitConfig[cent]['SgnFunc']
if not isinstance(SgnFuncStr, list):
    SgnFuncStr = [SgnFuncStr] * nPtBins

BkgFuncStr = fitConfig[cent]['BkgFunc']
if not isinstance(BkgFuncStr, list):
    BkgFuncStr = [BkgFuncStr] * nPtBins

SgnFunc, BkgFunc, degPol = [], [], []
for iPt, (bkgStr, sgnStr) in enumerate(zip(BkgFuncStr, SgnFuncStr)):
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
hMass, hMassForFit = [], []
inclSecPeak = [inclSecPeak] * len(ptMins) if not isinstance(inclSecPeak, list) else inclSecPeak
for iPt, (ptMin, ptMax, secPeak) in enumerate(zip(ptMins, ptMaxs, inclSecPeak)):
    if not args.isMC:
        hMass.append(infile.Get(f'hMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
        hMass[iPt].SetDirectory(0)
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

hEv = infile.Get('hEvForNorm')
hEv.SetDirectory(0)
hEv.Sumw2()
SetObjectStyle(hEv, color=kBlack, markerstyle=kFullCircle)
infile.Close()

hSigmaToFix = None
if sum(fixSigma) > 0:
    infileSigma = TFile.Open(fitConfig[cent]['SigmaFile'])
    if not infileSigma:
        print(f'ERROR: file "{infileSigma}" cannot be opened! Exit!')
        sys.exit()
    hSigmaToFix = infileSigma.Get('hRawYieldsSigma')
    hSigmaToFix.SetDirectory(0)
    if hSigmaToFix.GetNbinsX() != nPtBins:
        print('WARNING: Different number of bins for this analysis and histo for fix sigma')
    infileSigma.Close()

if fitConfig[cent]['FixSigmaRatio']:
    # load sigma of first gaussian
    infileSigma = TFile.Open(fitConfig[cent]['SigmaRatioFile'])
    if not infileSigma:
        print(f'ERROR: file "{infileSigma}" cannot be opened! Exit!')
        sys.exit()
    hSigmaToFix = infileSigma.Get('hRawYieldsSigma')
    hSigmaToFix.SetDirectory(0)
    if hSigmaToFix.GetNbinsX() != nPtBins:
        print('WARNING: Different number of bins for this analysis and histo for fix sigma')
    infileSigma.Close()
    # load sigma of second gaussian
    infileSigma2 = TFile.Open(fitConfig[cent]['SigmaRatioFile'])
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
    infileMean = TFile.Open(fitConfig[cent]['MeanFile'])
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
if fitConfig[cent]['SigmaFileSecPeak']:
    infileSigmaSecPeak = TFile.Open(fitConfig[cent]['SigmaFileSecPeak'])
if fitConfig[cent]['FixSigmaToFirstPeak'] and not infileSigmaSecPeak:
    print(f'ERROR: file "{fitConfig[cent]["SigmaFileSecPeak"]}" cannot be opened! Exit!')
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

hRawYields = TH1D('hRawYields', f';{ptTit};raw yield', nPtBins, ptBinsArr)
hRawYieldsSigma = TH1D('hRawYieldsSigma', f';{ptTit};width (GeV/#it{{c}}^{{2}})', nPtBins, ptBinsArr)
hRawYieldsSigma2 = TH1D('hRawYieldsSigma2', f';{ptTit};width (GeV/#it{{c}}^{{2}})', nPtBins, ptBinsArr)
hRawYieldsMean = TH1D('hRawYieldsMean', f';{ptTit};mean (GeV/#it{{c}}^{{2}})', nPtBins, ptBinsArr)
hRawYieldsFracGaus2 = TH1D('hRawYieldsFracGaus2', f';{ptTit};second-gaussian fraction', nPtBins, ptBinsArr)
hRawYieldsSignificance = TH1D('hRawYieldsSignificance', f';{ptTit};significance (3#sigma)', nPtBins, ptBinsArr)
hRawYieldsSoverB = TH1D('hRawYieldsSoverB', f';{ptTit};S/B (3#sigma)', nPtBins, ptBinsArr)
hRawYieldsSignal = TH1D('hRawYieldsSignal', f';{ptTit};Signal (3#sigma)', nPtBins, ptBinsArr)
hRawYieldsBkg = TH1D('hRawYieldsBkg', f';{ptTit};Background (3#sigma)', nPtBins, ptBinsArr)
hRawYieldsChiSquare = TH1D('hRawYieldsChiSquare', f';{ptTit};#chi^{{2}}/#it{{ndf}}', nPtBins, ptBinsArr)
hRawYieldsSecPeak = TH1D('hRawYieldsSecondPeak', f';{ptTit};raw yield second peak', nPtBins, ptBinsArr)
hRawYieldsMeanSecPeak = TH1D('hRawYieldsMeanSecondPeak', f';{ptTit};mean second peak (GeV/#it{{c}}^{{2}})',
                             nPtBins, ptBinsArr)
hRawYieldsSigmaSecPeak = TH1D('hRawYieldsSigmaSecondPeak', f';{ptTit};width second peak (GeV/#it{{c}}^{{2}})',
                              nPtBins, ptBinsArr)
hRawYieldsSignificanceSecPeak = TH1D('hRawYieldsSignificanceSecondPeak',
                                     f';{ptTit};signficance second peak (3#sigma)', nPtBins, ptBinsArr)
hRawYieldsSigmaRatioSecondFirstPeak = TH1D('hRawYieldsSigmaRatioSecondFirstPeak',
                                           f';{ptTit};width second peak / width first peak', nPtBins, ptBinsArr)
hRawYieldsSoverBSecPeak = TH1D('hRawYieldsSoverBSecondPeak', f';{ptTit};S/B second peak (3#sigma)',
                               nPtBins, ptBinsArr)
hRawYieldsSignalSecPeak = TH1D('hRawYieldsSignalSecondPeak', f';{ptTit};Signal second peak (3#sigma)',
                               nPtBins, ptBinsArr)
hRawYieldsBkgSecPeak = TH1D('hRawYieldsBkgSecondPeak', f';{ptTit};Background second peak (3#sigma)',
                            nPtBins, ptBinsArr)
hRawYieldsTrue = TH1D('hRawYieldsTrue', f';{ptTit};true signal', nPtBins, ptBinsArr)
hRawYieldsSecPeakTrue = TH1D('hRawYieldsSecondPeakTrue', f';{ptTit};true signal second peak',
                             nPtBins, ptBinsArr)
hRelDiffRawYieldsFitTrue = TH1D('hRelDiffRawYieldsFitTrue', f';{ptTit}; (Y_{{fit}} - Y_{{true}}) / Y_{{true}}',
                                nPtBins, ptBinsArr)
hRelDiffRawYieldsSecPeakFitTrue = TH1D('hRelDiffRawYieldsSecondPeakFitTrue',
                                       f';{ptTit};(Y_{{fit}} - Y_{{true}}) / Y_{{true}} second peak',
                                       nPtBins, ptBinsArr)

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

nMaxCanvases = 20 # do not put more than 20 bins per canvas to make them visible
nCanvases = int(np.ceil(nPtBins / nMaxCanvases))
cMass, cResiduals = [], []
for iCanv in range(nCanvases):
    nPads = nPtBins if nCanvases == 1 else nMaxCanvases
    cMass.append(TCanvas(f'cMass{iCanv}', f'cMass{iCanv}', canvSizes[0], canvSizes[1]))
    DivideCanvas(cMass[iCanv], nPads)
    cResiduals.append(TCanvas(f'cResiduals{iCanv}', f'cResiduals{iCanv}', canvSizes[0], canvSizes[1]))
    DivideCanvas(cResiduals[iCanv], nPads)

massFitter = []

rebins = fitConfig[cent]['Rebin']
if not isinstance(rebins, list):
    rebins = [rebins] * len(ptMins)

massMins = fitConfig[cent]['MassMin']
if not isinstance(massMins, list):
    massMins = [massMins] * len(ptMins)

massMaxs = fitConfig[cent]['MassMax']
if not isinstance(massMaxs, list):
    massMaxs = [massMaxs] * len(ptMins)

for iPt, (hM, ptMin, ptMax, reb, sgnEnum, bkgEnum, secPeak, massMin, massMax) in enumerate(
        zip(hMass, ptMins, ptMaxs, rebins, SgnFunc, BkgFunc, inclSecPeak, massMins, massMaxs)):
    iCanv = int(np.floor(iPt / nMaxCanvases))
    hMassForFit.append(TH1F())
    AliVertexingHFUtils.RebinHisto(hM, reb).Copy(hMassForFit[iPt]) #to cast TH1D to TH1F
    hMassForFit[iPt].SetDirectory(0)
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

    # MC
    if args.isMC:
        parRawYield, parMean, parSigma1 = 0, 1, 2 # always the same
        parSigma2, parFrac2Gaus, parRawYieldSecPeak, parMeanSecPeak, parSigmaSecPeak = (-1 for _ in range(5))

        if sgnEnum == AliHFInvMassFitter.kGaus:
            if not (secPeak and particleName == 'Ds'):
                massFunc = TF1(f'massFunc{iPt}', SingleGaus, massMin, massMax, 3)
                massFunc.SetParameters(hMassForFit[iPt].Integral() * binWidth, massForFit, 0.010)
            else:
                massFunc = TF1(f'massFunc{iPt}', DoublePeakSingleGaus, massMin, massMax, 6)
                massFunc.SetParameters(hMassForFit[iPt].Integral() * binWidth, massForFit, 0.010,
                                       hMassForFit[iPt].Integral() * binWidth, massDplus, 0.010)
                parRawYieldSecPeak = 3
                parMeanSecPeak = 4
                parSigmaSecPeak = 5

        elif sgnEnum == AliHFInvMassFitter.k2Gaus:
            parSigma2 = 3
            parFrac2Gaus = 4
            if not (secPeak and particleName == 'Ds'):
                massFunc = TF1(f'massFunc{iPt}', DoubleGaus, massMin, massMax, 5)
                if not particleName == 'Dstar':
                    massFunc.SetParameters(hMassForFit[iPt].Integral() * binWidth, massForFit, 0.010, 0.030, 0.9)
                else:
                    massFunc.SetParameters(hMassForFit[iPt].Integral() * binWidth, massForFit, 0.0010, 0.0030, 0.9)
            else:
                massFunc = TF1(f'massFunc{iPt}', DoublePeakDoubleGaus, massMin, massMax, 8)
                massFunc.SetParameters(hMassForFit[iPt].Integral() * binWidth, massForFit, 0.010, 0.030, 0.9,
                                       hMassForFit[iPt].Integral() * binWidth, massDplus, 0.010)
                parRawYieldSecPeak = 5
                parMeanSecPeak = 6
                parSigmaSecPeak = 7
        else:
            print("ERROR: Only kGaus and k2Gaus are supported for MC. Exit!") #TODO: add support for k2GausSigmaRatioPar
            sys.exit()

        if nPtBins > 1:
            cMass[iCanv].cd(iPt-nMaxCanvases*iCanv+1)
        else:
            cMass[iCanv].cd()
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

        if secPeak and particleName == 'Ds':
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

        if sgnEnum == AliHFInvMassFitter.k2Gaus:
            sigma2 = massFunc.GetParameter(parSigma2)
            sigma2err = massFunc.GetParError(parSigma2)
            frac2gaus = massFunc.GetParameter(parFrac2Gaus)
            frac2gauserr = massFunc.GetParError(parFrac2Gaus)
            hRawYieldsSigma2.SetBinContent(iPt+1, sigma2)
            hRawYieldsSigma2.SetBinError(iPt+1, sigma2err)
            hRawYieldsFracGaus2.SetBinContent(iPt+1, frac2gaus)
            hRawYieldsFracGaus2.SetBinError(iPt+1, frac2gauserr)

    else:  # data
        massFitter.append(AliHFInvMassFitter(hMassForFit[iPt], massMin, massMax, bkgEnum, sgnEnum))
        if degPol[iPt] > 0:
            massFitter[iPt].SetPolDegreeForBackgroundFit(degPol[iPt])

        if fitConfig[cent]['UseLikelihood']:
            massFitter[iPt].SetUseLikelihoodFit()
        if fixMean[iPt]:
            massFitter[iPt].SetFixGaussianMean(hMeanToFix.GetBinContent(iPt+1))
        if fitConfig[cent]['BoundMean']:
            massFitter[iPt].SetBoundGaussianMean(massForFit, massMin, massMax)
        else:
            massFitter[iPt].SetInitialGaussianMean(massForFit)
        if fitConfig[cent]['FixSigmaRatio']:
            massFitter[iPt].SetFixRatio2GausSigma(
                hSigmaToFix.GetBinContent(iPt+1)/hSigmaToFix2.GetBinContent(iPt+1))

        if fixSigma[iPt]:
            if isinstance(fitConfig[cent]['SigmaMultFactor'], (float, int)):
                massFitter[iPt].SetFixGaussianSigma(
                    hSigmaToFix.GetBinContent(iPt+1)*fitConfig[cent]['SigmaMultFactor'])
            else:
                if fitConfig[cent]['SigmaMultFactor'] == 'MinusUnc':
                    massFitter[iPt].SetFixGaussianSigma(
                        hSigmaToFix.GetBinContent(iPt+1)-hSigmaToFix.GetBinError(iPt+1))
                elif fitConfig[cent]['SigmaMultFactor'] == 'PlusUnc':
                    massFitter[iPt].SetFixGaussianSigma(
                        hSigmaToFix.GetBinContent(iPt+1)+hSigmaToFix.GetBinError(iPt+1))
                else:
                    print('WARNING: impossible to fix sigma! Wrong mult factor set in config file!')
        else:
            if hSigmaToFix:
                massFitter[iPt].SetInitialGaussianSigma(
                    hSigmaToFix.GetBinContent(iPt+1)*fitConfig[cent]['SigmaMultFactor'])
            else:
                if particleName == 'Dstar':
                    massFitter[iPt].SetInitialGaussianSigma(0.001)
                else:
                    massFitter[iPt].SetInitialGaussianSigma(0.008)

        if secPeak and particleName == 'Ds':
            if hSigmaToFixSecPeak:
                sigmaToFix = hSigmaToFixSecPeak.GetBinContent(iPt+1) * fitConfig[cent]['SigmaMultFactorSecPeak']
                massFitter[iPt].IncludeSecondGausPeak(massDplus, False, sigmaToFix, True)
                if fitConfig[cent]['FixSigmaToFirstPeak']:
                    # fix D+ peak to sigmaMC(D+)/sigmaMC(Ds+)*sigmaData(Ds+)
                    massFitter[iPt].MassFitter(False)
                    sigmaFirstPeak = massFitter[iPt].GetSigma()
                    sigmaRatioMC = hSigmaToFixSecPeak.GetBinContent(iPt+1) / hSigmaFirstPeakMC.GetBinContent(iPt+1)
                    massFitter[iPt].IncludeSecondGausPeak(massDplus, False, sigmaRatioMC * sigmaFirstPeak, True)
            else:
                massFitter[iPt].IncludeSecondGausPeak(massDplus, False, fitConfig[cent]['SigmaSecPeak'][iPt], True)
        if enableRef:
            rOverS = hMassForSig[iPt].Integral(
                hMassForSig[iPt].FindBin(massMin * 1.0001),
                hMassForSig[iPt].FindBin(massMax * 0.999)
            )
            rOverS /= hMassForRel[iPt].Integral(
                hMassForRel[iPt].FindBin(massMin * 1.0001),
                hMassForRel[iPt].FindBin(massMax * 0.999)
            )
            massFitter[iPt].SetFixReflOverS(rOverS)
            massFitter[iPt].SetTemplateReflections(hRel[iPt], "2gaus", massMin, massMax)
        massFitter[iPt].MassFitter(False)

        rawyield = massFitter[iPt].GetRawYield()
        rawyielderr = massFitter[iPt].GetRawYieldError()
        sigma = massFitter[iPt].GetSigma()
        sigmaerr = massFitter[iPt].GetSigmaUncertainty()
        mean = massFitter[iPt].GetMean()
        meanerr = massFitter[iPt].GetMeanUncertainty()
        redchi2 = massFitter[iPt].GetReducedChiSquare()
        signif, signiferr = ctypes.c_double(), ctypes.c_double()
        sgn, sgnerr = ctypes.c_double(), ctypes.c_double()
        bkg, bkgerr = ctypes.c_double(), ctypes.c_double()
        massFitter[iPt].Significance(3, signif, signiferr)
        massFitter[iPt].Signal(3, sgn, sgnerr)
        massFitter[iPt].Background(3, bkg, bkgerr)

        hRawYields.SetBinContent(iPt+1, rawyield)
        hRawYields.SetBinError(iPt+1, rawyielderr)
        hRawYieldsSigma.SetBinContent(iPt+1, sigma)
        hRawYieldsSigma.SetBinError(iPt+1, sigmaerr)
        hRawYieldsMean.SetBinContent(iPt+1, mean)
        hRawYieldsMean.SetBinError(iPt+1, meanerr)
        hRawYieldsSignificance.SetBinContent(iPt+1, signif.value)
        hRawYieldsSignificance.SetBinError(iPt+1, signiferr.value)
        hRawYieldsSoverB.SetBinContent(iPt+1, sgn.value/bkg.value)
        hRawYieldsSoverB.SetBinError(iPt+1, sgn.value/bkg.value*np.sqrt(
            sgnerr.value**2/sgn.value**2+bkgerr.value**2/bkg.value**2))
        hRawYieldsSignal.SetBinContent(iPt+1, sgn.value)
        hRawYieldsSignal.SetBinError(iPt+1, sgnerr.value)
        hRawYieldsBkg.SetBinContent(iPt+1, bkg.value)
        hRawYieldsBkg.SetBinError(iPt+1, bkgerr.value)
        hRawYieldsChiSquare.SetBinContent(iPt+1, redchi2)
        hRawYieldsChiSquare.SetBinError(iPt+1, 1.e-20)

        for iS, _ in enumerate(nSigma4SandB):
            massFitter[iPt].Significance(nSigma4SandB[iS],signif,signiferr)
            massFitter[iPt].Signal(nSigma4SandB[iS],sgn,sgnerr)
            massFitter[iPt].Background(nSigma4SandB[iS],bkg,bkgerr)

            hRawYieldsSignalDiffSigma[iS].SetBinContent(iPt+1,sgn.value)
            hRawYieldsSignalDiffSigma[iS].SetBinError(iPt+1,sgnerr.value)
            hRawYieldsBkgDiffSigma[iS].SetBinContent(iPt+1,bkg.value)
            hRawYieldsBkgDiffSigma[iS].SetBinError(iPt+1,bkgerr.value)
            hRawYieldsSoverBDiffSigma[iS].SetBinContent(iPt+1,sgn.value/bkg.value)
            hRawYieldsSoverBDiffSigma[iS].SetBinError(
                iPt+1,
                sgn.value/bkg.value*np.sqrt(sgnerr.value**2/sgn.value**2+bkgerr.value**2/bkg.value**2)
            )
            hRawYieldsSignifDiffSigma[iS].SetBinContent(iPt+1,signif.value)
            hRawYieldsSignifDiffSigma[iS].SetBinError(iPt+1,signiferr.value)

        fTotFunc = massFitter[iPt].GetMassFunc()
        fBkgFunc = massFitter[iPt].GetBackgroundRecalcFunc()

        parFrac2Gaus, parsecondsigma = -1, -1
        if sgnEnum == AliHFInvMassFitter.k2Gaus:
            if not (inclSecPeak and particleName == 'Ds'):
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

        if inclSecPeak and particleName == 'Ds':
            paryieldSecPeak = fTotFunc.GetNpar()-3
            parMeanSecPeak = fTotFunc.GetNpar()-2
            parSigmaSecPeak = fTotFunc.GetNpar()-1

            rawyieldSecPeak = fTotFunc.GetParameter(paryieldSecPeak) / binWidth
            rawyieldSecPeakerr = fTotFunc.GetParError(paryieldSecPeak) / binWidth
            meanSecPeak = fTotFunc.GetParameter(parMeanSecPeak)
            meanSecPeakerr = fTotFunc.GetParError(parMeanSecPeak)
            sigmaSecPeak = fTotFunc.GetParameter(parSigmaSecPeak)
            sigmaSecPeakerr = fTotFunc.GetParError(parSigmaSecPeak)

            bkgSecPeak = fBkgFunc.Integral(meanSecPeak - 3*sigmaSecPeak, meanSecPeak + 3*sigmaSecPeak) / binWidth
            bkgSecPeakerr = np.sqrt(bkgSecPeak)
            signalSecPeak = fTotFunc.Integral(meanSecPeak - 3*sigmaSecPeak, meanSecPeak + 3*sigmaSecPeak) / binWidth \
                            - bkgSecPeak
            signalSecPeakerr = np.sqrt(signalSecPeak+bkgSecPeak)
            signifSecPeak, signifSecPeakerr = ctypes.c_double(), ctypes.c_double()
            AliVertexingHFUtils.ComputeSignificance(signalSecPeak, signalSecPeakerr, bkgSecPeak, bkgSecPeakerr,
                                                    signifSecPeak, signifSecPeakerr)

            hRawYieldsSecPeak.SetBinContent(iPt+1, rawyieldSecPeak)
            hRawYieldsSecPeak.SetBinError(iPt+1, rawyieldSecPeakerr)
            hRawYieldsMeanSecPeak.SetBinContent(iPt+1, meanSecPeak)
            hRawYieldsMeanSecPeak.SetBinError(iPt+1, meanSecPeakerr)
            hRawYieldsSigmaSecPeak.SetBinContent(iPt+1, sigmaSecPeak)
            hRawYieldsSigmaSecPeak.SetBinError(iPt+1, sigmaSecPeakerr)
            hRawYieldsSignificanceSecPeak.SetBinContent(iPt+1, signifSecPeak.value)
            hRawYieldsSignificanceSecPeak.SetBinError(iPt+1, signifSecPeakerr.value)
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
            cMass[iCanv].cd(iPt-nMaxCanvases*iCanv+1)
        else:
            cMass[iCanv].cd()

        hMassForFit[iPt].GetYaxis().SetRangeUser(hMassForFit[iPt].GetMinimum()*0.95, hMassForFit[iPt].GetMaximum()*1.2)
        massFitter[iPt].GetSignalFunc().SetNpx(500)
        massFitter[iPt].GetMassFunc().SetNpx(500)

        massFitter[iPt].DrawHere(gPad)

        # residuals
        if nPtBins > 1:
            cResiduals[iCanv].cd(iPt-nMaxCanvases*iCanv+1)
        else:
            cResiduals[iCanv].cd()
        massFitter[iPt].DrawHistoMinusFit(gPad)

    cMass[iCanv].Modified()
    cMass[iCanv].Update()
    cResiduals[iCanv].Modified()
    cResiduals[iCanv].Update()

#save output histos
outFile = TFile(args.outFileName, 'recreate')
for canv in cMass:
    canv.Write()
if not args.isMC:
    for canv in cResiduals:
        canv.Write()
for hist in hMass:
    hist.Write()
for fitter, ptLow, ptHigh in zip(massFitter, ptMins, ptMaxs):
    fitter.GetMassFunc().SetName(f'fTot_{ptLow}_{ptHigh}')
    fitter.GetSignalFunc().SetName(f'fSgn_{ptLow}_{ptHigh}')
    fitter.GetBackgroundRecalcFunc().SetName(f'fBkg_{ptLow}_{ptHigh}')
    fitter.GetMassFunc().Write()
    fitter.GetSignalFunc().Write()
    fitter.GetBackgroundRecalcFunc().Write()
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
if not args.isMC:
    dirSB = TDirectoryFile('SandBDiffNsigma', 'SandBDiffNsigma')
    dirSB.Write()
    dirSB.cd()
    for iS, _ in enumerate(nSigma4SandB):
        hRawYieldsSignalDiffSigma[iS].Write()
        hRawYieldsBkgDiffSigma[iS].Write()
        hRawYieldsSoverBDiffSigma[iS].Write()
        hRawYieldsSignifDiffSigma[iS].Write()
    dirSB.Close()
outFile.Close()

outFileNamePDF = args.outFileName.replace('.root', '.pdf')
outFileNameResPDF = outFileNamePDF.replace('.pdf', '_Residuals.pdf')
for iCanv, (cM, cR) in enumerate(zip(cMass, cResiduals)):
    if iCanv == 0 and nCanvases > 1:
        cM.SaveAs(f'{outFileNamePDF}[')
    cM.SaveAs(outFileNamePDF)
    if iCanv == nCanvases-1 and nCanvases > 1:
        cM.SaveAs(f'{outFileNamePDF}]')
    if not args.isMC:
        if iCanv == 0 and nCanvases > 1:
            cR.SaveAs(f'{outFileNameResPDF}[')
        cR.SaveAs(outFileNameResPDF)
        if iCanv == nCanvases-1 and nCanvases > 1:
            cR.SaveAs(f'{outFileNameResPDF}]')

if not args.batch:
    input('Press enter to exit')
