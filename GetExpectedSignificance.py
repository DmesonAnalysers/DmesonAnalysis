#*******************************************************************************************#
# python script for the expected significance of D+ and Ds+ mesons using THnSparses         #
# run: python GetExpectedSignificance.py cfgFileName.yml cutSetFile.yml outFileName         #
# author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino                            #
#*******************************************************************************************#

from ROOT import TFile, TCanvas, TH1F, TSpline3, TF1, TGraphAsymmErrors, TLine, TLegend
from ROOT import gROOT, gStyle
from ROOT import kWhite, kBlack, kOrange, kRed, kGreen, kBlue, kAzure, kFullCircle, kFullSquare, kFullDiamond, kFullTriangleUp, kFullTriangleDown
import yaml, sys, array, math
from ReadModel import ReadFONLL, ReadTAMU, ReadPHSD, ReadGossiaux, ReadCatania
from ReadHepData import ReadHepDataROOT

def GetExpectedBackground(hSB, mean, sigma, Nexp, Nanal) :
  fMassBkg = TF1('fMassBkg','expo',1.8,2.14)
  hSB.Fit('fMassBkg','QR0')
  B = fMassBkg.Integral(mean-3*sigma,mean+3*sigma) / hSB.GetBinWidth(1) * Nexp / Nanal
  return B

def GetExpectedSignal(DeltaPt, sigma, Raa, Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp) :
  corrYield = sigma * fractoD * Raa * Taa / DeltaPt
  rawYield = 2 * DeltaPt * corrYield * effPrompt * Acc * BR * Nexp / fprompt
  return rawYield, corrYield

def GetSideBandHisto(hMassData, mean, sigma, name):
  hMassDataReb = hMassData.Clone('hSB')
  hMassDataReb.Rebin(5)
  hSB = hMassDataReb.Clone(name)
  massbinmin = hSB.GetXaxis().FindBin(mean-4*sigma)
  massbinmax = hSB.GetXaxis().FindBin(mean+4*sigma)
  for iMassBin in range(1,hSB.GetNbinsX()) :
    if iMassBin >= massbinmin and iMassBin <= massbinmax :
      hSB.SetBinContent(iMassBin,0.)
      hSB.SetBinError(iMassBin,0.)
    else :
      hSB.SetBinContent(iMassBin,hMassDataReb.GetBinContent(iMassBin))
      hSB.SetBinError(iMassBin,hMassDataReb.GetBinError(iMassBin))
  return hSB

def SetGraphStyle(graph, color, colorfill, marker, linewidth=2):
  if colorfill == kWhite:
    graph.SetFillStyle(0)
  else:
    graph.SetFillColorAlpha(colorfill,0.25)
  graph.SetMarkerColor(color)
  graph.SetLineColor(color)
  graph.SetLineWidth(linewidth)
  graph.SetMarkerStyle(marker)

def SetSplineStyle(spline, color):
  spline.SetLineColor(color)
  spline.SetLineWidth(2)

colors = {'TAMU' : kOrange+10, 'PHSD' : kGreen+3, 'Gossiaux' : kAzure+4, 'Catania' : kRed+2}
colorsFill = {'TAMU' : kOrange+7, 'PHSD' : kGreen+2, 'Gossiaux' : kAzure+2, 'Catania' : kRed+1}
markers = {'TAMU' : kFullSquare, 'PHSD' : kFullCircle, 'Gossiaux' : kFullDiamond, 'Catania' : kFullTriangleUp}

######################################################################################################
#main function

inputCfg = sys.argv[1]
cutSetFile = sys.argv[2]
outFileName = sys.argv[3]

with open(inputCfg, 'r') as ymlinputCfg:
    inputCfg = yaml.load(ymlinputCfg)

infileDataLowPt = TFile(inputCfg['fileDataLowPt']['filename'])
indirDataLowPt = infileDataLowPt.Get(inputCfg['fileDataLowPt']['dirname'])
inlistDataLowPt = indirDataLowPt.Get(inputCfg['fileDataLowPt']['listname'])
sMassPtCutVarsLowPt = inlistDataLowPt.FindObject(inputCfg['fileDataLowPt']['sparsenameAll'])
hEvLowPt = inlistDataLowPt.FindObject(inputCfg['fileDataLowPt']['histoevname'])

infileDataHighPt = TFile(inputCfg['fileDataHighPt']['filename'])
indirDataHighPt = infileDataHighPt.Get(inputCfg['fileDataHighPt']['dirname'])
inlistDataHighPt = indirDataHighPt.Get(inputCfg['fileDataHighPt']['listname'])
sMassPtCutVarsHighPt = inlistDataHighPt.FindObject(inputCfg['fileDataHighPt']['sparsenameAll'])
hEvHighPt = inlistDataHighPt.FindObject(inputCfg['fileDataHighPt']['histoevname'])

infileMCLowPt = TFile(inputCfg['fileMCLowPt']['filename'])
indirMCLowPt = infileMCLowPt.Get(inputCfg['fileMCLowPt']['dirname'])
inlistMCLowPt = indirMCLowPt.Get(inputCfg['fileMCLowPt']['listname'])
sMassPtCutVarsPromptLowPt = inlistMCLowPt.FindObject(inputCfg['fileMCLowPt']['sparsenamePrompt'])
sMassPtCutVarsFDLowPt = inlistMCLowPt.FindObject(inputCfg['fileMCLowPt']['sparsenameFD'])
sGenPromptLowPt = inlistMCLowPt.FindObject(inputCfg['fileMCLowPt']['sparsenameGenPrompt'])
sGenFDLowPt = inlistMCLowPt.FindObject(inputCfg['fileMCLowPt']['sparsenameGenFD'])

infileMCHighPt = TFile(inputCfg['fileMCHighPt']['filename'])
indirMCHighPt = infileMCHighPt.Get(inputCfg['fileMCHighPt']['dirname'])
inlistMCHighPt = indirMCHighPt.Get(inputCfg['fileMCHighPt']['listname'])
sMassPtCutVarsPromptHighPt = inlistMCHighPt.FindObject(inputCfg['fileMCHighPt']['sparsenamePrompt'])
sMassPtCutVarsFDHighPt = inlistMCHighPt.FindObject(inputCfg['fileMCHighPt']['sparsenameFD'])
sGenPromptHighPt = inlistMCHighPt.FindObject(inputCfg['fileMCHighPt']['sparsenameGenPrompt'])
sGenFDHighPt = inlistMCLowPt.FindObject(inputCfg['fileMCLowPt']['sparsenameGenFD'])

PtThreshold = inputCfg['PtThreshold']

nEvExp = inputCfg['nExpectedEvents']
Taa = inputCfg['Taa']
fprompt = inputCfg['fprompt']
BR = inputCfg['BR']
fractoD = inputCfg['fractoD']

infileAcc = TFile(inputCfg['filenameAcc'])
hPtGenAcc = infileAcc.Get('hPtGenAcc')
hPtGenLimAcc = infileAcc.Get('hPtGenLimAcc')

FONLL = ReadFONLL(inputCfg['filenameFONLL'])
TAMU = ReadTAMU(inputCfg['filenameRaaPredTAMU'])
PHSD = ReadPHSD(inputCfg['filenameRaaPredPHSD'])
Gossiaux = ReadGossiaux(inputCfg['filenameRaaPredGossiaux'])
Catania = ReadCatania(inputCfg['filenameRaaPredCatania'])
  
outfile = TFile('%s.root' % outFileName,'recreate')

sTAMU, sGossiaux, sPHSD, sCatania = ({} for iDic in range(4))
for iRaaTAMU in TAMU:
  if iRaaTAMU != 'PtCent':
    sTAMU[iRaaTAMU] = TSpline3('sTAMU%s' % iRaaTAMU,array.array('d',TAMU['PtCent']),array.array('d',TAMU[iRaaTAMU]),len(TAMU['PtCent']))
for iRaaPHSD in PHSD:
  if iRaaPHSD != 'PtCent':
    sPHSD[iRaaPHSD] = TSpline3('sPHSD%s' % iRaaPHSD,array.array('d',PHSD['PtCent']),array.array('d',PHSD[iRaaPHSD]),len(PHSD['PtCent']))
for iRaaGossiaux in Gossiaux:
  if iRaaGossiaux != 'PtCent':
    sGossiaux[iRaaGossiaux] = TSpline3('sGossiaux%s' % iRaaGossiaux,array.array('d',Gossiaux['PtCent']),array.array('d',Gossiaux[iRaaGossiaux]),len(Gossiaux['PtCent']))
for iRaaCatania in Catania:
  if iRaaCatania != 'PtCent':
    sCatania[iRaaCatania] = TSpline3('sCatania%s' % iRaaCatania,array.array('d',Catania['PtCent']),array.array('d',Catania[iRaaCatania]),len(Catania['PtCent']))

with open(cutSetFile, 'r') as ymlCutSetFile:
    cutSetFile = yaml.load(ymlCutSetFile)
cutVars = cutSetFile['cutvars']

gSignifTAMU, gSignifPHSD, gSignifGossiaux, gSignifCatania = (TGraphAsymmErrors(0) for iGraph in range(4))
gCorrYieldTAMU, gCorrYieldPHSD, gCorrYieldGossiaux, gCorrYieldCatania = (TGraphAsymmErrors(0) for iGraph in range(4))

gSignifTAMU.SetTitle(';#it{p}_{T} (GeV/#it{c});significance')
gSignifPHSD.SetTitle(';#it{p}_{T} (GeV/#it{c});significance')
gSignifGossiaux.SetTitle(';#it{p}_{T} (GeV/#it{c});significance')
gSignifCatania.SetTitle(';#it{p}_{T} (GeV/#it{c});significance')

gCorrYieldTAMU.SetTitle(';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}')
gCorrYieldPHSD.SetTitle(';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}')
gCorrYieldGossiaux.SetTitle(';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}')
gCorrYieldCatania.SetTitle(';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T}')

cutVars['Pt']['limits'] = list(cutVars['Pt']['min'])
cutVars['Pt']['limits'].append(cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1])

hAccEffPrompt = TH1F('hAccEffPrompt',';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon) #times 2#it{y}_{fid}',len(cutVars['Pt']['limits'])-1,array.array('f',cutVars['Pt']['limits']))
hAccEffFD = TH1F('hAccEffFD',';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon) #times 2#it{y}_{fid}',len(cutVars['Pt']['limits'])-1,array.array('f',cutVars['Pt']['limits']))
SetGraphStyle(hAccEffPrompt,kRed+1,kWhite,kFullSquare)
SetGraphStyle(hAccEffFD,kBlue+1,kWhite,kFullCircle)

for iPt in range(len(cutVars['Pt']['min'])):
  #check if low or high pt
  if cutVars['Pt']['max'][iPt] <= PtThreshold:
    nEvBkg = hEvLowPt.GetBinContent(5)
    sMassPtCutVars = sMassPtCutVarsLowPt.Clone('sMassPtCutVars')
    sMassPtCutVarsPrompt = sMassPtCutVarsPromptLowPt.Clone('sMassPtCutVars')
    sMassPtCutVarsFD = sMassPtCutVarsFDLowPt.Clone('sMassPtCutVars')
    sGenPrompt = sGenPromptLowPt.Clone('sMassPtCutVars')
    sGenFD = sGenPromptLowPt.Clone('sMassPtCutVars')
  else:
    nEvBkg = hEvHighPt.GetBinContent(5)
    sMassPtCutVars = sMassPtCutVarsHighPt.Clone('sMassPtCutVars')
    sMassPtCutVarsPrompt = sMassPtCutVarsPromptHighPt.Clone('sMassPtCutVars')
    sMassPtCutVarsFD = sMassPtCutVarsFDHighPt.Clone('sMassPtCutVars')
    sGenPrompt = sGenPromptHighPt.Clone('sMassPtCutVars')
    sGenFD = sGenPromptHighPt.Clone('sMassPtCutVars')

  #gen for efficiency
  binGenPtMin = sGenPrompt.GetAxis(0).FindBin(cutVars['Pt']['min'][iPt]*1.0001)
  binGenPtMax = sGenPrompt.GetAxis(0).FindBin(cutVars['Pt']['max'][iPt]*0.9999)
  hPtGenPrompt = sGenPrompt.GetAxis(0).SetRange(binGenPtMin,binGenPtMax)
  hPtGenPrompt = sGenPrompt.Projection(0)
  hPtGenPrompt.SetName('hPtPrompt_pT_%0.f_%0.f' % (cutVars['Pt']['min'][iPt],cutVars['Pt']['max'][iPt]))
  nGenPrompt = hPtGenPrompt.Integral()
  hPtGenFD = sGenFD.GetAxis(0).SetRange(binGenPtMin,binGenPtMax)
  hPtGenFD = sGenFD.Projection(0)
  hPtGenFD.SetName('hPtFD_pT_%0.f_%0.f' % (cutVars['Pt']['min'][iPt],cutVars['Pt']['max'][iPt]))
  nGenFD = hPtGenFD.Integral()
  
  #acceptance from TOY MC
  binAccPtMin = hPtGenAcc.GetXaxis().FindBin(cutVars['Pt']['min'][iPt]*1.0001)
  binAccPtMax = hPtGenAcc.GetXaxis().FindBin(cutVars['Pt']['max'][iPt]*0.9999)
  Acc = hPtGenAcc.Integral(binAccPtMin, binAccPtMax) / hPtGenLimAcc.Integral(binAccPtMin, binAccPtMax)

  #FONLL
  binFONLL = FONLL['PtMin'].index(cutVars['Pt']['min'][iPt])
  sigmaFONLL = []
  sigmaFONLL.append(FONLL['Cent'][binFONLL])
  sigmaFONLL.append(FONLL['Min'][binFONLL])
  sigmaFONLL.append(FONLL['Max'][binFONLL])

  #RAA
  RaaTAMU, RaaPHSD, RaaGossiaux, RaaCatania = ({} for iDic in range(4))
  for iRaaTAMU in sTAMU:
    RaaTAMU[iRaaTAMU] = sTAMU[iRaaTAMU].Eval((cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2)
  for iRaaPHSD in sPHSD:
    RaaPHSD[iRaaPHSD] = sPHSD[iRaaPHSD].Eval((cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2)
  for iRaaGossiaux in sGossiaux:
    RaaGossiaux[iRaaGossiaux] = sGossiaux[iRaaGossiaux].Eval((cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2)
  for iRaaCatania in sCatania:
    RaaCatania[iRaaCatania] = sCatania[iRaaCatania].Eval((cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2)

  for iVar in cutVars:
    if iVar == 'InvMass':
      continue
    binMin = sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min'][iPt]*1.0001)
    binMax = sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max'][iPt]*0.9999)
    sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin,binMax)
    sMassPtCutVarsPrompt.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin,binMax)
    sMassPtCutVarsFD.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin,binMax)
  
  hMassData = sMassPtCutVars.Projection(0)
  hMassData.SetName('hMassData_pT_%0.f_%0.f' % (cutVars['Pt']['min'][iPt],cutVars['Pt']['max'][iPt]))
  hMassPrompt = sMassPtCutVarsPrompt.Projection(0)
  hMassPrompt.SetName('hMassPrompt_pT_%0.f_%0.f' % (cutVars['Pt']['min'][iPt],cutVars['Pt']['max'][iPt]))
  hMassFD = sMassPtCutVarsFD.Projection(0)
  hMassFD.SetName('hMassFD_pT_%0.f_%0.f' % (cutVars['Pt']['min'][iPt],cutVars['Pt']['max'][iPt]))
  hMassSgn = hMassPrompt.Clone('hMassSgn_pT_%0.f_%0.f' % (cutVars['Pt']['min'][iPt],cutVars['Pt']['max'][iPt]))  
  hMassSgn.Add(hMassFD)
  
  fMassSgn = TF1('fMassSgn','gaus',1.7,2.15)
  hMassSgn.Fit('fMassSgn','Q0')
  mean = fMassSgn.GetParameter(1)
  sigma = fMassSgn.GetParameter(2)

  hMassSB = GetSideBandHisto(hMassData,mean,sigma,'hMassSB_pT_%0.f_%0.f' % (cutVars['Pt']['min'][iPt],cutVars['Pt']['max'][iPt]))
  B = GetExpectedBackground(hMassSB,mean,sigma,nEvExp,nEvBkg)

  nRecoPrompt = hMassPrompt.Integral()
  nRecoFD = hMassFD.Integral()
  effPrompt = nRecoPrompt / nGenPrompt
  effFD = nRecoFD / nGenFD
  effPromptUnc = math.sqrt(effPrompt*(1-effPrompt)/nGenPrompt)
  effFDUnc = math.sqrt(effFD*(1-effFD)/nGenFD)

  hAccEffPrompt.SetBinContent(iPt+1,effPrompt)
  hAccEffPrompt.SetBinError(iPt+1,effPromptUnc)
  hAccEffFD.SetBinContent(iPt+1,effFD)
  hAccEffFD.SetBinError(iPt+1,effFDUnc)

  SignifTAMU, SignifPHSD, SignifGossiaux, SignifCatania  = ({} for iDic in range(4))
  CorrYieldTAMU, CorrYieldPHSD, CorrYieldGossiaux, CorrYieldCatania = ({} for iDic in range(4))
  for iFONLL in FONLL:
    if iFONLL == 'Cent' or iFONLL == 'Min' or iFONLL == 'Max':

      SignifTAMU[iFONLL], CorrYieldTAMU[iFONLL] = (dict(RaaTAMU) for iDic in range(2)) 
      SignifPHSD[iFONLL], CorrYieldPHSD[iFONLL] = (dict(RaaPHSD) for iDic in range(2)) 
      SignifGossiaux[iFONLL], CorrYieldGossiaux[iFONLL] = (dict(RaaGossiaux) for iDic in range(2)) 
      SignifCatania[iFONLL], CorrYieldCatania[iFONLL] = (dict(RaaCatania) for iDic in range(2)) 

      for iRaaTAMU in RaaTAMU:
        rawyield, corryield = GetExpectedSignal(cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt],FONLL[iFONLL][binFONLL],RaaTAMU[iRaaTAMU],Taa,effPrompt,Acc,fprompt,BR,fractoD,nEvExp)
        SignifTAMU[iFONLL][iRaaTAMU] = rawyield/math.sqrt(rawyield+B)
        CorrYieldTAMU[iFONLL][iRaaTAMU] = corryield
      for iRaaPHSD in RaaPHSD:
        rawyield, corryield = GetExpectedSignal(cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt],FONLL[iFONLL][binFONLL],RaaPHSD[iRaaPHSD],Taa,effPrompt,Acc,fprompt,BR,fractoD,nEvExp)
        SignifPHSD[iFONLL][iRaaPHSD] = rawyield/math.sqrt(rawyield+B)
        CorrYieldPHSD[iFONLL][iRaaPHSD] = corryield
      for iRaaGossiaux in RaaGossiaux:
        rawyield, corryield = GetExpectedSignal(cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt],FONLL[iFONLL][binFONLL],RaaGossiaux[iRaaGossiaux],Taa,effPrompt,Acc,fprompt,BR,fractoD,nEvExp)
        SignifGossiaux[iFONLL][iRaaGossiaux] = rawyield/math.sqrt(rawyield+B)
        CorrYieldGossiaux[iFONLL][iRaaGossiaux] = corryield
      for iRaaCatania in RaaCatania:
        rawyield, corryield = GetExpectedSignal(cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt],FONLL[iFONLL][binFONLL],RaaCatania[iRaaCatania],Taa,effPrompt,Acc,fprompt,BR,fractoD,nEvExp)
        SignifCatania[iFONLL][iRaaCatania] = rawyield/math.sqrt(rawyield+B)
        CorrYieldCatania[iFONLL][iRaaCatania] = corryield

  PtCent = (cutVars['Pt']['min'][iPt]+cutVars['Pt']['max'][iPt])/2
  PtUnc = (cutVars['Pt']['max'][iPt]-cutVars['Pt']['min'][iPt])/2

  gSignifTAMU.SetPoint(iPt,PtCent,SignifTAMU['Max']['Min'])
  gSignifTAMU.SetPointError(iPt,PtUnc,PtUnc,SignifTAMU['Max']['Min']-SignifTAMU['Min'][min(SignifTAMU['Min'],key=SignifTAMU['Min'].get)],SignifTAMU['Max'][max(SignifTAMU['Max'],key=SignifTAMU['Max'].get)]-SignifTAMU['Max']['Min'])
  gSignifPHSD.SetPoint(iPt,PtCent,SignifPHSD['Max']['Cent'])
  gSignifPHSD.SetPointError(iPt,PtUnc,PtUnc,SignifPHSD['Max']['Cent']-SignifPHSD['Min'][min(SignifPHSD['Min'],key=SignifPHSD['Min'].get)],SignifPHSD['Max'][max(SignifPHSD['Max'],key=SignifPHSD['Max'].get)]-SignifPHSD['Max']['Cent'])
  gSignifGossiaux.SetPoint(iPt,PtCent,SignifGossiaux['Max']['ColRad'])
  gSignifGossiaux.SetPointError(iPt,PtUnc,PtUnc,SignifGossiaux['Max']['ColRad']-SignifGossiaux['Min'][min(SignifGossiaux['Min'],key=SignifGossiaux['Min'].get)],SignifGossiaux['Max'][max(SignifGossiaux['Max'],key=SignifGossiaux['Max'].get)]-SignifGossiaux['Max']['ColRad'])
  gSignifCatania.SetPoint(iPt,PtCent,SignifCatania['Max']['Cent'])
  gSignifCatania.SetPointError(iPt,PtUnc,PtUnc,SignifCatania['Max']['Cent']-SignifCatania['Min'][min(SignifCatania['Min'],key=SignifCatania['Min'].get)],SignifCatania['Max'][max(SignifCatania['Max'],key=SignifCatania['Max'].get)]-SignifCatania['Max']['Cent'])

  gCorrYieldTAMU.SetPoint(iPt,PtCent,CorrYieldTAMU['Max']['Min'])
  gCorrYieldTAMU.SetPointError(iPt,PtUnc,PtUnc,CorrYieldTAMU['Max']['Min']-CorrYieldTAMU['Min'][min(CorrYieldTAMU['Min'],key=CorrYieldTAMU['Min'].get)],CorrYieldTAMU['Max'][max(CorrYieldTAMU['Max'],key=CorrYieldTAMU['Max'].get)]-CorrYieldTAMU['Max']['Min'])
  gCorrYieldPHSD.SetPoint(iPt,PtCent,CorrYieldPHSD['Max']['Cent'])
  gCorrYieldPHSD.SetPointError(iPt,PtUnc,PtUnc,CorrYieldPHSD['Max']['Cent']-CorrYieldPHSD['Min'][min(CorrYieldPHSD['Min'],key=CorrYieldPHSD['Min'].get)],CorrYieldPHSD['Max'][max(CorrYieldPHSD['Max'],key=CorrYieldPHSD['Max'].get)]-CorrYieldPHSD['Max']['Cent'])
  gCorrYieldGossiaux.SetPoint(iPt,PtCent,CorrYieldGossiaux['Max']['ColRad'])
  gCorrYieldGossiaux.SetPointError(iPt,PtUnc,PtUnc,CorrYieldGossiaux['Max']['ColRad']-CorrYieldGossiaux['Min'][min(CorrYieldGossiaux['Min'],key=CorrYieldGossiaux['Min'].get)],CorrYieldGossiaux['Max'][max(CorrYieldGossiaux['Max'],key=CorrYieldGossiaux['Max'].get)]-CorrYieldGossiaux['Max']['ColRad'])
  gCorrYieldCatania.SetPoint(iPt,PtCent,CorrYieldCatania['Max']['Cent'])
  gCorrYieldCatania.SetPointError(iPt,PtUnc,PtUnc,CorrYieldCatania['Max']['Cent']-CorrYieldCatania['Min'][min(CorrYieldCatania['Min'],key=CorrYieldCatania['Min'].get)],CorrYieldCatania['Max'][max(CorrYieldCatania['Max'],key=CorrYieldCatania['Max'].get)]-CorrYieldCatania['Max']['Cent'])

  hMassData.Write()
  hMassPrompt.Write()
  hMassFD.Write()
  hMassSB.Write()

#load published result
hDsRaaStat, gDsRaaSyst = ReadHepDataROOT(inputCfg['PublishedResult']['filenameHepData'],inputCfg['PublishedResult']['tableRaa'])
hDsCorrYieldStat, gDsCorrYieldSyst = ReadHepDataROOT(inputCfg['PublishedResult']['filenameHepData'],inputCfg['PublishedResult']['tableCorrYield'])
SetGraphStyle(gDsRaaSyst,kBlack,kWhite,kFullSquare,1)
SetGraphStyle(gDsCorrYieldSyst,kBlack,kWhite,kFullSquare,1)
SetGraphStyle(hDsRaaStat,kBlack,kWhite,kFullSquare,1)
SetGraphStyle(hDsCorrYieldStat,kBlack,kWhite,kFullSquare,1)

#plots 
gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadTopMargin(0.035)
gStyle.SetTitleSize(0.045,'xy')
gStyle.SetLabelSize(0.040,'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)

legTAMU, legPHSD, legGossiaux, legCatania = ({'CorrYield' : TLegend(0.45,0.7,0.65,0.85), 'Raa' : TLegend(0.45,0.7,0.65,0.85), 'Signif' : TLegend(0.45,0.75,0.65,0.85)} for iDic in range(4))

for iLeg in legTAMU:
  legTAMU[iLeg].SetFillStyle(0)
  legTAMU[iLeg].SetTextSize(0.045)
  legPHSD[iLeg].SetFillStyle(0)
  legPHSD[iLeg].SetTextSize(0.045)
  legGossiaux[iLeg].SetFillStyle(0)
  legGossiaux[iLeg].SetTextSize(0.045)
  legCatania[iLeg].SetFillStyle(0)
  legCatania[iLeg].SetTextSize(0.045)
  if iLeg != 'Signif':
    legTAMU[iLeg].AddEntry(hDsCorrYieldStat,'ALICE - JHEP 10 (2018) 174','pl')
    legPHSD[iLeg].AddEntry(hDsCorrYieldStat,'ALICE - JHEP 10 (2018) 174','pl')
    legGossiaux[iLeg].AddEntry(hDsCorrYieldStat,'ALICE - JHEP 10 (2018) 174','pl')
    legCatania[iLeg].AddEntry(hDsCorrYieldStat,'ALICE - JHEP 10 (2018) 174','pl')
  
legTAMU['CorrYield'].AddEntry(gCorrYieldTAMU,'FONLL x TAMU','fp')
legTAMU['Raa'].AddEntry(sTAMU['Min'],'TAMU','l')
legTAMU['Signif'].AddEntry(gSignifTAMU,'FONLL x TAMU','fp')
legPHSD['CorrYield'].AddEntry(gCorrYieldPHSD,'FONLL x PHSD','fp')
legPHSD['Raa'].AddEntry(sPHSD['Cent'],'PHSD','l')
legPHSD['Signif'].AddEntry(gSignifPHSD,'FONLL x PHSD','fp')
legGossiaux['CorrYield'].AddEntry(gCorrYieldGossiaux,'FONLL x MC@sHQ+EPOS2','fp')
legGossiaux['Raa'].AddEntry(sGossiaux['Col'],'MC@sHQ+EPOS2','l')
legGossiaux['Signif'].AddEntry(gSignifGossiaux,'FONLL x MC@sHQ+EPOS2','fp')
legCatania['CorrYield'].AddEntry(gCorrYieldCatania,'FONLL x Catania (frag+coal)','fp')
legCatania['Raa'].AddEntry(sCatania['Cent'],'Catania (frag+coal)','l')
legCatania['Signif'].AddEntry(gSignifCatania,'FONLL x Catania (frag+coal)','fp')

lineAtOne = TLine(cutVars['Pt']['min'][0]-1,1.,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,1.)
lineAtOne.SetLineColor(kBlack)
lineAtOne.SetLineStyle(9)

outfile.cd()

cRaa = TCanvas('cRaa','',1000,1000)
cRaa.Divide(2,2)
cRaa.cd(1).DrawFrame(cutVars['Pt']['min'][0]-1,0,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,2,';#it{p}_{T} (GeV/#it{c});#it{R}_{AA}')
lineAtOne.Draw('same')
for iRaaTAMU in sTAMU:
  SetSplineStyle(sTAMU[iRaaTAMU],colorsFill['TAMU'])
  sTAMU[iRaaTAMU].Draw('same')
  sTAMU[iRaaTAMU].Write('sTAMU%s' % iRaaTAMU)
hDsRaaStat.Draw('same')
gDsRaaSyst.Draw('2')
legTAMU['Raa'].Draw('same')
cRaa.cd(2).DrawFrame(cutVars['Pt']['min'][0]-1,0,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,2,';#it{p}_{T} (GeV/#it{c});#it{R}_{AA}')
lineAtOne.Draw('same')
for iRaaPHSD in sPHSD:
  SetSplineStyle(sPHSD[iRaaPHSD],colorsFill['PHSD'])
  sPHSD[iRaaPHSD].Draw('same')
  sPHSD[iRaaPHSD].Write('sPHSD%s' % iRaaPHSD)
hDsRaaStat.Draw('same')
gDsRaaSyst.Draw('2')
legPHSD['Raa'].Draw('same')
cRaa.cd(3).DrawFrame(cutVars['Pt']['min'][0]-1,0,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,2,';#it{p}_{T} (GeV/#it{c});#it{R}_{AA}')
lineAtOne.Draw('same')
for iRaaGossiaux in sGossiaux:
  SetSplineStyle(sGossiaux[iRaaGossiaux],colorsFill['Gossiaux'])
  sGossiaux[iRaaGossiaux].Draw('same')
  sGossiaux[iRaaGossiaux].Write('sGossiaux%s' % iRaaGossiaux)
hDsRaaStat.Draw('same')
gDsRaaSyst.Draw('2')
legGossiaux['Raa'].Draw('same')
cRaa.cd(4).DrawFrame(cutVars['Pt']['min'][0]-1,0,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,2,';#it{p}_{T} (GeV/#it{c});#it{R}_{AA}')
lineAtOne.Draw('same')
for iRaaCatania in sCatania:
  SetSplineStyle(sCatania[iRaaCatania],colorsFill['Catania'])
  sCatania[iRaaCatania].Draw('same')
  sCatania[iRaaCatania].Write('sCatania%s' % iRaaCatania)
hDsRaaStat.Draw('same')
gDsRaaSyst.Draw('2')
legCatania['Raa'].Draw('same')
cRaa.Write()
cRaa.SaveAs('%s_Raa.pdf' % outFileName)

cCorrYield = TCanvas('cCorrYield','',1000,1000)
cCorrYield.Divide(2,2)
cCorrYield.cd(1).DrawFrame(cutVars['Pt']['min'][0]-1,3.e-6,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,2.,';#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}')
cCorrYield.cd(1).SetLogy()
SetGraphStyle(gCorrYieldTAMU, colors['TAMU'], colorsFill['TAMU'], markers['TAMU'],0)
gCorrYieldTAMU.Draw('P2')
hDsCorrYieldStat.Draw('same')
gDsCorrYieldSyst.Draw('2')
legTAMU['CorrYield'].Draw('same')
gCorrYieldTAMU.Write('gCorrYieldTAMU')
cCorrYield.cd(2).SetLogy()
cCorrYield.cd(2).DrawFrame(cutVars['Pt']['min'][0]-1,3.e-6,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,2.,';#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}')
SetGraphStyle(gCorrYieldPHSD, colors['PHSD'], colorsFill['PHSD'], markers['PHSD'],0)
gCorrYieldPHSD.Draw('P2')
hDsCorrYieldStat.Draw('same')
gDsCorrYieldSyst.Draw('2')
legPHSD['CorrYield'].Draw('same')
gCorrYieldPHSD.Write('gCorrYieldPHSD')
cCorrYield.cd(3).SetLogy()
cCorrYield.cd(3).DrawFrame(cutVars['Pt']['min'][0]-1,3.e-6,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,2.,';#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}')
SetGraphStyle(gCorrYieldGossiaux, colors['Gossiaux'], colorsFill['Gossiaux'], markers['Gossiaux'],0)
gCorrYieldGossiaux.Draw('P2')
hDsCorrYieldStat.Draw('same')
gDsCorrYieldSyst.Draw('2')
legGossiaux['CorrYield'].Draw('same')
gCorrYieldGossiaux.Write('gCorrYieldGossiaux')
cCorrYield.cd(4).SetLogy()
cCorrYield.cd(4).DrawFrame(cutVars['Pt']['min'][0]-1,3.e-6,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,2.,';#it{p}_{T} (GeV/#it{c});d#it{N}/d#it{p}_{T}')
SetGraphStyle(gCorrYieldCatania, colors['Catania'], colorsFill['Catania'], markers['Catania'],0)
gCorrYieldCatania.Draw('P2')
hDsCorrYieldStat.Draw('same')
gDsCorrYieldSyst.Draw('2')
gCorrYieldCatania.Write('gCorrYieldCatania')
legCatania['CorrYield'].Draw('same')
cCorrYield.Write()
cCorrYield.SaveAs('%s_CorrYield.pdf' % outFileName)

cSignif = TCanvas('cSignif','',1000,1000)
cSignif.Divide(2,2)
cSignif.cd(1).DrawFrame(cutVars['Pt']['min'][0]-1,0,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,25,';#it{p}_{T} (GeV/#it{c});significance')
SetGraphStyle(gSignifTAMU, colors['TAMU'], colorsFill['TAMU'], markers['TAMU'],0)
gSignifTAMU.Draw('P2')
legTAMU['Signif'].Draw('same')
gSignifTAMU.Write('gSignifTAMU')
cSignif.cd(2).DrawFrame(cutVars['Pt']['min'][0]-1,0,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,25,';#it{p}_{T} (GeV/#it{c});significance')
SetGraphStyle(gSignifPHSD, colors['PHSD'], colorsFill['PHSD'], markers['PHSD'],0)
gSignifPHSD.Draw('P2')
legPHSD['Signif'].Draw('same')
gSignifPHSD.Write('gSignifPHSD')
cSignif.cd(3).DrawFrame(cutVars['Pt']['min'][0]-1,0,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,25,';#it{p}_{T} (GeV/#it{c});significance')
SetGraphStyle(gSignifGossiaux, colors['Gossiaux'], colorsFill['Gossiaux'], markers['Gossiaux'],0)
gSignifGossiaux.Draw('P2')
legGossiaux['Signif'].Draw('same')
gSignifGossiaux.Write('gSignifGossiaux')
cSignif.cd(4).DrawFrame(cutVars['Pt']['min'][0]-1,0,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1]+1,25,';#it{p}_{T} (GeV/#it{c});significance')
SetGraphStyle(gSignifCatania, colors['Catania'], colorsFill['Catania'], markers['Catania'],0)
gSignifCatania.Draw('P2')
legCatania['Signif'].Draw('same')
gSignifCatania.Write('gSignifCatania')
cSignif.Write()
cSignif.SaveAs('%s_Significance.pdf' % outFileName)

cEff = TCanvas('cEff','',500,500)
cEff.DrawFrame(cutVars['Pt']['min'][0],1.e-4,cutVars['Pt']['max'][len(cutVars['Pt']['max'])-1],1.,';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon) #times 2#it{y}_{fid}')
cEff.SetLogy()
hAccEffPrompt.Draw('same')
hAccEffFD.Draw('same')
hAccEffPrompt.Write()
hAccEffFD.Write()
cEff.Write()
cEff.SaveAs('%s_Efficiency.pdf' % outFileName)