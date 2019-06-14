#********************************************************************************#
# python script for the significance scan of D+ and Ds+ mesons with THnSparses   #
# run: python ScanSignificanceSparse.py cfgFileName.yml outFileName.root         #
# author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino                 #
#********************************************************************************#

from ROOT import TFile, TCanvas, TH1F, TF1, TNtuple, TGraph, TSpline3, gROOT # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kGreen, kOrange # pylint: disable=import-error,no-name-in-module
import yaml, sys, itertools, math, array
from ReadModel import ReadFONLL, ReadTAMU

def ApplyCuts(sparse,bins,axesnum,upperlowercuts,name) :
  index = ''
  cutvalues = []
  for iVar in range(0,len(axesnum)) :
    if upperlowercuts[iVar] == 'Upper' :
      sparse.GetAxis(axesnum[iVar]).SetRange(0,bins[iVar])
      cutvalues.append(sparse.GetAxis(axesnum[iVar]).GetBinLowEdge(bins[iVar])+sparse.GetAxis(axesnum[iVar]).GetBinWidth(bins[iVar]))
    elif upperlowercuts[iVar] == 'Lower' :
      sparse.GetAxis(axesnum[iVar]).SetRange(bins[iVar],sparse.GetAxis(axesnum[iVar]).GetNbins()+1) 
      cutvalues.append(sparse.GetAxis(axesnum[iVar]).GetBinLowEdge(bins[iVar]))
    index += '%d' % bins[iVar]
  hProj = sparse.Projection(0)
  hProj.SetName(name)
  for iAxis in axesnum :
    sparse.GetAxis(iAxis).SetRange(-1,-1)
  return hProj, index, cutvalues

def GetSideBandHisto(hMassData, mean, sigma):
  hMassDataReb = hMassData.Clone('hSB')
  hMassDataReb.Rebin(5)
  hSB = hMassDataReb.Clone('hSB')
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

def GetExpectedBackground(hSB, mean, sigma, Nexp, Nanal) :
  fMassBkg = TF1('fMassBkg','expo',1.8,2.12)
  hSB.Fit('fMassBkg','QR')
  B = fMassBkg.Integral(mean-3*sigma,mean+3*sigma) / hSB.GetBinWidth(1) * Nexp / Nanal
  return B

def GetExpectedSignal(DeltaPt, sigma, Raa, Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp) :
  corrYield = sigma * fractoD * Raa * Taa / DeltaPt
  rawYield = 2 * DeltaPt * corrYield * effPrompt * Acc * BR * Nexp / fprompt
  return rawYield

#***************************************************************************
# Main function

cfgFileName = sys.argv[1]
outFileName = sys.argv[2]

with open(cfgFileName, 'r') as ymlCfgFile:
  inputCfg = yaml.load(ymlCfgFile)

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

PtThreshold = inputCfg['PtThreshold']

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

Nexp = inputCfg['nExpectedEvents']
Taa = inputCfg['Taa']
fprompt = inputCfg['fprompt']
BR = inputCfg['BR']
fractoD = inputCfg['fractoD']

FONLL = ReadFONLL(inputCfg['filenameFONLL'])
TAMU = ReadTAMU(inputCfg['filenameRaaPredTAMU']) #optimisation performed using TAMU
sTAMUmin = TSpline3('sTAMUmin',array.array('d',TAMU['PtCent']),array.array('d',TAMU['Min']),len(TAMU['PtCent']))
sTAMUmax = TSpline3('sTAMUmax',array.array('d',TAMU['PtCent']),array.array('d',TAMU['Max']),len(TAMU['PtCent']))

infileAcc = TFile(inputCfg['filenameAcc'])
hPtGenAcc = infileAcc.Get('hPtGenAcc')
hPtGenLimAcc = infileAcc.Get('hPtGenLimAcc')

cutVars = inputCfg['cutvars']
PtMin = inputCfg['PtMin']
PtMax = inputCfg['PtMax']

ranges = []
axesnum = []
upperlowercuts = []
steps = []

for varnum, iVar in enumerate(cutVars):
  if cutVars[iVar]['upperlowercut'] == 'Lower' :
    cutVars[iVar]['binmin'] = sMassPtCutVarsLowPt.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min']*1.0001)
    cutVars[iVar]['binmax'] = sMassPtCutVarsLowPt.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max']*1.0001)
  elif cutVars[iVar]['upperlowercut'] == 'Upper' :
    cutVars[iVar]['binmin'] = sMassPtCutVarsLowPt.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min']*0.9999)
    cutVars[iVar]['binmax'] = sMassPtCutVarsLowPt.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max']*0.9999)
  steps.append(int(cutVars[iVar]['step'] / sMassPtCutVarsLowPt.GetAxis(cutVars[iVar]['axisnum']).GetBinWidth(1)))
  ranges.append(range(cutVars[iVar]['binmin'],cutVars[iVar]['binmax']+1,steps[varnum]))
  axesnum.append(cutVars[iVar]['axisnum'])
  upperlowercuts.append(cutVars[iVar]['upperlowercut'])

gROOT.SetBatch(True)
gROOT.ProcessLine("gErrorIgnoreLevel = kFatal;")
outfile = TFile(outFileName,'recreate')

varsName4Tuple = ':'.join(cutVars) + ':PtMin:PtMax:Signif:SoverB:EffPrompt' 
tSignif = TNtuple('tSignif','tSignif',varsName4Tuple)

totSets = 1
for varnum, iVar in enumerate(cutVars):
  totSets *= int((cutVars[iVar]['binmax']-cutVars[iVar]['binmin']+1)/steps[varnum])

print('Total number of sets per pT bin: %d' % totSets)

for iPt in range(0,len(PtMin)) :
  
  #check if low or high pt
  if PtMax[iPt] <= PtThreshold:
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
  binGenPtMin = sGenPrompt.GetAxis(0).FindBin(PtMin[iPt]*1.0001)
  binGenPtMax = sGenPrompt.GetAxis(0).FindBin(PtMax[iPt]*0.9999)
  hPtGenPrompt = sGenPrompt.GetAxis(0).SetRange(binGenPtMin,binGenPtMax)
  hPtGenPrompt = sGenPrompt.Projection(0)
  nGenPrompt = hPtGenPrompt.Integral()
  
  #acceptance from TOY MC
  binAccPtMin = hPtGenAcc.GetXaxis().FindBin(PtMin[iPt]*1.0001)
  binAccPtMax = hPtGenAcc.GetXaxis().FindBin(PtMax[iPt]*0.9999)
  Acc = hPtGenAcc.Integral(binAccPtMin, binAccPtMax) / hPtGenLimAcc.Integral(binAccPtMin, binAccPtMax)

  #FONLL
  binFONLL = FONLL['PtMin'].index(PtMin[iPt])
  sigmaFONLL = FONLL['Max'][binFONLL]

  #RAA
  Raa = sTAMUmin.Eval((PtMin[iPt]+PtMax[iPt])/2)

  #reco 
  binRecoPtMin = sMassPtCutVars.GetAxis(1).FindBin(PtMin[iPt]*1.0001)
  binRecoPtMax = sMassPtCutVars.GetAxis(1).FindBin(PtMax[iPt]*0.9999)
  sMassPtCutVars.GetAxis(1).SetRange(binRecoPtMin,binRecoPtMax)
  sMassPtCutVarsPrompt.GetAxis(1).SetRange(binRecoPtMin,binRecoPtMax)
  sMassPtCutVarsFD.GetAxis(1).SetRange(binRecoPtMin,binRecoPtMax)

  cutSetCount = 0 
  for iBins in itertools.product(*ranges) :

    cutSetCount += 1
    if cutSetCount % 100 == 0:
      print('tested cut set number %d' % cutSetCount)

    hMassData, index, array4Ntuple = ApplyCuts(sMassPtCutVars,iBins,axesnum,upperlowercuts,'hMassData')  
    hMassPrompt, index, array4Ntuple = ApplyCuts(sMassPtCutVarsPrompt,iBins,axesnum,upperlowercuts,'hMassPrompt')  
    hMassFD, index, array4Ntuple = ApplyCuts(sMassPtCutVarsFD,iBins,axesnum,upperlowercuts,'hMassFD')  
    
    hMassSgn = hMassPrompt.Clone('hMassSgn') 
    hMassSgn.Add(hMassFD)
    fMassSgn = TF1('fMassSgn','gaus',1.7,2.15)
    hMassSgn.Fit('fMassSgn','Q')
    mean = fMassSgn.GetParameter(1)
    sigma = fMassSgn.GetParameter(2)

    hMassSB = GetSideBandHisto(hMassData,mean,sigma)
    B = GetExpectedBackground(hMassSB,mean,sigma,Nexp,nEvBkg)

    nRecoPrompt = hMassPrompt.Integral()
    effPrompt = nRecoPrompt / nGenPrompt
    S = GetExpectedSignal(PtMin[iPt]-PtMax[iPt],sigmaFONLL,Raa,Taa,effPrompt,Acc,fprompt,BR,fractoD,Nexp)

    array4Ntuple.append(PtMin[iPt])
    array4Ntuple.append(PtMax[iPt])
    array4Ntuple.append(S/math.sqrt(S+B))
    array4Ntuple.append(S/B)
    array4Ntuple.append(effPrompt)
    tSignif.Fill(array.array("f",array4Ntuple))

tSignif.Write()
outfile.Close()
