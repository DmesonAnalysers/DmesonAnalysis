#********************************************************************************#
# python script for the significance scan of Ds+ mesons with THnSparses          #
# run: python ScanSignificanceSparse.py cfgFileName.yml outFileName.root         #
# author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino                 #
#********************************************************************************#

import sys
import itertools
import math
import array
import time
import yaml
from ReadModel import ReadFONLL, ReadTAMU
from ROOT import TFile, TF1, TNtuple, TSpline3, gROOT # pylint: disable=import-error,no-name-in-module

# pylint: disable=redefined-outer-name,invalid-name
def ApplyCuts(sparse, bins, axesnum, upperlowercuts, name):
    index = ''
    cutvalues = []
    for iVar, _ in enumerate(axesnum):
        if upperlowercuts[iVar] == 'Upper':
            sparse.GetAxis(axesnum[iVar]).SetRange(0, bins[iVar])
            cutvalues.append(sparse.GetAxis(axesnum[iVar]).GetBinLowEdge(bins[iVar])+sparse.GetAxis(axesnum[iVar]).GetBinWidth(bins[iVar]))
        elif upperlowercuts[iVar] == 'Lower':
            sparse.GetAxis(axesnum[iVar]).SetRange(bins[iVar], sparse.GetAxis(axesnum[iVar]).GetNbins()+1)
            cutvalues.append(sparse.GetAxis(axesnum[iVar]).GetBinLowEdge(bins[iVar]))
        index += '%d' % bins[iVar]
    hProj = sparse.Projection(0)
    hProj.SetName(name)
    for iAxis in axesnum:
        sparse.GetAxis(iAxis).SetRange(-1, -1)
    return hProj, index, cutvalues

def GetSideBandHisto(hMassData, mean, sigma):
    hMassDataReb = hMassData.Clone('hSB')
    hMassDataReb.Rebin(5)
    hSB = hMassDataReb.Clone('hSB')
    massbinmin = hSB.GetXaxis().FindBin(mean-4*sigma)
    massbinmax = hSB.GetXaxis().FindBin(mean+4*sigma)
    for iMassBin in range(1, hSB.GetNbinsX()):
        if massbinmin <= iMassBin <= massbinmax:
            hSB.SetBinContent(iMassBin, 0.)
            hSB.SetBinError(iMassBin, 0.)
        else:
            hSB.SetBinContent(iMassBin, hMassDataReb.GetBinContent(iMassBin))
            hSB.SetBinError(iMassBin, hMassDataReb.GetBinError(iMassBin))
    return hSB

def GetExpectedBackgroundFromSB(hSB, mean, sigma, Nexp, Nanal):
    fMassBkg = TF1('fMassBkg', 'expo', 1.8, 2.12)
    hSB.Fit('fMassBkg', 'QR')
    B = fMassBkg.Integral(mean-3*sigma, mean+3*sigma) / hSB.GetBinWidth(1) * Nexp / Nanal
    return B

def GetExpectedBackgroundFromMC(hBkgMC, mean, sigma, Nexp, Nanal):
    binmin = hBkgMC.GetXaxis().FindBin(mean-3*sigma)
    binmax = hBkgMC.GetXaxis().FindBin(mean+3*sigma)
    B = hBkgMC.Integral(binmin, binmax) / hBkgMC.GetBinWidth(1) * Nexp / Nanal
    return B

def GetExpectedSignal(DeltaPt, sigma, Raa, Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp):
    corrYield = sigma * fractoD * Raa * Taa / DeltaPt
    rawYield = 2 * DeltaPt * corrYield * effPrompt * Acc * BR * Nexp / fprompt
    return rawYield

def ComputeExpectedFprompt(ptmin, ptmax, effPrompt, hPredPrompt, effFD, hPredFD, RatioRaaFDPrompt):
    ptbinmin = hPredPrompt.GetXaxis().FindBin(ptmin*1.001)
    ptbinmax = hPredFD.GetXaxis().FindBin(ptmax*0.999)
    predPrompt, predFD = (0 for iPred in range(2))
    for iPt in range(ptbinmin, ptbinmax+1):
        predPrompt += hPredPrompt.GetBinContent(iPt) * hPredPrompt.GetBinWidth(iPt)
        predFD += hPredFD.GetBinContent(iPt) * hPredFD.GetBinWidth(iPt)
    predPrompt /= (ptmax-ptmin)
    predFD /= (ptmax-ptmin)
    PromptSignal = predPrompt * effPrompt
    FDSignal = predFD * effFD
    fprompt = 1. / (1 + RatioRaaFDPrompt * FDSignal / PromptSignal)
    return fprompt

#***************************************************************************
# Main function

cfgFileName = sys.argv[1]
outFileName = sys.argv[2]

with open(cfgFileName, 'r') as ymlCfgFile:
    if six.PY2:
        inputCfg = yaml.load(ymlCfgFile)
    else:
        inputCfg = yaml.safe_load(ymlCfgFile)

if inputCfg['getbkgfromMC'] is False:
    infileDataLowPt = TFile.Open(inputCfg['fileDataLowPt']['filename'])
    indirDataLowPt = infileDataLowPt.Get(inputCfg['fileDataLowPt']['dirname'])
    inlistDataLowPt = indirDataLowPt.Get(inputCfg['fileDataLowPt']['listname'])
    sMassPtCutVarsLowPt = inlistDataLowPt.FindObject(inputCfg['fileDataLowPt']['sparsenameAll'])
    hEvLowPt = inlistDataLowPt.FindObject(inputCfg['fileDataLowPt']['histoevname'])

    infileDataHighPt = TFile.Open(inputCfg['fileDataHighPt']['filename'])
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
if inputCfg['getbkgfromMC']:
    sMassPtCutVarsBkgLowPt = inlistMCLowPt.FindObject(inputCfg['fileMCLowPt']['sparsenameBkg'])
    hEvLowPt = inlistMCLowPt.FindObject(inputCfg['fileMCLowPt']['histoevname'])

infileMCHighPt = TFile(inputCfg['fileMCHighPt']['filename'])
indirMCHighPt = infileMCHighPt.Get(inputCfg['fileMCHighPt']['dirname'])
inlistMCHighPt = indirMCHighPt.Get(inputCfg['fileMCHighPt']['listname'])
sMassPtCutVarsPromptHighPt = inlistMCHighPt.FindObject(inputCfg['fileMCHighPt']['sparsenamePrompt'])
sMassPtCutVarsFDHighPt = inlistMCHighPt.FindObject(inputCfg['fileMCHighPt']['sparsenameFD'])
sGenPromptHighPt = inlistMCHighPt.FindObject(inputCfg['fileMCHighPt']['sparsenameGenPrompt'])
sGenFDHighPt = inlistMCLowPt.FindObject(inputCfg['fileMCHighPt']['sparsenameGenFD'])
if inputCfg['getbkgfromMC']:
    sMassPtCutVarsBkgHighPt = inlistMCLowPt.FindObject(inputCfg['fileMCHighPt']['sparsenameBkg'])
    hEvHighPt = inlistMCLowPt.FindObject(inputCfg['fileMCHighPt']['histoevname'])

PtThreshold = inputCfg['PtThreshold']

Nexp = inputCfg['nExpectedEvents']
Taa = inputCfg['Taa']
BR = inputCfg['BR']
fractoD = inputCfg['fractoD']

if inputCfg['PredForFprompt']['estimateFprompt']:
    infilePred = TFile.Open(inputCfg['PredForFprompt']['filename'])
    hPredPrompt = infilePred.Get(inputCfg['PredForFprompt']['histonamePrompt'])
    hPredFD = infilePred.Get(inputCfg['PredForFprompt']['histonameFD'])
    RatioRaaFDPrompt = inputCfg['PredForFprompt']['RatioRaaFDPrompt']

FONLL = ReadFONLL(inputCfg['filenameFONLL'])
TAMU = ReadTAMU(inputCfg['filenameRaaPredTAMU']) #optimisation performed using TAMU
sTAMUmin = TSpline3('sTAMUmin', array.array('d', TAMU['PtCent']), array.array('d', TAMU['Min']), len(TAMU['PtCent']))
sTAMUmax = TSpline3('sTAMUmax', array.array('d', TAMU['PtCent']), array.array('d', TAMU['Max']), len(TAMU['PtCent']))

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
    if cutVars[iVar]['upperlowercut'] == 'Lower':
        cutVars[iVar]['binmin'] = sMassPtCutVarsPromptLowPt.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min']*1.0001)
        cutVars[iVar]['binmax'] = sMassPtCutVarsPromptLowPt.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max']*1.0001)
    elif cutVars[iVar]['upperlowercut'] == 'Upper':
        cutVars[iVar]['binmin'] = sMassPtCutVarsPromptLowPt.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min']*0.9999)
        cutVars[iVar]['binmax'] = sMassPtCutVarsPromptLowPt.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max']*0.9999)
    steps.append(int(cutVars[iVar]['step'] / sMassPtCutVarsPromptLowPt.GetAxis(cutVars[iVar]['axisnum']).GetBinWidth(1)))
    ranges.append(range(cutVars[iVar]['binmin'], cutVars[iVar]['binmax']+1, steps[varnum]))
    axesnum.append(cutVars[iVar]['axisnum'])
    upperlowercuts.append(cutVars[iVar]['upperlowercut'])

gROOT.SetBatch(True)
gROOT.ProcessLine("gErrorIgnoreLevel = kFatal;")
outfile = TFile(outFileName, 'recreate')

varsName4Tuple = ':'.join(cutVars) + ':PtMin:PtMax:S:B:Signif:SoverB:EffPrompt:EffFD:fprompt'
tSignif = TNtuple('tSignif', 'tSignif', varsName4Tuple)

totSets = 1
for varnum, iVar in enumerate(cutVars):
    totSets *= int((cutVars[iVar]['binmax']-cutVars[iVar]['binmin'])/steps[varnum])+1

print('Total number of sets per pT bin: %d' % totSets)

for iPt, _ in enumerate(PtMin):
    #check if low or high pt
    if PtMax[iPt] <= PtThreshold:
        nEvBkg = hEvLowPt.GetBinContent(5)
        sMassPtCutVarsPrompt = sMassPtCutVarsPromptLowPt.Clone('sMassPtCutVarsPrompt')
        sMassPtCutVarsFD = sMassPtCutVarsFDLowPt.Clone('sMassPtCutVarsFD')
        sGenPrompt = sGenPromptLowPt.Clone('sGenPrompt')
        sGenFD = sGenFDLowPt.Clone('sGenFD')
        if inputCfg['getbkgfromMC']:
            sMassPtCutVarsBkg = sMassPtCutVarsBkgLowPt.Clone('sMassPtCutVarsBkg')
        else:
            sMassPtCutVars = sMassPtCutVarsLowPt.Clone('sMassPtCutVars')
    else:
        nEvBkg = hEvHighPt.GetBinContent(5)
        sMassPtCutVarsPrompt = sMassPtCutVarsPromptHighPt.Clone('sMassPtCutVarsPrompt')
        sMassPtCutVarsFD = sMassPtCutVarsFDHighPt.Clone('sMassPtCutVarsFD')
        sGenPrompt = sGenPromptHighPt.Clone('sGenPrompt')
        sGenFD = sGenFDHighPt.Clone('sGenFD')
        if inputCfg['getbkgfromMC']:
            sMassPtCutVarsBkg = sMassPtCutVarsBkgHighPt.Clone('sMassPtCutVarsBkg')
        else:
            sMassPtCutVars = sMassPtCutVarsHighPt.Clone('sMassPtCutVars')

    #gen for efficiency
    binGenPtMin = sGenPrompt.GetAxis(0).FindBin(PtMin[iPt]*1.0001)
    binGenPtMax = sGenPrompt.GetAxis(0).FindBin(PtMax[iPt]*0.9999)
    hPtGenPrompt = sGenPrompt.GetAxis(0).SetRange(binGenPtMin, binGenPtMax)
    hPtGenPrompt = sGenPrompt.Projection(0)
    hPtGenFD = sGenFD.GetAxis(0).SetRange(binGenPtMin, binGenPtMax)
    hPtGenFD = sGenFD.Projection(0)
    nGenPrompt = hPtGenPrompt.Integral()
    nGenFD = hPtGenPrompt.Integral()

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
    binRecoPtMin = sMassPtCutVarsPrompt.GetAxis(1).FindBin(PtMin[iPt]*1.0001)
    binRecoPtMax = sMassPtCutVarsPrompt.GetAxis(1).FindBin(PtMax[iPt]*0.9999)
    sMassPtCutVarsPrompt.GetAxis(1).SetRange(binRecoPtMin, binRecoPtMax)
    sMassPtCutVarsFD.GetAxis(1).SetRange(binRecoPtMin, binRecoPtMax)
    if inputCfg['getbkgfromMC'] is False:
        sMassPtCutVars.GetAxis(1).SetRange(binRecoPtMin, binRecoPtMax)

    cutSetCount = 0
    start_time = time.time()
    for iBins in itertools.product(*ranges):
        cutSetCount += 1
        if cutSetCount % 100 == 0:
            elapsed_time = time.time() - start_time
            print('tested cut set number %d, elapsed time: %f s' % (cutSetCount, elapsed_time))

        hMassPrompt, index, array4Ntuple = ApplyCuts(sMassPtCutVarsPrompt, iBins, axesnum, upperlowercuts, 'hMassPrompt')
        hMassFD, index, array4Ntuple = ApplyCuts(sMassPtCutVarsFD, iBins, axesnum, upperlowercuts, 'hMassFD')

        if inputCfg['getbkgfromMC']:
            hMassBkg, index, array4Ntuple = ApplyCuts(sMassPtCutVarsBkg, iBins, axesnum, upperlowercuts, 'hMassBkg')
        else:
            hMassData, index, array4Ntuple = ApplyCuts(sMassPtCutVars, iBins, axesnum, upperlowercuts, 'hMassData')

        hMassSgn = hMassPrompt.Clone('hMassSgn')
        hMassSgn.Add(hMassFD)
        fMassSgn = TF1('fMassSgn', 'gaus', 1.7, 2.15)
        hMassSgn.Fit('fMassSgn', 'Q')
        mean = fMassSgn.GetParameter(1)
        sigma = fMassSgn.GetParameter(2)

        nRecoPrompt = hMassPrompt.Integral()
        nRecoFD = hMassFD.Integral()
        effPrompt = nRecoPrompt / nGenPrompt
        effFD = nRecoFD / nGenFD

        if inputCfg['getbkgfromMC']:
            B = GetExpectedBackgroundFromMC(hMassBkg, mean, sigma, Nexp, nEvBkg)        
        else:
            hMassSB = GetSideBandHisto(hMassData, mean, sigma)
            B = GetExpectedBackgroundFromSB(hMassSB, mean, sigma, Nexp, nEvBkg)

        if inputCfg['PredForFprompt']['estimateFprompt']:
            fprompt = ComputeExpectedFprompt(PtMin[iPt], PtMax[iPt], effPrompt, hPredPrompt, effFD, hPredFD, RatioRaaFDPrompt)
        else:
            fprompt = inputCfg['fprompt']

        S = GetExpectedSignal(PtMin[iPt]-PtMax[iPt], sigmaFONLL, Raa, Taa, effPrompt, Acc, fprompt, BR, fractoD, Nexp)

        array4Ntuple.append(PtMin[iPt])
        array4Ntuple.append(PtMax[iPt])
        array4Ntuple.append(S)
        array4Ntuple.append(B)
        array4Ntuple.append(S/math.sqrt(S+B))
        array4Ntuple.append(S/B)
        array4Ntuple.append(effPrompt)
        array4Ntuple.append(effFD)
        array4Ntuple.append(fprompt)
        tSignif.Fill(array.array("f", array4Ntuple))

tSignif.Write()
outfile.Close()
