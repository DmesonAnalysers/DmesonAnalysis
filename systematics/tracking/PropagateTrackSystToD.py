'''
python script for the propagation of the single-track uncertainty to D+ and Ds+ mesons using TTrees
run: python PropagateTrackSystToD.py cfgFileName.yml cutSetFileName.yml outFileName.root period [--Ds] [--Dplus]
'''

import sys
import argparse
import yaml
import pandas as pd
import numpy as np
from ROOT import TFile, TCanvas, TH2F, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kAzure, kRainBow, kFullCircle, kFullSquare, kFullDiamond # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.DfUtils import FilterBitDf, LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ApplyHistoEntriesToColumn #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

SetGlobalStyle(palette=kRainBow, padleftmargin=0.14, padrightmargin=0.14, padbottommargin=0.14, titleoffsety=1.4)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml', help='input file with cut set')
parser.add_argument('outFileName', metavar='text', default='outFileName.root', help='output root file name')
parser.add_argument('period', metavar='text', default='LHC17pq', help='data period for systematic evaluation')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('--Dplus', action='store_true', default=False, help='enable comparison for D+')
group.add_argument('--Ds', action='store_true', default=False, help='enable comparison for Ds')
group.add_argument('--Lc2pK0s', action='store_true', default=False, help='enable comparison for Lc2pK0s')
args = parser.parse_args()

nDau = 3
if args.Dplus:
    particle = 'Dplus'
elif args.Ds:
    particle = 'Ds'
else:
    particle = 'Lc2pK0s'
    nDau = 1

# input systematic uncertainties
inFileSystUnc = TFile.Open(f'singletracksyst/traking_ME_piK_syst_{particle}_{args.period}.root')
hTrkEff = inFileSystUnc.Get('hTrEff')
hME = inFileSystUnc.Get('h')
hTrkEff.SetDirectory(0)
hME.SetDirectory(0)
inFileSystUnc.Close()

# config with input file details
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
inFileNames = inputCfg['filename']
if not isinstance(inFileNames, list):
    inFileNames = [inFileNames]
isMC = inputCfg['isMC']
if not isMC:
    if args.ptweights:
        print('ERROR: systematic uncertainty cannot be computed since it is not MC! Exit')
        sys.exit()

# selections to be applied
with open(args.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']
selToApply = []
for iPt, _ in enumerate(cutVars['Pt']['min']):
    selToApply.append('')
    for varName in cutVars:
        if varName == 'InvMass':
            continue
        if selToApply[iPt] != '':
            selToApply[iPt] += ' & '
        selToApply[iPt] += f"{cutVars[varName]['min'][iPt]}<{cutVars[varName]['name']}<{cutVars[varName]['max'][iPt]}"

# define filter bits
bitSignal = 0
bitPrompt = 2
bitFD = 3
bitRefl = 4

dataFramePrompt = LoadDfFromRootOrParquet(inputCfg['tree']['filenamePrompt'], inputCfg['tree']['dirname'],
                                          inputCfg['tree']['treename'])
        
if 'cand_type' in dataFramePrompt.columns: #if not filtered tree, select only prompt and not reflected
    dataFramePrompt = FilterBitDf(dataFramePrompt, 'cand_type', [bitSignal, bitPrompt], 'and')
    dataFramePrompt = FilterBitDf(dataFramePrompt, 'cand_type', [bitRefl], 'not')
dataFramePrompt.reset_index(inplace=True)

dataFrameFD = LoadDfFromRootOrParquet(inputCfg['tree']['filenameFD'], inputCfg['tree']['dirname'],
                                      inputCfg['tree']['treename'])
if 'cand_type' in dataFrameFD.columns: #if not filtered tree, select only FD and not reflected
    dataFrameFD = FilterBitDf(dataFrameFD, 'cand_type', [bitSignal, bitFD], 'and')
    dataFrameFD = FilterBitDf(dataFrameFD, 'cand_type', [bitRefl], 'not')
dataFrameFD.reset_index(inplace=True)

if 'pt_prong0' not in dataFramePrompt.columns or 'pt_prong0' not in dataFrameFD.columns:
    print('ERROR: input dataframe does not contain daughter-track pt, propagation not possible! Exit')
    sys.exit()

for iPt, (cuts, ptMin, ptMax) in enumerate(zip(selToApply, cutVars['Pt']['min'], cutVars['Pt']['max'])):
    if iPt == 0:
        dataFramePromptSel = dataFramePrompt.astype(float).query(cuts)
        dataFrameFDSel = dataFrameFD.astype(float).query(cuts)
    else:
        dataFramePromptSel = dataFramePromptSel.append(dataFramePrompt.astype(float).query(cuts))
        dataFrameFDSel = dataFrameFDSel.append(dataFrameFD.astype(float).query(cuts))

dataFramePromptSel.reset_index(inplace=True)
dataFrameFDSel.reset_index(inplace=True)

for iDau in range(nDau):
    dataFramePromptSel[f'ME_dau{iDau}'] = ApplyHistoEntriesToColumn(dataFramePromptSel, f'pt_prong{iDau}', hME)
    dataFrameFDSel[f'ME_dau{iDau}'] = ApplyHistoEntriesToColumn(dataFrameFDSel, f'pt_prong{iDau}', hME)

if args.Lc2pK0s: # Use only the track of the proton
    dataFramePromptSel['ME_tot'] = dataFramePromptSel.apply(
        lambda row: row['ME_dau0'], axis=1)
    dataFrameFDSel['ME_tot'] = dataFrameFDSel.apply(
        lambda row: row['ME_dau0'], axis=1)
else:
    dataFramePromptSel['ME_tot'] = dataFramePromptSel.apply(
        lambda row: row['ME_dau0']+row['ME_dau1']+row['ME_dau2'], axis=1)
    dataFrameFDSel['ME_tot'] = dataFrameFDSel.apply(
        lambda row: row['ME_dau0']+row['ME_dau1']+row['ME_dau2'], axis=1)

# dataFramePromptSel=dataFramePromptSel[dataFramePromptSel['pt_prong0']==dataFramePromptSel['pt_prong0']]
# dataFrameFDSel=dataFrameFDSel[dataFrameFDSel['pt_prong0']==dataFrameFDSel['pt_prong0']]

nPtBins = hTrkEff.GetNbinsX()
ptLims = cutVars['Pt']['min'].copy()
ptLims.append(cutVars['Pt']['max'][-1])

hSystVsPtPrompt = TH2F('hSystVsPtPrompt', ';#it{p}_{T} (GeV/#it{c});ME relative systematic uncertainty',
                       nPtBins, np.array(ptLims, 'd'), 30, 0., 0.15)
hSystVsPtFD = TH2F('hSystVsPtFD', ';#it{p}_{T} (GeV/#it{c});ME relative systematic uncertainty',
                   nPtBins, np.array(ptLims, 'd'), 30, 0., 0.15)
hPtDauVsPtDPrompt = TH2F('hPtDauVsPtDPrompt', ';#it{p}_{T}^{ D} (GeV/#it{c});#it{p}_{T}^{ dau} (GeV/#it{c})',
                         int(ptLims[-1]-ptLims[0])*10, ptLims[0], ptLims[-1],
                         int(ptLims[-1]-ptLims[0])*10, ptLims[0], ptLims[-1])
hPtDauVsPtDFD = TH2F('hPtDauVsPtDFD', ';#it{p}_{T}^{ D} (GeV/#it{c});#it{p}_{T}^{ dau} (GeV/#it{c})',
                     int(ptLims[-1]-ptLims[0])*10, ptLims[0], ptLims[-1],
                     int(ptLims[-1]-ptLims[0])*10, ptLims[0], ptLims[-1])

for pt, meTot in zip(dataFramePromptSel['pt_cand'].to_numpy(), dataFramePromptSel['ME_tot'].to_numpy()):
    hSystVsPtPrompt.Fill(pt, meTot)
for pt, meTot in zip(dataFrameFDSel['pt_cand'].to_numpy(), dataFrameFDSel['ME_tot'].to_numpy()):
    hSystVsPtFD.Fill(pt, meTot)
hSystVsPtAll = hSystVsPtPrompt.Clone('hSystVsPtAll')
hSystVsPtAll.Add(hSystVsPtFD)

hSystMeanPrompt = hSystVsPtPrompt.ProfileX()
hSystMeanFD = hSystVsPtFD.ProfileX()
hSystMeanAll = hSystVsPtAll.ProfileX()
SetObjectStyle(hSystMeanPrompt, fillstyle=0, markersize=0.8)
SetObjectStyle(hSystMeanFD, fillstyle=0, markersize=0.8)
SetObjectStyle(hSystMeanAll, fillstyle=0, markersize=0.8)

for iDau in range(nDau):
    hTmp = hPtDauVsPtDPrompt.Clone('hTmp')
    hTmp.Reset()
    for pt, ptProng in zip(dataFramePromptSel['pt_cand'].to_numpy(), dataFramePromptSel[f'pt_prong{iDau}'].to_numpy()):
        hTmp.Fill(pt, ptProng)
    hPtDauVsPtDPrompt.Add(hTmp)
    hTmp = hPtDauVsPtDFD.Clone('hTmp')
    hTmp.Reset()
    for pt, ptProng in zip(dataFrameFDSel['pt_cand'].to_numpy(), dataFrameFDSel[f'pt_prong{iDau}'].to_numpy()):
        hTmp.Fill(pt, ptProng)
    hPtDauVsPtDFD.Add(hTmp)
hPtDauVsPtDAll = hPtDauVsPtDPrompt.Clone('hPtDauVsPtDAll')
hPtDauVsPtDAll.Add(hPtDauVsPtDFD)

hPtDauMeanPrompt = hPtDauVsPtDPrompt.ProfileX()
hPtDauMeanFD = hPtDauVsPtDFD.ProfileX()
hPtDauMeanAll = hPtDauVsPtDAll.ProfileX()
SetObjectStyle(hPtDauMeanPrompt, fillstyle=0, markersize=0.5)
SetObjectStyle(hPtDauMeanFD, fillstyle=0, markersize=0.5)
SetObjectStyle(hPtDauMeanAll, fillstyle=0, markersize=0.5)

hTotSystPrompt = hTrkEff.Clone('hTotSystPrompt')
hTotSystFD = hTrkEff.Clone('hTotSystFD')
hTotSystAll = hTrkEff.Clone('hTotSystAll')
hTotSystPrompt.GetYaxis().SetTitle('Relative systematic uncertainty')
hTotSystFD.GetYaxis().SetTitle('Relative systematic uncertainty')
hTotSystAll.GetYaxis().SetTitle('Relative systematic uncertainty')
SetObjectStyle(hTotSystPrompt, color=kRed+1, fillstyle=0, markerstyle=kFullCircle)
SetObjectStyle(hTotSystFD, color=kAzure+4, fillstyle=0, markerstyle=kFullSquare)
SetObjectStyle(hTotSystAll, fillstyle=0, markerstyle=kFullDiamond)
for iPt in range(1, hTotSystPrompt.GetNbinsX()+1):
    hTotSystPrompt.SetBinContent(iPt, np.sqrt(hTrkEff.GetBinContent(iPt)**2 + hSystMeanPrompt.GetBinContent(iPt)**2))
    hTotSystFD.SetBinContent(iPt, np.sqrt(hTrkEff.GetBinContent(iPt)**2 + hSystMeanFD.GetBinContent(iPt)**2))
    hTotSystAll.SetBinContent(iPt, np.sqrt(hTrkEff.GetBinContent(iPt)**2 + hSystMeanAll.GetBinContent(iPt)**2))
    hTotSystPrompt.SetBinError(iPt, 1.e-20)
    hTotSystFD.SetBinError(iPt, 1.e-20)
    hTotSystAll.SetBinError(iPt, 1.e-20)

legAverage = TLegend(0.2, 0.2, 0.4, 0.3)
legAverage.SetTextSize(0.045)
legAverage.SetFillStyle(0)
legAverage.AddEntry(hPtDauMeanPrompt, 'average', 'pl')

leg = TLegend(0.2, 0.2, 0.4, 0.5)
leg.SetTextSize(0.045)
leg.SetFillStyle(0)
leg.AddEntry(hTotSystPrompt, 'prompt', 'pl')
leg.AddEntry(hTotSystFD, 'FD', 'pl')
leg.AddEntry(hTotSystAll, 'prompt + FD', 'pl')

cPrompt = TCanvas('cPrompt', 'Prompt', 1000, 500)
cPrompt.Divide(2, 1)
cPrompt.cd(1).SetLogz()
hSystVsPtPrompt.Draw('colz')
hSystMeanPrompt.DrawCopy('PH][ E0 same')
legAverage.Draw()
cPrompt.cd(2).SetLogz()
hPtDauVsPtDPrompt.Draw('colz')
hPtDauMeanPrompt.DrawCopy('PH][ E0 same')
cPrompt.Modified()
cPrompt.Update()

cFD = TCanvas('cFD', 'FD', 1000, 500)
cFD.Divide(2, 1)
cFD.cd(1).SetLogz()
hSystVsPtFD.Draw('colz')
hSystMeanFD.DrawCopy('PH][ E0 same')
legAverage.Draw()
cFD.cd(2).SetLogz()
hPtDauVsPtDFD.Draw('colz')
hPtDauMeanFD.DrawCopy('PH][ E0 same')
cFD.Modified()
cFD.Update()

cAll = TCanvas('cAll', 'Prompt+FD', 1000, 500)
cAll.Divide(2, 1)
cAll.cd(1).SetLogz()
hSystVsPtAll.Draw('colz')
hSystMeanAll.DrawCopy('PH][ E0 same')
legAverage.Draw()
cAll.cd(2).SetLogz()
hPtDauVsPtDAll.Draw('colz')
hPtDauMeanAll.DrawCopy('PH][ E0 same')
cAll.Modified()
cAll.Update()

cFinalSyst = TCanvas('cFinalSyst', '', 500, 500)
cFinalSyst.DrawFrame(ptLims[0], 0., ptLims[-1], 0.15, ';#it{p}_{T} (GeV/#it{c});Relative systematic uncertainty')
hTotSystFD.DrawCopy('PH][ E0 same')
hTotSystPrompt.DrawCopy('PH][ E0 same')
hTotSystAll.DrawCopy('PH][ E0 same')
leg.Draw()
cFinalSyst.Modified()
cFinalSyst.Update()

outFile = TFile(args.outFileName, 'recreate')
cPrompt.Write()
cFD.Write()
cAll.Write()
cFinalSyst.Write()
hSystVsPtPrompt.Write()
hSystMeanPrompt.Write()
hPtDauVsPtDPrompt.Write()
hPtDauMeanPrompt.Write()
hSystVsPtFD.Write()
hSystMeanFD.Write()
hPtDauVsPtDFD.Write()
hPtDauMeanFD.Write()
hSystVsPtAll.Write()
hSystMeanAll.Write()
hPtDauVsPtDAll.Write()
hSystMeanAll.Write()
hTotSystPrompt.Write()
hTotSystFD.Write()
hTotSystAll.Write()
outFile.Close()

outFileNamePDF = args.outFileName
outFileNamePDF = outFileNamePDF.replace('.root', '')

cPrompt.SaveAs(f'{outFileNamePDF}_Prompt.pdf')
cFD.SaveAs(f'{outFileNamePDF}_FD.pdf')
cAll.SaveAs(f'{outFileNamePDF}_All.pdf')
cFinalSyst.SaveAs(f'{outFileNamePDF}_FinalSyst.pdf')

print('\nTotal tracking systematic uncertainty:')
print('\nPrompt')
for iPt in range(hTotSystPrompt.GetNbinsX()):
    print(
        f'\t\t pT (GeV/c): [{ptLims[iPt]:.1f} - {ptLims[iPt+1]:.1f}] syst: {hTotSystPrompt.GetBinContent(iPt+1):.3f}')
print('\nFD')
for iPt in range(hTotSystFD.GetNbinsX()):
    print(
        f'\t\t pT (GeV/c): [{ptLims[iPt]:.1f} - {ptLims[iPt+1]:.1f}] syst: {hTotSystFD.GetBinContent(iPt+1):.3f}')
print('\nAll')
for iPt in range(hTotSystAll.GetNbinsX()):
    print(
        f'\t\t pT (GeV/c): [{ptLims[iPt]:.1f} - {ptLims[iPt+1]:.1f}] syst: {hTotSystAll.GetBinContent(iPt+1):.3f}')

# print(hPtDauVsPtDPrompt.GetEntries())

input('\nPress enter to exit')
