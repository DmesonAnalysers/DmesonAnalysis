'''
Script for the combination of efficiency and acceptance factors
run: python ComputeEffAccWeightedAvg.py effFileAll.root effFileNonRes.root effFileKStar.root effFileDelta.root effFileLambda1520.root outFileName.root
'''

import argparse
import numpy as np
from ROOT import gROOT, TFile, TCanvas, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kFullDiamond, kRed, kAzure, kFullCircle, kOpenSquare  # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('effFileAll', metavar='text', default='')
parser.add_argument('effFileNonRes', metavar='text', default='')
parser.add_argument('effFileKStar', metavar='text', default='')
parser.add_argument('effFileDelta', metavar='text', default='')
parser.add_argument('effFileLambda1520', metavar='text', default='')
parser.add_argument('outFileName', metavar='text', default='')
parser.add_argument("--batch", help="suppress video output", action="store_true")
args = parser.parse_args()

gROOT.SetBatch(args.batch)
SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, titlesize=0.045, labelsize=0.04)

effFileAll = TFile.Open(args.effFileAll)
hEffAllPrompt = effFileAll.Get('hAccEffPrompt')
hEffAllFD = effFileAll.Get('hAccEffFD')
SetObjectStyle(hEffAllPrompt, color=kRed+1, markerstyle=kFullCircle)
SetObjectStyle(hEffAllFD, color=kAzure+4, markerstyle=kOpenSquare, markersize=1.5, linewidh=2, linestyle=7)

effFileNonRes = TFile.Open(args.effFileNonRes)
hEffNonResPrompt = effFileNonRes.Get('hAccEffPrompt')
hEffNonResFD = effFileNonRes.Get('hAccEffFD')
SetObjectStyle(hEffNonResPrompt, color=kRed+1, markerstyle=kFullCircle)
SetObjectStyle(hEffNonResFD, color=kAzure+4, markerstyle=kOpenSquare, markersize=1.5, linewidh=2, linestyle=7)

effFileKStar = TFile.Open(args.effFileKStar)
hEffKStarPrompt = effFileKStar.Get('hAccEffPrompt')
hEffKStarFD = effFileKStar.Get('hAccEffFD')
SetObjectStyle(hEffKStarPrompt, color=kRed+1, markerstyle=kFullCircle)
SetObjectStyle(hEffKStarFD, color=kAzure+4, markerstyle=kOpenSquare, markersize=1.5, linewidh=2, linestyle=7)

effFileDelta = TFile.Open(args.effFileDelta)
hEffDeltaPrompt = effFileDelta.Get('hAccEffPrompt')
hEffDeltaFD = effFileDelta.Get('hAccEffFD')
SetObjectStyle(hEffDeltaPrompt, color=kRed+1, markerstyle=kFullCircle)
SetObjectStyle(hEffDeltaFD, color=kAzure+4, markerstyle=kOpenSquare, markersize=1.5, linewidh=2, linestyle=7)

effFileLambda1520 = TFile.Open(args.effFileLambda1520)
hEffLambda1520Prompt = effFileLambda1520.Get('hAccEffPrompt')
hEffLambda1520FD = effFileLambda1520.Get('hAccEffFD')
SetObjectStyle(hEffLambda1520Prompt, color=kRed+1, markerstyle=kFullCircle)
SetObjectStyle(hEffLambda1520FD, color=kAzure+4, markerstyle=kOpenSquare, markersize=1.5, linewidh=2, linestyle=7)


hEffCw = hEffAllPrompt.Clone("hAccEffPrompt")
hEffBw = hEffAllFD.Clone("hAccEffFD")

BR = [3.5 * 1e-02, 1.96 * 0.667 * 1e-02, 1.08 * 1e-02, 2.2 * 0.225 * 1e-02]


nPtBins = hEffNonResPrompt.GetNbinsX()
for iPt in range(nPtBins):
    hEffCw.SetBinContent(iPt+1, ((hEffNonResPrompt.GetBinContent(iPt+1)*BR[0])+(hEffKStarPrompt.GetBinContent(iPt+1)*BR[1])+(hEffDeltaPrompt.GetBinContent(iPt+1)*BR[2])+(hEffLambda1520Prompt.GetBinContent(iPt+1)*BR[3]))/(BR[0]+BR[1]+BR[2]+BR[3]))
    hEffCw.SetBinError(iPt+1, hEffAllPrompt.GetBinError(iPt+1))
    hEffBw.SetBinContent(iPt+1, ((hEffNonResFD.GetBinContent(iPt+1)*BR[0])+(hEffKStarFD.GetBinContent(iPt+1)*BR[1])+(hEffDeltaFD.GetBinContent(iPt+1)*BR[2])+(hEffLambda1520FD.GetBinContent(iPt+1)*BR[3]))/(BR[0]+BR[1]+BR[2]+BR[3]))
    hEffBw.SetBinError(iPt+1, hEffAllFD.GetBinError(iPt+1))

outFile = TFile(args.outFileName, 'recreate')
hEffCw.Write()
hEffBw.Write()
outFile.Close()