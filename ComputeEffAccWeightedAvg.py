'''
Script for the combination of efficiency and acceptance factors
run: python ComputeEffAccWeightedAvg.py effFileAll.root effFileNonRes.root effFileKStar.root effFileDelta.root effFileLambda1520.root outFileName.root
'''

import sys
import argparse
import ctypes
import numpy as np
from ROOT import gROOT, TFile, TCanvas, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure  # pylint: disable=import-error,no-name-in-module
from ROOT import kFullDiamond, kOpenDiamond, kFullTriangleUp, kOpenTriangleUp, kFullStar, kOpenStar, kFullCircle, kOpenCircle, kFullSquare, kOpenSquare, kFullCross, kOpenCross # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('effFileNonRes', metavar='text', default='')
parser.add_argument('effFileKStar', metavar='text', default='')
parser.add_argument('effFileDelta', metavar='text', default='')
parser.add_argument('effFileLambda1520', metavar='text', default='')
parser.add_argument('outFileName', metavar='text', default='')
parser.add_argument("--batch", help="suppress video output", action="store_true")
args = parser.parse_args()

gROOT.SetBatch(args.batch)
SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, titlesize=0.045, labelsize=0.04)

inputFile = [args.effFileNonRes, args.effFileKStar, args.effFileDelta, args.effFileLambda1520]
colorMarkerPrompt = [kBlack, kRed, kRed, kRed, kRed]
colorMarkerFD = [kBlack, kAzure, kAzure, kAzure, kAzure]
styleMarkerPrompt = [kFullDiamond, kFullSquare, kFullStar, kFullTriangleUp, kFullCross]
styleMarkerFD = [kOpenDiamond, kOpenSquare, kOpenStar, kOpenTriangleUp, kOpenCross]
hEffPrompt, hEffFD = [], []

for iReso, fileName, in enumerate(inputFile):
    infile = TFile.Open(fileName)
    hEffPrompt.append(infile.Get('hAccEffPrompt'))
    hEffFD.append(infile.Get('hAccEffFD'))
    hEffPrompt[iReso].SetDirectory(0)
    hEffFD[iReso].SetDirectory(0)
    SetObjectStyle(
        hEffPrompt[iReso],
        color=colorMarkerPrompt[iReso],
        markerstyle=styleMarkerPrompt[iReso]
    )
    SetObjectStyle(
        hEffFD[iReso],
        color=colorMarkerFD[iReso],
        markerstyle=styleMarkerFD[iReso],
        markersize=1.5, linewidh=2, linestyle=7
    )
    
ptMin = hEffPrompt[0].GetBinLowEdge(1)
ptMax = hEffPrompt[0].GetBinLowEdge(hEffPrompt[0].GetNbinsX()) + hEffPrompt[0].GetBinWidth(hEffPrompt[0].GetNbinsX())

cEffPrompt = TCanvas('cEffPrompt', '', 800, 800)
cEffPrompt.DrawFrame(ptMin, 1.e-5, ptMax, 1.,
               ';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon);')
cEffPrompt.SetLogy()
for histo in hEffPrompt:
    histo.Draw('same')

cEffFD = TCanvas('cEffFD', '', 800, 800)
cEffFD.DrawFrame(ptMin, 1.e-5, ptMax, 1.,
               ';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon);')
cEffFD.SetLogy()
for histo in hEffFD:
    histo.Draw('same')

outFileNamePDF = args.outFileName.replace('.root', '')
cEffPrompt.SaveAs('%s_prompt.pdf' % (outFileNamePDF))
cEffFD.SaveAs('%s_FD.pdf' % (outFileNamePDF))

hEffCw = hEffPrompt[0].Clone("hAccEffPrompt")
hEffBw = hEffFD[0].Clone("hAccEffFD")

BR = [6.28 * 1e-02, 3.5 * 1e-02, 1.96 * 0.667 * 1e-02, 1.08 * 1e-02, 2.2 * 0.225 * 1e-02]
effC, wC = 0,0
nPtBins = hEffPrompt[0].GetNbinsX()
for iPt in range(nPtBins):
    effC, effB, uncEffC, uncEffB, sumOfW = (0. for _ in range(5))
    for histo, br, in zip(hEffPrompt, BR):
        effC += histo.GetBinContent(iPt+1) * br
        effB += histo.GetBinContent(iPt+1) * br
        uncEffC += histo.GetBinError(iPt+1)**2 * br**2
        uncEffB += histo.GetBinError(iPt+1)**2 * br**2
        sumOfW += br
    hEffCw.SetBinContent(iPt+1, (effC/sumOfW))
    hEffCw.SetBinError(iPt+1, np.sqrt(uncEffC)/sumOfW)
    hEffBw.SetBinContent(iPt+1, (effB/sumOfW))
    hEffBw.SetBinError(iPt+1, np.sqrt(uncEffB)/sumOfW)

outFile = TFile(args.outFileName, 'recreate')
hEffCw.Write()
hEffBw.Write()
outFile.Close()