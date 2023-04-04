'''
Script for the combination of efficiency and acceptance factors
run: python CombineAccTimesEff.py effFileName.root accFileName.root outFileName.root
'''

import argparse
import numpy as np
from ROOT import gROOT, TFile, TCanvas, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kFullDiamond, kRed, kAzure, kFullCircle, kOpenSquare  # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.AnalysisUtils import ComputeEfficiency

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('effFileName', metavar='text', default='')
parser.add_argument('accFileName', metavar='text', default='')
parser.add_argument('outFileName', metavar='text', default='')
parser.add_argument("--batch", help="suppress video output", action="store_true")
args = parser.parse_args()

gROOT.SetBatch(args.batch)
SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, titlesize=0.045, labelsize=0.04)

effFile = TFile.Open(args.effFileName)
hEffPrompt = effFile.Get('hEffPrompt')
hEffFD = effFile.Get('hEffFD')
SetObjectStyle(hEffPrompt, color=kRed+1, markerstyle=kFullCircle)
SetObjectStyle(hEffFD, color=kAzure+4, markerstyle=kOpenSquare, markersize=1.5, linewidh=2, linestyle=7)

accFile = TFile.Open(args.accFileName)
hAccNum = accFile.Get('hPtGenAcc')
hAccDen = accFile.Get('hPtGenLimAcc')
hAcc = hEffPrompt.Clone('hAcc')
hAcc.SetMarkerColor(kBlack)
hAcc.SetLineColor(kBlack)
hAcc.SetMarkerStyle(kFullDiamond)

nPtBins = hEffPrompt.GetNbinsX()
for iPt in range(nPtBins):
    PtMin = hEffPrompt.GetBinLowEdge(iPt+1)
    PtMax = PtMin+hEffPrompt.GetBinWidth(iPt+1)
    PtBinMin = hAccNum.GetXaxis().FindBin(PtMin*1.0001)
    PtBinMax = hAccNum.GetXaxis().FindBin(PtMax*0.9999)

    acc, accUnc = ComputeEfficiency(hAccNum.Integral(PtBinMin, PtBinMax), hAccDen.Integral(PtBinMin, PtBinMax),
                                    np.sqrt(hAccNum.Integral(PtBinMin, PtBinMax)),
                                    np.sqrt(hAccDen.Integral(PtBinMin, PtBinMax)))
    hAcc.SetBinContent(iPt+1, acc)
    hAcc.SetBinError(iPt+1, accUnc)

hAccEffPrompt = hEffPrompt.Clone('hAccEffPrompt')
hAccEffPrompt.Multiply(hAcc)
hAccEffFD = hEffFD.Clone('hAccEffFD')
hAccEffFD.Multiply(hAcc)

hAccEffPrompt.SetTitle(';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon);')
hAccEffFD.SetTitle(';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon);')

cAccEff = TCanvas('cAccEff', '', 800, 800)
cAccEff.DrawFrame(hEffPrompt.GetBinLowEdge(1), 1.e-5, hEffPrompt.GetBinLowEdge(nPtBins) +
                  hEffPrompt.GetBinWidth(nPtBins), 1., ';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon);')
cAccEff.SetLogy()
hAccEffPrompt.Draw('same')
hAccEffFD.Draw('same')

leg = TLegend(0.6, 0.2, 0.8, 0.4)
leg.SetTextSize(0.045)
leg.SetFillStyle(0)
leg.AddEntry(hAccEffPrompt, "Prompt", "p")
leg.AddEntry(hAccEffFD, "Feed-down", "p")
leg.Draw()

outFile = TFile(args.outFileName, 'recreate')
cAccEff.Write()
hEffPrompt.Write()
hEffFD.Write()
hAccEffPrompt.Write()
hAccEffFD.Write()
outFile.Close()

outFileNamePDF = args.outFileName.replace('.root', '.pdf')
cAccEff.SaveAs(outFileNamePDF)

if not args.batch:
    input('Press enter to exit')
