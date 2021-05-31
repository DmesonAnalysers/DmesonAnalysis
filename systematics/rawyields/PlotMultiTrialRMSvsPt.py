'''
Script for plotting the RMS and shift of the raw-yield multi-trial distributions vs pT
run: python PlotMultiTrialRMSvsPt MultiTrial.root RawYieldsDefault.root Output.pdf
'''

import sys
import argparse
from ROOT import TCanvas, TFile, TGaxis, gPad, TMath # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

# set global style
SetGlobalStyle(padleftmargin=0.14, padrightmargin=0.14, padbottommargin=0.14, titleoffsety=1.2, padticky=0)

# load inputs
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('inFileNameMultiTrial', metavar='text', default='MultiTrial.root')
parser.add_argument('inFileNameRawYields', metavar='text', default='RawYieldsDefault.root')
parser.add_argument('outFileName', metavar='text', default='Output.pdf')
args = parser.parse_args()

inFileRaw = TFile.Open(args.inFileNameRawYields)
hRawYields = inFileRaw.Get('hRawYields')
hSoverB = inFileRaw.Get('hRawYieldsSoverB')
hRawYields.SetDirectory(0)
hSoverB.SetDirectory(0)
SetObjectStyle(hRawYields, color=kBlack, fillstyle=0)
SetObjectStyle(hSoverB, color=kRed+1, fillstyle=0)
inFileRaw.Close()

nPtBins = hRawYields.GetNbinsX()

hRawYieldDistr = []
hRMS = hRawYields.Clone('hRMS')
hRMS.Reset()
hMeanShift = hRawYields.Clone('hMeanShift')
hMeanShift.Reset()
hSyst = hRawYields.Clone('hSyst')
hSyst.Reset()
inFileMT = TFile.Open(args.inFileNameMultiTrial)
for iPt in range(1, nPtBins+1):
    ptMin = hRawYields.GetBinLowEdge(iPt)
    ptMax = ptMin+hRawYields.GetBinWidth(iPt)
    hRawYieldDistr.append(inFileMT.Get(f'hRawYield_pT_{ptMin*10:.0f}-{ptMax*10:.0f}'))
    rms =  hRawYieldDistr[-1].GetRMS() / hRawYieldDistr[-1].GetMean() * 100
    shift = TMath.Abs(hRawYields.GetBinContent(iPt) - hRawYieldDistr[-1].GetMean()) / hRawYieldDistr[-1].GetMean() * 100
    syst = TMath.Sqrt(rms**2 + shift**2)
    hRMS.SetBinContent(iPt, rms)
    hMeanShift.SetBinContent(iPt, shift)
    hSyst.SetBinContent(iPt, syst)

cRMS = TCanvas('cRMS', '', 800, 800)
cRMS.DrawFrame(hRMS.GetBinLowEdge(1), 0.1, ptMax, 20., ';#it{p}_{T} (GeV/#it{c});RMS (%)')
hRMS.DrawCopy('same')
axisSoverB = TGaxis(gPad.GetUxmax(), gPad.GetUymin(), gPad.GetUxmax(), gPad.GetUymax(), 0.01, 20., 510, "+LG")
axisSoverB.SetLineColor(kRed+1)
axisSoverB.SetLabelColor(kRed+1)
axisSoverB.SetLabelFont(42)
axisSoverB.SetLabelSize(0.045)
axisSoverB.SetTitle('S/B (3#sigma)')
axisSoverB.SetTitleOffset(1.2)
axisSoverB.SetLabelOffset(0.012)
axisSoverB.SetTitleColor(kRed+1)
axisSoverB.SetTitleFont(42)
axisSoverB.SetTitleSize(0.05)
axisSoverB.SetMaxDigits(3)
axisSoverB.Draw()
hSoverB.DrawCopy('same')
cRMS.Update()

cShift = TCanvas('cShift', '', 800, 800)
cShift.DrawFrame(hMeanShift.GetBinLowEdge(1), 0.1, ptMax, 20., ';#it{p}_{T} (GeV/#it{c}); Mean Shift (%)')
hMeanShift.DrawCopy('same')
cShift.Update()

cSyst = TCanvas('cSyst', '', 800, 800)
cSyst.DrawFrame(hSyst.GetBinLowEdge(1), 0.1, ptMax, 20., ';#it{p}_{T} (GeV/#it{c}); Syst (%)')
hSyst.DrawCopy('same')
cSyst.Update()

outFile = TFile(f'{args.outFileName}.root', 'recreate')
cRMS.Write()
cShift.Write()
cSyst.Write()

input('Press Enter to exit')
