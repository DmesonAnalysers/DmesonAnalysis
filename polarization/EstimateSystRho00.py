'''
Script for the evaluation of systematic uncertainties on the rho_00 parameter
run: python EstimateSystRho00.py cfgFileName.yml
'''

import os
import sys
import argparse
import yaml
from ROOT import TFile, TLine, TCanvas, TLegend, gROOT
sys.path.append('../')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor

colors = [
    GetROOTColor("kBlack"),
    GetROOTColor("kAzure+4"),
    GetROOTColor("kRed+1"),
    GetROOTColor("kGreen+2"),
    GetROOTColor("kMagenta+1"),
    GetROOTColor("kOrange+7"),
    GetROOTColor("kAzure+2"),
    GetROOTColor("kRed-7"),
    GetROOTColor("kSpring-5"),
    GetROOTColor("kViolet+5"),
    GetROOTColor("kOrange-2"),
    GetROOTColor("kBlue+2"),
    GetROOTColor("kRed+3"),
    GetROOTColor("kGreen+3"),
    GetROOTColor("kMagenta+3")
]

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('--batch', action='store_true', default=False,
                    help='flag to run the script in batch mode')
args = parser.parse_args()

SetGlobalStyle(
    padleftmargin=0.16,
    padbottommargin=0.15,
    titleoffsety=1.5
)

if args.batch:
    gROOT.SetBatch(True)

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

cSyst, hRho00Charm, hRho00Beauty, hRho00CharmRatio, hRho00BeautyRatio, legSyst = ({} for _ in range(6))
for syst in inputCfg:
    inputDir = inputCfg[syst]['inputs']['directory']
    inputFiles = inputCfg[syst]['inputs']['files']
    labels = inputCfg[syst]['labels']

    if len(inputFiles) != len(labels):
        print(f'ERROR: inconsistent number of files and labels for {syst} syst please check! Exit')
        sys.exit()

    hRho00Charm[syst] = []
    hRho00Beauty[syst] = []
    hRho00CharmRatio[syst] = []
    hRho00BeautyRatio[syst] = []

    cSyst[syst] = TCanvas(f'cSyst_{syst}', '', 1000, 1000)
    cSyst[syst].Divide(2, 2)

    legSyst[syst] = TLegend(0.2, 0.2, 0.4, 0.5)
    legSyst[syst].SetTextSize(0.035)
    legSyst[syst].SetBorderSize(0)
    legSyst[syst].SetFillStyle(0)
    if len(inputFiles) > 5:
        legSyst[syst].SetNColumns(2)
        legSyst[syst].SetX2(0.8)

    for iFile, inFileName in enumerate(inputFiles):
        inFile = TFile.Open(os.path.join(inputDir, inFileName))
        hRho00Charm[syst].append(inFile.Get('hRho00Charm'))
        hRho00Charm[syst][iFile].SetDirectory(0)
        hRho00Charm[syst][iFile].SetName(f'hRho00Charm_{syst}')
        SetObjectStyle(hRho00Charm[syst][iFile], color=colors[iFile])
        hRho00Beauty[syst].append(inFile.Get('hRho00Beauty'))
        hRho00Beauty[syst][iFile].SetDirectory(0)
        hRho00Beauty[syst][iFile].SetName(f'hRho00Beauty_{syst}')
        SetObjectStyle(hRho00Beauty[syst][iFile], color=colors[iFile])
        legSyst[syst].AddEntry(hRho00Charm[syst][iFile], labels[iFile], 'p')
        nPtBins = hRho00Charm[syst][iFile].GetNbinsX()
        ptMin = hRho00Charm[syst][iFile].GetXaxis().GetBinLowEdge(1)
        ptMax = hRho00Charm[syst][iFile].GetXaxis().GetBinUpEdge(nPtBins)

        lineNoPolVsPt = TLine(ptMin, 1./3, ptMax, 1./3)
        lineNoPolVsPt.SetLineWidth(2)
        lineNoPolVsPt.SetLineStyle(9)
        lineNoPolVsPt.SetLineColor(GetROOTColor('kGrey+2'))

        if iFile == 0:
            hFrameCharm = cSyst[syst].cd(1).DrawFrame(
                ptMin,
                0.,
                ptMax,
                0.6,
                ';#it{p}_{T} (GeV/#it{c}; #it{#rho}_{00} (charm)'
            )
            hFrameCharm.GetYaxis().SetDecimals()
            hFrameBeauty = cSyst[syst].cd(2).DrawFrame(
                ptMin,
                0.,
                ptMax,
                0.6,
                ';#it{p}_{T} (GeV/#it{c}; #it{#rho}_{00} (beauty)'
            )
            hFrameBeauty.GetYaxis().SetDecimals()
            hFrameCharmRatio = cSyst[syst].cd(3).DrawFrame(
                ptMin,
                0.8,
                ptMax,
                1.2,
                ';#it{p}_{T} (GeV/#it{c}; #it{#rho}_{00} / #it{#rho}_{00}^{def} (charm)'
            )
            hFrameCharmRatio.GetYaxis().SetDecimals()
            hFrameBeautyRatio = cSyst[syst].cd(4).DrawFrame(
                ptMin,
                0.8,
                ptMax,
                1.2,
                ';#it{p}_{T} (GeV/#it{c}; #it{#rho}_{00} / #it{#rho}_{00}^{def} (beauty)'
            )
            hFrameBeautyRatio.GetYaxis().SetDecimals()

        cSyst[syst].cd(1)
        lineNoPolVsPt.Draw('same')
        hRho00Charm[syst][iFile].DrawCopy('same')
        cSyst[syst].cd(2)
        lineNoPolVsPt.Draw('same')
        hRho00Beauty[syst][iFile].DrawCopy('same')
        legSyst[syst].Draw()
        if iFile > 0:
            hRho00CharmRatio[syst].append(hRho00Charm[syst][iFile].Clone(f'hRho00CharmRatio_{syst}_{iFile}'))
            hRho00CharmRatio[syst][iFile].SetDirectory(0)
            hRho00BeautyRatio[syst].append(hRho00Beauty[syst][iFile].Clone(f'hRho00BeautyRatio_{syst}_{iFile}'))
            hRho00BeautyRatio[syst][iFile].SetDirectory(0)
            hRho00CharmRatio[syst][iFile].Divide(hRho00Charm[syst][iFile], hRho00Charm[syst][0])
            hRho00BeautyRatio[syst][iFile].Divide(hRho00Beauty[syst][iFile], hRho00Beauty[syst][0])
            for iPt in range(1, hRho00CharmRatio[syst][iFile].GetNbinsX()+1):
                hRho00CharmRatio[syst][iFile].SetBinError(iPt, 1.e-20)
                hRho00BeautyRatio[syst][iFile].SetBinError(iPt, 1.e-20)
            cSyst[syst].cd(3)
            hRho00CharmRatio[syst][iFile].DrawCopy('same')
            cSyst[syst].cd(4)
            hRho00BeautyRatio[syst][iFile].DrawCopy('same')
        else:
            hRho00CharmRatio[syst].append(None)
            hRho00BeautyRatio[syst].append(None)

    cSyst[syst].Modified()
    cSyst[syst].Update()
    cSyst[syst].SaveAs(os.path.join(inputDir, f'Rho00_syst_{syst}.pdf'))

outFile = TFile.Open(os.path.join(inputDir, 'Rho00_syst.root'), 'recreate')
for syst in inputCfg:
    cSyst[syst].Write()
    for iFile, _ in enumerate(inputCfg[syst]['inputs']['files']):
        hRho00Charm[syst][iFile].Write()
        hRho00Beauty[syst][iFile].Write()
        if iFile > 0:
            hRho00CharmRatio[syst][iFile].Write()
            hRho00BeautyRatio[syst][iFile].Write()
outFile.Close()
