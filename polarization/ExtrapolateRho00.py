'''
Script for the extrapolation of the rho_00 parameter to fprompt (f_nonprompt) = 0, 1 (1, 0)
run: python ExtrapolateRho00.py cfgFileName.yml
'''

import sys
import argparse
import numpy as np
import yaml
from ROOT import TFile, TLine, TGraphAsymmErrors, TH1D, TCanvas, TF1, TVirtualFitter, TLatex, TLegend, gROOT
sys.path.append('../')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, DivideCanvas, GetROOTColor, GetROOTMarker

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('--batch', action='store_true', default=False,
                    help='flag to run the script in batch mode')
args = parser.parse_args()

SetGlobalStyle(
    padleftmargin=0.15,
    padbottommargin=0.15,
)

if args.batch:
    gROOT.SetBatch(True)

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

inFilesRho00 = inputCfg['inputs']['rho00']['files']
inFilesFrac = inputCfg['inputs']['fractions']['files']
histoNameRho00 = inputCfg['inputs']['rho00']['histoname']
histoNamesFrac = inputCfg['inputs']['fractions']['histonames']
ptMins = inputCfg['inputs']['ptbins']['min']
ptMaxs = inputCfg['inputs']['ptbins']['max']

if len(inFilesRho00) != len(inFilesFrac):
    print('ERROR: inconsistent number of files, please check! Exit')
    sys.exit()

hRhoVsPt, hPromptFracVsPt, hNonPromptFracVsPt = ([] for _ in range(3))
for iFile, fileNameRho in enumerate(inFilesRho00):
    inFile = TFile.Open(fileNameRho)
    hRhoVsPt.append(inFile.Get(histoNameRho00))
    hRhoVsPt[iFile].SetName(f'hRhoVsPt_{iFile}')
    hRhoVsPt[iFile].SetDirectory(0)
    inFile.Close()

for iFile, fileNameFrac in enumerate(inFilesFrac):
    inFile = TFile.Open(fileNameFrac)
    hPromptFracVsPt.append(inFile.Get(histoNamesFrac['prompt']))
    hPromptFracVsPt[iFile].SetName(f'hPromptFracVsPt_{iFile}')
    hPromptFracVsPt[iFile].SetDirectory(0)
    hNonPromptFracVsPt.append(inFile.Get(histoNamesFrac['nonprompt']))
    hNonPromptFracVsPt[iFile].SetName(f'hPromptFracVsPt_{iFile}')
    hNonPromptFracVsPt[iFile].SetDirectory(0)
    inFile.Close()

lineNoPol = TLine(0., 1./3, 1., 1./3)
lineNoPol.SetLineWidth(2)
lineNoPol.SetLineStyle(9)
lineNoPol.SetLineColor(GetROOTColor('kGrey+2'))

lineNoPolVsPt = TLine(ptMins[0], 1./3, ptMaxs[-1], 1./3)
lineNoPolVsPt.SetLineWidth(2)
lineNoPolVsPt.SetLineStyle(9)
lineNoPolVsPt.SetLineColor(GetROOTColor('kGrey+2'))

nPtBins = len(ptMins)
cRho00VsPromptFrac = TCanvas('cRho00VsPromptFrac', '', 1200, 400)
cRho00VsNonPromptFrac = TCanvas('cRho00VsNonPromptFrac', '', 1200, 400)
DivideCanvas(cRho00VsPromptFrac, nPtBins)
DivideCanvas(cRho00VsNonPromptFrac, nPtBins)
gRho00VsPromptFrac, gRho00VsNonPromptFrac, fRho00VsPromptFrac, fRho00VsNonPromptFrac = ([] for _ in range(4))
ciRho00VsPromptFrac, ciRho00VsNonPromptFrac = [], []

ptLims = ptMins.copy()
ptLims.append(ptMaxs[-1])
ptLims = np.array(ptLims, 'f')

hRho00Charm = TH1D('hRho00Charm', ';#it{p}_{T} (GeV/#it{c});#it{#rho}_{00}', nPtBins, ptLims)
SetObjectStyle(hRho00Charm, color=GetROOTColor('kRed+1'), marker=GetROOTMarker('kFullSquare'))
hRho00Beauty = TH1D('hRho00Beauty', ';#it{p}_{T} (GeV/#it{c});#it{#rho}_{00}', nPtBins, ptLims)
SetObjectStyle(hRho00Beauty, color=GetROOTColor('kAzure+4'), marker=GetROOTMarker('kFullCircle'))

lat = TLatex()
lat.SetNDC()
lat.SetTextFont(42)
lat.SetTextSize(0.045)
lat.SetTextColor(GetROOTColor('kBlack'))

for iPt, (ptMin, ptMax) in enumerate(zip(ptMins, ptMaxs)):

    gRho00VsPromptFrac.append(TGraphAsymmErrors(1))
    gRho00VsNonPromptFrac.append(TGraphAsymmErrors(1))
    gRho00VsPromptFrac[iPt].SetNameTitle(
        f'gRho00VsPromptFrac_pt{ptMin}-{ptMax}',
        ';#it{f}_{prompt};#it{#rho}_{00}'
    )
    gRho00VsNonPromptFrac[iPt].SetNameTitle(
        f'gRho00VsPromptFrac_pt{ptMin}-{ptMax}',
        ';#it{f}_{non-prompt};#it{#rho}_{00}'
    )
    SetObjectStyle(gRho00VsPromptFrac[iPt])
    SetObjectStyle(gRho00VsNonPromptFrac[iPt])

    ciRho00VsPromptFrac.append(TH1D(f'ciRho00VsPromptFrac_pt{ptMin}-{ptMax}', '', 1000, 0., 1.))
    ciRho00VsNonPromptFrac.append(TH1D(f'ciRho00VsNonPromptFrac_pt{ptMin}-{ptMax}', '', 1000, 0., 1.))
    SetObjectStyle(ciRho00VsPromptFrac[iPt], alpha=0.3, color=GetROOTColor('kAzure+4'), markerstyle=0)
    SetObjectStyle(ciRho00VsNonPromptFrac[iPt], alpha=0.3, color=GetROOTColor('kAzure+4'), markerstyle=0)

    fRho00VsPromptFrac.append(TF1(f'fRho00VsPromptFrac_pt{ptMin}-{ptMax}', 'pol1', 0., 1.))
    fRho00VsNonPromptFrac.append(TF1(f'fRho00VsNonPromptFrac_pt{ptMin}-{ptMax}', 'pol1', 0., 1.))
    fRho00VsPromptFrac[iPt].SetLineWidth(2)
    fRho00VsPromptFrac[iPt].SetLineColor(GetROOTColor('kAzure+4'))
    fRho00VsNonPromptFrac[iPt].SetLineWidth(2)
    fRho00VsNonPromptFrac[iPt].SetLineColor(GetROOTColor('kAzure+4'))

    for iPoint, (hRho, hPFrac, hNPFrac) in enumerate(zip(hRhoVsPt, hPromptFracVsPt, hNonPromptFracVsPt)):
        ptCent = (ptMin + ptMax) / 2
        binRho = hRho.GetXaxis().FindBin(ptCent)
        binFrac = hPFrac.GetXaxis().FindBin(ptCent)
        rho00 = hRho.GetBinContent(binRho)
        uncRho00 = hRho.GetBinError(binRho)
        promptFrac = hPFrac.GetBinContent(binFrac)
        nonpromptFrac = hNPFrac.GetBinContent(binFrac)
        uncPromptFrac = hPFrac.GetBinError(binFrac)
        uncNonpromptFrac = hNPFrac.GetBinError(binFrac)

        gRho00VsPromptFrac[iPt].SetPoint(iPoint, promptFrac, rho00)
        gRho00VsPromptFrac[iPt].SetPointError(iPoint, uncPromptFrac, uncPromptFrac, uncRho00, uncRho00)

        gRho00VsNonPromptFrac[iPt].SetPoint(iPoint, nonpromptFrac, rho00)
        gRho00VsNonPromptFrac[iPt].SetPointError(iPoint, uncNonpromptFrac, uncNonpromptFrac, uncRho00, uncRho00)

    cRho00VsPromptFrac.cd(iPt+1).DrawFrame(
        0., 0., 1., 0.6,
        ';#it{f}_{prompt};#it{#rho}_{00}'
    )
    lineNoPol.Draw('same')
    gRho00VsPromptFrac[iPt].Fit(f'fRho00VsPromptFrac_pt{ptMin}-{ptMax}', '0')
    TVirtualFitter.GetFitter().GetConfidenceIntervals(ciRho00VsPromptFrac[iPt], 0.68)
    ciRho00VsPromptFrac[iPt].DrawCopy('same')
    fRho00VsPromptFrac[iPt].Draw('same')
    gRho00VsPromptFrac[iPt].Draw('PZ')
    lat.DrawLatex(0.5, 0.4, f'{ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}}')
    lat.DrawLatex(
        0.5,
        0.3,
        f'#chi^{{2}}/ndf = {fRho00VsPromptFrac[iPt].GetChisquare():0.2f}/{fRho00VsPromptFrac[iPt].GetNDF()}'
    )

    cRho00VsNonPromptFrac.cd(iPt+1).DrawFrame(
        0., 0., 1., 0.6,
        ';#it{f}_{non-prompt};#it{#rho}_{00}'
    )
    lineNoPol.Draw('same')
    gRho00VsNonPromptFrac[iPt].Fit(f'fRho00VsNonPromptFrac_pt{ptMin}-{ptMax}', '0')
    TVirtualFitter.GetFitter().GetConfidenceIntervals(ciRho00VsNonPromptFrac[iPt], 0.68)
    ciRho00VsNonPromptFrac[iPt].DrawCopy('same')
    fRho00VsNonPromptFrac[iPt].Draw('same')
    gRho00VsNonPromptFrac[iPt].Draw('PZ')
    lat.DrawLatex(0.5, 0.4, f'{ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}}')
    lat.DrawLatex(
        0.5,
        0.3,
        f'#chi^{{2}}/ndf = {fRho00VsNonPromptFrac[iPt].GetChisquare():0.2f}/{fRho00VsNonPromptFrac[iPt].GetNDF()}'
    )

    hRho00Charm.SetBinContent(iPt+1, ciRho00VsNonPromptFrac[iPt].GetBinContent(1))
    hRho00Charm.SetBinError(iPt+1, ciRho00VsNonPromptFrac[iPt].GetBinError(1))

    hRho00Beauty.SetBinContent(iPt+1, ciRho00VsNonPromptFrac[iPt].GetBinContent(1000))
    hRho00Beauty.SetBinError(iPt+1, ciRho00VsNonPromptFrac[iPt].GetBinError(1000))

cRho00VsPromptFrac.Modified()
cRho00VsPromptFrac.Update()

cRho00VsNonPromptFrac.Modified()
cRho00VsNonPromptFrac.Update()

cRho00vsPt = TCanvas('cRho00vsPt', '', 500, 500)
cRho00vsPt.DrawFrame(
    ptMins[0], 0., ptMaxs[-1], 0.6,
    ';#it{p}_{T} (GeV/#it{c});#it{#rho}_{00}'
)
lineNoPolVsPt.Draw('same')
hRho00Charm.DrawCopy('same')
hRho00Beauty.DrawCopy('same')

leg = TLegend(0.2, 0.2, 0.5, 0.35)
leg.SetTextSize(0.045)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.AddEntry(hRho00Charm, 'c #rightarrow D*^{+}', 'pl')
leg.AddEntry(hRho00Beauty, 'b #rightarrow D*^{+}', 'pl')
leg.Draw()

# output
outFileName = inputCfg['output']['file']

cRho00vsPt.SaveAs(f'{outFileName}.pdf')
cRho00VsPromptFrac.SaveAs(f'{outFileName}_vs_promptfrac.pdf')
cRho00VsNonPromptFrac.SaveAs(f'{outFileName}_vs_nonpromptfrac.pdf')

outFile = TFile.Open(f'{outFileName}.root', 'recreate')
cRho00vsPt.Write()
hRho00Charm.Write()
hRho00Beauty.Write()
cRho00VsPromptFrac.Write()
cRho00VsNonPromptFrac.Write()
for iPt, _ in enumerate(ptMins):
    gRho00VsPromptFrac[iPt].Write()
    gRho00VsNonPromptFrac[iPt].Write()
    fRho00VsPromptFrac[iPt].Write()
    fRho00VsNonPromptFrac[iPt].Write()
    ciRho00VsPromptFrac[iPt].Write()
    ciRho00VsNonPromptFrac[iPt].Write()

if not args.batch:
    input('Press enter to exit')
