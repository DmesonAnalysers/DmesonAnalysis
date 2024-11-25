"""
python script for the computation of the prompt or non-prompt v2 via extrapolation
run: python ComputeV2vsFDFrac.py cfgFileName.yml outFileName.root

"""

import sys
import argparse
import uproot
import os
import numpy as np
import yaml
from ROOT import (
    TFile,
    TH1F,
    TH2F,
    TCanvas,
    TLegend,
    TGraphAsymmErrors,
    TLatex,
    gRandom,
    TGraphErrors,
    TF1,
    TH1D,
    TVirtualFitter,
    Double_t,
)  # pylint: disable=import-error,no-name-in-module
from ROOT import (
    kBlack,
    kRed,
    kAzure,
    kGreen,
    kRainBow,
)  # pylint: disable=import-error,no-name-in-module
from ROOT import (
    kFullCircle,
    kFullSquare,
    kOpenSquare,
    kOpenCircle,
    kOpenCross,
    kOpenDiamond,
)  # pylint: disable=import-error,no-name-in-module

sys.path.append('/home/wuct/ALICE/local/DmesonAnalysis')
from utils.StyleFormatter import (
    SetGlobalStyle,
    SetObjectStyle,
    GetROOTColor,
)  # pylint: disable=wrong-import-position,import-error,no-name-in-module

parser = argparse.ArgumentParser(description="Arguments to pass")
parser.add_argument(
    "cfgFileName",
    metavar="text",
    default="cfgFileName.yml",
    help="config file name with root input files",
)
parser.add_argument(
    "outFileName", metavar="text", default="outFile.root", help="output root file name"
)
args = parser.parse_args()

outFileNameFracV2AllPt = args.outFileName.replace(".root", "_Frac_pt.pdf")
outFileNameNpV2PDFFinal = args.outFileName.replace(".root", "_Frac_total_np.pdf")
outFileNamePrV2PDFFinal = args.outFileName.replace(".root", "_Frac_total_pr.pdf")
outFileNameV2PDFFinal = args.outFileName.replace(".root", "_Frac_total_compare.pdf")

with open(args.cfgFileName, "r") as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)

inputFilesV2 = cutSetCfg["v2"]["inputfiles"]
inputFilesFrac = cutSetCfg["npfraction"]["inputfiles"]

inputDirV2 = cutSetCfg["v2"]["inputdir"]
inputDirFrac = cutSetCfg["npfraction"]["inputdir"]

histoNameV2 = cutSetCfg["v2"]["histoname"]
graphNameV2 = cutSetCfg["v2"]["graphname"]

histoNameFracFD = cutSetCfg["npfraction"]["histonamesFD"]
histoNameFracPr = cutSetCfg["npfraction"]["histonamesPrompt"]

nSets = len(cutSetCfg["v2"]["inputfiles"])

(
    hV2,
    hNpFrac,
    hFracvsV2,
    cFrac,
    hNpV2vsPt,
    hPrV2vsPt,
    hGraphV2,
    hPrFrac,
    AverageV2,
    FinalErroXl,
    FinalErroXh,
) = ([], [], [], [], [], [], [], [], [], [], [])

for (
    inFileNameV2,
    inFileNameFrac,
) in zip(inputFilesV2, inputFilesFrac):
    inFileNameV2 = os.path.join(inputDirV2, inFileNameV2)
    inFileV2 = TFile.Open(inFileNameV2)
    hV2.append(inFileV2.Get(histoNameV2))
    AverageV2.append(inFileV2.Get(graphNameV2))
    hV2[-1].SetDirectory(0)

    inFileNameFrac = os.path.join(inputDirFrac, inFileNameFrac)
    inFileFrac = TFile.Open(inFileNameFrac)
    hNpFrac.append(inFileFrac.Get(histoNameFracFD))
    hNpFrac[-1].SetDirectory(0)
    hPrFrac.append(inFileFrac.Get(histoNameFracPr))
    hPrFrac[-1].SetDirectory(0)

SetGlobalStyle(
    padbottommargin=0.13,
    padleftmargin=0.12,
    titleoffsety=1.0,
    titlesize=0.05,
    labelsize=0.05,
)

legDistr = TLegend(0.45, 0.69, 0.75, 0.89)
legDistr.SetFillStyle(0)
legDistr.SetBorderSize(0)
legDistr.SetTextSize(0.045)

cNonPromptV2 = TCanvas("cNonPromptV2", "non-prompt v2 versus pt")
cNonPromptV2.SetCanvasSize(800, 800)
cPromptV2 = TCanvas("cPromptV2", "prompt v2 versus pt")
cPromptV2.SetCanvasSize(800, 800)
cPromptAndNonPromptV2 = TCanvas(
    "cPromptAndNonPromptV2", "prompt and non-prompt v2 versus pt"
)
cPromptAndNonPromptV2.SetCanvasSize(800, 800)

leg = TLegend(0.55, 0.75, 0.88, 0.89)
leg.SetTextSize(0.045)
leg.SetBorderSize(0)
leg.SetFillStyle(0)

outFile = TFile(args.outFileName, "recreate")
hNpV2vsPt = hV2[0].Clone("hNpV2vsPt")
hPrV2vsPt = hV2[0].Clone("hPrV2vsPt")

for iPt in range(hNpFrac[0].GetNbinsX()):

    hGraphV2.append(TGraphErrors(-1))
    hFracvsV2.append(TH1D(f"hFracvsV2{iPt}", "", 1000, 0.0, 1.0))
    hFracvsV2[-1].SetDirectory(0)

    TotalV2Xl = Double_t()
    TotalV2Xh = Double_t()
    for i in range(4):
        TotalV2Xl += AverageV2[i].GetErrorXlow(iPt)
        TotalV2Xh += AverageV2[i].GetErrorXhigh(iPt)

    FinalErroXl.append(TotalV2Xl / 4)
    FinalErroXh.append(TotalV2Xh / 4)

gV2Np = TGraphAsymmErrors("gNpV2vsPt")
gV2Pr = TGraphAsymmErrors("gPrV2vsPt")

t = TLatex(8, 8, "")
t.SetNDC()
t.SetTextFont(42)
t.SetTextSize(0.047)
t.SetTextColor(kBlack)

for iPt in range(hNpFrac[0].GetNbinsX()):

    ptMin = hNpFrac[0].GetBinLowEdge(iPt + 1)
    ptMax = ptMin + hNpFrac[0].GetBinWidth(iPt + 1)
    ptCent = hNpFrac[0].GetBinCenter(iPt + 1)

    listV2 = [hV2Tmp.GetBinContent(iPt + 1) for hV2Tmp in hV2]
    listV2Unc = [hV2Tmp.GetBinError(iPt + 1) for hV2Tmp in hV2]
    listNpFrac = [hNpFracTmp.GetBinContent(iPt + 1) for hNpFracTmp in hNpFrac]
    listNpFracUnc = [hNpFracTmp.GetBinError(iPt + 1) for hNpFracTmp in hNpFrac]

    for ipoint, (v2, npfrac, v2unc, npfracunc) in enumerate(
        zip(listV2, listNpFrac, listV2Unc, listNpFracUnc)
    ):

        if iPt == 0 and ipoint == 3:  # only 3 points for the first pt bin
            continue
        hGraphV2[iPt].SetPoint(ipoint, npfrac, v2)
        hGraphV2[iPt].SetPointError(ipoint, npfracunc, v2unc)

    pol = TF1("pol", "pol1", 0, 1)
    hGraphV2[iPt].Fit("pol", "", "", 0, 1)
    chi2 = pol.GetChisquare()
    ndf = pol.GetNDF()
    TVirtualFitter.GetFitter().GetConfidenceIntervals(hFracvsV2[iPt], 0.683)
    hFracvsV2[iPt].SetLineColorAlpha(4, 0.15)
    hNpV2vsPt.SetBinContent(
        iPt + 1, hFracvsV2[iPt].GetBinContent(hFracvsV2[iPt].FindBin(1.0) - 1)
    )
    hNpV2vsPt.SetBinError(
        iPt + 1, hFracvsV2[iPt].GetBinError(hFracvsV2[iPt].FindBin(1.0) - 1)
    )
    hPrV2vsPt.SetBinContent(
        iPt + 1, hFracvsV2[iPt].GetBinContent(hFracvsV2[iPt].FindBin(0.0))
    )
    hPrV2vsPt.SetBinError(
        iPt + 1, hFracvsV2[iPt].GetBinError(hFracvsV2[iPt].FindBin(0.0))
    )

    gV2Np.SetPoint(
        iPt,
        hNpV2vsPt.GetBinLowEdge(iPt + 1) + FinalErroXl[iPt],
        hFracvsV2[iPt].GetBinContent(hFracvsV2[iPt].FindBin(1.0) - 1),
    )
    gV2Np.SetPointError(
        iPt,
        FinalErroXl[iPt],
        FinalErroXh[iPt],
        hFracvsV2[iPt].GetBinError(hFracvsV2[iPt].FindBin(1.0) - 1),
        hFracvsV2[iPt].GetBinError(hFracvsV2[iPt].FindBin(1.0) - 1),
    )
    gV2Pr.SetPoint(
        iPt,
        hNpV2vsPt.GetBinLowEdge(iPt + 1) + FinalErroXl[iPt],
        hFracvsV2[iPt].GetBinContent(hFracvsV2[iPt].FindBin(0.0)),
    )
    gV2Pr.SetPointError(
        iPt,
        FinalErroXl[iPt],
        FinalErroXh[iPt],
        hFracvsV2[iPt].GetBinError(hFracvsV2[iPt].FindBin(0.0)),
        hFracvsV2[iPt].GetBinError(hFracvsV2[iPt].FindBin(0.0)),
    )

    ptString = f"{ptMin:.0f} < #it{{p}}_{{T}} < {ptMax:.0f} GeV/#it{{c}}"
    cFrac.append(TCanvas(f"cFrac_{ptString}", "", 500, 500))
    hFrame = cFrac[iPt].DrawFrame(
        0.0,
        0.0001,
        1,
        0.35,
        f"{ptString};Non-prompt D^{{0}} fraction; #it{{v}}_{{2}}^{{#it{{obs}}}}",
    )

    hFrame.GetYaxis().SetDecimals()
    hFrame.GetYaxis().SetNoExponent()
    hFrame.GetXaxis().SetMoreLogLabels()
    hFrame.GetYaxis().SetTitleSize(0.04)
    hFrame.GetYaxis().SetTitleOffset(1.4)
    hFrame.GetYaxis().SetLabelSize(0.04)
    hFrame.GetXaxis().SetTitleSize(0.04)
    hFrame.GetXaxis().SetLabelSize(0.04)
    hFrame.GetXaxis().SetTitleOffset(1.4)
    hFrame.GetYaxis().SetNdivisions(505)

    chi2String = f"#chi^{{2}}/n.d.f = {chi2/ndf:.2f}"
    hGraphV2[iPt].Draw("pZ")
    t.SetTextSize(0.06)
    t.DrawLatex(0.18, 0.88, "ALICE Preliminary")
    t.SetTextSize(0.045)
    t.DrawLatex(
        0.18,
        0.81,
        "30#font[122]{-}50% Pb#font[122]{-}Pb, #sqrt{#it{s}_{NN}} = 5.02 TeV",
    )
    t.SetTextSize(0.04)
    t.DrawLatex(0.18, 0.74, "D^{0}#rightarrow K^{#minus}#pi^{+} and charge conj.")
    t.DrawLatex(0.18, 0.67, f"{ptString}")
    t.SetTextSize(0.035)
    # t.DrawLatex(0.250, 0.23, f'{chi2String}');

    hGraphV2[iPt].SetLineColor(1)
    hGraphV2[iPt].SetLineWidth(2)
    hGraphV2[iPt].SetMarkerStyle(20)
    hGraphV2[iPt].SetMarkerSize(1)
    hGraphV2[iPt].SetMarkerColor(1)

    hFracvsV2[iPt].DrawCopy("same")
    hFracvsV2[iPt].Write()
    hGraphV2[iPt].Write()
    if iPt == 0:
        cFrac[iPt].SaveAs(f"{outFileNameFracV2AllPt}[")
    cFrac[iPt].SaveAs(outFileNameFracV2AllPt)
    if iPt == hNpFrac[0].GetNbinsX() - 1:
        cFrac[iPt].SaveAs(f"{outFileNameFracV2AllPt}]")

SetObjectStyle(hNpV2vsPt, color=GetROOTColor("kAzure+4 "), fillstyle=1)
SetObjectStyle(hPrV2vsPt, color=GetROOTColor("kRed+1"), fillstyle=1)
cNonPromptV2.cd()
hNpV2vsPt.Draw("")
hNpV2vsPt.GetXaxis().SetTitle("#it{p}_{T} GeV/#it{c}")
hNpV2vsPt.GetYaxis().SetTitle("Non-prompt #it{v_{2}}")
hNpV2vsPt.GetYaxis().SetRangeUser(-0.05, 0.35)
hNpV2vsPt.SetMarkerStyle(20)
hNpV2vsPt.SetMarkerSize(2)
hNpV2vsPt.GetYaxis().SetNoExponent()

cPromptV2.cd()
hPrV2vsPt.Draw("")
hPrV2vsPt.GetXaxis().SetTitle("#it{p}_{T} GeV/#it{c}")
hPrV2vsPt.GetYaxis().SetTitle("Prompt #it{v_{2}}")
hPrV2vsPt.GetYaxis().SetRangeUser(-0.05, 0.35)
hPrV2vsPt.SetMarkerStyle(20)
hPrV2vsPt.SetMarkerSize(2)
hPrV2vsPt.GetYaxis().SetNoExponent()

cPromptAndNonPromptV2.cd()
hNpV2vsPt.GetYaxis().SetTitle("#it{v_{2}}")
hNpV2vsPt.Draw("")
hPrV2vsPt.Draw("same")

leg.AddEntry(hNpV2vsPt, "Non-prompt #it{v_{2}}", "lp")
leg.AddEntry(hPrV2vsPt, "Prompt #it{v_{2}}", "lp")
leg.Draw("same")

gV2Np.Write("gNpV2vsPt")
gV2Pr.Write("gPrV2vsPt")
hNpV2vsPt.Write()
hPrV2vsPt.Write()

cNonPromptV2.SaveAs(f"{outFileNameNpV2PDFFinal}]")
cPromptV2.SaveAs(f"{outFileNamePrV2PDFFinal}]")
cPromptAndNonPromptV2.SaveAs(f"{outFileNameV2PDFFinal}]")
