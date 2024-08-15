import ROOT
import numpy as np
import argparse
import yaml
import ctypes
from ROOT import TFile, TGraphAsymmErrors, TH1F, TCanvas, TLegend, TLine, gStyle, TMath

# Systematic uncertainties
absSystUncFit = [0.017, 0.004, 0.008, 0.004, 0.006, 0.007, 0.004, 0.010, 0.009, 0.018, 0.017, 0.036, 0.034]
relResolUnc = 0.01

def ComputePromptVnWithSyst(inFileNameVn, inFileNamefPrompt, outDir, suffix):
    SetStyle()

    inFileVn = TFile.Open(inFileNameVn)
    if not inFileVn or not inFileVn.IsOpen():
        return
    gVnObsStat = inFileVn.Get("gvnSimFit")
    if not gVnObsStat:
        print(f"ERROR: graph gvnSimFit does not exist! Exit")
        return
    gVnObsStat.SetName("gVnObsStat")
    nPtBins = gVnObsStat.GetN()
    SetGraphStyle(gVnObsStat, ROOT.kRed, ROOT.kOpenCircle)

    hfPromptCent = None
    hfPromptMin = None
    hfPromptMax = None
    gVnPromptStat = gVnObsStat.Clone("gVnPromptStat")
    gVnPromptSystFeedDown = gVnObsStat.Clone("gVnPromptSystFeedDown")
    gfPrompt = gVnObsStat.Clone("gfPrompt")
    inFilefPrompt = None
    if inFileNamefPrompt:
        inFilefPrompt = TFile.Open(inFileNamefPrompt)
        if not inFilefPrompt or not inFilefPrompt.IsOpen():
            return
        gfPromptCent = inFilefPrompt.Get("gprompt_frac")
        hfPromptCent = inFilefPrompt.Get("hRawYieldsSimFit")
        hfPromptMin = inFilefPrompt.Get("hRawYieldsSimFit")
        hfPromptMax = inFilefPrompt.Get("hRawYieldsSimFit")
        hfPromptCent.Reset()
        hfPromptMin.Reset()
        hfPromptMax.Reset()

        # convert TGRaphAsymmErrors to TH1F
        for iPt in range(gfPromptCent.GetN()):
            pt, promptfrac = ctypes.c_double(), ctypes.c_double()
            gfPromptCent.GetPoint(iPt, pt, promptfrac)
            hfPromptCent.SetBinContent(iPt+1, promptfrac.value)
            hfPromptMin.SetBinContent(iPt+1, hfPromptCent.GetBinContent(iPt+1) - gfPromptCent.GetErrorYlow(iPt))
            hfPromptMax.SetBinContent(iPt+1, hfPromptCent.GetBinContent(iPt+1) + gfPromptCent.GetErrorYhigh(iPt))


        #if gVnObsStat.GetN() != hfPromptCent.GetNbinsX():
        #    print("ERROR: different number of pT bins in input vn graph and fprompt histo! Exit")
        #    return
        gVnPromptStat = CorrectVnForFeedDown(gVnPromptSystFeedDown, gVnObsStat,
                                             hfPromptCent, hfPromptMin, hfPromptMax)
        gVnPromptStat.SetName("gVnPromptStat")
        gVnPromptSystFeedDown.SetName("gVnPromptSystFeedDown")

        gfPrompt = TGraphAsymmErrors(1)
        gfPrompt.SetName("gfPrompt")
        SetGraphStyle(gfPrompt, ROOT.kBlack, ROOT.kFullCircle)
        for iPt in range(nPtBins):
            gfPrompt.SetPoint(iPt, hfPromptCent.GetBinCenter(iPt+1), hfPromptCent.GetBinContent(iPt+1))
            gfPrompt.SetPointError(iPt, hfPromptCent.GetBinWidth(iPt+1)/2, hfPromptCent.GetBinWidth(iPt+1)/2, hfPromptCent.GetBinContent(iPt+1)-hfPromptMin.GetBinContent(iPt+1), hfPromptMax.GetBinContent(iPt+1)-hfPromptCent.GetBinContent(iPt+1))
    else:
        gVnPromptStat = gVnObsStat.Clone("gVnPromptStat")
        gVnPromptSystFeedDown = gVnObsStat.Clone("gVnPromptSystFeedDown")
        for iPt in range(nPtBins):
            gVnPromptSystFeedDown.SetPointError(iPt, 0.3, 0.3, 0., 0.)

    #SetGraphStyle(gVnPromptStat, ROOT.kBlack, ROOT.kFullCircle)
    #SetGraphStyle(gVnPromptSystFeedDown, ROOT.kGray+2, ROOT.kFullCircle, True)

    gVnPromptSystFit = gVnPromptStat.Clone("gVnPromptSystFit")
    gVnPromptSystResol = gVnPromptStat.Clone("gVnPromptSystResol")
    gVnPromptSystData = gVnPromptStat.Clone("gVnPromptSystData")
    gVnPromptSystTot = gVnPromptStat.Clone("gVnPromptSystTot")

    ptmin = -1
    ptmax = -1
    for iPt in range(nPtBins):
        pt, vnprompt = ctypes.c_double(), ctypes.c_double()
        gVnPromptStat.GetPoint(iPt, pt, vnprompt)
        vnobs = gVnObsStat.GetY()[iPt]
        vnprompt = float(vnprompt.value)
        pt = float(pt.value)
        scalefactor = vnprompt/vnobs
        fitsyst = absSystUncFit[iPt] * scalefactor
        resolsyst = relResolUnc #* vnprompt
        feeddownsystlow = gVnPromptSystFeedDown.GetErrorYlow(iPt)
        feeddownsysthigh = gVnPromptSystFeedDown.GetErrorYhigh(iPt)
        datasyst = TMath.Sqrt(fitsyst**2 + resolsyst**2)
        totsystlow = TMath.Sqrt(datasyst**2 + feeddownsystlow**2)
        totsysthigh = TMath.Sqrt(datasyst**2 + feeddownsysthigh**2)
        gVnPromptSystFit.SetPointError(iPt, 0., 0., fitsyst, fitsyst)
        gVnPromptSystResol.SetPointError(iPt, 0., 0., resolsyst, resolsyst)
        gVnPromptSystData.SetPointError(iPt, 0., 0., datasyst, datasyst)
        gVnPromptSystTot.SetPointError(iPt, 0., 0., totsystlow, totsysthigh)
        if iPt == 0:
            ptmin = pt - gVnObsStat.GetErrorXlow(iPt)
        elif iPt == nPtBins - 1:
            ptmax = pt + gVnObsStat.GetErrorXhigh(iPt)

    lineatzero = TLine(ptmin, 0., ptmax, 0.)
    lineatzero.SetLineWidth(2)
    lineatzero.SetLineColor(ROOT.kBlack)
    lineatzero.SetLineStyle(9)

    legPrompt = TLegend(0.5, 0.75, 0.85, 0.85)
    legPrompt.SetFillStyle(0)
    legPrompt.SetTextSize(0.04)
    legPrompt.SetBorderSize(0)
    legPrompt.AddEntry(gVnObsStat, "Observed #it{v}_{2}", "p")
    legPrompt.AddEntry(gVnPromptStat, "Prompt #it{v}_{2}", "p")

    leg = TLegend(0.5, 0.7, 0.85, 0.85)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.04)
    leg.SetBorderSize(0)
    leg.AddEntry(gVnPromptStat, "Prompt #it{v}_{2}", "p")
    leg.AddEntry(gVnPromptSystData, "Data syst. unc.", "f")
    leg.AddEntry(gVnPromptSystFeedDown, "B feed-down syst. unc.", "f")

    cVnPrompt = TCanvas("cVnPrompt", "", 800, 800)
    cVnPrompt.DrawFrame(ptmin, -0.2, ptmax, 0.4, ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}")
    lineatzero.Draw("same")
    gVnObsStat.Draw("P")
    gVnPromptStat.Draw("P")
    legPrompt.Draw()

    cVn = TCanvas("cVn", "", 800, 800)
    cVn.DrawFrame(ptmin, -0.2, ptmax, 0.4, ";#it{p}_{T} (GeV/#it{c});#it{v}_{2}")
    lineatzero.Draw("same")
    gVnPromptSystFeedDown.Draw("2")
    gVnPromptSystData.Draw("2")
    gVnPromptStat.Draw("P")
    leg.Draw()

    outFileName = f"{outDir}/promptvn_withsyst{suffix}.root"
    outFile = TFile(outFileName, "recreate")
    gVnObsStat.Write()
    gVnPromptStat.Write()
    gVnPromptSystFit.Write()
    gVnPromptSystResol.Write()
    gVnPromptSystData.Write()
    gVnPromptSystFeedDown.Write()
    gVnPromptSystTot.Write()
    if gfPrompt:
        gfPrompt.Write()
    outFile.Close()

def CorrectVnForFeedDown(gVnPromptSyst, gVnObsStat, hfPromptCent, hfPromptMin, hfPromptMax):
    gVnPromptSyst = gVnObsStat.Clone("gVnPromptSyst")
    for iPt in range(gVnObsStat.GetN()):
        vnobs = gVnObsStat.GetY()[iPt]
        vnobserr = max(gVnObsStat.GetErrorYlow(iPt), gVnObsStat.GetErrorYhigh(iPt))
        pt = gVnObsStat.GetX()[iPt]
        ptbin = hfPromptCent.GetXaxis().FindBin(pt)
        bwidth = hfPromptCent.GetXaxis().GetBinWidth(iPt+1)
        fPrompt = hfPromptCent.GetBinContent(iPt+1)
        fPromptMin = hfPromptMin.GetBinContent(ptbin)
        fPromptMax = hfPromptMax.GetBinContent(ptbin)
        fPromptErrLow = fPrompt - fPromptMin
        fPromptErrHigh = fPromptMax - fPrompt
        fPromptErrLow = max(fPromptErrLow, 0.01)
        fPromptErrHigh = max(fPromptErrHigh, 0.01)
        vnPrompt = vnobs / fPrompt
        vnPromptErr = TMath.Sqrt((vnobserr/fPrompt)**2 + (vnPrompt*fPromptErrLow/fPrompt)**2)
        gVnPromptSyst.SetPoint(iPt, pt, vnPrompt)
        gVnPromptSyst.SetPointError(iPt, bwidth/2., bwidth/2., vnPromptErr, vnPromptErr)
    gVnPromptSyst.Draw()
    input("Press Enter to continue...")
    return gVnPromptSyst

def SetStyle():
    gStyle.SetOptStat(0)
    gStyle.SetOptTitle(0)
    gStyle.SetTitleSize(0.045, "xy")
    gStyle.SetLabelSize(0.045, "xy")
    gStyle.SetTitleOffset(1.2, "x")
    gStyle.SetTitleOffset(1.4, "y")
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetPadTopMargin(0.03)
    gStyle.SetPadBottomMargin(0.12)
    gStyle.SetPadLeftMargin(0.12)
    gStyle.SetPadRightMargin(0.04)
    gStyle.SetPadGridX(0)
    gStyle.SetPadGridY(0)
    gStyle.SetFrameFillStyle(0)
    gStyle.SetCanvasColor(10)
    gStyle.SetPadColor(10)
    gStyle.SetTitleColor(1, "XYZ")
    gStyle.SetLabelColor(1, "XYZ")
    gStyle.SetAxisColor(1, "XYZ")
    gStyle.SetLegendBorderSize(0)

def SetGraphStyle(g, col, msty, isBand=False):
    g.SetLineColor(col)
    g.SetMarkerColor(col)
    g.SetMarkerStyle(msty)
    g.SetMarkerSize(1.3)
    if isBand:
        g.SetLineWidth(0)
        g.SetFillColor(col)
        g.SetFillStyle(3001)
    else:
        g.SetLineWidth(2)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("vn_file", metavar="text",
                        default="config.yaml", help="vn file")
    parser.add_argument("frac_file", metavar="text",
                        default="config.yaml", help="frac file")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    ComputePromptVnWithSyst(args.vn_file,
                            args.frac_file,
                            args.outputdir,
                            args.suffix)
