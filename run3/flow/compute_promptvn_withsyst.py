import ROOT
import ctypes
import argparse
import yaml
from ROOT import TFile, TGraphAsymmErrors, TH1F, TCanvas, TLine, TLegend, gStyle, TMath

# Systematic uncertainties
#absSystUncFit = [0.017, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.01, 0.01, 0.01, 0.015, 0.03, 0.03]  # v3 010
absSystUncFit = [0.007, 0.01, 0.004, 0.01, 0.006, 0.015, 0.093, 0.14]  # v3 010
relResolUnc = 0.015  # v2 0-10

# Method to set graph style
def SetGraphStyle(graph, color, markerstyle, dofill=False, fillstyle=1000, markersize=1.5, linewidth=2):
    graph.SetLineColor(color)
    graph.SetMarkerColor(color)
    if dofill:
        graph.SetFillColor(color)
        graph.SetFillStyle(fillstyle)
    else:
        graph.SetFillStyle(0)
    graph.SetMarkerStyle(markerstyle)
    graph.SetMarkerSize(markersize)
    graph.SetLineWidth(linewidth)

# Method to set plots style
def SetStyle():
    gStyle.SetPadRightMargin(0.035)
    gStyle.SetPadLeftMargin(0.15)
    gStyle.SetPadBottomMargin(0.1)
    gStyle.SetPadTopMargin(0.07)
    gStyle.SetTitleSize(0.045, "xy")
    gStyle.SetLabelSize(0.04, "xy")
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetHistLineWidth(2)
    gStyle.SetOptStat(0)

# Method to compute prompt vn and associated systematic uncertainties
def CorrectVnForFeedDown(gVnPromptSyst, gVnObsStat, hfPromptCent, hfPromptMin, hfPromptMax):
    gVnPromptStat = TGraphAsymmErrors(1)

    for iPt in range(gVnObsStat.GetN()):
        pt, vnobs = ctypes.c_double(), ctypes.c_double()
        gVnObsStat.GetPoint(iPt, pt, vnobs)
        vnobs = float(vnobs.value)
        pt = float(pt.value)
        vnobsstatunc = gVnObsStat.GetErrorYlow(iPt)
        ptunclow = gVnObsStat.GetErrorXlow(iPt)
        ptunchigh = gVnObsStat.GetErrorXhigh(iPt)

        fpromptcent = hfPromptCent.GetBinContent(iPt + 1)
        fpromptmin = hfPromptMin.GetBinContent(iPt + 1)
        fpromptmax = hfPromptMax.GetBinContent(iPt + 1)

        vnprompt = 2 * vnobs / (1 + fpromptcent)
        vnpromptstatunc = 2 * vnobsstatunc / (1 + fpromptcent)
        vnpromptmin = 2 * TMath.Sqrt(3) * vnobs / ((TMath.Sqrt(3) + 1) + fpromptmax * (TMath.Sqrt(3) - 1))
        vnpromptmax = 2 * TMath.Sqrt(3) * vnobs / ((TMath.Sqrt(3) - 1) + fpromptmin * (TMath.Sqrt(3) + 1))

        gVnPromptStat.SetPoint(iPt, pt, vnprompt)
        gVnPromptStat.SetPointError(iPt, ptunclow, ptunchigh, vnpromptstatunc, vnpromptstatunc)
        gVnPromptSyst.SetPoint(iPt, pt, vnprompt)
        gVnPromptSyst.SetPointError(iPt, 0.3, 0.3, vnprompt - vnpromptmin, vnpromptmax - vnprompt)

    return gVnPromptStat

# Main function to get prompt vn with all systematic uncertainties
def ComputePromptVnWithSyst(outDir, inFileNameVn, inFileNamefPrompt, suffix):
    SetStyle()

    nPtBins = len(absSystUncFit)

    inFileVn = TFile.Open(inFileNameVn)
    if not inFileVn or not inFileVn.IsOpen():
        return
    gVnObsStat = inFileVn.Get("gvnSimFit")
    if not gVnObsStat:
        print(f"ERROR: graph gvnSimFit does not exist! Exit")
        return
    gVnObsStat.SetName("gVnObsStat")
    SetGraphStyle(gVnObsStat, ROOT.kRed, ROOT.kOpenCircle)
    if nPtBins != gVnObsStat.GetN():
        print(f"ERROR: different number of pT bins in systematic-uncertainty array and input vn graph {nPtBins, gVnObsStat.GetN()}! Exit")
        return

    hfPromptCent, hfPromptMin, hfPromptMax = TGraphAsymmErrors(1), TGraphAsymmErrors(1), TGraphAsymmErrors(1)
    gVnPromptStat, gVnPromptSystFeedDown, gfPrompt = TGraphAsymmErrors(1), TGraphAsymmErrors(1), TGraphAsymmErrors(1)
    inFilefPrompt = TGraphAsymmErrors(1)
    if inFileNamefPrompt:
        inFilefPrompt = TFile.Open(inFileNamefPrompt)
        if not inFilefPrompt or not inFilefPrompt.IsOpen():
            return
        hfPromptCent = inFilefPrompt.Get("hprompt_frac")
        hfPromptMin = inFilefPrompt.Get("hprompt_frac_min")
        hfPromptMax = inFilefPrompt.Get("hprompt_frac_max")
        hfPromptCent.SetDirectory(0)
        hfPromptMin.SetDirectory(0)
        hfPromptMax.SetDirectory(0)
        if gVnObsStat.GetN() != hfPromptCent.GetNbinsX():
            print("ERROR: different number of pT bins in input vn graph and fprompt histo! Exit")
            return
        gVnPromptStat = CorrectVnForFeedDown(gVnPromptSystFeedDown, gVnObsStat, hfPromptCent, hfPromptMin, hfPromptMax)
        gVnPromptStat.SetName("gVnPromptStat")
        gVnPromptSystFeedDown.SetName("gVnPromptSystFeedDown")

        gfPrompt = TGraphAsymmErrors(1)
        gfPrompt.SetName("gfPrompt")
        SetGraphStyle(gfPrompt, ROOT.kBlack, ROOT.kFullCircle)
        for iPt in range(nPtBins):
            gfPrompt.SetPoint(iPt, hfPromptCent.GetBinCenter(iPt + 1), hfPromptCent.GetBinContent(iPt + 1))
            gfPrompt.SetPointError(iPt, hfPromptCent.GetBinWidth(iPt + 1) / 2, hfPromptCent.GetBinWidth(iPt + 1) / 2,
                                   hfPromptCent.GetBinContent(iPt + 1) - hfPromptMin.GetBinContent(iPt + 1),
                                   hfPromptMax.GetBinContent(iPt + 1) - hfPromptCent.GetBinContent(iPt + 1))
    else:
        gVnPromptStat = gVnObsStat.Clone("gVnPromptStat")
        gVnPromptSystFeedDown = gVnObsStat.Clone("gVnPromptSystFeedDown")
        for iPt in range(nPtBins):
            gVnPromptSystFeedDown.SetPointError(iPt, 0.3, 0.3, 0., 0.)

    SetGraphStyle(gVnPromptStat, ROOT.kBlack, ROOT.kFullCircle)
    SetGraphStyle(gVnPromptSystFeedDown, ROOT.kGray + 2, ROOT.kFullCircle, True, 3154)
    

    gVnPromptSystFit = gVnPromptStat.Clone("gVnPromptSystFit")
    gVnPromptSystResol = gVnPromptStat.Clone("gVnPromptSystResol")
    gVnPromptSystData = gVnPromptStat.Clone("gVnPromptSystData")
    SetGraphStyle(gVnPromptSystData, ROOT.kAzure+4, ROOT.kFullCircle, True, 3145)
    gVnPromptSystTot = gVnPromptStat.Clone("gVnPromptSystTot")

    ptmin, ptmax = -1., -1.
    for iPt in range(nPtBins):
        pt, vnobs, vnprompt = ctypes.c_double(), ctypes.c_double(), ctypes.c_double()
        gVnPromptStat.GetPoint(iPt, pt, vnprompt)
        gVnObsStat.GetPoint(iPt, pt, vnobs)
        pt = float(pt.value)
        vnprompt = float(vnprompt.value)
        vnobs = float(vnobs.value)
        scalefactor = vnprompt / vnobs
        fitsyst = absSystUncFit[iPt] * scalefactor
        resolsyst = relResolUnc * vnprompt
        feeddownsystlow = gVnPromptSystFeedDown.GetErrorYlow(iPt)
        feeddownsysthigh = gVnPromptSystFeedDown.GetErrorYhigh(iPt)
        datasyst = TMath.Sqrt(fitsyst * fitsyst + resolsyst * resolsyst)
        totsystlow = TMath.Sqrt(datasyst * datasyst + feeddownsystlow * feeddownsystlow)
        totsysthigh = TMath.Sqrt(datasyst * datasyst + feeddownsysthigh * feeddownsysthigh)
        gVnPromptSystFit.SetPointError(iPt, 0.5, 0.5, fitsyst, fitsyst)
        gVnPromptSystResol.SetPointError(iPt, 0.5, 0.5, resolsyst, resolsyst)
        gVnPromptSystData.SetPointError(iPt, 0.5, 0.5, datasyst, datasyst)
        gVnPromptSystTot.SetPointError(iPt, 0.5, 0.5, totsystlow, totsysthigh)
        if iPt == 0:
            ptmin = pt - gVnObsStat.GetErrorXlow(iPt)
        if iPt == nPtBins - 1:
            ptmax = pt + gVnObsStat.GetErrorXhigh(iPt)

    outFileName = f"{outDir}/promptvn_withsyst{suffix}.root"
    outFile = TFile(outFileName, "recreate")
    outFile.cd()

    if inFileNamefPrompt:
        gVnObsStat.Write()
        gVnPromptStat.Write()
        gVnPromptSystFeedDown.Write()
        gfPrompt.Write()

    gVnPromptSystFit.Write()
    gVnPromptSystResol.Write()
    gVnPromptSystData.Write()
    gVnPromptSystTot.Write()

    if inFileNamefPrompt:
        c = TCanvas("c", "c", 800, 600)
        c.cd()
        haxis = TH1F("haxis", "; #it{p}_{T} (GeV/#it{c});", 100, ptmin, ptmax)
        haxis.GetYaxis().SetRangeUser(-0.002, 0.3)
        haxis.Draw("axis")
        gVnPromptSystFeedDown.Draw("e2 same")
        gVnPromptSystData.Draw("e2 same")
        gVnPromptSystTot.Draw("e2 same")
        gVnPromptStat.Draw("ep same")
        gVnObsStat.Draw("ep same")
        line = TLine(ptmin, 0, ptmax, 0)
        line.SetLineStyle(2)
        line.SetLineWidth(2)
        line.Draw("same")
        legend = TLegend(0.53, 0.63, 0.83, 0.83)
        legend.SetTextSize(0.04)
        legend.SetTextFont(42)
        legend.AddEntry(gVnObsStat, "observed", "ep")
        legend.AddEntry(gVnPromptStat, "prompt", "ep")
        legend.AddEntry(gVnPromptSystFeedDown, "syst. unc. (FD corr.)", "f")
        legend.AddEntry(gVnPromptSystData, "syst. unc. (fit and reso. corr.)", "f")
        legend.AddEntry(gVnPromptSystTot, "syst. unc. (total)", "f")
        legend.Draw("same")
        c.Update()
        input()
        c.Write()

    outFile.Close()
    inFileVn.Close()
    if inFilefPrompt:
        inFilefPrompt.Close()

    print(f"Output written to: {outFileName}")
    input("Press Enter to exit.")

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

    ComputePromptVnWithSyst(args.outputdir,
                            args.vn_file,
                            args.frac_file,
                            args.suffix)
