import ROOT
import ctypes
import argparse
import yaml
from ROOT import TFile, TGraphAsymmErrors, TH1F, TCanvas, TLine, TLegend, gStyle, TMath

class Debug:
    enable = False

# Systematic uncertainties
absSystUncFit = [0.017, 0.006, 0.006, 0.006, 0.006, 0.006, 0.006, 0.01, 0.01, 0.01, 0.015, 0.03, 0.03]  # v3 010
relSystUncFit_prompt = [0.197166927240137, 0.05198153724713105, 0.009238430426062955, 0.003807074793593355, 0.00573406806087924, 0.00420063043530908, 0.011980737414319336, 0.014549226030331559, 0.04149920978209411, 0.03639053094367019, 0.04822245036276569, 0.07488006355516015, 0.08312806950388721, 0.20016242362285763, 0.3536755872283735]
relSystUncFit_fd = [0.5983447664451705, 0.32959917437091374, 0.4414108325857921, 0.4623692737694891, 0.06848127616231166, 0.5470543835149732, 0.1955550028330369, 1.078900463525467, 0.19517561852344925, 0.5951551943767623, 0.74038858864252, 0.43742966818775986, 0.6450615616944664, 0.708602993677684, 0.93417399175986]
relResolUnc = 0.015  # v2 0-10

def ReadSysUncFile(filenamesyst, nonprompt):
    gsys_prompt = TGraphAsymmErrors(1)
    gsys_fd = TGraphAsymmErrors(1)

    file = TFile.Open(filenamesyst)

    if nonprompt:
        gsys_prompt = file.Get("gsyst_frac_prompt")
        gsys_fd = file.Get("gsyst_frac_fd")
    else:
        gsys_prompt = file.Get("gsyst_frac_prompt")
        gsys_fd = None

    file.Close()

    return gsys_prompt, gsys_fd

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
        print("ERROR: graph gvnSimFit does not exist! Exit")
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
    
# Method to compute prompt vn with all systematic uncertainties vn, fprompt, and fit
def ComputePromptVnWithSyst_bdt(outDir, inFileNameVn, inFileNamefPrompt, nonprompt, inFileNameFit, suffix):
    SetStyle()

    # Read v2 histograms, and convert them to TGraphAsymmErrors
    inFileVn = TFile.Open(inFileNameVn)
    hV2PromptStat  = inFileVn.Get("hV2VsPtPrompt")
    hV2PromptStat.SetDirectory(0)
    gV2PromptStat = TGraphAsymmErrors(hV2PromptStat)
    gV2PromptStat.SetName("gV2PromptStat")
    SetGraphStyle(gV2PromptStat, ROOT.kRed, ROOT.kFullCircle, markersize=1)
    if nonprompt:
        hV2FDStat = inFileVn.Get("hV2VsPtFD")
        hV2FDStat.SetDirectory(0)
        gV2FDStat = TGraphAsymmErrors(hV2FDStat)
        gV2FDStat.SetName("gV2FDStat")
        SetGraphStyle(gV2FDStat, ROOT.kRed, ROOT.kFullCircle, markersize=1)
    inFileVn.Close()
    nPtBins = gV2PromptStat.GetN()
    ptmin, ptmax = hV2PromptStat.GetXaxis().GetXmin(), hV2PromptStat.GetXaxis().GetXmax()
    if Debug.enable:
        print(f"Number of pT bins: {nPtBins}, pT range: {ptmin} - {ptmax}")

    if len(relSystUncFit_prompt) != nPtBins:
        print(f"ERROR: different number of pT bins in systematic-uncertainty array and input v2 graph {nPtBins, len(relSystUncFit_prompt)}! Exit")
        return

    # Read frac systematic uncertainties
    gsys_prompt, gsys_fd = ReadSysUncFile(inFileNamefPrompt, nonprompt)
    
    #______________________________________________________________________________________________________________________________________________________
    # graph with frac systematic uncertainties
    gfPrompt, gfFD = gV2PromptStat.Clone("gfPrompt"), gV2FDStat.Clone("gfFD")
    if Debug.enable:
        print('Debug: frac systematic uncertainties')
    for iPt in range(1, nPtBins):
        gfPrompt.SetPoint(iPt, gV2PromptStat.GetX()[iPt], gV2PromptStat.GetY()[iPt])
        #! absolute systematic uncertainty
        gfPrompt.SetPointError(iPt, gV2PromptStat.GetErrorXlow(iPt)*0.8, gV2PromptStat.GetErrorXhigh(iPt)*0.8,
                               gsys_prompt.GetErrorYlow(iPt), gsys_prompt.GetErrorYhigh(iPt))
        gfPrompt.SetName("gfPrompt")
        SetGraphStyle(gfPrompt, ROOT.kAzure+4, ROOT.kFullCircle, True, 3145)
        if Debug.enable:
            print(f"pt: {gfPrompt.GetX()[iPt]}, prompt v2: {gfPrompt.GetY()[iPt]}, error: {gfPrompt.GetErrorYlow(iPt)}")
        if nonprompt:
            gfFD.SetPoint(iPt, gV2FDStat.GetX()[iPt], gV2FDStat.GetY()[iPt])
            gfFD.SetPointError(iPt, gV2FDStat.GetErrorXlow(iPt)*0.8, gV2FDStat.GetErrorXhigh(iPt)*1.2,
                               gsys_fd.GetErrorYlow(iPt), gsys_fd.GetErrorYhigh(iPt))
            gfFD.SetName("gfFD")
            SetGraphStyle(gfFD, ROOT.kAzure+4, ROOT.kFullCircle, True, 3145)
            if Debug.enable:
                print(f"pt: {gfFD.GetX()[iPt]}, fd v2: {gfFD.GetY()[iPt]}, error: {gfFD.GetErrorYlow(iPt)}")
    if Debug.enable:
        input()

    #______________________________________________________________________________________________________________________________________________________
    # graph of fit systematic uncertainties
    gFitOptPrompt, gFitOptFD = gV2PromptStat.Clone("gFitOptPrompt"), gV2FDStat.Clone("gFitOptFD")
    if Debug.enable:
        print('Debug: fit systematic uncertainties')
    for iPt in range(nPtBins):
        gFitOptPrompt.SetPoint(iPt, gV2PromptStat.GetX()[iPt], gV2PromptStat.GetY()[iPt])
        #! relative systematic uncertainty
        gFitOptPrompt.SetPointError(iPt, gV2PromptStat.GetErrorXlow(iPt), gV2PromptStat.GetErrorXhigh(iPt),
                                relSystUncFit_prompt[iPt]*gV2PromptStat.GetY()[iPt], relSystUncFit_prompt[iPt]*gV2PromptStat.GetY()[iPt])
        gFitOptPrompt.SetName("gFitOptPrompt")
        # SetGraphStyle(gFitOptPrompt, ROOT.kGray + 2, ROOT.kFullCircle, True, 3154)
        if Debug.enable:
            print(f"pt: {gFitOptPrompt.GetX()[iPt]}, prompt v2: {gFitOptPrompt.GetY()[iPt]}, error: {gFitOptPrompt.GetErrorYlow(iPt)}")
        if nonprompt:
            gFitOptFD.SetPoint(iPt, gV2FDStat.GetX()[iPt], gV2FDStat.GetY()[iPt])
            gFitOptFD.SetPointError(iPt, gV2FDStat.GetErrorXlow(iPt), gV2FDStat.GetErrorXhigh(iPt),
                                relSystUncFit_fd[iPt]*gV2FDStat.GetY()[iPt], relSystUncFit_fd[iPt]*gV2FDStat.GetY()[iPt])
            gFitOptFD.SetName("gFitOptFD")
            # SetGraphStyle(gFitOptFD, ROOT.kGray + 2, ROOT.kFullCircle, True, 3154)
            if Debug.enable:
                print(f"pt: {gFitOptFD.GetX()[iPt]}, fd v2: {gFitOptFD.GetY()[iPt]}, error: {gFitOptFD.GetErrorYlow(iPt)}")
    if Debug.enable:
        input()
    
    # graph of resolution systematic uncertainties
    gResoPrompt, gResoFD = gV2PromptStat.Clone("gResoPrompt"), gV2FDStat.Clone("gResoFD")
    if Debug.enable:
        print('Debug: resolution systematic uncertainties')
    for iPt in range(nPtBins):
        gResoPrompt.SetPoint(iPt, gV2PromptStat.GetX()[iPt], gV2PromptStat.GetY()[iPt])
        #! relative systematic uncertainty
        gResoPrompt.SetPointError(iPt, gV2PromptStat.GetErrorXlow(iPt), gV2PromptStat.GetErrorXhigh(iPt),
                                relResolUnc*gV2PromptStat.GetY()[iPt], relResolUnc*gV2PromptStat.GetY()[iPt])
        gResoPrompt.SetName("gReso")
        # SetGraphStyle(gResoPrompt, ROOT.kBlack, ROOT.kFullCircle, True, 3154) #TODO: change color
        if Debug.enable:
            print(f"pt: {gResoPrompt.GetX()[iPt]}, prompt v2: {gResoPrompt.GetY()[iPt]}, error: {gResoPrompt.GetErrorYlow(iPt)}")
        if nonprompt:
            gResoFD.SetPoint(iPt, gV2FDStat.GetX()[iPt], gV2FDStat.GetY()[iPt])
            gResoFD.SetPointError(iPt, gV2FDStat.GetErrorXlow(iPt), gV2FDStat.GetErrorXhigh(iPt),
                                relResolUnc*gV2FDStat.GetY()[iPt], relResolUnc*gV2FDStat.GetY()[iPt])
            gResoFD.SetName("gResoFD")
            # SetGraphStyle(gResoFD, ROOT.kBlack, ROOT.kFullCircle, True, 3154)
            if Debug.enable:
                print(f"pt: {gResoFD.GetX()[iPt]}, fd v2: {gResoFD.GetY()[iPt]}, error: {gResoFD.GetErrorYlow(iPt)}")
    if Debug.enable:
        input()
            
    # graph of data systematic uncertainties
    gDataPrompt, gDataFD = gV2PromptStat.Clone("gDataPrompt"), gV2FDStat.Clone("gDataFD")
    if Debug.enable:
        print('Debug: data systematic uncertainties: resolution and fit')
    for iPt in range(nPtBins):
        gDataPrompt.SetPoint(iPt, gV2PromptStat.GetX()[iPt], gV2PromptStat.GetY()[iPt])
        #! data systematic uncertainty (fit and reso. corr.) prompt
        gDataPrompt.SetPointError(iPt, gV2PromptStat.GetErrorXlow(iPt)*0.5, gV2PromptStat.GetErrorXhigh(iPt)*0.5,
                              TMath.Sqrt(gFitOptPrompt.GetErrorYlow(iPt)*gFitOptPrompt.GetErrorYlow(iPt) + gResoPrompt.GetErrorYlow(iPt)*gResoPrompt.GetErrorYlow(iPt)),
                              TMath.Sqrt(gFitOptPrompt.GetErrorYhigh(iPt)*gFitOptPrompt.GetErrorYhigh(iPt) + gResoPrompt.GetErrorYhigh(iPt)*gResoPrompt.GetErrorYhigh(iPt)))
        gDataPrompt.SetName("gDataPrompt")
        SetGraphStyle(gDataPrompt, ROOT.kGray + 2, ROOT.kFullCircle, True, 3154)
        if Debug.enable:
            print(f"pt: {gDataPrompt.GetX()[iPt]}, prompt v2: {gDataPrompt.GetY()[iPt]}, error: {gDataPrompt.GetErrorYlow(iPt)}")
        if nonprompt:
            gDataFD.SetPoint(iPt, gV2FDStat.GetX()[iPt], gV2FDStat.GetY()[iPt])
            #! data systematic uncertainty (fit and reso. corr.) fd
            gDataFD.SetPointError(iPt, gV2FDStat.GetErrorXlow(iPt)*0.9, gV2FDStat.GetErrorXhigh(iPt)*0.9,
                              TMath.Sqrt(gFitOptFD.GetErrorYlow(iPt)*gFitOptFD.GetErrorYlow(iPt) + gResoFD.GetErrorYlow(iPt)*gResoFD.GetErrorYlow(iPt)),
                              TMath.Sqrt(gFitOptFD.GetErrorYhigh(iPt)*gFitOptFD.GetErrorYhigh(iPt) + gResoFD.GetErrorYhigh(iPt)*gResoFD.GetErrorYhigh(iPt)))
            gDataFD.SetName("gDataFD")
            SetGraphStyle(gDataFD, ROOT.kGray + 2, ROOT.kFullCircle, True, 3154)
            if Debug.enable:
                print(f"pt: {gDataFD.GetX()[iPt]}, fd v2: {gDataFD.GetY()[iPt]}, error: {gDataFD.GetErrorYlow(iPt)}")
    if Debug.enable:
        input()

    #______________________________________________________________________________________________________________________________________________________
    # graph of total systematic uncertainties
    gSystTotPrompt, gSystTotFD = gV2PromptStat.Clone("gSystTotPrompt"), gV2FDStat.Clone("gSystTotFD")
    for iPt in range(nPtBins):
        gSystTotPrompt.SetPoint(iPt, gV2PromptStat.GetX()[iPt], gV2PromptStat.GetY()[iPt])
        #! total systematic uncertainty prompt
        gSystTotPrompt.SetPointError(iPt, gV2PromptStat.GetErrorXlow(iPt)*0.6, gV2PromptStat.GetErrorXhigh(iPt)*0.6,
                                TMath.Sqrt(gDataPrompt.GetErrorYlow(iPt)*gDataPrompt.GetErrorYlow(iPt) + gfPrompt.GetErrorYlow(iPt)*gfPrompt.GetErrorYlow(iPt)),
                                TMath.Sqrt(gDataPrompt.GetErrorYhigh(iPt)*gDataPrompt.GetErrorYhigh(iPt) + gfPrompt.GetErrorYhigh(iPt)*gfPrompt.GetErrorYhigh(iPt)))
        gSystTotPrompt.SetName("gSystTotPrompt")
        SetGraphStyle(gSystTotPrompt, ROOT.kBlack, ROOT.kFullCircle) #TODO: change color
        if nonprompt:
            gSystTotFD.SetPoint(iPt, gV2FDStat.GetX()[iPt], gV2FDStat.GetY()[iPt])
            #! total systematic uncertainty fd
            gSystTotFD.SetPointError(iPt, gV2FDStat.GetErrorXlow(iPt)*0.6, gV2FDStat.GetErrorXhigh(iPt)*0.6,
                                TMath.Sqrt(gDataPrompt.GetErrorYlow(iPt)*gDataPrompt.GetErrorYlow(iPt) + gfFD.GetErrorYlow(iPt)*gfFD.GetErrorYlow(iPt)),
                                TMath.Sqrt(gDataPrompt.GetErrorYhigh(iPt)*gDataPrompt.GetErrorYhigh(iPt) + gfFD.GetErrorYhigh(iPt)*gfFD.GetErrorYhigh(iPt)))
            gSystTotFD.SetName("gSystTotFD")
            SetGraphStyle(gSystTotFD, ROOT.kBlack, ROOT.kFullCircle) #TODO: change color

    outFileName = f"{outDir}/promptvn_withsyst{suffix}.root"
    outFile = TFile(outFileName, "recreate")
    outFile.cd()
    
    gV2PromptStat.Write()
    gfPrompt.Write()
    gFitOptPrompt.Write()
    gResoPrompt.Write()
    gDataPrompt.Write()
    gSystTotPrompt.Write()
    if nonprompt:
        gV2FDStat.Write()
        gfFD.Write()
        gFitOptFD.Write()
        gResoFD.Write()
        gDataFD.Write()
        gSystTotFD.Write()

    c = TCanvas("c", "c", 800, 600)
    c.cd()
    haxis = TH1F("", "; #it{p}_{T} (GeV/#it{c});Prompt #it{v}_{2}", 100, ptmin, ptmax)
    haxis.GetYaxis().SetRangeUser(-0.002, 0.2)
    haxis.Draw("axis")
    gfPrompt.Draw("e2 same")
    gDataPrompt.Draw("e2 same")
    gSystTotPrompt.Draw("e2 same")
    gV2PromptStat.Draw("ep same")
    line = TLine(ptmin, 0, ptmax, 0)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("same")
    legend = TLegend(0.53, 0.63, 0.83, 0.83)
    legend.SetTextSize(0.04)
    legend.SetTextFont(42)
    legend.AddEntry(gV2PromptStat, "prompt", "ep")
    legend.AddEntry(gfPrompt, "fraction", "f")
    legend.AddEntry(gDataPrompt, "syst. unc. (fit and reso. corr.)", "f")
    legend.AddEntry(gSystTotPrompt, "syst. unc. (total)", "f")
    legend.Draw("same")
    c.Update()
    input()
    c.Write()
    
    if nonprompt:
        cFD = TCanvas("cFD", "cFD", 800, 600)
        cFD.cd()
        haxisFD = TH1F("haxisFD", "; #it{p}_{T} (GeV/#it{c});Non-prompt #it{v}_{2}", 100, ptmin, ptmax)
        haxisFD.GetYaxis().SetRangeUser(-0.2, 0.3)
        haxisFD.Draw("axis")
        gfFD.Draw("e2 same")
        gDataFD.Draw("e2 same")
        gSystTotFD.Draw("e2 same")
        gV2FDStat.Draw("ep same")
        lineFD = TLine(ptmin, 0, ptmax, 0)
        lineFD.SetLineStyle(2)
        lineFD.SetLineWidth(2)
        lineFD.Draw("same")
        legendFD = TLegend(0.53, 0.63, 0.83, 0.83)
        legendFD.SetTextSize(0.04)
        legendFD.SetTextFont(42)
        legendFD.AddEntry(gV2FDStat, "non-prompt", "ep")
        legendFD.AddEntry(gfFD, "fraction", "f")
        legendFD.AddEntry(gDataFD, "syst. unc. (fit and reso. corr.)", "f")
        legendFD.AddEntry(gSystTotFD, "syst. unc. (total)", "f")
        legendFD.Draw("same")
        cFD.Update()
        input()
        cFD.Write()

    outFile.Close()

    print(f"Output written to: {outFileName}")
    input("Press Enter to exit.")

    pass
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("vn_file", metavar="text",
                        default="path/to/v2VSFrac.root", help="central vn file")
    parser.add_argument("frac_file", metavar="text",
                        default="path/to/syst_fFD_suffix.root", help="frac file")
    parser.add_argument("fit_file", metavar="text",
                        default="path/to/syst_prompt.root", help="fit file")
    parser.add_argument("--nonprompt", "-np", action="store_true",
                        help="compute non-prompt vn")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    Debug.enable = False

    if not args.fit_file:
        ComputePromptVnWithSyst(args.outputdir,
                                args.vn_file,
                                args.frac_file,
                                args.suffix)
    else:
        ComputePromptVnWithSyst_bdt(args.outputdir,
                                    args.vn_file,
                                    args.frac_file,
                                    args.nonprompt,
                                    args.fit_file,
                                    args.suffix)
                                