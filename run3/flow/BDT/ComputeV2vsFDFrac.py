"""
python script for the computation of the prompt or non-prompt v2 via extrapolation
run: python ComputeV2vsFDFrac.py config.yaml --inputdir path/to/input --outputdir path/to/output --suffix text

"""
import argparse
import os
import yaml
import sys
import ROOT
from ROOT import TFile, TCanvas, TLegend, TLatex, TGraphErrors, TF1, TH1D, TVirtualFitter, Double_t
from ROOT import kBlack, kAzure, kCyan, kOrange
from ROOT import kFullCircle, kOpenCircle
sys.path.append('../../../')
sys.path.append('../')
from flow_analysis_utils import get_particle_info, get_cut_sets_config
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor

def load_v2_files(inputdir, suffix):
    if os.path.exists(f'{inputdir}/ry'):
        v2Files = [f'{inputdir}/ry/{file}'
                    for file in os.listdir(f'{inputdir}/ry') if file.endswith('.root') and suffix in file]
        v2Files.sort()
    else:
        raise ValueError(f'No ry folder found in {inputdir}')
    return v2Files

def load_frac_files(inputdir, suffix):
    if os.path.exists(f'{inputdir}/DataDrivenFrac'):
        fracFiles = [f'{inputdir}/DataDrivenFrac/{file}'
                        for file in os.listdir(f'{inputdir}/DataDrivenFrac') if file.endswith('.root') and suffix in file]
        fracFiles.sort()
    else:
        raise ValueError(f'No DataDrivenFrac folder found in {inputdir}')
    return fracFiles

def set_frame_style(canv, Title, particleTit):
    canv.SetLeftMargin(0.15)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.15)
    canv.SetTopMargin(0.05)
    hFrame = canv.DrawFrame(0.0, -0.2, 1, 0.35, f";Non-prompt {particleTit} fraction; #it{{v}}_{{2}}^{{#it{{obs}}}}")
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

def set_frame_margin(canv):
    canv.SetLeftMargin(0.15)
    canv.SetRightMargin(0.05)
    canv.SetBottomMargin(0.15)
    canv.SetTopMargin(0.05)

def v2_vs_frac(config_flow, inputdir, outputdir, suffix, fracFiles, v2Files):

    CutSets, _, _, _, _ = get_cut_sets_config(config_flow)
    nCutSets = max(CutSets)
    with open(config_flow, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)
        
    ptmins = config['ptmins']
    ptmaxs = config['ptmaxs']

    particleName = config['Dmeson']

    particleTit, _, decay, _ = get_particle_info(particleName)
    
    CutSets, _, _, _, _ = get_cut_sets_config(config_flow)

    if len(fracFiles) != len(v2Files):
        raise ValueError('Number of eff and frac files do not match')

    hV2, gV2, hFracFD, hFracPrompt = [], [], [], []
    avrV2XErrL, avrV2XErrH = [], []

    for fracFile, v2File in zip(fracFiles, v2Files):
        inV2File = TFile.Open(v2File)
        hV2.append(inV2File.Get('hvnSimFit'))
        gV2.append(inV2File.Get('gvnSimFit'))
        hV2[-1].SetDirectory(0)

        inFracFile = TFile.Open(fracFile)
        hFracFD.append(inFracFile.Get('hFDFrac'))
        hFracPrompt.append(inFracFile.Get('hPromptFrac'))
        hFracFD[-1].SetDirectory(0)
        hFracPrompt[-1].SetDirectory(0)

    gFracVsV2, hV2VsFrac = [], [] # gFracVsV2 used for fitting, hV2VsFrac used for plotting
    hV2VsPtFD = hV2[0].Clone("hV2VsPtFD")
    hV2VsPtPrompt = hV2[0].Clone("hV2VsPtPrompt")

    cFrac, ptStrings, chi2Strings = [], [], []
    
    nPtBins = len(ptmins)
    for iPt, (ptMin, ptMax) in enumerate(zip(ptmins, ptmaxs)):
        ptCent = (ptMin + ptMax) / 2
        nSets = CutSets[iPt]

        gFracVsV2.append(TGraphErrors(-1))
        hV2VsFrac.append(TH1D(f"hV2VsFrac_{iPt}", "", 1000, 0.0, 1.0))
        hV2VsFrac[-1].SetDirectory(0)
        SetObjectStyle(hV2VsFrac[-1], markerstyle=kFullCircle, markersize=0)
        SetObjectStyle(gFracVsV2[-1], linecolor=kAzure+4, linewidth=2, markerstyle=kFullCircle, markersize=1, markercolor=kAzure+4)

        avrV2XErrL.append(Double_t(sum(gV2[i].GetErrorXlow(iPt) for i in range(nSets)) / nSets))
        avrV2XErrH.append(Double_t(sum(gV2[i].GetErrorXhigh(iPt) for i in range(nSets)) / nSets))
        

        v2Values = [hV2[i].GetBinContent(iPt + 1) for i in range(nSets)]
        v2Unc = [hV2[i].GetBinError(iPt + 1) for i in range(nSets)]
        fracFDValues = [hFracFD[i].GetBinContent(iPt + 1) for i in range(nSets)]
        fracFDUnc = [hFracFD[i].GetBinError(iPt + 1) for i in range(nSets)]

        for iSet, (v2, fracFD, v2Unc, fracFDUnc) in enumerate(zip(v2Values, fracFDValues, v2Unc, fracFDUnc)):
            print(f"pt: {ptCent:.4f}, v2: {v2:.4f}, fracFD: {fracFD:.4f}")
            gFracVsV2[iPt].SetPoint(iSet, fracFD, v2)
            gFracVsV2[iPt].SetPointError(iSet, fracFDUnc, v2Unc)
        
        # gFracVsV2Fit = TGraphErrors(gFracVsV2[-1])
        linFunc = TF1("linear", "pol1", 0, 1)
        SetObjectStyle(linFunc, color=kOrange+1, linestyle=9, linewidth=2)
        gFracVsV2[-1].Fit("linear", "", "", 0, 1)
        chi2 = linFunc.GetChisquare()
        ndf = linFunc.GetNDF()

        # get the confidence intervals 0.683
        fitter = TVirtualFitter.GetFitter()
        fitter.GetConfidenceIntervals(hV2VsFrac[-1], 0.683)
        hV2VsFrac[-1].SetLineColorAlpha(kAzure+5, 0.15)

        # get the v2 value at the FD fraction = 1, and it is not the last bin?
        hV2VsPtFD.SetBinContent(iPt + 1, 
                                hV2VsFrac[-1].GetBinContent(hV2VsFrac[-1].GetNbinsX()))
        hV2VsPtFD.SetBinError(iPt + 1,
                                hV2VsFrac[-1].GetBinError(hV2VsFrac[-1].GetNbinsX()))
        
        # get the v2 value at the FD fraction = 0, and it is the first bin?
        hV2VsPtPrompt.SetBinContent(iPt + 1, 
                                    hV2VsFrac[-1].GetBinContent(hV2VsFrac[-1].GetBin(1)))
        hV2VsPtPrompt.SetBinError(iPt + 1,
                                    hV2VsFrac[-1].GetBinError(hV2VsFrac[-1].GetBin(1)))
        
        #TODO: plot the v2 vs pt, and the center of the pt bin is calculate by the average of pT

        ptStrings.append(f"{ptMin:.1f} < #it{{p}}_{{T}} < {ptMax:.1f} GeV/#it{{c}}")
        chi2Strings.append(f"#chi^{{2}}/n.d.f = {chi2:.2f}/{ndf:.2f}")


    # save the results
    os.makedirs(outputdir + f'/V2VsFrac', exist_ok=True)
    outFile = TFile(f'{outputdir}/V2VsFrac/V2VsFrac_{suffix}.root', "recreate")
    
    t = TLatex(8, 8, "")
    t.SetNDC()
    t.SetTextFont(42)
    t.SetTextColor(kBlack)

    for iPt, (ptMin, ptMax) in enumerate(zip(ptmins, ptmaxs)):
        if iPt == 0:
            suffix_pdf = '('
        elif iPt == nPtBins-1:
            suffix_pdf = ')'
        else:
            suffix_pdf = ''
        if nPtBins == 1:
            suffix_pdf = ''

        cFrac.append(TCanvas(f"cFrac_{ptMin}_{ptMax}", "", 1200, 1200))
        set_frame_style(cFrac[-1], ptStrings[iPt], particleTit)

        t.SetTextSize(0.04)
        t.DrawLatex(0.25, 0.85, decay)
        t.DrawLatex(0.25, 0.78, f"{ptStrings[iPt]}")
        t.SetTextSize(0.035)
        t.DrawLatex(0.250, 0.23, f'{chi2Strings[iPt]}')

        hV2VsFrac[iPt].Draw("same pZ")
        gFracVsV2[iPt].Draw("same pZ")

        cFrac[-1].Update()
        cFrac[-1].Write()

        cFrac[iPt].SaveAs(f"{outputdir}/V2VsFrac/FracV2_{suffix}.pdf{suffix_pdf}")

        outFile.mkdir(f"pt_{int(ptMin*10)}_{int(ptMax*10)}")
        outFile.cd(f"pt_{int(ptMin*10)}_{int(ptMax*10)}")
        gFracVsV2[iPt].Write('gV2VsFrac')
        hV2VsFrac[iPt].Write('hV2VsFrac')

    outFile.cd()
    PtTit = "#it{p}_{T} GeV/#it{c}"
    leg = TLegend(0.55, 0.75, 0.88, 0.89)
    leg.SetTextSize(0.045)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    SetObjectStyle(hV2VsPtFD, color=GetROOTColor("kAzure+4"), fillstyle=1)
    SetObjectStyle(hV2VsPtPrompt, color=GetROOTColor("kRed+1"), fillstyle=1)

    cV2VsPtFD = TCanvas("cV2VsPtFD", "non-prompt v2 versus pt", 800, 800)
    set_frame_margin(cV2VsPtFD)    
    cV2VsPtFD.cd()
    hV2VsPtFD.Draw("")
    hV2VsPtFD.GetXaxis().SetTitle(PtTit)
    hV2VsPtFD.GetYaxis().SetTitle("Non-prompt #it{v_{2}}")
    hV2VsPtFD.GetYaxis().SetRangeUser(-0.05, 0.35)
    hV2VsPtFD.SetMarkerStyle(20)
    hV2VsPtFD.SetMarkerSize(2)
    hV2VsPtFD.GetYaxis().SetNoExponent()

    cV2VsPtPrompt = TCanvas("cV2VsPtPrompt", "prompt v2 versus pt", 800, 800)
    set_frame_margin(cV2VsPtPrompt)
    cV2VsPtPrompt.cd()
    hV2VsPtPrompt.Draw("")
    hV2VsPtPrompt.GetXaxis().SetTitle(PtTit)
    hV2VsPtPrompt.GetYaxis().SetTitle("Prompt #it{v_{2}}")
    hV2VsPtPrompt.GetYaxis().SetRangeUser(-0.05, 0.35)
    hV2VsPtPrompt.SetMarkerStyle(20)
    hV2VsPtPrompt.SetMarkerSize(2)
    hV2VsPtPrompt.GetYaxis().SetNoExponent()

    cPromptAndFDV2 = TCanvas("cPromptAndFDV2", "prompt and non-prompt v2 versus pt", 800, 800)
    set_frame_margin(cPromptAndFDV2)
    cPromptAndFDV2.cd()
    hV2VsPtFD.GetYaxis().SetTitle("#it{v_{2}}")
    hV2VsPtFD.Draw("")
    hV2VsPtPrompt.Draw("same")

    leg.AddEntry(hV2VsPtFD, "Non-prompt #it{v_{2}}", "lp")
    leg.AddEntry(hV2VsPtPrompt, "Prompt #it{v_{2}}", "lp")
    leg.Draw("same")

    hV2VsPtFD.Write()
    hV2VsPtPrompt.Write()
    cV2VsPtFD.SaveAs(f"{outputdir}/V2VsFrac/V2VsPtFD_{suffix}.pdf")
    cV2VsPtPrompt.SaveAs(f"{outputdir}/V2VsFrac/V2VsPtPrompt_{suffix}.pdf")
    cPromptAndFDV2.SaveAs(f"{outputdir}/V2VsFrac/V2VsPtPromptAndFD_{suffix}.pdf")

def main_v2_vs_frac(config, inputdir, outputdir, suffix, combined=False, inputdir_combined='', outputdir_combined=''):

    v2_vs_frac(
        config,
        inputdir,
        outputdir,
        suffix,
        load_frac_files(inputdir, suffix),
        load_v2_files(inputdir, suffix)
    )
    if combined:
        v2_vs_frac(
            config,
            inputdir_combined,
            outputdir_combined,
            suffix,
            load_frac_files(inputdir_combined, suffix),
            load_v2_files(inputdir, suffix)
        )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="flow configuration file")
    parser.add_argument('--inputdir', '-i', metavar='text',
                        default='.', help='input directory containing rawyields and frac files')
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--combined", '-comb', metavar='bool', required=False,
                        action='store_true', help="combined method")
    parser.add_argument("--inputdir_combined", "-ic", metavar="text", required=False,
                        default="", help="input directory containing the frac files for the combined method")
    parser.add_argument("--output_combined", "-oc", metavar="text", required=False,
                        default="", help="output directory for the combined method")
    args = parser.parse_args()

    main_v2_vs_frac(
        args.config,
        args.inputdir,
        args.outputdir,
        args.suffix,
        combined=args.combined if args.combined else False,
        inputdir_combined=args.inputdir_combined if args.combined else '',
        outputdir_combined=args.output_combined if args.combined else ''
    )