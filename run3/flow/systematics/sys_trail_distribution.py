import os
import sys
import yaml
import argparse
sys.path.append('/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT')
sys.path.append('/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/systematics')
import ROOT
from ROOT import TFile, TCanvas, TLegend, TLatex, TGraphErrors, TF1, TH1D, TVirtualFitter, Double_t
from ROOT import kBlack, kAzure, kCyan, kOrange
from ROOT import kFullCircle, kOpenCircle
sys.path.append('../../../')
sys.path.append('../')
from flow_analysis_utils import get_particle_info, get_cut_sets_config
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor

#def v2_vs_frac(config_flow, inputdir, outputdir, suffix, fracFiles, v2Files, systematics=False):

ROOT.gROOT.SetBatch(True)

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

def get_mean_vn_rms_vn(hvn, hsyst):
    mean_vn = 0
    rms_vn = 0
    for i in range(hvn.GetNbinsX()):
        mean_vn = hvn.GetBinContent(i) / hvn.GetNbinsX() + mean_vn
    rms_vn = hsyst.GetMean()
    return mean_vn, rms_vn
config_flow_path = '/home/wuct/ALICE/local/Results/test/flow/systematics/trails/trails_pt_30_35/cutvar_config_pt3_3.5_Sigma-0.018930816044528502_BkgFuncVn-kLin_BkgFunc-kExpo_Rebin-1_MassMin-1.7_MassMax-2.0/config_flow/config_flow_config_pt3_3.5_Sigma-0.018930816044528502_BkgFuncVn-kLin_BkgFunc-kExpo_Rebin-1_MassMin-1.7_MassMax-2.0_0.yml'

ry_template_path = '/home/wuct/ALICE/local/Results/test/flow/systematics/trails/trails_pt_30_35/cutvar_config_pt3_3.5_Sigma-0.018930816044528502_BkgFuncVn-kLin_BkgFunc-kExpo_Rebin-1_MassMin-1.7_MassMax-2.0/ry/raw_yields_config_pt3_3.5_Sigma-0.018930816044528502_BkgFuncVn-kLin_BkgFunc-kExpo_Rebin-1_MassMin-1.7_MassMax-2.0_00.root'

ry_template_file = TFile.Open(ry_template_path)

gry = ry_template_file.Get('gvnSimFit')
hry = ry_template_file.Get('hvnSimFit')
hry.SetDirectory(0)
gry.Set(0)
# hry.Set(0)
ry_template_file.Close()

syst_mult_path = '/home/wuct/ALICE/local/Results/test/flow/systematics/sys/pt_30_35/syst_multitrial'
outputdir = '/home/wuct/ALICE/local/Results/test/flow/systematics/sys/pt_30_35/syst_multitrial'
syst_file_names = os.listdir(syst_mult_path)
for file in syst_file_names:
    if 'syst_multitrial' not in file or 'ry' not in file:
        syst_file_names.remove(file)

print(syst_file_names)
syst_file_names.sort()
syst_files = []
mean_v2s = []
rms_v2s = []
syst_finals = []

for syst_file_name in syst_file_names:
    syst_file = TFile.Open(os.path.join(syst_mult_path, syst_file_name))
    print(f'Loading {syst_file_name}')
    syst_files.append(syst_file)
    hvn = syst_file.Get('hvn')
    hsyst = syst_file.Get('hsyst')
    hsyst_final = syst_file.Get('hsyst_final')
    hvn.SetDirectory(0)
    hsyst.SetDirectory(0)
    hsyst_final.SetDirectory(0)
    mean_v2, rms_v2 = get_mean_vn_rms_vn(hvn, hsyst)
    print(f'mean_v2: {mean_v2}, rms_v2: {rms_v2}')
    mean_v2s.append(mean_v2)
    rms_v2s.append(rms_v2)
    syst_final = hsyst_final.GetBinContent(1)
    syst_finals.append(syst_final)
    print(f'syst_final: {syst_final}')
    syst_file.Close()

from itertools import product

n_files = len(mean_v2s)

combinations = product([0, 1, 2], repeat=n_files)

all_gry_configs = []
all_hry_configs = []

for icombo, combo in enumerate(combinations):
    gry_sys = []
    hry_sys = []
    
    for iFile, choice in enumerate(combo):
        gry_sy = gry.Clone(f"gvnSimFit_combo_{iFile}_choice_{choice}")
        hry_sy = hry.Clone(f"hvnSimFit_combo_{iFile}_choice_{choice}")
        
        if choice == 0:
            value = mean_v2s[iFile]
        elif choice == 1:
            value = mean_v2s[iFile] - rms_v2s[iFile]
        elif choice == 2:
            value = mean_v2s[iFile] + rms_v2s[iFile]
        
        gry_sy.SetPoint(0, 0, value)
        gry_sy.SetPointError(0, 0.25, 0.25, syst_finals[iFile], syst_finals[iFile])

        hry_sy.SetBinContent(1, value)
        hry_sy.SetBinError(1, syst_finals[iFile])
        hry_sy.SetDirectory(0)
        
        gry_sys.append(gry_sy)
        hry_sys.append(hry_sy)

# for ll in range(3):
#     gry_sys, hry_sys = [], []
#     for iFile in range(len(mean_v2s)):
#         gry_sy = gry.Clone('gvnSimFit')
#         if ll == 0:
#             gry_sy.SetPoint(0, 0, mean_v2s[iFile])
#         elif ll == 1:
#             gry_sy.SetPoint(0, 0, mean_v2s[iFile] - rms_v2s[iFile])
#         elif ll == 2:
#             gry_sy.SetPoint(0, 0, mean_v2s[iFile] + rms_v2s[iFile])
#         gry_sy.SetPointError(0, 0.25, 0.25, syst_finals[iFile], syst_finals[iFile])
#         hry_sy = hry.Clone('hvnSimFit')
#         hry_sy.SetBinContent(1, mean_v2s[iFile])
#         hry_sy.SetBinError(1, syst_finals[iFile])
#         hry_sy.SetDirectory(0)
#         gry_sys.append(gry_sy)
#         hry_sys.append(hry_sy)

    CutSets = len(syst_files)
    with open(config_flow_path, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)
        
    ptmins = config['ptmins']
    ptmaxs = config['ptmaxs']

    particleName = config['Dmeson']

    particleTit, _, decay, _ = get_particle_info(particleName)

    fracFiles = load_frac_files('/home/wuct/ALICE/local/Results/test/flow/systematics/trails/trails_pt_30_35/cutvar_config_pt3_3.5_Sigma-0.018930816044528502_BkgFuncVn-kLin_BkgFunc-kExpo_Rebin-1_MassMin-1.7_MassMax-2.0', 'DataDrivenFrac')


    hV2, gV2, hFracFD, hFracPrompt = [], [], [], []
    avrV2XErrL, avrV2XErrH = [], []

    for iFile in range(len(fracFiles)):
        
        fracFile = fracFiles[iFile]
        inFracFile = TFile.Open(fracFile)
        hFracFD.append(inFracFile.Get('hFDFrac'))
        hFracPrompt.append(inFracFile.Get('hPromptFrac'))
        hFracFD[-1].SetDirectory(0)
        hFracPrompt[-1].SetDirectory(0)

    gFracVsV2, hV2VsFrac = [], [] # gFracVsV2 used for fitting, hV2VsFrac used for plotting
    hV2VsPtFD = hry.Clone("hV2VsPtFD")
    hV2VsPtPrompt = hry.Clone("hV2VsPtPrompt")

    cFrac, ptStrings, chi2Strings = [], [], []

    nPtBins = len(ptmins)
    for iPt, (ptMin, ptMax) in enumerate(zip(ptmins, ptmaxs)):
        ptCent = (ptMin + ptMax) / 2
        nSets = CutSets

        gFracVsV2.append(TGraphErrors(-1))
        hV2VsFrac.append(TH1D(f"hV2VsFrac_{iPt}", "", 1000, 0.0, 1.0))
        hV2VsFrac[-1].SetDirectory(0)
        SetObjectStyle(hV2VsFrac[-1], markerstyle=kFullCircle, markersize=0)
        SetObjectStyle(gFracVsV2[-1], linecolor=kAzure+4, linewidth=2, markerstyle=kFullCircle, markersize=1, markercolor=kAzure+4)

        avrV2XErrL.append(Double_t(sum(gry_sys[i].GetErrorXlow(iPt) for i in range(nSets)) / nSets))
        avrV2XErrH.append(Double_t(sum(gry_sys[i].GetErrorXhigh(iPt) for i in range(nSets)) / nSets))
        

        v2Values = [hry_sys[i].GetBinContent(iPt + 1) for i in range(nSets)]
        v2Unc = [hry_sys[i].GetBinError(iPt + 1) for i in range(nSets)]
        fracBins = [hFracFD[i].GetXaxis().FindBin(ptCent) for i in range(nSets)]
        fracFDValues = [hFracFD[i].GetBinContent(fracBins[i]) for i in range(nSets)]
        fracFDUnc = [hFracFD[i].GetBinError(fracBins[i]) for i in range(nSets)]

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
    outFile = TFile(f'{outputdir}/V2VsFrac/V2VsFrac_sys_trail_distribution_{icombo}.root', "recreate")

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

        cFrac[iPt].SaveAs(f"{outputdir}/V2VsFrac/FracV2{icombo}.pdf{suffix_pdf}")

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
    cV2VsPtFD.SaveAs(f"{outputdir}/V2VsFrac/V2VsPtFD{icombo}.pdf")
    cV2VsPtPrompt.SaveAs(f"{outputdir}/V2VsFrac/V2VsPtPrompt{icombo}.pdf")
    cPromptAndFDV2.SaveAs(f"{outputdir}/V2VsFrac/V2VsPtPromptAndFD{icombo}.pdf")
    print(syst_files)