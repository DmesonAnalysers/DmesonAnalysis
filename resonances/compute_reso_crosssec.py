'''
python script to compute the integrated reco efficiency for resonances
'''

import os
import argparse
import sys
import yaml
import numpy as np
from ROOT import TFile, TH1F, TCanvas, kRed, kAzure, TLegend, TGraphAsymmErrors, gROOT  # pylint: disable=import-error,no-name-in-module
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ComputeCrossSection #pylint: disable=wrong-import-position,import-error

def compute_crosssec(file_rawy, file_eff, file_frac, outputdir, suffix):
    """
    Function to compute cross section and corrected yield
    """

    SetGlobalStyle(padtopmargin=0.05, padleftmargin=0.18, padbottommargin=0.15, palette=55,
                   labelsize=0.04, titlesize=0.05, labeloffset=0.008, titleoffsety=1.7,
                   titleoffsetx=1.2, titleoffsetz=1., opttitle=0, optstat=0)

    gROOT.SetBatch(True)

    delta_y = 1.
    sigmaMB = 57.95e+3 # in mb
    eff_trigger = 0.92 # to be considered when computing the yields

    if 'Ds1plus' in file_rawy:
        reso = 10433
        reso_label = 'Ds1plus'
        d_meson = 'Dstar'
        v0 = 'K0S'
        reso_label_plot = 'D_{s1}^{+}'
    elif 'Ds2starplus' in file_rawy:
        reso = 435
        reso_label = 'Ds2starplus'
        d_meson = 'Dplus'
        v0 = 'K0S'
        reso_label_plot = 'D_{s2}^{*+}'
    else:
        raise ValueError('Resonance not supported')

    ptshape = "Dsptshape"
    if "Lcptshape" in file_eff:
        ptshape = "Lcptshape"
    elif "DsHarderptshape" in file_eff:
        ptshape = "DsHarderptshape"
    elif "DsVeryHarderptshape" in file_eff:
        ptshape = "DsVeryHarderptshape"
    elif "DsSofterptshape" in file_eff:
        ptshape = "DsSofterptshape"
    elif "DsVerySofterptshape" in file_eff:
        ptshape = "DsVerySofterptshape"

    mult_weights_suffix = ""
    if "_multweights_all" in file_eff:
        mult_weights_suffix = "_multweights_all"
    elif "_multweights_candinmass" in file_eff:
        mult_weights_suffix = "_multweights_candinmass"
    elif "_multweights_cand" in file_eff:
        mult_weights_suffix = "_multweights_cand"

    trigger = ""
    if "MB" in file_rawy:
        trigger = "MB"
    elif "HM" in file_rawy:
        trigger = "HM"

    #_____________________________________________________________________________
    # load BR and syst information
    with open('config_syst_list.yml', 'r', encoding='utf-8') as ymlCfgFile:
        inputSyst = yaml.load(ymlCfgFile, yaml.FullLoader)

    BR_K0s = inputSyst['K0s']['BR']
    BR_K0s_unc = inputSyst['K0s']['BR_unc']
    BR_meson = inputSyst[d_meson]['BR']
    BR_meson_unc = inputSyst[d_meson]['BR_unc']

    BR = BR_K0s * BR_meson
    BR_unc = np.sqrt((BR_K0s_unc / BR_K0s)**2 + (BR_meson_unc / BR_meson)**2) * BR

    rawYieldSyst = inputSyst[reso][trigger]['raw_yield']
    effSyst = inputSyst[reso][trigger]['eff']
    ptShapeSyst = inputSyst[reso][trigger]['pT_shape']
    multSyst = inputSyst[reso][trigger]['mult']
    trackSyst = inputSyst[reso][trigger]['tracking']

    totSimmSystUnc = np.sqrt(rawYieldSyst**2 + effSyst**2 + ptShapeSyst**2 + multSyst**2 + trackSyst**2)

    #_____________________________________________________________________________
    # load input files
    # raw yields
    rawYieldFile = TFile.Open(file_rawy)
    hRawYields = rawYieldFile.Get('h_rawyields')
    rawYield = hRawYields.GetBinContent(1)
    rawYieldUnc = hRawYields.GetBinError(1)
    delta_pt = hRawYields.GetBinWidth(1)
    pt_min = hRawYields.GetXaxis().GetBinLowEdge(1)
    pt_max = hRawYields.GetXaxis().GetBinUpEdge(1)
    pt_cent = pt_min + delta_pt / 2

    # nevents
    input_dir = file_rawy.split("mass_")[0]
    file_norm = os.path.join(input_dir, f"normalisation_{d_meson}_{v0}_{trigger}.root")
    normFile = TFile.Open(file_norm)
    hnev = normFile.Get('hist_events')
    nev_sigma = hnev.GetBinContent(1)
    nev_yield = hnev.GetBinContent(2)

    # eff. acc.
    effAccFile = TFile.Open(file_eff)
    hEffAccPrompt = effAccFile.Get('heff_prompt_ptint')
    hEffAccNonPrompt = effAccFile.Get('heff_nonprompt_ptint')
    effAcc_prompt = hEffAccPrompt.GetBinContent(1)
    effAcc_nonprompt = hEffAccNonPrompt.GetBinContent(1)
    effAccUncPrompt = hEffAccPrompt.GetBinError(1)
    effAccUncNonPrompt = hEffAccNonPrompt.GetBinError(1)

    # fraction of prompt and FD
    fracFile = TFile.Open(file_frac)
    gFracPrompt = fracFile.Get('graph_pfrac_dreso_stat_ptint')
    gFracNonPrompt = fracFile.Get('graph_npfrac_dreso_stat_ptint')
    gFracPromptSyst = fracFile.Get('graph_pfrac_dreso_hypo_ptint')
    gFracNonPromptSyst = fracFile.Get('graph_npfrac_dreso_hypo_ptint')

    frac_prompt = gFracPrompt.GetY()[0]
    frac_nonprompt = gFracNonPrompt.GetY()[0]
    frac_prompt_statunc = gFracPrompt.GetErrorY(0)
    frac_nonprompt_statunc = gFracNonPrompt.GetErrorY(0)
    frac_prompt_systunc_low= gFracPromptSyst.GetErrorYlow(0)
    frac_prompt_systunc_high = gFracPromptSyst.GetErrorYhigh(0)
    frac_nonprompt_systunc_low = gFracNonPromptSyst.GetErrorYlow(0)
    frac_nonprompt_systunc_high = gFracNonPromptSyst.GetErrorYhigh(0)

    # sigma_MB
    hNormIngr = TH1F('hNormIngredients', 'hNormIngredients', 10, 0, 10)
    hNormIngr.SetBinContent(1, sigmaMB)
    hNormIngr.GetXaxis().SetBinLabel(1, 'sigma_MB')
    hNormIngr.SetBinContent(2, eff_trigger)
    hNormIngr.GetXaxis().SetBinLabel(2, 'eff_trigger (num)')
    hNormIngr.SetBinContent(3, BR)
    hNormIngr.GetXaxis().SetBinLabel(3, 'BR')
    hNormIngr.SetBinContent(4, BR_unc)
    hNormIngr.GetXaxis().SetBinLabel(4, 'BR_unc')
    hNormIngr.SetBinContent(5, nev_sigma)
    hNormIngr.GetXaxis().SetBinLabel(5, 'nev_sigma')
    hNormIngr.SetBinContent(6, nev_yield)
    hNormIngr.GetXaxis().SetBinLabel(6, 'nev_yield')

    print('\033[1m\033[92mIngredients\033[0m')
    print(f'resonance: {reso}')
    print(f'trigger: {trigger}')
    print(f'delta pt: {delta_pt}')
    print(f'delta y: {delta_y}')
    print(f'rawyield: {rawYield:.2f}, rawyiledunc: {rawYieldUnc:.2f}')
    print(f'nev_sigma: {nev_sigma}, nev_yield: {nev_yield}')
    print(f'effAcc_prompt: {effAcc_prompt:.2e}, effAccUncPrompt: {effAccUncPrompt:.2e}')
    print(f'effAcc_nonprompt: {effAcc_nonprompt:.2e}, effAccUncNonPrompt: {effAccUncNonPrompt:.2e}')
    print(f'eff_trigger: {eff_trigger}')
    print(f'frac_prompt: {frac_prompt:.2e} +- {frac_prompt_statunc:.2e} (stat) - {frac_prompt_systunc_low:.2e} (syst low) + {frac_prompt_systunc_high:.2e} (syst high)')
    print(f'frac_nonprompt: {frac_nonprompt:.2e} +- {frac_nonprompt_statunc:.2e} (stat) - {frac_nonprompt_systunc_low:.2e} (syst low) + {frac_nonprompt_systunc_high:.2e} (syst high)')
    print(f'BR: {BR:.2e} +- {BR_unc:.2e}')
    print(f'sigmaMB (mb): {sigmaMB:.2e}')
    print('____________________________________________________')

    #_____________________________________________________________________________
    # compute yield and cross section
    print('\033[1m\033[92mCompute yield\033[0m')
    hYieldPrompt = TH1F('hYieldPrompt', 'hYield', 1, pt_min, pt_max)
    gYieldPromptSystotal = TGraphAsymmErrors()
    gYieldPromptSystotal.SetName(f'gYieldPromptSystotal')
    hYieldNonPrompt = TH1F('hYieldNonPrompt', 'hYield', 1, pt_min, pt_max)
    gYieldNonPromptSystotal = TGraphAsymmErrors()
    gYieldNonPromptSystotal.SetName(f'gYieldNonPromptSystotal')

    # compute yield
    yield_prompt, stat_prompt = ComputeCrossSection(rawYield, rawYieldUnc, frac_prompt, frac_prompt_statunc,
                                                    effAcc_prompt, delta_pt, delta_y, 1., nev_yield, BR, 'uncorr')
    yield_nonprompt, stat_nonprompt = ComputeCrossSection(rawYield, rawYieldUnc, frac_nonprompt,
                                                          frac_nonprompt_statunc, effAcc_nonprompt, delta_pt, delta_y,
                                                          1., nev_yield, BR, 'uncorr')
    # account for trigger efficiency
    yield_prompt *= eff_trigger
    stat_prompt *= eff_trigger
    yield_nonprompt *= eff_trigger
    stat_nonprompt *= eff_trigger

    totSimmSystUncPrompt = np.sqrt(totSimmSystUnc**2 + (effAccUncPrompt / effAcc_prompt)**2)
    totSimmSystUncNonPrompt = np.sqrt(totSimmSystUnc**2 + (effAccUncNonPrompt / effAcc_nonprompt)**2)
    totSystUncPrompt_high = np.sqrt(totSimmSystUncPrompt**2 + (frac_prompt_systunc_high / frac_prompt)**2)
    totSystUncPrompt_low = np.sqrt(totSimmSystUncPrompt**2 + (frac_prompt_systunc_low / frac_prompt)**2)
    totSystUncNonPrompt_high = np.sqrt(totSimmSystUncNonPrompt**2 + (frac_nonprompt_systunc_high / frac_nonprompt)**2)
    totSystUncNonPrompt_low = np.sqrt(totSimmSystUncNonPrompt**2 + (frac_nonprompt_systunc_low / frac_nonprompt)**2)

    syst_prompt_high = yield_prompt * totSystUncPrompt_high
    syst_prompt_low = yield_prompt * totSystUncPrompt_low
    syst_nonprompt_high = yield_nonprompt * totSystUncNonPrompt_high
    syst_nonprompt_low = yield_nonprompt * totSystUncNonPrompt_low

    print(f'yield_prompt: {yield_prompt:.2e}')
    print(f'yield_nonprompt: {yield_nonprompt:.2e}')
    print(f'stat_prompt: {stat_prompt:.2e}, stat_nonprompt: {stat_nonprompt:.2e}')
    print(f'syst_prompt_low: {syst_prompt_low:.2e}, syst_prompt_high: {syst_prompt_high:.2e}')
    print(f'syst_nonprompt_low: {syst_nonprompt_low:.2e}, syst_nonprompt_high: {syst_nonprompt_high:.2e}')
    hYieldPrompt.SetBinContent(1, yield_prompt)
    hYieldPrompt.SetBinError(1, stat_prompt)
    hYieldNonPrompt.SetBinContent(1, yield_nonprompt)
    hYieldNonPrompt.SetBinError(1, stat_nonprompt)
    gYieldPromptSystotal.SetPoint(0, pt_cent,  yield_prompt)
    gYieldNonPromptSystotal.SetPoint(0, pt_cent,  yield_nonprompt)
    gYieldPromptSystotal.SetPointError(0, delta_pt * 0.25, delta_pt * 0.25, syst_prompt_low, syst_prompt_high)
    gYieldNonPromptSystotal.SetPointError(0, delta_pt * 0.25, delta_pt * 0.25, syst_nonprompt_low, syst_nonprompt_high)

    # compute cross section
    print('\033[1m\033[92mCompute cross section (mb)\033[0m')
    hCrossSecPrompt = TH1F('hCrossSecPrompt', 'hCrossSec', 1, pt_min, pt_max)
    hCrossSecNonPrompt = TH1F('hCrossSecNonPrompt', 'hCrossSec', 1, pt_min, pt_max)
    gCrossSecPromptSystotal = TGraphAsymmErrors()
    gCrossSecPromptSystotal.SetName(f'gCrossSecPromptSystotal')
    gCrossSecNonPromptSystotal = TGraphAsymmErrors()
    gCrossSecNonPromptSystotal.SetName(f'gCrossSecNonPromptSystotal')

    cross_prompt, stat_prompt = ComputeCrossSection(rawYield, rawYieldUnc, frac_prompt, frac_prompt_statunc,
                                                    effAcc_prompt, delta_pt, delta_y, sigmaMB, nev_sigma, BR, 'uncorr')
    cross_nonprompt, stat_nonprompt = ComputeCrossSection(rawYield, rawYieldUnc, frac_nonprompt,
                                                          frac_nonprompt_statunc, effAcc_nonprompt, delta_pt, delta_y,
                                                          sigmaMB, nev_sigma, BR, 'uncorr')

    syst_prompt_high = cross_prompt * totSystUncPrompt_high
    syst_prompt_low = cross_prompt * totSystUncPrompt_low
    syst_nonprompt_high = cross_nonprompt * totSystUncNonPrompt_high
    syst_nonprompt_low = cross_nonprompt *totSystUncNonPrompt_low

    print(f'cross_prompt: {cross_prompt:.2e}')
    print(f'cross_nonprompt: {cross_nonprompt:.2e}')
    print(f'stat_prompt: {stat_prompt:.2e}, stat_nonprompt: {stat_nonprompt:.2e}')
    print(f'syst_prompt_low: {syst_prompt_low:.2e}, syst_prompt_high: {syst_prompt_high:.2e}')
    print(f'syst_nonprompt_low: {syst_nonprompt_low:.2e}, syst_nonprompt_high: {syst_nonprompt_high:.2e}\n')
    hCrossSecPrompt.SetBinContent(1, cross_prompt)
    hCrossSecPrompt.SetBinError(1, stat_prompt)
    hCrossSecNonPrompt.SetBinContent(1, cross_nonprompt)
    hCrossSecNonPrompt.SetBinError(1, stat_nonprompt)
    gCrossSecPromptSystotal.SetPoint(0, pt_cent,  cross_prompt)
    gCrossSecNonPromptSystotal.SetPoint(0, pt_cent,  cross_nonprompt)
    gCrossSecPromptSystotal.SetPointError(0, delta_pt * 0.25, delta_pt * 0.25, syst_prompt_low, syst_prompt_high)
    gCrossSecNonPromptSystotal.SetPointError(0, delta_pt * 0.25, delta_pt * 0.25, syst_nonprompt_low,
                                             syst_nonprompt_high)

    #_____________________________________________________________________________
    # output
    SetObjectStyle(hYieldPrompt, color=kRed+1, markerstyle=20, markersize=1.2, linecolor=kRed+1, markercolor=kRed+1,
                   fillstyle=0, linewidth=2)
    SetObjectStyle(hCrossSecPrompt,  color=kRed+1, markerstyle=20, markersize=1.2, linecolor=kRed+1, markercolor=kRed+1,
                   fillstyle=0, linewidth=2)
    SetObjectStyle(gYieldPromptSystotal, color=kRed+1, markerstyle=20, markersize=1.2, linecolor=kRed+1,
                   markercolor=kRed+1, fillstyle=0, linewidth=2)
    SetObjectStyle(gCrossSecPromptSystotal, color=kRed+1, markerstyle=20, markersize=1.2, linecolor=kRed+1,
                   markercolor=kRed+1, fillstyle=0, linewidth=2)

    SetObjectStyle(hYieldNonPrompt, color=kAzure+4, markerstyle=20, markersize=1.2, linecolor=kAzure+4,
                   markercolor=kAzure+4, fillstyle=0, linewidth=2)
    SetObjectStyle(hCrossSecNonPrompt, color=kAzure+4, markerstyle=20, markersize=1.2, linecolor=kAzure+4,
                   markercolor=kAzure+4, fillstyle=0, linewidth=2)
    SetObjectStyle(gYieldNonPromptSystotal, color=kAzure+4, markerstyle=20, markersize=1.2, linecolor=kAzure+4,
                   markercolor=kAzure+4, fillstyle=0, linewidth=2)
    SetObjectStyle(gCrossSecNonPromptSystotal, color=kAzure+4, markerstyle=20, markersize=1.2, linecolor=kAzure+4,
                   markercolor=kAzure+4, fillstyle=0, linewidth=2)

    outFileName = f'integrated_crosssec_pt{pt_min}-{pt_max}_{reso}_{trigger}_{ptshape}{mult_weights_suffix}{suffix}.root'
    output = TFile.Open(os.path.join(outputdir, outFileName), 'recreate')

    leg = TLegend(0.2, 0.75, 0.7, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetHeader(f'{reso_label} - {trigger}')
    leg.AddEntry(hYieldPrompt, 'Prompt', 'p')
    leg.AddEntry(hYieldNonPrompt, 'Non-prompt', 'p')

    canvas = TCanvas('canvas', 'canvas', 1600, 800)
    canvas.Divide(2)

    yieldTitle = (';#it{p}_{T} (GeV/#it{c}); #frac{d^{2}#it{N}}{d#it{p}_{T}d#it{y}} (GeV^{-1} #it{c})'
                  f' #times BR({reso_label_plot})')
    canvas.cd(1).SetLogy()
    hFrame = canvas.cd(1).DrawFrame(0., yield_nonprompt * 1.e-2, pt_max + 2, yield_prompt * 1.e+2, yieldTitle)
    hFrame.GetXaxis().SetNdivisions(505)
    hYieldNonPrompt.Draw('same')
    hYieldPrompt.Draw('same')
    gYieldPromptSystotal.Draw('2 same')
    gYieldNonPromptSystotal.Draw('2 same')
    leg.Draw()

    crossTitle = (';#it{p}_{T} (GeV/#it{c}); #frac{d^{2}#sigma}{d#it{p}_{T}d#it{y}} (mb GeV^{-1} #it{c})'
                  f' #times BR({reso_label_plot})')
    canvas.cd(2).SetLogy()
    hFrame2 = canvas.cd(2).DrawFrame(0., cross_nonprompt * 1.e-2, pt_max + 2, cross_prompt * 1.e+2, crossTitle)
    hFrame2.GetXaxis().SetNdivisions(505)
    hCrossSecNonPrompt.Draw('pz same')
    hCrossSecPrompt.Draw('pz same')
    gCrossSecPromptSystotal.Draw('2 same')
    gCrossSecNonPromptSystotal.Draw('2 same')
    leg.Draw()
    canvas.SaveAs(os.path.join(outputdir, outFileName.replace('.root', '.pdf')))

    canvas.Write()
    hCrossSecPrompt.Write()
    gCrossSecPromptSystotal.Write()
    hYieldPrompt.Write()
    gYieldPromptSystotal.Write()
    hCrossSecNonPrompt.Write()
    gCrossSecNonPromptSystotal.Write()
    hYieldNonPrompt.Write()
    gYieldNonPromptSystotal.Write()

    # save ingredients
    hRawYields.Write()
    hnev.Write()
    hEffAccPrompt.Write()
    hEffAccNonPrompt.Write()
    gFracPrompt.Write()
    gFracNonPrompt.Write()
    gFracPromptSyst.Write()
    gFracNonPromptSyst.Write()
    hNormIngr.Write()
    output.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("file_rawy", metavar="text",
                        default="rawyields.root", help="input ROOT file with raw yields")
    parser.add_argument("file_eff", metavar="text",
                        default="effacc.root", help="input ROOT file with eff x acc")
    parser.add_argument("file_frac", metavar="text",
                        default="frac.root", help="input ROOT file with fraction")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    compute_crosssec(
        args.file_rawy,
        args.file_eff,
        args.file_frac,
        args.outputdir,
        args.suffix
    )
