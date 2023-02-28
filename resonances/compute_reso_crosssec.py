'''
python script to compute the integrated reco efficiency for resonances
'''

import os
import argparse
import sys
import numpy as np
from ROOT import TFile, TH1F, TCanvas, kRed, kAzure, TLegend, TGraphAsymmErrors, gROOT
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

def compute_crosssec(file_rawy, file_eff, file_frac, outputdir, suffix):
    """
    """

    SetGlobalStyle(padtopmargin=0.05, padleftmargin=0.18,padbottommargin=0.15, palette=55,
                   labelsize=0.04, titlesize=0.05, labeloffset=0.008, titleoffsety=1.7,
                   titleoffsetx=1.2, titleoffsetz=1., opttitle=0, optstat=0)

    gROOT.SetBatch(True)

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
    elif "_multweights_cand" in file_eff:
        mult_weights_suffix = "_multweights_cand"
    elif "_multweights_candinmass" in file_eff:
        mult_weights_suffix = "_multweights_candinmass"

    trigger = ""
    if "MB" in file_rawy:
        trigger = "MB"
    elif "HM" in file_rawy:
        trigger = "HM"

    #_____________________________________________________________________________
    # load input files
    # raw yields
    rawYieldFile = TFile.Open(file_rawy)
    hRawYields = rawYieldFile.Get('h_rawyields')
    rawyield = hRawYields.GetBinContent(1)
    rawyiledunc = hRawYields.GetBinError(1)
    delta_pt = hRawYields.GetBinWidth(1)
    pt_min = hRawYields.GetXaxis().GetBinLowEdge(1)
    pt_max = hRawYields.GetXaxis().GetBinUpEdge(1)
    delta_y = 1. # hard coded

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
    # TODO: missing systematic uncertainties
    effAccUncPrompt = hEffAccPrompt.GetBinError(1)
    effAccUncNonPrompt = hEffAccNonPrompt.GetBinError(1)
    eff_trigger = 0.92 # to be considered when computing the yields

    # fraction of prompt and FD
    fracFile = TFile.Open(file_frac)
    gFracPrompt = fracFile.Get('graph_pfrac_dreso_stat_ptint')
    gFracNonPrompt = fracFile.Get('graph_npfrac_dreso_stat_ptint')
    gFracPromptSyst = fracFile.Get('graph_pfrac_dreso_hypo_ptint')
    gFracNonPromptSyst = fracFile.Get('graph_npfrac_dreso_hypo_ptint')
    gFracNonPromptvspt = fracFile.Get('graph_npfrac_dreso_vspt_stat')
    gFracPromptvspt = gFracNonPromptvspt.Clone()
    for i in range(gFracPromptvspt.GetN()):
        gFracPromptvspt.SetPoint(i, gFracNonPromptvspt.GetX()[i],
                                 1. - gFracNonPromptvspt.GetY()[i])

    frac_prompt = gFracPrompt.GetY()[0]
    frac_nonprompt = gFracNonPrompt.GetY()[0]
    frac_prompt_statunc = gFracPrompt.GetErrorY(0)
    frac_nonprompt_statunc = gFracNonPrompt.GetErrorY(0)
    frac_prompt_systunc_low= gFracPromptSyst.GetErrorYlow(0)
    frac_prompt_systunc_high = gFracPromptSyst.GetErrorYhigh(0)
    frac_nonprompt_systunc_low = gFracNonPromptSyst.GetErrorYlow(0)
    frac_nonprompt_systunc_high = gFracNonPromptSyst.GetErrorYhigh(0)

    # BR: not available for resonances
    BR_Dplus = 9.38e-2 # https://pdglive.lbl.gov/Particle.action?init=0&node=M062&home=MXXX035
    BR_Dplus_unc = 0.16e-2
    BR_Dstar = 2.67e-2 # https://pdglive.lbl.gov/Particle.action?init=0&node=S031&home=MXXX035#decayclump_B
    BR_Dstar_unc = 0.03e-2
    BR_K0s = 69.20e-2 # https://pdglive.lbl.gov/Particle.action?init=0&node=S012&home=MXXX020#decayclump_A
    if d_meson == 'Dplus':
        BR = BR_Dplus
        BR_unc = BR_Dplus_unc
    else:
        BR = BR_Dstar
        BR_unc = BR_Dstar_unc
    BR *= BR_K0s

    # sigma_MB
    sigmaMB = 57.95e+3 # in mb, taken from
    hIngredients = TH1F('hIngredients', 'hIngredients', 10, 0, 10)
    hIngredients.SetBinContent(1, sigmaMB)
    hIngredients.GetXaxis().SetBinLabel(1, 'sigma_MB')
    hIngredients.SetBinContent(2, eff_trigger)
    hIngredients.GetXaxis().SetBinLabel(2, 'eff_trigger (num)')
    hIngredients.SetBinContent(3, BR)
    hIngredients.GetXaxis().SetBinLabel(3, 'BR')
    hIngredients.SetBinContent(4, BR_unc)
    hIngredients.GetXaxis().SetBinLabel(4, 'BR_unc')
    hIngredients.SetBinContent(5, nev_sigma)
    hIngredients.GetXaxis().SetBinLabel(5, 'nev_sigma')
    hIngredients.SetBinContent(6, nev_yield)
    hIngredients.GetXaxis().SetBinLabel(6, 'nev_yield')

    print('\033[1m\033[92mIngredients\033[0m')
    print(f'resonance: {reso}')
    print(f'trigger: {trigger}')
    print(f'delta pt: {delta_pt}')
    print(f'delta y: {delta_y}')
    print(f'rawyield: {rawyield}, rawyiledunc: {rawyiledunc}')
    print(f'nev_sigma: {nev_sigma}, nev_yield: {nev_yield}')
    print(f'effAcc_prompt: {effAcc_prompt}, effAccUncPrompt: {effAccUncPrompt}, effAccUncPrompt: {effAccUncPrompt}')
    print(f'effAcc_nonprompt: {effAcc_nonprompt}, effAccUncNonPrompt: {effAccUncNonPrompt}, effAccUncNonPrompt: {effAccUncNonPrompt}')
    print(f'eff_trigger: {eff_trigger}')
    print(f'frac_prompt: {frac_prompt} +- {frac_prompt_statunc} (stat) - {frac_prompt_systunc_low} (syst low) + {frac_prompt_systunc_high} (syst high)')
    print(f'frac_nonprompt: {frac_nonprompt} +- {frac_nonprompt_statunc} (stat) - {frac_nonprompt_systunc_low} (syst low) + {frac_nonprompt_systunc_high} (syst high)')
    print(f'BR ({d_meson}): {BR} +- {BR_unc}')
    print(f'sigmaMB (mb): {sigmaMB}')
    print('____________________________________________________')

    #_____________________________________________________________________________
    # compute yield and cross section
    print('\033[1m\033[92mCompute yield\033[0m')
    hYieldPrompt = TH1F('hYieldPrompt', 'hYield', 1, -0.5, 0.5)
    gYieldPrompt = TGraphAsymmErrors()
    gYieldPrompt.SetName(f'gYieldPrompt_{reso}_{trigger}')
    gYieldPromptSystotal = TGraphAsymmErrors()
    gYieldPromptSystotal.SetName(f'gYieldPromptSystotal_{reso}_{trigger}')
    gYieldNonPrompt = TGraphAsymmErrors()
    hYieldNonPrompt = TH1F('hYieldNonPrompt', 'hYield', 1, -0.5, 0.5)
    gYieldNonPrompt.SetName(f'gYieldNonPrompt_{reso}_{trigger}')
    gYieldNonPromptSystotal = TGraphAsymmErrors()
    gYieldNonPromptSystotal.SetName(f'gYieldNonPromptSystotal_{reso}_{trigger}')

    # compute yield
    num_prompt = rawyield * frac_prompt * eff_trigger
    num_nonprompt = rawyield * frac_nonprompt * eff_trigger
    den_yield_prompt = 2 * delta_pt * delta_y * effAcc_prompt * nev_yield * BR
    den_yield_nonprompt = 2 * delta_pt * delta_y * effAcc_nonprompt * nev_yield * BR
    yield_prompt = num_prompt/den_yield_prompt
    yield_nonprompt = num_nonprompt/den_yield_nonprompt
    stat_prompt = yield_prompt * np.sqrt((frac_prompt_statunc/frac_prompt)**2 + (rawyiledunc/rawyield)**2 + (effAccUncPrompt/effAcc_prompt)**2)
    stat_nonprompt = yield_nonprompt * np.sqrt((frac_nonprompt_statunc/frac_nonprompt)**2 + (rawyiledunc/rawyield)**2 + (effAccUncNonPrompt/effAcc_nonprompt)**2)

    # TODO: add other systematic uncertainties
    syst_prompt_high = yield_prompt * np.sqrt((frac_prompt_systunc_high/frac_prompt)**2)
    syst_prompt_low = yield_prompt * np.sqrt((frac_prompt_systunc_low/frac_prompt)**2)
    syst_nonprompt_high = yield_nonprompt * np.sqrt((frac_nonprompt_systunc_high/frac_nonprompt)**2)
    syst_nonprompt_low = yield_nonprompt * np.sqrt((frac_nonprompt_systunc_low/frac_nonprompt))

    print(f'num_prompt: {num_prompt:.5f}, den_prompt: {den_yield_prompt:.5f}, yield_prompt: {yield_prompt}')
    print(f'num_nonprompt: {num_nonprompt:.5f}, den_nonprompt: {den_yield_nonprompt:.5f}, yield_nonprompt: {yield_nonprompt}')
    print(f'stat_prompt: {stat_prompt}, stat_nonprompt: {stat_nonprompt}')
    print(f'syst_prompt_low: {syst_prompt_low}, syst_prompt_high: {syst_prompt_high}')
    print(f'syst_nonprompt_low: {syst_nonprompt_low}, syst_nonprompt_high: {syst_nonprompt_high}')
    hYieldPrompt.SetBinContent(1, yield_prompt)
    hYieldPrompt.SetBinError(1, stat_prompt)
    gYieldPrompt.SetPoint(0, 1,  yield_prompt)
    gYieldNonPrompt.SetPoint(0, 1,  yield_nonprompt)
    hYieldNonPrompt.SetBinContent(1, yield_nonprompt)
    hYieldNonPrompt.SetBinError(1, stat_nonprompt)
    gYieldPrompt.SetPointError(0, 0.5, 0.5, stat_prompt, stat_prompt)
    gYieldNonPrompt.SetPointError(0, 0.5, 0.5, stat_nonprompt, stat_nonprompt)
    gYieldPromptSystotal.SetPoint(0, 1,  yield_prompt)
    gYieldNonPromptSystotal.SetPoint(0, 1,  yield_nonprompt)
    gYieldPromptSystotal.SetPointError(0, 0.25, 0.25, syst_prompt_low, syst_prompt_high)
    gYieldNonPromptSystotal.SetPointError(0, 0.25, 0.25, syst_nonprompt_low, syst_nonprompt_high)

    # compute cross section
    print('\033[1m\033[92mCompute cross section (mb)\033[0m')
    hCrossSecPrompt = TH1F('hCrossSecPrompt', 'hCrossSec', 1, -0.5, 0.5)
    gCrossSecPrompt = TGraphAsymmErrors()
    gCrossSecPrompt.SetName(f'gCrossSecPrompt_{reso}_{trigger}')
    hCrossSecNonPrompt = TH1F('hCrossSecNonPrompt', 'hCrossSec', 1, -0.5, 0.5)
    gCrossSecPromptSystotal = TGraphAsymmErrors()
    gCrossSecPromptSystotal.SetName(f'gCrossSecPromptSystotal_{reso}_{trigger}')
    gCrossSecNonPrompt = TGraphAsymmErrors()
    gCrossSecNonPrompt.SetName(f'gCrossSecNonPrompt_{reso}_{trigger}')
    gCrossSecNonPromptSystotal = TGraphAsymmErrors()
    gCrossSecNonPromptSystotal.SetName(f'gCrossSecNonPromptSystotal_{reso}_{trigger}')
    num_prompt = rawyield * frac_prompt * sigmaMB
    num_nonprompt = rawyield * frac_nonprompt * sigmaMB
    den_crosssec_prompt = 2 * delta_pt * delta_y * effAcc_prompt *  nev_sigma
    den_crosssec_nonprompt = 2 * delta_pt * delta_y * effAcc_nonprompt * nev_sigma
    cross_prompt = num_prompt/den_crosssec_prompt
    cross_nonprompt = num_nonprompt/den_crosssec_nonprompt

    stat_prompt = cross_prompt * np.sqrt((frac_prompt_statunc/frac_prompt)**2 +
                                         (rawyiledunc/rawyield)**2 +
                                         (effAccUncPrompt/effAcc_prompt)**2)
    stat_nonprompt = cross_nonprompt * np.sqrt((frac_nonprompt_statunc/frac_nonprompt)**2 +
                                               (rawyiledunc/rawyield)**2 +
                                               (effAccUncNonPrompt/effAcc_nonprompt)**2)

    # TODO: add other systematic uncertainties
    syst_prompt_high = cross_prompt * np.sqrt((frac_prompt_systunc_high/frac_prompt)**2)
    syst_prompt_low = cross_prompt * np.sqrt((frac_prompt_systunc_low/frac_prompt)**2)
    syst_nonprompt_high = cross_nonprompt * np.sqrt((frac_nonprompt_systunc_high/frac_nonprompt)**2)
    syst_nonprompt_low = cross_nonprompt * np.sqrt((frac_nonprompt_systunc_low/frac_nonprompt))

    print(f'num_prompt: {num_prompt:.5f}, den_prompt: {den_crosssec_prompt:.5f}, cross_prompt: {cross_prompt}')
    print(f'num_nonprompt: {num_nonprompt:.5f}, den_nonprompt: {den_crosssec_nonprompt:.5f}, cross_nonprompt: {cross_nonprompt}')
    print(f'stat_prompt: {stat_prompt}, stat_nonprompt: {stat_nonprompt}')
    print(f'syst_prompt_low: {syst_prompt_low}, syst_prompt_high: {syst_prompt_high}')
    print(f'syst_nonprompt_low: {syst_nonprompt_low}, syst_nonprompt_high: {syst_nonprompt_high}\n')
    hCrossSecPrompt.SetBinContent(1, cross_prompt)
    hCrossSecPrompt.SetBinError(1, stat_prompt)
    gCrossSecPrompt.SetPoint(0, 1,  cross_prompt)
    hCrossSecNonPrompt.SetBinContent(1, cross_nonprompt)
    hCrossSecNonPrompt.SetBinError(1, stat_nonprompt)
    gCrossSecNonPrompt.SetPoint(0, 1,  cross_nonprompt)
    gCrossSecPrompt.SetPointError(0, 0.5, 0.5, stat_prompt, stat_prompt)
    gCrossSecNonPrompt.SetPointError(0, 0.5, 0.5, stat_nonprompt, stat_nonprompt)
    gCrossSecPromptSystotal.SetPoint(0, 1,  cross_prompt)
    gCrossSecNonPromptSystotal.SetPoint(0, 1,  cross_nonprompt)
    gCrossSecPromptSystotal.SetPointError(0, 0.25, 0.25, syst_prompt_low, syst_prompt_high)
    gCrossSecNonPromptSystotal.SetPointError(0, 0.25, 0.25, syst_nonprompt_low, syst_nonprompt_high)

    #_____________________________________________________________________________
    # output
    SetObjectStyle(gYieldPrompt, color=kRed+1, markerstyle=20, markersize=1.2,
                   linecolor=kRed+1, markercolor=kRed+1, fillalpha=0.2, fillstyle=0,
                   linewidth=2)
    SetObjectStyle(gCrossSecPrompt, color=kRed+1, markerstyle=20, markersize=1.2,
                   linecolor=kRed+1, markercolor=kRed+1, fillalpha=0.2, fillstyle=0,
                   linewidth=2)
    SetObjectStyle(gYieldPromptSystotal, color=kRed+1, markerstyle=20, markersize=1.2,
                   linecolor=kRed+1, markercolor=kRed+1, fillalpha=0.2, fillstyle=0,
                   linewidth=2)
    SetObjectStyle(gCrossSecPromptSystotal, color=kRed+1, markerstyle=20, markersize=1.2,
                   linecolor=kRed+1, markercolor=kRed+1, fillalpha=0.2, fillstyle=0,
                   linewidth=2)

    SetObjectStyle(gYieldNonPrompt, color=kAzure+4, markerstyle=20, markersize=1.2,
                   linecolor=kAzure+4, markercolor=kAzure+4, fillalpha=0.2, fillstyle=0,
                   linewidth=2)
    SetObjectStyle(gCrossSecNonPrompt, color=kAzure+4, markerstyle=20, markersize=1.2,
                   linecolor=kAzure+4, markercolor=kAzure+4, fillalpha=0.2, fillstyle=0,
                   linewidth=2)
    SetObjectStyle(gYieldNonPromptSystotal, color=kAzure+4, markerstyle=20, markersize=1.2,
                   linecolor=kAzure+4, markercolor=kAzure+4, fillalpha=0.2, fillstyle=0,
                   linewidth=2)
    SetObjectStyle(gCrossSecNonPromptSystotal, color=kAzure+4, markerstyle=20, markersize=1.2,
                   linecolor=kAzure+4, markercolor=kAzure+4, fillalpha=0.2, fillstyle=0,
                   linewidth=2)

    output = TFile.Open(
        os.path.join(outputdir, f'integrated_crosssec_pt{pt_min}-{pt_max}_{reso}_{trigger}_{ptshape}{mult_weights_suffix}{suffix}.root'),
        'recreate'
    )
    canvas = TCanvas('canvas', 'canvas', 1600, 900)
    canvas.Divide(2)
    leg = TLegend(0.2, 0.75, 0.7, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetHeader(f'{reso_label} - {trigger}')
    leg2 = TLegend(0.2, 0.75, 0.7, 0.9)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextSize(0.03)
    leg2.SetHeader(f'{reso_label} - {trigger}')
    leg.AddEntry(gYieldPrompt, 'Prompt', 'p')
    leg.AddEntry(gYieldNonPrompt, 'Non-prompt', 'p')
    leg2.AddEntry(gCrossSecPrompt, 'Prompt', 'p')
    leg2.AddEntry(gCrossSecNonPrompt, 'Non-prompt', 'p')

    ymin = yield_nonprompt * 1.e-2
    ymax = yield_prompt * 1.e+2
    canvas.cd(1).SetLogy()
    hFrame = canvas.cd(1).DrawFrame(0., ymin, 2, ymax, f';; #frac{{d#it{{N}}}}{{d#it{{y}}}}#timesBR({reso_label_plot})')
    hFrame.GetXaxis().SetNdivisions(505)
    gYieldNonPrompt.Draw('pz same')
    gYieldPrompt.Draw('pz same')
    gYieldPromptSystotal.Draw('2 same')
    gYieldNonPromptSystotal.Draw('2 same')
    leg.Draw()

    ymin = cross_nonprompt * 1.e-2
    ymax = cross_prompt * 1.e+2
    canvas.cd(2).SetLogy()
    hFrame2 = canvas.cd(2).DrawFrame(0., ymin, 2, ymax, ';; integrated cross section (mb)')
    hFrame2.GetXaxis().SetNdivisions(505)
    gCrossSecNonPrompt.Draw('pz same')
    gCrossSecPrompt.Draw('pz same')
    gCrossSecPromptSystotal.Draw('2 same')
    gCrossSecNonPromptSystotal.Draw('2 same')
    leg2.Draw()
    canvas.SaveAs(
        os.path.join(
            outputdir, f'integrated_crosssec_pt{pt_min}-{pt_max}_{reso}_{trigger}_{ptshape}{mult_weights_suffix}{suffix}.pdf'
        )
    )

    canvas.Write()
    hCrossSecPrompt.Write()
    hYieldPrompt.Write()
    gYieldPrompt.Write()
    hCrossSecNonPrompt.Write()
    hYieldNonPrompt.Write()
    gYieldNonPrompt.Write()
    gYieldPromptSystotal.Write()
    gYieldNonPromptSystotal.Write()
    gCrossSecNonPrompt.Write()
    gCrossSecPrompt.Write()
    gCrossSecPromptSystotal.Write()
    gCrossSecNonPromptSystotal.Write()

    # save ingredients
    hRawYields.Write()
    hnev.Write()
    hEffAccPrompt.Write()
    hEffAccNonPrompt.Write()
    gFracPrompt.Write()
    gFracNonPrompt.Write()
    gFracPromptSyst.Write()
    gFracNonPromptSyst.Write()
    hIngredients.Write()
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
