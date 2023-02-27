'''
python script to compute the integrated reco efficiency for resonances
run: python compute_reso_crossec.py
'''
import sys
import argparse
import yaml
import numpy as np
import pandas as pd
import uproot
from  multiprocessing import Process
from alive_progress import alive_bar
from particle import Particle
from ROOT import TFile, TH2F, TH1F, TCanvas, TLatex, kBlack, kRed, kAzure, TLegend, TGraphAsymmErrors
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.AnalysisUtils import ComputeCrossSection
SetGlobalStyle(padtopmargin=0.05, padleftmargin=0.18,padbottommargin=0.15, palette=55, labelsize=0.04, titlesize=0.05,
               labeloffset=0.008, titleoffsety=1.7, titleoffsetx=1.2, titleoffsetz=1.,
               opttitle=0, optstat=0)

resos = ['10433', '10433']#, '10433', '435', '10433']
mults = ['MB', 'HM']#, 'HM', 'MB', 'MB']
pt_min = 2
pt_max = 24
delta_pt = pt_max - pt_min
delta_y = 1

for i, (reso, mult) in enumerate(zip(resos, mults)):
    if reso == '435':
        reso_label = 'Ds2starplus'
        reso_name = Particle.from_pdgid(reso).name
        d_meson = 'Dplus'
        v0 = 'K0S'
        path = f'./435_Ds2_Dplus_K0S/{mult}'
        reso_label_plot = 'D_{s2}^{*+}'
    elif reso == '10433':
        reso_label = 'Ds1plus'
        reso_name = Particle.from_pdgid(reso).name
        d_meson = 'Dstar'
        v0 = 'K0S'
        path = f'./10433_Ds1_Dstar_K0S/{mult}'
        reso_label_plot = 'D_{s1}^{+}'
    else:
        raise ValueError('Resonance not supported')

    #_____________________________________________________________________________
    # load input files
    # raw yields
    rawYieldFile = TFile.Open(f'./input/rawyields/{reso_label}/mass_{reso_label}_pt2.0-24.0_{mult}.root')
    hRawYields = rawYieldFile.Get('h_rawyields')
    rawyield = hRawYields.GetBinContent(1)
    rawyiledunc = hRawYields.GetBinError(1)

    # nevents
    norm_file = TFile.Open(f'./input/rawyields/{reso_label}/normalisation_{d_meson}_{v0}_{mult}.root')
    hnev = norm_file.Get('hist_events')
    nev_sigma = hnev.GetBinContent(1)
    nev_yield = hnev.GetBinContent(2)

    # eff. acc.
    effAccFile = TFile.Open(f'{path}/integrated_recoeff_reso{reso}_{mult}_Dsptshape.root')
    gEffAccPrompt = effAccFile.Get(f'geff_int_prompt_{reso}_{mult}_Dsptshape')
    gEffAccFD = effAccFile.Get(f'geff_int_nonprompt_{reso}_{mult}_Dsptshape')
    hEffAccPromptvspt = effAccFile.Get('hreco_pt') # to check pT-differential spectra
    hEffAccFDvspt = effAccFile.Get('hreco_pt_np') # to check pT-differential spectra
    effAcc_prompt = gEffAccPrompt.GetY()[0]
    effAcc_nonprompt = gEffAccFD.GetY()[0]
    # TODO: missing systematic uncertainties
    effAccUncLow_prompt = gEffAccPrompt.GetErrorYlow(0)
    effAccUncHigh_prompt = gEffAccPrompt.GetErrorYhigh(0)
    effAccUncLow_nonprompt = gEffAccFD.GetErrorYlow(0)
    effAccUncHigh_nonprompt = gEffAccFD.GetErrorYhigh(0)
    eff_trigger = 0.92 # to be considered when computing the yields

    # fraction of prompt and FD
    fracFile = TFile.Open(f'{path}/frac_{reso}_{mult}_Feb9.root') if reso == '10433' else TFile.Open(f'{path}/frac_{reso}_{mult}.root') 
    gFracPrompt = fracFile.Get(f'graph_pfrac_dreso_stat_pt2.0-24.0')
    gFracNonPrompt = fracFile.Get(f'graph_npfrac_dreso_stat_pt2.0-24.0')
    gFracPromptSyst = fracFile.Get(f'graph_pfrac_dreso_hypo_pt2.0-24.0')
    gFracNonPromptSyst = fracFile.Get(f'graph_npfrac_dreso_hypo_pt2.0-24.0')
    gFracNonPromptvspt = fracFile.Get(f'graph_npfrac_dreso_vspt_stat')
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
    print(f'multiplicity: {mult}')
    print(f'delta pt: {delta_pt}')
    print(f'delta y: {delta_y}')
    print(f'rawyield: {rawyield}, rawyiledunc: {rawyiledunc}')
    print(f'nev_sigma: {nev_sigma}, nev_yield: {nev_yield}')
    print(f'effAcc_prompt: {effAcc_prompt}, effAccUncLow_prompt: {effAccUncLow_prompt}, effAccUncHigh_prompt: {effAccUncHigh_prompt}')
    print(f'effAcc_nonprompt: {effAcc_nonprompt}, effAccUncLow_nonprompt: {effAccUncLow_nonprompt}, effAccUncHigh_nonprompt: {effAccUncHigh_nonprompt}')
    print(f'eff_trigger: {eff_trigger}')
    print(f'frac_prompt: {frac_prompt} +- {frac_prompt_statunc} (stat) - {frac_prompt_systunc_low} (syst low) + {frac_prompt_systunc_high} (syst high)')
    print(f'frac_nonprompt: {frac_nonprompt} +- {frac_nonprompt_statunc} (stat) - {frac_nonprompt_systunc_low} (syst low) + {frac_nonprompt_systunc_high} (syst high)')
    print(f'BR ({d_meson}): {BR} +- {BR_unc}')
    print(f'sigmaMB (mb): {sigmaMB}')
    print('____________________________________________________')
    input('Press enter to continue')

    #_____________________________________________________________________________
    # compute yield and cross section
    print('\033[1m\033[92mCompute yield\033[0m')
    hYieldPrompt = TH1F('hYieldPrompt', 'hYield', 1, -0.5, 0.5)
    gYieldPrompt = TGraphAsymmErrors()
    gYieldPrompt.SetName(f'gYieldPrompt_{reso}_{mult}')
    gYieldPromptSystotal = TGraphAsymmErrors()
    gYieldPromptSystotal.SetName(f'gYieldPromptSystotal_{reso}_{mult}')
    gYieldNonPrompt = TGraphAsymmErrors()
    hYieldNonPrompt = TH1F('hYieldNonPrompt', 'hYield', 1, -0.5, 0.5)
    gYieldNonPrompt.SetName(f'gYieldNonPrompt_{reso}_{mult}')
    gYieldNonPromptSystotal = TGraphAsymmErrors()
    gYieldNonPromptSystotal.SetName(f'gYieldNonPromptSystotal_{reso}_{mult}')

    # compute yield
    num_prompt = rawyield * frac_prompt * eff_trigger
    num_nonprompt = rawyield * frac_nonprompt * eff_trigger
    den_yield_prompt = 2 * delta_pt * delta_y * effAcc_prompt * nev_yield * BR
    den_yield_nonprompt = 2 * delta_pt * delta_y * effAcc_nonprompt * nev_yield * BR
    yield_prompt = num_prompt/den_yield_prompt
    yield_nonprompt = num_nonprompt/den_yield_nonprompt
    stat_prompt = yield_prompt * np.sqrt((frac_prompt_statunc/frac_prompt)**2 + (rawyiledunc/rawyield)**2 + (effAccUncLow_prompt/effAcc_prompt)**2)
    stat_nonprompt = yield_nonprompt * np.sqrt((frac_nonprompt_statunc/frac_nonprompt)**2 + (rawyiledunc/rawyield)**2 + (effAccUncLow_nonprompt/effAcc_nonprompt)**2)

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
    gCrossSecPrompt.SetName(f'gCrossSecPrompt_{reso}_{mult}')
    hCrossSecNonPrompt = TH1F('hCrossSecNonPrompt', 'hCrossSec', 1, -0.5, 0.5)
    gCrossSecPromptSystotal = TGraphAsymmErrors()
    gCrossSecPromptSystotal.SetName(f'gCrossSecPromptSystotal_{reso}_{mult}')
    gCrossSecNonPrompt = TGraphAsymmErrors()
    gCrossSecNonPrompt.SetName(f'gCrossSecNonPrompt_{reso}_{mult}')
    gCrossSecNonPromptSystotal = TGraphAsymmErrors()
    gCrossSecNonPromptSystotal.SetName(f'gCrossSecNonPromptSystotal_{reso}_{mult}')
    num_prompt = rawyield * frac_prompt * sigmaMB
    num_nonprompt = rawyield * frac_nonprompt * sigmaMB
    den_crosssec_prompt = 2 * delta_pt * delta_y * effAcc_prompt *  nev_sigma
    den_crosssec_nonprompt = 2 * delta_pt * delta_y * effAcc_nonprompt * nev_sigma
    cross_prompt = num_prompt/den_crosssec_prompt
    cross_nonprompt = num_nonprompt/den_crosssec_nonprompt

    stat_prompt = cross_prompt * np.sqrt((frac_prompt_statunc/frac_prompt)**2 +
                                         (rawyiledunc/rawyield)**2 +
                                         (effAccUncLow_prompt/effAcc_prompt)**2)
    stat_nonprompt = cross_nonprompt * np.sqrt((frac_nonprompt_statunc/frac_nonprompt)**2 +
                                               (rawyiledunc/rawyield)**2 +
                                               (effAccUncLow_nonprompt/effAcc_nonprompt)**2)

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

    output = TFile.Open(f'{path}/integrated_crosssec_{reso}_{mult}_pt{pt_min}-{pt_max}_Feb17.root', 'recreate')
    canvas = TCanvas('canvas', 'canvas', 1600, 900)
    canvas.Divide(2)
    leg = TLegend(0.2, 0.75, 0.7, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetHeader(f'{reso_label} - {mult}')
    leg2 = TLegend(0.2, 0.75, 0.7, 0.9)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextSize(0.03)
    leg2.SetHeader(f'{reso_label} - {mult}')
    leg.AddEntry(gYieldPrompt, f'Prompt', 'p')
    leg.AddEntry(gYieldNonPrompt, 'Non-prompt', 'p')
    leg2.AddEntry(gCrossSecPrompt, f'Prompt', 'p')
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
    canvas.SaveAs(f'{path}/integrated_crosssec_{reso}_{mult}_pt{pt_min}-{pt_max}_Feb17.pdf')

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
    gEffAccPrompt.Write()
    gEffAccFD.Write()
    gFracPrompt.Write()
    gFracNonPrompt.Write()
    gFracPromptSyst.Write()
    gFracNonPromptSyst.Write()
    hIngredients.Write()
    output.Close()

input('Press enter to exit')
