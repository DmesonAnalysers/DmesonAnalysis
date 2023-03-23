'''
python script to compute the ratio of the integrated yields of the resonance to the ground state
run: python compute_resotogroundstate_ratio.py
'''
import sys
import numpy as np
from ROOT import TFile, TColor, TCanvas, TLatex, TLegend, TGraphAsymmErrors # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure, kGreen, kOrange, kFullCross, kFullSquare # pylint: disable=import-error,no-name-in-module
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.16, titleoffsety=1.55)

reso_list = [435, 10433]
base_path = '/home/fcatalan/Desktop/Analyses/Ds_resonances'

delta_pt = 22
avg_mults = [6.9, 31.5] # MB, HM

for reso in reso_list:

    if reso == 10433:
        reso_label = 'Ds1plus'
        meson = 'Dstar'
        meson_label_plot = 'D^{*+}'
        reso_label_plot = 'D_{s1}^{+}'
    elif reso == 435:
        reso_label = 'Ds2starplus'
        meson = 'Dplus'
        meson_label_plot = 'D^{+}'
        reso_label_plot = 'D_{s2}^{*+}'
    else:
        sys.exit()

    print(f'\033[1m\033[4mRESONANCE: {reso}\033[0m')
    #_____________________________________________________________
    # Collect input files
    # Resonances
    infile_reso_HM = TFile.Open(f'{base_path}/corr_yields/integrated_crosssec_pt2.0-24.0_{reso}_HM_Dsptshape_multweights_candinmass.root')
    infile_reso_MB = TFile.Open(f'{base_path}/corr_yields/integrated_crosssec_pt2.0-24.0_{reso}_MB_Dsptshape_multweights_candinmass.root')
    hreso_HM = infile_reso_HM.Get('hYieldPrompt')
    greso_HM_syst = infile_reso_HM.Get('gYieldPromptSystotal')
    hreso_MB = infile_reso_MB.Get('hYieldPrompt')
    greso_MB_syst = infile_reso_MB.Get('gYieldPromptSystotal')

    # Ground state
    infile_gs = TFile.Open(f'{base_path}/Ds_pp13TeV_vsMult/V0M_ptintegratedyields_D0DsLc_pt2-24.root') # from Luuk
    hgs_MB = infile_gs.Get('h_vis_dNdy_1_0')
    ggs_MB = infile_gs.Get('gr_vis_dNdy_TotSyst_1_0')
    hgs_HM = infile_gs.Get('h_vis_dNdy_1_4')
    ggs_HM = infile_gs.Get('gr_vis_dNdy_TotSyst_1_4')

    #_____________________________________________________________
    # Scale for deltapt
    hgs_HM.Scale(1./delta_pt)
    ggs_HM.SetPoint(0, ggs_HM.GetPointX(0), ggs_HM.GetPointY(0) / delta_pt)
    ggs_HM.SetPointError(0, ggs_HM.GetErrorXlow(0), ggs_HM.GetErrorXhigh(0),
                        ggs_HM.GetErrorYlow(0) / delta_pt, ggs_HM.GetErrorYhigh(0) / delta_pt)
    hgs_MB.Scale(1./delta_pt)
    ggs_MB.SetPoint(0, ggs_MB.GetPointX(0), ggs_MB.GetPointY(0) / delta_pt)
    ggs_MB.SetPointError(0, ggs_MB.GetErrorXlow(0), ggs_MB.GetErrorXhigh(0),
                        ggs_MB.GetErrorYlow(0) / delta_pt, ggs_MB.GetErrorYhigh(0) / delta_pt)

    #_____________________________________________________________
    # Compute ratio
    gratio = TGraphAsymmErrors()
    gratio.SetName('gratio')
    gratio_sys = TGraphAsymmErrors()
    gratio_sys.SetName('gratio_sys')

    # MB
    print('\033[1m\033[4mMB\033[0m')
    num = hreso_MB.GetBinContent(1)
    num_unc = hreso_MB.GetBinError(1)
    den = hgs_MB.GetBinContent(1)
    den_unc = hgs_MB.GetBinError(1)
    ratio_MB = num / den
    ratio_MB_stat = ratio_MB * np.sqrt((num_unc / num)**2 + (den_unc / den)**2)
    ratio_MB_sys_low = ratio_MB * np.sqrt((greso_MB_syst.GetErrorYlow(0) / num)**2 + (ggs_MB.GetErrorYlow(0) / den)**2)
    ratio_MB_sys_high = ratio_MB * np.sqrt((greso_MB_syst.GetErrorYhigh(0) / num)**2 + (ggs_MB.GetErrorYhigh(0) / den)**2)
    print(f'Reso integrated yield = {num:.2e} +/- {num_unc:.2e}')
    print(f'Ds integrated yield = {den:.2e} +/- {den_unc:.2e}')
    print(f'ratio_MB = {ratio_MB:.3f} +/- {ratio_MB_stat:.3f}')
    print(f'ratio_MB_sys_low = {ratio_MB_sys_low:.3f}, ratio_MB_sys_high = {ratio_MB_sys_high:.3f}')

    # HM
    print('\033[1m\033[4mHM\033[0m')
    num = hreso_HM.GetBinContent(1)
    num_unc = hreso_HM.GetBinError(1)
    den = hgs_HM.GetBinContent(1)
    den_unc = hgs_HM.GetBinError(1)
    ratio_HM = num / den
    ratio_HM_stat = ratio_HM * np.sqrt((num_unc / num)**2 + (den_unc / den)**2)
    ratio_HM_sys_low = ratio_HM * np.sqrt((greso_HM_syst.GetErrorYlow(0) / num)**2 + (ggs_HM.GetErrorYlow(0) / den)**2)
    ratio_HM_sys_high = ratio_HM * np.sqrt((greso_HM_syst.GetErrorYhigh(0) / num)**2 + (ggs_HM.GetErrorYhigh(0) / den)**2)
    print(f'Reso integrated yield = {num:.2e} +/- {num_unc:.2e}')
    print(f'Ds integrated yield = {den:.2e} +/- {den_unc:.2e}')
    print(f'ratio_HM = {ratio_HM:.3f} +/- {ratio_HM_stat:.3f}')
    print(f'ratio_MB_sys_low = {ratio_HM_sys_low:.3f}, ratio_MB_sys_high = {ratio_HM_sys_high:.3f}')

    gratio.SetPoint(0, avg_mults[0], ratio_MB)
    gratio.SetPointError(0, 0, 0, ratio_MB_stat, ratio_MB_stat)
    gratio_sys.SetPoint(0, avg_mults[0], ratio_MB)
    gratio_sys.SetPointError(0, 1., 1., ratio_MB_sys_low, ratio_MB_sys_high)
    gratio.SetPoint(1, avg_mults[1], ratio_HM)
    gratio.SetPointError(1, 0, 0, ratio_HM_stat, ratio_HM_stat)
    gratio_sys.SetPoint(1, avg_mults[1], ratio_HM)
    gratio_sys.SetPointError(1, 1., 1., ratio_HM_sys_low, ratio_HM_sys_high)
    color = kGreen+3 if reso == 435 else kRed+2
    marker = kFullCross if reso == 435 else kFullSquare
    markersize = 1.4 if reso == 435 else 1.2
    SetObjectStyle(gratio, color=color, markerstyle=marker, markersize=markersize, fillstyle=0, linewidth=2)
    SetObjectStyle(gratio_sys, color=color, markerstyle=marker, markersize=markersize, fillstyle=0, linewidth=2)

    #_____________________________________________________________
    # Collect theoretical prediction
    br_reso = 23.35e-2 if reso == 435 else 22.e-2 # taken from AN (see also https://indico.cern.ch/event/1253172/)
    br_reso_unc = 0. if reso == 435 else 2.e-2
    if reso == 435:
        pred_shm = 0.075 # taken from ./predictions/MinHe/
        pred_shmc = 0.034 / 0.357 # taken from Ds2*+/D0 = 0.034 and Ds/D0 = 0.357
        pred_shm_cent = pred_shm * br_reso
        pred_shm_low = pred_shm * (br_reso - br_reso_unc)
        pred_shm_up = pred_shm * (br_reso + br_reso_unc)
        pred_shmc_cent = pred_shmc * br_reso
        pred_shmc_low = pred_shmc * (br_reso - br_reso_unc)
        pred_shmc_up = pred_shmc * (br_reso + br_reso_unc)
        br_str = f'BR = {br_reso*100:.2f}%'
        paper_str = 'PRD 93 (2016) 034035'
        model_leg_opt = 'l'
        linewidth = 3

    else:
        pred_shm = 0.054 # taken from ./predictions/MinHe/
        pred_shmc = 0.026 / 0.357 # taken from Ds1+/D0 = 0.026 and Ds/D0 = 0.357
        pred_shm_cent = pred_shm * br_reso
        pred_shm_low = pred_shm * (br_reso - br_reso_unc)
        pred_shm_up = pred_shm * (br_reso + br_reso_unc)
        pred_shmc_cent = pred_shmc * br_reso
        pred_shmc_low = pred_shmc * (br_reso - br_reso_unc)
        pred_shmc_up = pred_shmc * (br_reso + br_reso_unc)
        br_str = f'BR = ({br_reso*100:.0f} #pm {br_reso_unc*100:.0f})%'
        paper_str = 'PRD 93 (2016) 034035 and PTEP 2022 083C01'
        model_leg_opt = 'f'
        linewidth = 0

    gtheo_shm, gtheo_shmc = TGraphAsymmErrors(), TGraphAsymmErrors()
     # no explicit multiplicity dependence indicated from theoreticians
    gtheo_shm.SetPoint(0, 4, pred_shm_cent)
    gtheo_shm.SetPointError(0, 0., 0., pred_shm_cent - pred_shm_low, pred_shm_up - pred_shm_cent)
    gtheo_shm.SetPoint(1, 36, pred_shm_cent)
    gtheo_shm.SetPointError(1, 0., 0., pred_shm_cent - pred_shm_low, pred_shm_up - pred_shm_cent)
    gtheo_shmc.SetPoint(0, 4, pred_shmc_cent)
    gtheo_shmc.SetPointError(0, 0., 0., pred_shmc_cent - pred_shmc_low, pred_shmc_up - pred_shmc_cent)
    gtheo_shmc.SetPoint(1, 36, pred_shmc_cent)
    gtheo_shmc.SetPointError(1, 0., 0., pred_shmc_cent - pred_shmc_low, pred_shmc_up - pred_shmc_cent)

    # define custom colors to mimic transparency
    kAzureMy = TColor.GetFreeColorIndex()
    cAzureMy = TColor(kAzureMy, 194./255, 213./255, 237./255, 'kAzureMy', 1.0)
    kOrangeMy = TColor.GetFreeColorIndex()
    cOrangeMy = TColor(kOrangeMy, 255./255, 190./255, 92./255, 'kOrangeMy', 1.0)
    SetObjectStyle(gtheo_shm, linecolor=kAzure+4, linewidth=linewidth, fillcolor=kAzureMy, linestyle=5)
    SetObjectStyle(gtheo_shmc, linecolor=kOrange-3, linewidth=linewidth, fillcolor=kOrangeMy, linestyle=2)

    #_____________________________________________________________
    # nSigma calculation
    # use lower syst. unc. since model predictions are lower than data
    ratio_tot_unc_MB = np.sqrt(ratio_MB_stat**2 + ratio_MB_sys_low**2)
    ratio_tot_unc_HM = np.sqrt(ratio_HM_stat**2 + ratio_HM_sys_low**2)
    nsigma_mb_SHM = (ratio_MB - gtheo_shmc.GetY()[0]) / ratio_tot_unc_MB
    nsigma_mb_He = (ratio_MB - gtheo_shm.GetY()[0]) / ratio_tot_unc_MB
    nsigma_hm_SHM = (ratio_HM - gtheo_shmc.GetY()[0]) / ratio_tot_unc_HM
    nsigma_hm_He = (ratio_HM - gtheo_shm.GetY()[0]) / ratio_tot_unc_HM

    print('\n\033[1m\033[4mnSigma (only stat)\033[0m')
    print(f'nsigma_mb_SHM = {nsigma_mb_SHM:.2f}')
    print(f'nsigma_mb_He = {nsigma_mb_He:.2f}')
    print(f'nsigma_hm_SHM = {nsigma_hm_SHM:.2f}')
    print(f'nsigma_hm_He = {nsigma_hm_He:.2f}\n\n')

    #_____________________________________________________________
    # Plot and save output

    leg = TLegend(0.4, 0.65, 0.9, 0.8)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.AddEntry(gratio, 'Data #scale[0.8]{(2 < #it{p}_{T} < 24 GeV/#it{c})}', 'p')
    leg.AddEntry(gtheo_shm, 'SHM M. He, R. Rapp #scale[0.8]{(#it{p}_{T} > 0)}', model_leg_opt)
    leg.AddEntry(gtheo_shmc, 'SHMc GSI#minusHeidelberg #scale[0.8]{(#it{p}_{T} > 0)}', model_leg_opt)

    latALICE = TLatex()
    latALICE.SetNDC()
    latALICE.SetTextSize(0.06)
    latALICE.SetTextFont(42)
    latALICE.SetTextColor(kBlack)

    latLabel = TLatex()
    latLabel.SetNDC()
    latLabel.SetTextSize(0.05)
    latLabel.SetTextFont(42)
    latLabel.SetTextColor(kBlack)

    latNote = TLatex()
    latNote.SetNDC()
    latNote.SetTextSize(0.035)
    latNote.SetTextFont(42)
    latNote.SetTextColor(kBlack)

    canvas = TCanvas('canvas', 'canvas', 800, 800)
    decay_tag = f'{reso_label_plot} #times BR({reso_label_plot} #rightarrow {meson_label_plot} K^{{0}}_{{S}})'
    ytitle = f';#LTd#it{{N}}_{{ch}}/d#it{{#eta}}#GT_{{|#it{{#eta}}| < 0.5}}; {decay_tag} / D_{{s}}^{{+}}'
    ymax = 0.109 if reso == 435 else 0.054
    hFrame = canvas.cd().DrawFrame(0., 0., 39, ymax, ytitle)
    latALICE.DrawLatex(0.19, 0.89, 'ALICE Preliminary')
    latLabel.DrawLatex(0.19, 0.83, 'pp, #sqrt{#it{s}} = 13 TeV, |#it{y}| < 0.5')
    latNote.DrawLatex(0.19, 0.21, 'Model predictions only')
    latNote.DrawLatex(0.19, 0.175, f'{br_str} #scale[0.8]{{{paper_str}}}')
    gtheo_shm.Draw('c3 same')
    gtheo_shmc.Draw('c3 same')
    gratio_sys.Draw('2 same')
    gratio.Draw('same pz')
    leg.Draw()

    canvas.SaveAs(f'{base_path}/prel_figures/Ds{reso}_over_Ds_ratio_vs_mult.pdf')
    canvas.SaveAs(f'{base_path}/prel_figures/Ds{reso}_over_Ds_ratio_vs_mult.eps')
    outfile = TFile(f'{base_path}/prel_figures/Ds{reso}_over_Ds_ratio_vs_mult.root', 'recreate')
    canvas.Write()
    gratio.Write()
    gratio_sys.Write()
    gtheo_shm.Write()
    gtheo_shmc.Write()
    outfile.Close()
