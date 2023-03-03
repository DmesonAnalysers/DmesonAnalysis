'''
python script to compute the ratio of the integrated yields of the resonance to the ground state
run: python compute_integrated_recoeff.py
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
from ROOT import TFile, TH2F, TH1F, TCanvas, TLatex, kBlack, kRed, kAzure, kOrange, kOpenCross, kOpenDiamond, TLegend, TGraphAsymmErrors, TGraph
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.AnalysisUtils import ComputeCrossSection
SetGlobalStyle(padtopmargin=0.05,padbottommargin=0.15, padleftmargin=0.2, padrightmargin=0.05,
               palette=55, labelsize=0.04, titlesize=0.05,
               labeloffset=0.008, titleoffsety=1.4, titleoffsetx=1.2, titleoffsetz=1.,
               opttitle=0, optstat=0)

reso = '10433'
use_mb_cross = False # if True, pt-integrated cross section adopted as ground state, otherwise vsmult paper yield
outlabel = '_030323_new'

if reso == '435':
    meson = 'Dplus'
    path = '435_Ds2_Dplus_K0S'
else:
    meson = 'Dstar'
    path = '10433_Ds1_Dstar_K0S'
delta_pt = 22
    
print(f'\033[1m\033[4mRESONANCE: {reso}\033[0m')
#_____________________________________________________________
# Collect input files
# Resonances
infile_reso_HM = TFile.Open(f'/home/stefano/Desktop/cernbox/Ds_resonances/corr_yields/integrated_crosssec_pt2.0-24.0_{reso}_HM_Dsptshape_multweights_candinmass.root')
infile_reso_MB = TFile.Open(f'/home/stefano/Desktop/cernbox/Ds_resonances/corr_yields/integrated_crosssec_pt2.0-24.0_{reso}_MB_Dsptshape_multweights_candinmass.root')
hreso_HM = infile_reso_HM.Get('hYieldPrompt')
greso_HM_syst = infile_reso_HM.Get(f'gYieldPromptSystotal_{reso}_HM')
hreso_MB = infile_reso_MB.Get('hYieldPrompt')
greso_MB_syst = infile_reso_MB.Get(f'gYieldPromptSystotal_{reso}_MB')

# Ground state
if use_mb_cross:
    infile_gs_MB = TFile.Open('input/cross_section_and_yields/MB/yieldDs-PtIntegr2-24_pp13TeV_MB.root') # from Stefano, using the same macro to integrate the cross section scaling the Ds cross section by (sigmaV0M*BR)
    hgs_MB = infile_gs_MB.Get('hPtIntCrossSecStat')
    ggs_MB = infile_gs_MB.Get('gPtIntCrossSecDataSyst')
else:
    infile_gs_MB = TFile.Open('input/cross_section_and_yields/HM/V0M_ptintegratedyields_D0DsLc_pt2-24.root') # from Stefano, using the same macro to integrate the cross section scaling the Ds cross section by (sigmaV0M*BR)
    hgs_MB = infile_gs_MB.Get('h_vis_dNdy_1_0')
    ggs_MB = infile_gs_MB.Get('gr_vis_dNdy_TotSyst_1_0')
infile_gs_HM = TFile.Open('input/cross_section_and_yields/HM/V0M_ptintegratedyields_D0DsLc_pt2-24.root') # from Luuk
hgs_HM = infile_gs_HM.Get('h_vis_dNdy_1_4')
ggs_HM = infile_gs_HM.Get('gr_vis_dNdy_TotSyst_1_4')

#_____________________________________________________________
# Scale for deltapt
hgs_HM.Scale(1./delta_pt)
ggs_HM.SetPoint(0, ggs_HM.GetX()[0], ggs_HM.GetY()[0]*1/delta_pt)
ggs_HM.SetPointError(0, ggs_HM.GetErrorXlow(0),
                     ggs_HM.GetErrorXhigh(0),
                     ggs_HM.GetErrorYlow(0)*1/delta_pt,
                     ggs_HM.GetErrorYhigh(0)*1/delta_pt)
hgs_MB.Scale(1./delta_pt)
ggs_MB.SetPoint(0, ggs_MB.GetX()[0], ggs_MB.GetY()[0]*1/delta_pt)
ggs_MB.SetPointError(0, ggs_MB.GetErrorXlow(0),
                     ggs_MB.GetErrorXhigh(0),
                     ggs_MB.GetErrorYlow(0)*1/delta_pt,
                     ggs_MB.GetErrorYhigh(0)*1/delta_pt)

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
ratio_MB = num/den
ratio_MB_stat = ratio_MB * np.sqrt((num_unc/num)**2 + (den_unc/den)**2)
ratio_MB_sys_low = ratio_MB * np.sqrt((greso_MB_syst.GetErrorYlow(0)/num)**2 + (ggs_MB.GetErrorYlow(0)/den)**2)
ratio_MB_sys_high = ratio_MB * np.sqrt((greso_MB_syst.GetErrorYhigh(0)/num)**2 + (ggs_MB.GetErrorYhigh(0)/den)**2)
print(f'Reso integrated yield = {num} +/- {num_unc}')
print(f'Ds integrated yield = {den} +/- {den_unc}')
print(f'ratio_MB = {ratio_MB} +/- {ratio_MB_stat}')
print(f'ratio_MB_sys_low = {ratio_MB_sys_low} + {ratio_MB_sys_high}')

# HM
print('\033[1m\033[4mHM\033[0m')
num = hreso_HM.GetBinContent(1)
num_unc = hreso_HM.GetBinError(1)
den = hgs_HM.GetBinContent(1)
den_unc = hgs_HM.GetBinError(1)
ratio_HM = num/den
ratio_HM_stat = ratio_HM * np.sqrt((num_unc/num)**2 + (den_unc/den)**2)
ratio_HM_sys_low = ratio_HM * np.sqrt((greso_HM_syst.GetErrorYlow(0)/num)**2 + (ggs_HM.GetErrorYlow(0)/den)**2)
ratio_HM_sys_high = ratio_HM * np.sqrt((greso_HM_syst.GetErrorYhigh(0)/num)**2 + (ggs_HM.GetErrorYhigh(0)/den)**2)
print(f'Reso integrated yield = {num} +/- {num_unc}')
print(f'Ds integrated yield = {den} +/- {den_unc}')
print(f'ratio_HM = {ratio_HM} +/- {ratio_HM_stat}')
print(f'ratio_HM_sys_low = {ratio_HM_sys_low} + {ratio_HM_sys_high}')

gratio.SetPoint(0, 7, ratio_MB)
gratio.SetPointError(0, 0, 0, ratio_MB_stat, ratio_MB_stat)
gratio_sys.SetPoint(0, 7, ratio_MB)
gratio_sys.SetPointError(0, 1., 1., ratio_MB_sys_low, ratio_MB_sys_high)
gratio.SetPoint(1, 31, ratio_HM)
gratio.SetPointError(1, 0, 0, ratio_HM_stat, ratio_HM_stat)
gratio_sys.SetPoint(1, 31, ratio_HM)
gratio_sys.SetPointError(1, 1., 1., ratio_HM_sys_low, ratio_HM_sys_high)
SetObjectStyle(gratio, color=kRed+1, markerstyle=20, markersize=1.2,
               linecolor=kRed+1, markercolor=kRed+1, fillalpha=0.2, fillstyle=0,
               linewidth=2)
SetObjectStyle(gratio_sys, color=kRed+1, markerstyle=20, markersize=1.2,
               linecolor=kRed+1, markercolor=kRed+1, fillalpha=0.2, fillstyle=0,
               linewidth=2)

#_____________________________________________________________
# Collect theoretical prediction
br_reso = 23.35e-2 if reso == '435' else 22.e-2 # taken from AN (see also https://indico.cern.ch/event/1253172/)
if reso == '435':
    pred_shm = 0.075 # taken from ./predictions/MinHe/
    pred_shmc = 0.034/0.357 # taken from Ds2*+/D0 = 0.034 and Ds / D0 = 0.357
    pred_shm *= br_reso
    pred_shmc *= br_reso
else:
    pred_shm = 0.054 # taken from ./predictions/MinHe/
    pred_shmc = 0.026/0.357 # taken from Ds1+/D0 = 0.026 and Ds / D0 = 0.357
    pred_shm *= br_reso
    pred_shmc *= br_reso

gtheo_shm, gtheo_shmc = TGraph(), TGraph()
gtheo_shm.SetPoint(0, 7, pred_shm)
gtheo_shm.SetPoint(1, 31, pred_shm) # no explicit multiplicity dependence indicated from theoreticians
gtheo_shmc.SetPoint(0, 7, pred_shmc)
gtheo_shmc.SetPoint(1, 31, pred_shmc) # no explicit multiplicity dependence indicated from theoreticians
SetObjectStyle(gtheo_shm, color=kAzure+4, markerstyle=kOpenCross, markersize=1.5,
               linewidth=2, fillalpha=0.5, fillcolor=kAzure+4, linestyle=7)
SetObjectStyle(gtheo_shmc, color=kOrange+1, markerstyle=kOpenDiamond, markersize=1.5,
                linewidth=2, fillalpha=0.5, fillcolor=kOrange+1, linestyle=3)

#_____________________________________________________________
# nSigma calculation
nsigma_mb_SHM = (ratio_MB - gtheo_shmc.GetY()[0])/ratio_MB_stat
nsigma_mb_SHM_sys = (ratio_MB - gtheo_shmc.GetY()[0])/np.sqrt(ratio_MB_stat**2 + ratio_MB_sys_low**2)
nsigma_mb_He = (ratio_MB - gtheo_shm.GetY()[0])/ratio_MB_stat
nsigma_mb_He_sys = (ratio_MB - gtheo_shm.GetY()[0])/np.sqrt(ratio_MB_stat**2 + ratio_MB_sys_low**2)
nsigma_hm_SHM = (ratio_HM - gtheo_shmc.GetY()[0])/ratio_HM_stat
nsigma_hm_SHM_sys = (ratio_HM - gtheo_shmc.GetY()[0])/np.sqrt(ratio_HM_stat**2 + ratio_HM_sys_low**2)
nsigma_hm_He = (ratio_HM - gtheo_shm.GetY()[0])/ratio_HM_stat
nsigma_hm_He_sys = (ratio_HM - gtheo_shm.GetY()[0])/np.sqrt(ratio_HM_stat**2 + ratio_HM_sys_low**2)

print('\n\033[1m\033[4mnSigma (only stat)\033[0m')
print(f'nsigma_mb_SHM = {nsigma_mb_SHM} +/- {nsigma_mb_SHM_sys}')
print(f'nsigma_mb_He = {nsigma_mb_He} +/- {nsigma_mb_He_sys}')
print(f'nsigma_hm_SHM = {nsigma_hm_SHM} +/- {nsigma_hm_SHM_sys}')
print(f'nsigma_hm_He = {nsigma_hm_He} +/- {nsigma_hm_He_sys}')


#_____________________________________________________________
# Plot and save output
canvas = TCanvas('canvas', 'canvas', 850, 800)
#canvas.cd().SetLogy()
ytitle = ';#LTd#it{N}_{ch}/d#eta#GT_{|#eta| < 0.5}; D_{s2}^{*+}#timesBR(D_{s2}^{*+} #rightarrow D^{+}K^{0}_{S}) / D_{s}^{+}' if reso == '435' else ';#LTd#it{N}_{ch}/d#eta#GT_{|#eta| < 0.5}; D_{s1}^{+}#timesBR(D_{s1}^{+} #rightarrow D^{+}K^{0}_{S}) / D_{s}^{+}'
ymin = 0. if reso == '435' else 0.
ymax = 0.18 if reso == '435' else 0.08
hFrame = canvas.cd().DrawFrame(0., ymin, 40, ymax, ytitle)
gtheo_shm.Draw('csame')
gtheo_shmc.Draw('csame')
gratio_sys.Draw('5 samez')
gratio.Draw('same pz')
leg = TLegend(0.45, 0.5, 0.65, 0.75)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.SetTextFont(42)
leg.SetTextSize(0.035)
leg.AddEntry(gratio, 'Data (2.0 < #it{p}_{T} < 24.0 GeV/#it{c})', 'p')
leg.AddEntry(gtheo_shm, 'SHM M. He, R. Rapp (#it{p}_{T} > 0)', 'l')
leg.AddEntry(gtheo_shmc, 'SHMc GSI-Heidelberg (#it{p}_{T} > 0)', 'l')
leg.Draw()
lat = TLatex()
lat.SetTextSize(0.04)
lat.SetTextFont(42)
lat.DrawLatexNDC(0.25, 0.88, 'ALICE Preliminary')
lat.DrawLatexNDC(0.25, 0.8, 'pp, #sqrt{#it{s}} = 13 TeV')

canvas.SaveAs(f'Ds{reso}_over_Ds_ratio_vs_mult{outlabel}.pdf')
outfile = TFile(f'Ds{reso}_over_Ds_ratio_vs_mult{outlabel}.root', 'recreate')
canvas.Write()
gratio.Write()
gratio_sys.Write()
gtheo_shm.Write()
gtheo_shmc.Write()
outfile.Close()
input()
