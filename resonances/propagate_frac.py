"""
Script for fraction propagation
"""

import os
import sys
import argparse
import ctypes
import numpy as np
import uproot
import ROOT

sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.AnalysisUtils import GetPromptFDFractionCutSet

# pylint: disable=too-many-arguments,too-many-locals,too-many-statements
def propagate(cutvarfile_name, efffile_name, kinefile_name, beautyhypofile_name,
              output_dir, pt_min_reso, pt_max_reso, suffix, seed):
    """
    function for the propagation of the non-prompt fraction
    from the D daughter to the D resonance mother

    Parameters
    ----------

    - cutvarfile_name (str): input ROOT file with the cut-variation output
    - efffile_name (str): input ROOT file for the eff x acc factors
    - kinefile_name (str): input ROOT file with D* -> D kinematics
    - beautyhypofile_name (str): input ROOT file with Ds/D Np/p ratio for beauty hypo
    - output_dir (str): output directory
    - pt_min_reso (float): min pT for pT integrated Nnp/Np
    - pt_max_reso (float): max pT for pT integrated Nnp/Np
    - suffix (str): suffix for output file
    - seed (int): seed for random sampling
    """

    SetGlobalStyle(padleftmargin=0.18, padrightmargin=0.1, padbottommargin=0.14)

    ROOT.gROOT.SetBatch(True)
    ROOT.gRandom.SetSeed(seed)

    if 'Ds2starplus' in kinefile_name:
        reso_name = 'Ds2starplus'
        d_name = "Dplus"
    elif 'Ds1plus' in kinefile_name:
        reso_name = 'Ds1plus'
        d_name = "Dstar"
    else:
        print(f'\033[1m\033[91mResonance not supported\033[0m')
        sys.exit()

    ptshape = "Dsptshape"
    if "Lcptshape" in kinefile_name:
        ptshape = "Lcptshape"
    elif "DsHarderptshape" in kinefile_name:
        ptshape = "DsHarderptshape"
    elif "DsVeryHarderptshape" in kinefile_name:
        ptshape = "DsVeryHarderptshape"
    elif "DsSofterptshape" in kinefile_name:
        ptshape = "DsSofterptshape"
    elif "DsVerySofterptshape" in kinefile_name:
        ptshape = "DsVerySofterptshape"

    trigger = ""
    if "MB" in kinefile_name:
        trigger = "MB"
    elif "HM" in kinefile_name:
        trigger = "HM"

    mult_weights_suffix = ""
    if "_multweights_all" in efffile_name:
        mult_weights_suffix = "_multweights_all"
    elif "_multweights_cand" in efffile_name:
        mult_weights_suffix = "_multweights_cand"
    elif "_multweights_candinmass" in efffile_name:
        mult_weights_suffix = "_multweights_candinmass"

    infile_cutvar = ROOT.TFile.Open(cutvarfile_name)
    hist_corryield_prompt = infile_cutvar.Get('hCorrYieldPrompt')
    hist_corryield_nonprompt = infile_cutvar.Get('hCorrYieldFD')
    hist_cov_pp = infile_cutvar.Get('hCovPromptPrompt')
    hist_cov_pnp = infile_cutvar.Get('hCovPromptFD')
    hist_cov_npnp = infile_cutvar.Get('hCovFDFD')
    hist_corryield_prompt.SetDirectory(0)
    hist_corryield_nonprompt.SetDirectory(0)
    hist_cov_pp.SetDirectory(0)
    hist_cov_pnp.SetDirectory(0)
    hist_cov_npnp.SetDirectory(0)
    infile_cutvar.Close()

    infile_eff = ROOT.TFile.Open(efffile_name)
    hist_reco_prompt = infile_eff.Get(f"{d_name}_Prompt/h_pt_reco_{d_name}_Prompt")
    hist_gen_prompt = infile_eff.Get(f"{d_name}_Prompt/h_pt_gen_{d_name}_Prompt")
    hist_reco_nonprompt = infile_eff.Get(f"{d_name}_NonPrompt/h_pt_reco_{d_name}_NonPrompt")
    hist_gen_nonprompt = infile_eff.Get(f"{d_name}_NonPrompt/h_pt_gen_{d_name}_NonPrompt")
    nptbins = hist_corryield_prompt.GetNbinsX()
    ptlims = []
    for ipt in range(1, nptbins + 1):
        ptlims.append(hist_corryield_prompt.GetBinLowEdge(ipt))
    ptlims.append(hist_corryield_prompt.GetXaxis().GetBinUpEdge(nptbins))
    ptlims = np.array(ptlims)

    hist_reco_prompt_reb = hist_reco_prompt.Rebin(nptbins, "hist_reco_prompt_reb", ptlims)
    hist_gen_prompt_reb = hist_gen_prompt.Rebin(nptbins, "hist_gen_prompt_reb", ptlims)
    hist_reco_nonprompt_reb = hist_reco_nonprompt.Rebin(nptbins, "hist_reco_nonprompt_reb", ptlims)
    hist_gen_nonprompt_reb = hist_gen_nonprompt.Rebin(nptbins, "hist_gen_nonprompt_reb", ptlims)

    hist_acceff_prompt = hist_reco_prompt_reb.Clone("hist_acceff_prompt")
    hist_acceff_prompt.Divide(hist_reco_prompt_reb, hist_gen_prompt_reb, 1., 1., "B")
    hist_acceff_nonprompt = hist_reco_nonprompt_reb.Clone("hist_acceff_prompt")
    hist_acceff_nonprompt.Divide(hist_reco_nonprompt_reb, hist_gen_nonprompt_reb, 1., 1., "B")
    hist_acceff_prompt.SetDirectory(0)
    hist_acceff_nonprompt.SetDirectory(0)
    infile_eff.Close()

    hist_npfrac_d = hist_corryield_nonprompt.Clone("hist_npfrac_d")
    hist_npfrac_d.SetDirectory(0)
    SetObjectStyle(hist_npfrac_d, color=ROOT.kAzure+4)

    for ipt in range(1, hist_npfrac_d.GetNbinsX()+1):
        frac, unc_frac = GetPromptFDFractionCutSet(
            hist_acceff_prompt.GetBinContent(ipt),
            hist_acceff_nonprompt.GetBinContent(ipt),
            hist_corryield_prompt.GetBinContent(ipt),
            hist_corryield_nonprompt.GetBinContent(ipt),
            hist_cov_pp.GetBinContent(ipt),
            hist_cov_npnp.GetBinContent(ipt),
            hist_cov_pnp.GetBinContent(ipt)
        )
        hist_npfrac_d.SetBinContent(ipt, frac[1])
        hist_npfrac_d.SetBinError(ipt, unc_frac[1])

    pt_mins, pt_maxs = [], []
    for ipt in range(1, hist_npfrac_d.GetNbinsX()+1):
        pt_mins.append(hist_npfrac_d.GetXaxis().GetBinLowEdge(ipt))
        pt_maxs.append(hist_npfrac_d.GetXaxis().GetBinUpEdge(ipt))
    pt_lims = pt_mins.copy()
    pt_lims.append(pt_maxs[-1])
    nptbins = len(pt_mins)

    frac_hypo_cent = 2.
    frac_hypo_min = 1.
    frac_hypo_max = 3.
    if beautyhypofile_name != "":
        infile_beautyhypo = ROOT.TFile.Open(beautyhypofile_name)
        graph_double_ratio_dsd_npp = infile_beautyhypo.Get("double_ratio_all")
        x, y = ctypes.c_double(), ctypes.c_double()
        graph_double_ratio_dsd_npp.GetPoint(0, x, y)
        frac_hypo_cent = y.value
        frac_hypo_min = frac_hypo_cent - graph_double_ratio_dsd_npp.GetErrorYlow(0)
        frac_hypo_max = frac_hypo_cent + graph_double_ratio_dsd_npp.GetErrorYhigh(0)

    hist_npfrac_dreso_vspt_2d = {}
    hist_npfrac_dreso_ptint_distr = {}
    for hist_suffix in ["cent", "hypo_min", "hypo_max", "stat_min", "stat_max"]:
        hist_npfrac_dreso_vspt_2d[hist_suffix] = ROOT.TH2F(
            f"hist_npfrac_dreso_vspt_2d_{hist_suffix}",
            ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}",
            nptbins, np.array(pt_lims), 1000, 0., 1.)
        hist_npfrac_dreso_ptint_distr[hist_suffix] = ROOT.TH1F(
            f"hist_npfrac_dreso_ptint_distr{hist_suffix}",
            ";#it{f}_{non-prompt};entries", 1000, 0., 1.)

    hist_pt_corr = ROOT.TH2F("hist_pt_corr",
                             ";#it{p}_{T}^{reso} (GeV/#it{c});#it{p}_{T}^{D} (GeV/#it{c})",
                             int(pt_max_reso) * 2, 0., pt_max_reso,
                             int(pt_max_reso) * 2, 0., pt_max_reso)

    # prompt and non-prompt should be very similar in terms of pT shape
    df_kine = uproot.open(kinefile_name)["recoTreePrompt"].arrays(library="pd")

    for pt_reso, pt_d in zip(df_kine["pt_reso"].to_numpy(), df_kine["pt_d"].to_numpy()):
        hist_pt_corr.Fill(pt_reso, pt_d)
        pt_bin = hist_npfrac_d.GetXaxis().FindBin(pt_d)

        frac_d_reso_cent, frac_d_reso_centunc = GetPromptFDFractionCutSet(
            hist_acceff_prompt.GetBinContent(pt_bin),
            hist_acceff_nonprompt.GetBinContent(pt_bin),
            hist_corryield_prompt.GetBinContent(pt_bin),
            hist_corryield_nonprompt.GetBinContent(pt_bin) * frac_hypo_cent, # multiplication factor enters here
            hist_cov_pp.GetBinContent(pt_bin),
            hist_cov_npnp.GetBinContent(pt_bin) * frac_hypo_cent**2,
            hist_cov_pnp.GetBinContent(pt_bin) * frac_hypo_cent
        )
        frac_d_reso_cent = frac_d_reso_cent[1]
        frac_d_reso_stat_min = frac_d_reso_cent - frac_d_reso_centunc[1]
        frac_d_reso_stat_max = frac_d_reso_cent + frac_d_reso_centunc[1]

        frac_d_reso_hypo_min, _ = GetPromptFDFractionCutSet(
            hist_acceff_prompt.GetBinContent(pt_bin),
            hist_acceff_nonprompt.GetBinContent(pt_bin),
            hist_corryield_prompt.GetBinContent(pt_bin),
            hist_corryield_nonprompt.GetBinContent(pt_bin) * frac_hypo_min, # multiplication factor enters here
            hist_cov_pp.GetBinContent(pt_bin),
            hist_cov_npnp.GetBinContent(pt_bin) * frac_hypo_cent**2,
            hist_cov_pnp.GetBinContent(pt_bin) * frac_hypo_cent
        )
        frac_d_reso_hypo_min = frac_d_reso_hypo_min[1]

        frac_d_reso_hypo_max, _ = GetPromptFDFractionCutSet(
            hist_acceff_prompt.GetBinContent(pt_bin),
            hist_acceff_nonprompt.GetBinContent(pt_bin),
            hist_corryield_prompt.GetBinContent(pt_bin),
            hist_corryield_nonprompt.GetBinContent(pt_bin) * frac_hypo_max, # multiplication factor enters here
            hist_cov_pp.GetBinContent(pt_bin),
            hist_cov_npnp.GetBinContent(pt_bin) * frac_hypo_cent**2,
            hist_cov_pnp.GetBinContent(pt_bin) * frac_hypo_cent
        )
        frac_d_reso_hypo_max = frac_d_reso_hypo_max[1]

        hist_npfrac_dreso_vspt_2d["cent"].Fill(pt_reso, frac_d_reso_cent)
        hist_npfrac_dreso_vspt_2d["hypo_min"].Fill(pt_reso, frac_d_reso_hypo_min)
        hist_npfrac_dreso_vspt_2d["hypo_max"].Fill(pt_reso, frac_d_reso_hypo_max)
        hist_npfrac_dreso_vspt_2d["stat_min"].Fill(pt_reso, frac_d_reso_stat_min)
        hist_npfrac_dreso_vspt_2d["stat_max"].Fill(pt_reso, frac_d_reso_stat_max)
        if pt_min_reso < pt_reso < pt_max_reso:
            hist_npfrac_dreso_ptint_distr["cent"].Fill(frac_d_reso_cent)
            hist_npfrac_dreso_ptint_distr["hypo_min"].Fill(frac_d_reso_hypo_min)
            hist_npfrac_dreso_ptint_distr["hypo_max"].Fill(frac_d_reso_hypo_max)
            hist_npfrac_dreso_ptint_distr["stat_min"].Fill(frac_d_reso_stat_min)
            hist_npfrac_dreso_ptint_distr["stat_max"].Fill(frac_d_reso_stat_max)

    graph_npfrac_dreso_stat = ROOT.TGraphAsymmErrors(1)
    graph_npfrac_dreso_vspt_stat = ROOT.TGraphAsymmErrors(nptbins)
    graph_npfrac_dreso_hypo = ROOT.TGraphAsymmErrors(1)
    graph_npfrac_dreso_vspt_hypo = ROOT.TGraphAsymmErrors(nptbins)
    graph_pfrac_dreso_stat = ROOT.TGraphAsymmErrors(1)
    graph_pfrac_dreso_vspt_stat = ROOT.TGraphAsymmErrors(nptbins)
    graph_pfrac_dreso_hypo = ROOT.TGraphAsymmErrors(1)
    graph_pfrac_dreso_vspt_hypo = ROOT.TGraphAsymmErrors(nptbins)
    SetObjectStyle(graph_npfrac_dreso_stat, color=ROOT.kRed+1)
    SetObjectStyle(graph_npfrac_dreso_vspt_stat, color=ROOT.kRed+1)
    SetObjectStyle(graph_npfrac_dreso_hypo, color=ROOT.kRed+1, fillstyle=0)
    SetObjectStyle(graph_npfrac_dreso_vspt_hypo, color=ROOT.kRed+1, fillstyle=0)
    SetObjectStyle(graph_pfrac_dreso_stat, color=ROOT.kRed+1)
    SetObjectStyle(graph_pfrac_dreso_vspt_stat, color=ROOT.kRed+1)
    SetObjectStyle(graph_pfrac_dreso_hypo, color=ROOT.kRed+1, fillstyle=0)
    SetObjectStyle(graph_pfrac_dreso_vspt_hypo, color=ROOT.kRed+1, fillstyle=0)
    graph_npfrac_dreso_stat.SetNameTitle("graph_npfrac_dreso_stat_ptint",
                                         ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}")
    graph_npfrac_dreso_vspt_stat.SetNameTitle("graph_npfrac_dreso_vspt_stat",
                                              ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}")
    graph_npfrac_dreso_hypo.SetNameTitle("graph_npfrac_dreso_hypo_ptint",
                                         ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}")
    graph_npfrac_dreso_vspt_hypo.SetNameTitle("graph_npfrac_dreso_vspt_hypo",
                                              ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}")
    graph_pfrac_dreso_stat.SetNameTitle("graph_pfrac_dreso_stat_ptint",
                                        ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{prompt}")
    graph_pfrac_dreso_vspt_stat.SetNameTitle("graph_pfrac_dreso_vspt_stat",
                                             ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{prompt}")
    graph_pfrac_dreso_hypo.SetNameTitle("graph_pfrac_dreso_hypo_ptint",
                                        ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{prompt}")
    graph_pfrac_dreso_vspt_hypo.SetNameTitle("graph_pfrac_dreso_vspt_hypo",
                                             ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{prompt}")

    npfrac_cent = hist_npfrac_dreso_ptint_distr["cent"].GetMean()
    graph_npfrac_dreso_stat.SetPoint(0, (pt_min_reso + pt_max_reso) / 2, npfrac_cent)
    graph_npfrac_dreso_hypo.SetPoint(0, (pt_min_reso + pt_max_reso) / 2, npfrac_cent)
    graph_pfrac_dreso_stat.SetPoint(0, (pt_min_reso + pt_max_reso) / 2, 1 - npfrac_cent)
    graph_pfrac_dreso_hypo.SetPoint(0, (pt_min_reso + pt_max_reso) / 2, 1 - npfrac_cent)
    pt_unc = (pt_max_reso - pt_min_reso) / 2
    frac_unc_stat_high = hist_npfrac_dreso_ptint_distr["stat_max"].GetMean() - npfrac_cent
    frac_unc_stat_low = npfrac_cent - hist_npfrac_dreso_ptint_distr["stat_min"].GetMean()
    frac_unc_hypo_high = hist_npfrac_dreso_ptint_distr["hypo_max"].GetMean() - npfrac_cent
    frac_unc_hypo_low = npfrac_cent - hist_npfrac_dreso_ptint_distr["hypo_min"].GetMean()
    graph_npfrac_dreso_stat.SetPointError(0, pt_unc, pt_unc, frac_unc_stat_low, frac_unc_stat_high)
    graph_npfrac_dreso_hypo.SetPointError(0, pt_unc, pt_unc, frac_unc_hypo_low, frac_unc_hypo_high)
    graph_pfrac_dreso_stat.SetPointError(0, pt_unc, pt_unc, frac_unc_stat_high, frac_unc_stat_low)
    graph_pfrac_dreso_hypo.SetPointError(0, pt_unc, pt_unc, frac_unc_hypo_high, frac_unc_hypo_low)

    profile_npfrac_dreso_vspt = {}
    for var in hist_npfrac_dreso_vspt_2d:
        profile_npfrac_dreso_vspt[var] = hist_npfrac_dreso_vspt_2d[var].ProfileX()

    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        npfrac_cent = profile_npfrac_dreso_vspt["cent"].GetBinContent(ipt+1)
        graph_npfrac_dreso_vspt_stat.SetPoint(ipt, (pt_min + pt_max) / 2, npfrac_cent)
        graph_npfrac_dreso_vspt_hypo.SetPoint(ipt, (pt_min + pt_max) / 2, npfrac_cent)
        graph_pfrac_dreso_vspt_stat.SetPoint(ipt, (pt_min + pt_max) / 2, 1 - npfrac_cent)
        graph_pfrac_dreso_vspt_hypo.SetPoint(ipt, (pt_min + pt_max) / 2, 1 - npfrac_cent)
        pt_unc = (pt_max - pt_min) / 2
        frac_unc_stat_high = profile_npfrac_dreso_vspt["stat_max"].GetBinContent(ipt+1) - npfrac_cent
        frac_unc_stat_low = npfrac_cent - profile_npfrac_dreso_vspt["stat_min"].GetBinContent(ipt+1)
        frac_unc_hypo_high = profile_npfrac_dreso_vspt["hypo_max"].GetBinContent(ipt+1) - npfrac_cent
        frac_unc_hypo_low = npfrac_cent - profile_npfrac_dreso_vspt["hypo_min"].GetBinContent(ipt+1)
        graph_npfrac_dreso_vspt_stat.SetPointError(ipt, pt_unc, pt_unc, frac_unc_stat_low, frac_unc_stat_high)
        graph_npfrac_dreso_vspt_hypo.SetPointError(ipt, pt_unc, pt_unc, frac_unc_hypo_low, frac_unc_hypo_high)
        graph_pfrac_dreso_vspt_stat.SetPointError(ipt, pt_unc, pt_unc, frac_unc_stat_high, frac_unc_stat_low)
        graph_pfrac_dreso_vspt_hypo.SetPointError(ipt, pt_unc, pt_unc, frac_unc_hypo_high, frac_unc_hypo_low)

    canv_npfrac_vspt = ROOT.TCanvas("canv_npfrac_vspt", "", 1000, 500)
    leg = ROOT.TLegend(0.2, 0.75, 0.6, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.045)
    leg.AddEntry(hist_npfrac_d, "D daughter", "pl")
    leg.AddEntry(graph_npfrac_dreso_vspt_stat, "D resonance", "pl")
    canv_npfrac_vspt.Divide(2, 1)
    canv_npfrac_vspt.cd(1).SetLogz()
    hist_pt_corr.Draw("colz")
    canv_npfrac_vspt.cd(2).DrawFrame(pt_mins[0], 0., pt_maxs[-1], 1.,
                                     ";#it{p}_{T} (GeV/#it{c});#it{f}_{non-prompt}")
    graph_npfrac_dreso_vspt_hypo.Draw("2")
    graph_npfrac_dreso_vspt_stat.Draw("pz")
    hist_npfrac_d.Draw("same")
    leg.Draw()

    canv_frac_ptint = ROOT.TCanvas("canv_frac_ptint", "", 1000, 500)
    canv_frac_ptint.Divide(2, 1)
    canv_frac_ptint.cd(1).DrawFrame(pt_min_reso, 0., pt_max_reso, 1.,
                                    ";#it{p}_{T} (GeV/#it{c});#it{f}_{non-prompt}")
    graph_npfrac_dreso_hypo.Draw("2")
    graph_npfrac_dreso_stat.Draw("pz")
    canv_frac_ptint.cd(2).DrawFrame(pt_min_reso, 0., pt_max_reso, 1.,
                                    ";#it{p}_{T} (GeV/#it{c});#it{f}_{prompt}")
    graph_pfrac_dreso_hypo.Draw("2")
    graph_pfrac_dreso_stat.Draw("pz")

    outfile_name = os.path.join(
        output_dir,
        f"fraction_{reso_name}_{trigger}_{ptshape}{mult_weights_suffix}_propagated{suffix}.root"
    )
    canv_npfrac_vspt.SaveAs(outfile_name.replace(".root", "_vspt.pdf"))
    canv_frac_ptint.SaveAs(outfile_name.replace(".root", "_ptint.pdf"))

    outfile = ROOT.TFile(outfile_name, "recreate")
    canv_npfrac_vspt.Write()
    canv_frac_ptint.Write()
    hist_pt_corr.Write()
    hist_npfrac_d.SetName("hist_npfrac_d")
    hist_npfrac_d.Write()
    for var in hist_npfrac_dreso_vspt_2d:
        hist_npfrac_dreso_vspt_2d[var].Write()
        hist_npfrac_dreso_ptint_distr[var].Write()
        profile_npfrac_dreso_vspt[var].Write()
    graph_npfrac_dreso_stat.Write()
    graph_npfrac_dreso_vspt_stat.Write()
    graph_npfrac_dreso_hypo.Write()
    graph_npfrac_dreso_vspt_hypo.Write()
    graph_pfrac_dreso_stat.Write()
    graph_pfrac_dreso_vspt_stat.Write()
    graph_pfrac_dreso_hypo.Write()
    graph_pfrac_dreso_vspt_hypo.Write()
    outfile.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("cutvarfile", metavar="text",
                        default="cutvar.root", help="input ROOT file for cut-variation")
    parser.add_argument("efffile", metavar="text",
                        default="eff.root", help="input ROOT file for eff x acc")
    parser.add_argument("kinefile", metavar="text",
                        default="kine.root", help="ROOT file with kinematics")
    parser.add_argument("--file_beauty_hypo", "-b", metavar="text",
                        default="",
                        help="file with hypothesis for strange/nonstrange Np/p ratio")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--pt_min_reso", type=float,
                        default=2., help="min pT for pT integrated fraction")
    parser.add_argument("--pt_max_reso", type=float,
                        default=24., help="max pT for pT integrated fraction")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--seed", type=int,
                        default=42, help="seed for random sampling")
    args = parser.parse_args()

    propagate(args.cutvarfile, args.efffile, args.kinefile,
              args.file_beauty_hypo, args.outputdir,
              args.pt_min_reso, args.pt_max_reso, args.suffix, args.seed)
