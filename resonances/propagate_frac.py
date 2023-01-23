import sys
import argparse
import numpy as np
import uproot
import ROOT

sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

def propagate(infile_name, outfile_name, kinefile_name,
              frac_hypo_cent, frac_hypo_min, frac_hypo_max,
              pt_min_reso, pt_max_reso):
    """
    function for the propagation of the non-prompt fraction
    from the D daughter to the D resonance mother

    Parameters
    ----------

    - infile_name (str): input ROOT file
    - outfile_name (str): output ROOT file
    - kinefile_name (str): file with D* -> D kinematics
    - frac_hypo_cent (float): central hypothesis for Dreso/D fraction ratio
    - frac_hypo_min (float): min hypothesis for Dreso/D fraction ratio
    - frac_hypo_max (float): max hypothesis for Dreso/D fraction ratio
    - pt_min_reso (float): min pT for pT integrated fraction
    - pt_max_reso (float): max pT for pT integrated fraction
    """

    ROOT.gROOT.SetBatch()
    SetGlobalStyle(padleftmargin=0.18, padrightmargin=0.1, padbottommargin=0.14)

    infile_npfrac = ROOT.TFile.Open(infile_name)
    hist_npfrac_d = infile_npfrac.Get("hFDFrac")
    hist_npfrac_d.SetDirectory(0)
    SetObjectStyle(hist_npfrac_d, color=ROOT.kAzure+4, markerstyle=ROOT.kFullSquare)
    spl_npfrac_d = ROOT.TSpline3(hist_npfrac_d)
    pt_mins, pt_maxs = [], []
    for ipt in range(1, hist_npfrac_d.GetNbinsX()+1):
        pt_mins.append(hist_npfrac_d.GetXaxis().GetBinLowEdge(ipt))
        pt_maxs.append(hist_npfrac_d.GetXaxis().GetBinUpEdge(ipt))
    pt_lims = pt_mins.copy()
    pt_lims.append(pt_maxs[-1])
    nptbins = len(pt_mins)
    infile_npfrac.Close()

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

    df_kine = uproot.open(kinefile_name)["tupleDreso"].arrays(library="pd")
    for pt_reso, pt_d in zip(df_kine["pt_reso"].to_numpy(), df_kine["pt_d"].to_numpy()):
        hist_pt_corr.Fill(pt_reso, pt_d)
        pt_bin = hist_npfrac_d.GetXaxis().FindBin(pt_d)
        frac_d = spl_npfrac_d.Eval(pt_d) #hist_npfrac_d.GetBinContent(pt_bin)
        frac_d_statunc = hist_npfrac_d.GetBinError(pt_bin)
        frac_d_reso_cent = frac_d * frac_hypo_cent
        frac_d_reso_hypo_min = frac_d * frac_hypo_min
        frac_d_reso_hypo_max = frac_d * frac_hypo_max
        frac_d_reso_stat_min = (frac_d - frac_d_statunc) * frac_hypo_cent
        frac_d_reso_stat_max = (frac_d + frac_d_statunc) * frac_hypo_cent
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
    graph_npfrac_dreso_stat.SetNameTitle(f"graph_npfrac_dreso_stat_pt{pt_min_reso}-{pt_max_reso}",
                                       ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}")
    graph_npfrac_dreso_vspt_stat.SetNameTitle("graph_npfrac_dreso_vspt_stat",
                                            ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}")
    graph_npfrac_dreso_hypo.SetNameTitle(f"graph_npfrac_dreso_hypo_pt{pt_min_reso}-{pt_max_reso}",
                                       ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}")
    graph_npfrac_dreso_vspt_hypo.SetNameTitle("graph_npfrac_dreso_vspt_hypo",
                                            ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{non-prompt}")
    graph_pfrac_dreso_stat.SetNameTitle(f"graph_pfrac_dreso_stat_pt{pt_min_reso}-{pt_max_reso}",
                                        ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{prompt}")
    graph_pfrac_dreso_vspt_stat.SetNameTitle("graph_pfrac_dreso_vspt_stat",
                                             ";#it{p}_{T}^{reso} (GeV/#it{c}); #it{f}_{prompt}")
    graph_pfrac_dreso_hypo.SetNameTitle(f"graph_pfrac_dreso_hypo_pt{pt_min_reso}-{pt_max_reso}",
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
    parser.add_argument("inputfile", metavar="text",
                        default="fraction_D.root", help="input ROOT file")
    parser.add_argument("kinefile", metavar="text",
                        default="kine.root", help="ROOT file with kinematics")
    parser.add_argument("outputfile", metavar="text",
                        default="fraction_reso.root", help="output ROOT file")
    parser.add_argument("--f_cent", type=float,
                        default=2., help="central hypothesis for Dreso/D fraction ratio")
    parser.add_argument("--f_min", type=float,
                        default=1., help="min hypothesis for Dreso/D fraction ratio")
    parser.add_argument("--f_max", type=float,
                        default=3., help="max hypothesis for Dreso/D fraction ratio")
    parser.add_argument("--pt_min_reso", type=float,
                        default=2., help="min pT for pT integrated fraction")
    parser.add_argument("--pt_max_reso", type=float,
                        default=24., help="max pT for pT integrated fraction")
    args = parser.parse_args()

    propagate(args.inputfile, args.outputfile, args.kinefile,
              args.f_cent, args.f_min, args.f_max,
              args.pt_min_reso, args.pt_max_reso)
