import os
import argparse
import yaml
import matplotlib.pyplot as plt
import uproot
from hist import Hist
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from particle import Particle

def project(config, trigger, pdg_d, pdg_v0):
    """
    function for raw yield computation

    Parameters
    ----------

    - config (str): path of config file 
    - trigger (str): trigger class (MB or HM)
    - pdg_d (int): PDG code of the D meson
    - pdg_v0 (int): PDG code of the V0

    """

    with open(config, "r") as yml_cfg:  # pylint: disable=bad-option-value
        cfg = yaml.load(yml_cfg, yaml.FullLoader)

    # define titles and labels
    if pdg_d == 411:
        name_d = "Dplus"
        name_decay = "DplustoKpipi"
        label_mass_d = r"M(K$\pi\pi$) (GeV/$c^2$)"
        mass_d = Particle.from_pdgid(pdg_d).mass*1e-3
    elif pdg_d == 413:
        name_d = "Dstar"
        name_decay = "DstartoD0pi"
        label_mass_d = r"M(K$\pi\pi$)-M(K$\pi$) (GeV/$c^2$)"
        mass_d = Particle.from_pdgid(pdg_d).mass*1e-3-Particle.from_pdgid(421).mass*1e-3
    else:
        print(f"ERROR: D pdg code {pdg_d} not supported")

    if pdg_v0 == 310:
        name_v0 = "K0S"
        label_mass_v0 = r"M($\pi\pi$) (GeV/$c^2$)"
    elif pdg_v0 == 3122:
        name_v0 = "Lambda"
        label_mass_v0 = r"M(p$\pi$) (GeV/$c^2$)"
    else:
        print(f"ERROR: V0 pdg code {pdg_v0} not supported")
    mass_v0 = Particle.from_pdgid(pdg_v0).mass*1e-3

    pdg_reso = -1
    if pdg_d == 411:
        if pdg_v0 == 310:
            pdg_reso = 435
            label_mass_reso = r"M(D$^+$K$_\mathrm{S}^0$)$-$M(D$^+$) (GeV/$c^2$)"
        else:
            print(f"ERROR: combination of D and V0 pdg codes {pdg_d}-{pdg_v0} not supported")
    elif pdg_d == 413:
        if pdg_v0 == 310:
            pdg_reso = 10433
            label_mass_reso = r"M(D$^{*+}$K$_\mathrm{S}^0$)$-$M(D$^{*+}$) (GeV/$c^2$)"
        else:
            print(f"ERROR: combination of D and V0 pdg codes {pdg_d}-{pdg_v0} not supported")

    # get inputs
    df_all = uproot.concatenate(
        os.path.join(f"{cfg['input_dirs']}",
                     "AnalysisResults*.root:"
                     f"PWGHF_D2H_HFResoBuilder{name_decay}{name_d}K0S_{trigger}/fNtupleCharmReso"),
        library="pd"
    )

    # apply selections to D and V0
    pt_min_d = cfg['selections'][pdg_d][trigger]['pt_min']
    mass_min_d = cfg['selections'][pdg_d][trigger]['mass_min']
    mass_max_d = cfg['selections'][pdg_d][trigger]['mass_max']
    sel_string_d = f"pt_D > {pt_min_d} and {mass_min_d} < inv_mass_D < {mass_max_d} and"
    pt_mins_d = cfg["selections"][pdg_d][trigger]["BDT"]["pt_mins"]
    pt_maxs_d = cfg["selections"][pdg_d][trigger]["BDT"]["pt_maxs"]
    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins_d, pt_maxs_d)):
        bdt_bkg = cfg["selections"][pdg_d][trigger]["BDT"]["BDT_bkg"][ipt]
        bdt_prompt = cfg["selections"][pdg_d][trigger]["BDT"]["BDT_prompt"][ipt]
        sel_string_pt = f"({pt_min} < pt_D < {pt_max} and "
        sel_string_pt += f"outputscore_bkg_D < {bdt_bkg} and outputscore_prompt_D > {bdt_prompt})"
        if ipt == 0:
            sel_string_d += f" ({sel_string_pt}"
        elif ipt < len(pt_mins_d) - 1:
            sel_string_d += f" or {sel_string_pt}"
        else:
            sel_string_d += f" or {sel_string_pt})"

    pt_v0_min = cfg['selections'][pdg_v0][trigger]['pt_min']
    cosp_v0_min = cfg['selections'][pdg_v0][trigger]['cosp_min']
    declen_xy_v0_min = cfg['selections'][pdg_v0][trigger]['declen_xy_min']
    dca_dau_v0_min = cfg['selections'][pdg_v0][trigger]['dca_dau_min']
    sel_string_v0 = f"pt_v0 > {pt_v0_min} and cosp_v0 > {cosp_v0_min} and" \
        + f" declen_xy_v0 > {declen_xy_v0_min} and dca_dau_min_v0 > {dca_dau_v0_min}"

    sel_string = sel_string_d + " and " + sel_string_v0

    df_sel = df_all[df_all["id_v0"] == pdg_v0]
    df_sel.query(sel_string, inplace=True)

    # fit D and v0 mass (TODO: to be improved vs pT)
    data_hdl_v0 = DataHandler(df_sel, var_name="inv_mass_v0", limits=[mass_v0-0.03, mass_v0+0.03])
    if pdg_v0 == 310:
        fitter_v0 = F2MassFitter(data_hdl_v0, ["doublegaus"], ["expo"], "K0S_ptint")
        fitter_v0.set_signal_initpar(0, "sigma1", 0.0001)
    elif pdg_v0 == 3122:
        fitter_v0 = F2MassFitter(data_hdl_v0, ["gaussian"], ["chebpol1"], "Lambda_ptint")
        fitter_v0.set_signal_initpar(0, "sigma", 0.0001)

    fitter_v0.set_particle_mass(0, pdg_id=pdg_v0)
    fitter_v0.mass_zfit()
    fig_v0 = fitter_v0.plot_mass_fit(style="ATLAS", bins=200, axis_title=label_mass_v0)
    mean_v0, _ = fitter_v0.get_mass(0)
    frac_gaus_v0 = fitter_v0.get_signal_parameter(0, "frac1")[0]
    if pdg_v0 == 310:
        width_v0 = frac_gaus_v0*fitter_v0.get_signal_parameter(0, "sigma1")[0] + \
            (1-frac_gaus_v0)*fitter_v0.get_signal_parameter(0, "sigma2")[0]
    elif pdg_v0 == 3122:
        width_v0, _ = fitter_v0.get_signal_parameter(0, "sigma")

    fig_v0.savefig(os.path.join(cfg["output_dir"], f"mass_{name_v0}_{trigger}.pdf"))

    if pdg_d == 411:
        data_hdl_D = DataHandler(df_sel, var_name="inv_mass_D", limits=[1.75, 1.98])
        fitter_d = F2MassFitter(data_hdl_D, ["gaussian"], ["expo"], "Dplus_ptint")
        fitter_d.set_signal_initpar(0, "sigma", 0.015)
    elif pdg_d == 413:
        data_hdl_D = DataHandler(df_sel, var_name="inv_mass_D", limits=[Particle.from_pdgid(211).mass*1e-3, 0.16])
        fitter_d = F2MassFitter(data_hdl_D, ["gaussian"], ["expopow"], "Dstar_ptint")
        fitter_d.set_signal_initpar(0, "sigma", 0.0008)
    fitter_d.set_particle_mass(0, mass=mass_d)
    fitter_d.mass_zfit()
    fig_d = fitter_d.plot_mass_fit(style="ATLAS", bins=200, axis_title=label_mass_d)
    mean_d, _ = fitter_d.get_mass(0)
    width_d, _ = fitter_d.get_sigma(0)
    fig_d.savefig(os.path.join(cfg["output_dir"], f"mass_{name_d}_{trigger}.pdf"))

    # plot mass D vs. mass V0
    inv_mass_d = df_sel["inv_mass_D"].to_numpy()
    inv_mass_v0 = df_sel["inv_mass_v0"].to_numpy()
    h_massd_vs_massv0 = (
        Hist.new
        .Reg(200, min(inv_mass_d), max(inv_mass_d), name="x", label=label_mass_d)
        .Reg(200, min(inv_mass_v0), max(inv_mass_v0), name="y", label=label_mass_v0)
        .Double()
    )
    h_massd_vs_massv0.fill(x=inv_mass_d, y=inv_mass_v0)

    fig_masscorr = plt.figure(figsize=(10, 8))
    grid_masscorr = fig_masscorr.add_gridspec(
        2, 2, hspace=0, wspace=0, width_ratios=[4, 1], height_ratios=[1, 4]
    )
    ax_masscorr = {}
    ax_masscorr["main_ax"] = fig_masscorr.add_subplot(grid_masscorr[1, 0])
    ax_masscorr["top_ax"] = fig_masscorr.add_subplot(
        grid_masscorr[0, 0], sharex=ax_masscorr["main_ax"])
    ax_masscorr["top_ax"].yaxis.set_ticks_position('left')
    ax_masscorr["side_ax"] = fig_masscorr.add_subplot(
        grid_masscorr[1, 1], sharey=ax_masscorr["main_ax"])
    ax_masscorr["side_ax"].xaxis.set_ticks_position('bottom')
    h_massd_vs_massv0.plot2d_full(main_cmap="turbo",
                                  top_yerr=True, top_color='black', top_histtype='errorbar',
                                  side_yerr=True, side_color='black', side_histtype='errorbar',
                                  ax_dict=ax_masscorr)

    ax_masscorr["main_ax"].plot([min(inv_mass_d), max(inv_mass_d)],
                                [mean_v0-3*width_v0, mean_v0-3*width_v0], color="darkred", ls="--")
    ax_masscorr["main_ax"].plot([min(inv_mass_d), max(inv_mass_d)],
                                [mean_v0+3*width_v0, mean_v0+3*width_v0], color="darkred", ls="--")
    ax_masscorr["main_ax"].plot([mean_d-3*width_d, mean_d-3*width_d],
                                [min(inv_mass_v0), max(inv_mass_v0)], color="darkred", ls="--")
    ax_masscorr["main_ax"].plot([mean_d+3*width_d, mean_d+3*width_d],
                                [min(inv_mass_v0), max(inv_mass_v0)], color="darkred", ls="--")

    fig_masscorr.savefig(os.path.join(cfg["output_dir"], f"mass{name_d}_vs_mass{name_v0}_{trigger}.pdf"))

    # select signal region
    n_sigma_d = cfg["selections"][pdg_reso][trigger]["delta_mass_D"]
    n_sigma_v0 = cfg["selections"][pdg_reso][trigger]["delta_mass_D"]
    inv_mass_d_sel = f"{mean_d-n_sigma_d*width_d} < inv_mass_D < {mean_d+n_sigma_d*width_d}"
    inv_mass_v0_sel = f"{mean_v0-n_sigma_v0*width_v0} < inv_mass_v0 < {mean_v0+n_sigma_v0*width_v0}"
    mass_min_reso = max(Particle.from_pdgid(pdg_v0).mass*1e-3, cfg["selections"][pdg_reso][trigger]["mass_min"])
    mass_max_reso = cfg["selections"][pdg_reso][trigger]["mass_max"]
    inv_mass_reso_sel = f"{mass_min_reso} < delta_inv_mass_reso < {mass_max_reso}"
    df_sel.query(f"{inv_mass_d_sel} and {inv_mass_v0_sel} and {inv_mass_reso_sel}", inplace=True)
    df_sel.to_parquet(os.path.join(cfg["output_dir"],
                                   f"{name_d}_{name_v0}_{trigger}_centsel.parquet.gzip"), compression="gzip")

    # plot mass reso vs pT
    h_massreso_vs_pt = (
        Hist.new
        .Reg(62, 0., 31., name="x", label=r"$p_\mathrm{T}$ (GeV/$c$)")
        .Reg(100, mass_min_reso, mass_max_reso, name="y", label=label_mass_reso)
        .Double()
    )
    h_massreso_vs_pt.fill(x=df_sel["pt_reso"].to_numpy(), y=df_sel["delta_inv_mass_reso"].to_numpy())

    fig_massvspt = plt.figure(figsize=(10, 8))
    grid_massvspt = fig_massvspt.add_gridspec(
        2, 2, hspace=0, wspace=0, width_ratios=[4, 1], height_ratios=[1, 4]
    )
    ax_massvspt = {}
    ax_massvspt["main_ax"] = fig_massvspt.add_subplot(grid_massvspt[1, 0])
    ax_massvspt["top_ax"] = fig_massvspt.add_subplot(
        grid_massvspt[0, 0], sharex=ax_massvspt["main_ax"])
    ax_massvspt["top_ax"].yaxis.set_ticks_position('left')
    ax_massvspt["side_ax"] = fig_massvspt.add_subplot(
        grid_massvspt[1, 1], sharey=ax_massvspt["main_ax"])
    ax_massvspt["side_ax"].xaxis.set_ticks_position('bottom')
    h_massreso_vs_pt.plot2d_full(main_cmap="turbo",
                                 top_yerr=True, top_color='black', top_histtype='errorbar',
                                 side_yerr=True, side_color='black', side_histtype='errorbar',
                                 ax_dict=ax_massvspt)

    fig_massvspt.savefig(os.path.join(cfg["output_dir"], f"mass{name_d}{name_v0}_vs_pt_{trigger}.pdf"))

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config_proj.yml", help="input config file")
    parser.add_argument("--trigger", "-t", metavar="text",
                        default="HM", required=True, help="trigger class (options: [HM, MB]")
    parser.add_argument("--pdg_D", "-d", type=int, required=True, default=413,
                        help="pdg code of the D meson (options: [411, 413])")
    parser.add_argument("--pdg_V0", "-v0", type=int, required=True, default=310,
                        help="pdg code of the V0 (options: [310, 3122])")
    args = parser.parse_args()

    project(args.config, args.trigger, args.pdg_D, args.pdg_V0)
