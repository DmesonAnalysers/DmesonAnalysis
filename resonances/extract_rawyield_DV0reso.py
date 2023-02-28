"""
Script for raw yield extraction
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
from hist import Hist
import uproot
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from particle import Particle

def create_hist(pt_lims, contents, errors, label_pt=r"$p_\mathrm{T}~(\mathrm{GeV}/c)$"):
    """
    Helper method to create histogram

    Parameters
    ----------

    - pt_lims (list): list of pt limits
    - contents (list): histogram contents
    - errors (list): histogram errors
    - label_pt (str): label for x axis

    Returns
    ----------
    - histogram (hist.Hist)

    """

    pt_cent = [(pt_min+pt_max)/2 for pt_min, pt_max in zip(pt_lims[:-1], pt_lims[1:])]
    histo = Hist.new.Var(pt_lims, name="x", label=label_pt).Weight()
    histo.fill(pt_cent, weight=contents)
    histo.view(flow=False).variance = np.array(errors)**2

    return histo


def extract_rawyield(input_file, config, output_dir, suffix):
    """
    function for raw yield computation

    Parameters
    ----------

    - input_file (str): path of input file
    - config (str): path of config file
    - output_dir (str): output directory
    - suffix (str): suffix to append to output files
    """

    with open(config, "r") as yml_cfg:  # pylint: disable=bad-option-value
        cfg = yaml.load(yml_cfg, yaml.FullLoader)

    pdg_reso = cfg["pdg_reso"]

    # define titles and labels
    extra_info_loc = None
    if pdg_reso == 10433:
        name_reso = "Ds1plus"
        label_mass = r"M($\mathrm{D}^{*+}\mathrm{K_S^0}$)$-$M($\mathrm{D}^{*+}$)"
        pdg_v0 = 310
        pdg_d = 413
        extra_info_loc = ["lower right", "upper center"]
    elif pdg_reso == 435:
        name_reso = "Ds2starplus"
        label_mass = r"M($\mathrm{D}^+\mathrm{K_S^0}$)$-$M($\mathrm{D}^+$)"
        pdg_v0 = 310
        pdg_d = 411
        extra_info_loc = ["upper left", "lower right"]
    else:
        print(f"ERROR: pdg code {pdg_reso} not supported")
        sys.exit()

    trigger = ""
    if "MB" in input_file:
        trigger = "MB"
    elif "HM" in input_file:
        trigger = "HM"

    df_reso = pd.read_parquet(input_file)

    pt_mins_reso = cfg["fit_config"]["pt_mins"]
    pt_maxs_reso = cfg["fit_config"]["pt_maxs"]
    mass_mins_reso = cfg["fit_config"]["mass_mins"]
    mass_maxs_reso = cfg["fit_config"]["mass_maxs"]
    signal_pdf_reso = cfg["fit_config"]["signal_pdf"]
    bkg_pdf_reso = cfg["fit_config"]["bkg_pdf"]
    pt_lims_reso = pt_mins_reso.copy()
    pt_lims_reso.append(pt_maxs_reso[-1])

    # define output file
    outfile_name = os.path.join(
        output_dir, f"mass_{name_reso}_pt{pt_mins_reso[0]}-{pt_maxs_reso[-1]}_{trigger}{suffix}.root")
    file_root = uproot.recreate(outfile_name)
    file_root.close()

    # fit resonance
    signif, signif_unc = [], []
    s_over_b, s_over_b_unc = [], []
    raw_yields, raw_yields_unc = [], []
    means, means_unc = [], []
    sigmas, sigmas_unc = [], []
    for pt_min, pt_max, mass_min, mass_max, spdf, bpdf in zip(
        pt_mins_reso, pt_maxs_reso, mass_mins_reso, mass_maxs_reso, signal_pdf_reso, bkg_pdf_reso):
        
        mass_min = max([mass_min, Particle.from_pdgid(pdg_v0).mass*1e-3])
        df_sel_pt = df_reso.query(f"{pt_min} < pt_reso < {pt_max}")
        data_hdl = DataHandler(df_sel_pt, var_name="delta_inv_mass_reso", limits=[mass_min, mass_max], nbins=100)

        fitter_reso = F2MassFitter(data_hdl, [spdf], [bpdf], name=f"{name_reso}_pt{pt_min}_{pt_max}")
        delta_mass = Particle.from_pdgid(pdg_reso).mass*1e-3 - Particle.from_pdgid(pdg_d).mass*1e-3
        width = Particle.from_pdgid(pdg_reso).width*1.e-3
        fitter_reso.set_signal_initpar(0, "gamma", width, fix=True)
        fitter_reso.set_signal_initpar(0, "sigma", 0.0008, limits=[0.0006, 0.0012])
        fitter_reso.set_particle_mass(0, mass=delta_mass, limits=[delta_mass-0.1, delta_mass+0.1])
        fitter_reso.set_signal_initpar(0, "frac", 0.01, limits=[0., 1.])
        fitter_reso.set_background_initpar(0, "mass", Particle.from_pdgid(pdg_v0).mass*1e-3)
        fitter_reso.set_background_initpar(0, "c0", 100)
        fitter_reso.set_background_initpar(0, "c1", 1)
        fitter_reso.set_background_initpar(0, "c2", 1)
        fitter_reso.mass_zfit()
        fig_reso = fitter_reso.plot_mass_fit(style="ATLAS",
                                             figsize=(8, 8),
                                             axis_title=rf"{label_mass} (GeV/$c^2$)",
                                             show_extra_info=True,
                                             extra_info_massnhwhm=3.,
                                             extra_info_loc=extra_info_loc)
        fig_reso_res = fitter_reso.plot_raw_residuals(style="ATLAS", figsize=(8, 8), axis_title=rf"{label_mass} (GeV/$c^2$)")
        fig_reso.savefig(os.path.join(output_dir, f"mass_{name_reso}_pt{pt_min}-{pt_max}_{trigger}{suffix}.pdf"))
        fig_reso_res.savefig(os.path.join(output_dir,
                                          f"mass_residuals_{name_reso}_pt{pt_min}-{pt_max}_{trigger}{suffix}.pdf"))
        fitter_reso.dump_to_root(outfile_name, option="update")

        rawy, rawy_unc = fitter_reso.get_raw_yield(0)

        sign, sign_unc = fitter_reso.get_significance(0, nhwhm=3.)
        soverb, soverb_unc = fitter_reso.get_signal_over_background(0, nhwhm=3.)
        mean, mean_unc = fitter_reso.get_signal_parameter(0, "mu")
        sigma, sigma_unc = fitter_reso.get_signal_parameter(0, "sigma")

        raw_yields.append(rawy)
        raw_yields_unc.append(rawy_unc)
        signif.append(sign)
        signif_unc.append(sign_unc)
        s_over_b.append(soverb)
        s_over_b_unc.append(soverb_unc)
        means.append(mean)
        means_unc.append(mean_unc)
        sigmas.append(sigma)
        sigmas_unc.append(sigma_unc)

    file_root = uproot.update(outfile_name)
    file_root["h_rawyields"] = create_hist(pt_lims_reso, raw_yields, raw_yields_unc)
    file_root["h_significance"] = create_hist(pt_lims_reso, signif, signif_unc)
    file_root["h_soverb"] = create_hist(pt_lims_reso, s_over_b, s_over_b_unc)
    file_root["h_means"] = create_hist(pt_lims_reso, means, means_unc)
    file_root["h_sigmas"] = create_hist(pt_lims_reso, sigmas, sigmas_unc)
    file_root.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("input_file", metavar="text",
                        default="projected_data.parquet.gzip",
                        help="input parquet file")
    parser.add_argument("config", metavar="text",
                        default="config_fit.yml", help="input config file")
    parser.add_argument("--output_dir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    extract_rawyield(args.input_file, args.config, args.output_dir, args.suffix)
