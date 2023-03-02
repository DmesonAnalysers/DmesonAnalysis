"""
Script for raw yield systematics
run: python rawyield_syst_Dsreso.py
"""

import os
import sys
import argparse
import yaml
import pandas as pd
import numpy as np
import uproot
import hist
from hist import Hist
from alive_progress import alive_bar
import itertools
from flarefly.data_handler import DataHandler
from flarefly.fitter import F2MassFitter
from particle import Particle
from extract_rawyield_DV0reso import create_hist

def blockPrint():
    sys.stdout = open(os.devnull, 'w')

def enablePrint():
    sys.stdout = sys.__stdout__

def main():
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config_fit.yml", help="input config file")
    parser.add_argument("inputfile", metavar="text",
                        default="input.parquet.gzip", help="input file")
    parser.add_argument('--histoname', metavar=('text'), nargs=1, required=False,
                        help='Histogram name')
    args = parser.parse_args()

    if args.inputfile.endswith("root") and not args.histoname:
        print("ERROR: Please provide a histogram name")
        sys.exit()
    isparquet = True if args.inputfile.endswith("gzip") else False
    config = args.config
    with open(config, "r") as yml_cfg:  # pylint: disable=bad-option-value
        cfg = yaml.load(yml_cfg, yaml.FullLoader)

    #______________________________________________________________________________
    # Collect multi-trial parameters
    if "Dstar" in args.inputfile or "Ds1" in args.inputfile:
        name_reso = "Ds1plus"
        pdg_reso = 10433
        pdg_d = 413
        pdg_v0 = 310
    elif "Dplus" in args.inputfile or "Ds2" in args.inputfile:
        name_reso = "Ds2starplus"
        pdg_reso = 435
        pdg_d = 411
        pdg_v0 = 310
    else:
        print(f"ERROR: Resonance {cfg['input']} not supported")
        sys.exit()

    trigger = ""
    if "MB" in args.inputfile:
        trigger = "MB"
    elif "HM" in args.inputfile:
        trigger = "HM"

    print(f'\033[1m\033[92m Loading input file {args.inputfile} \033[0m')
    if isparquet:
        df_reso = pd.read_parquet(args.inputfile)
    pt_mins_reso = cfg["fit_config"][trigger]["pt_mins"]
    pt_maxs_reso = cfg["fit_config"][trigger]["pt_maxs"]
    mass_mins_reso = cfg["fit_config"][trigger]["mass_mins"]
    mass_maxs_reso = cfg["fit_config"][trigger]["mass_maxs"]
    signal_pdf_reso = cfg["fit_config"][trigger]["signal_pdf"]
    bkg_pdf_reso = cfg["fit_config"][trigger]["bkg_pdf"]
    gamma_reso = cfg["fit_config"][trigger]["gamma"]
    quality_criteria = cfg["quality_criteria"]
    tot_combinations = list(itertools.product(pt_mins_reso, mass_mins_reso, mass_maxs_reso, signal_pdf_reso, bkg_pdf_reso, gamma_reso))
    print(f"Total number of combinations: {len(tot_combinations)}")

    #______________________________________________________________________________
    # Define output file
    outlabel = cfg["outlabel"]
    if args.histoname:
        outlabel += f"_{args.histoname[0]}"
    file_root = uproot.recreate(os.path.join(cfg["output_dir"], f"Multitrial_{name_reso}_pt{pt_mins_reso[0]}-{pt_maxs_reso[-1]}_{trigger}{outlabel}.root"))


    #______________________________________________________________________________
    # Multitrial loop
    signif, signif_unc = [], []
    raw_yields, raw_yields_unc = [], []
    means, means_unc = [], []
    sigmas, sigmas_unc = [], []
    chi2s, chi2s_unc = [], []
    gamma_fits, gamma_fits_unc = [], []
    delta_mass = Particle.from_pdgid(pdg_reso).mass*1e-3 - Particle.from_pdgid(pdg_d).mass*1e-3

    with alive_bar(len(tot_combinations), bar='smooth',
                   title='Collecting multitrials') as bar:
        for pt_min, pt_max in zip(pt_mins_reso, pt_maxs_reso): # loop over pt bins
            if isparquet:
                df_sel_pt = df_reso.query(f"{pt_min} < pt_reso < {pt_max}")
            else:
                print("WARNING: root input not queried yet")
            for spdf in signal_pdf_reso: # loop over signal pdfs (tipically only one)
                for gamma in gamma_reso: # loop over gamma
                    if gamma == "width":
                        width = Particle.from_pdgid(pdg_reso).width*1.e-3
                    elif gamma == "width_upper":
                        width = (Particle.from_pdgid(pdg_reso).width_upper + Particle.from_pdgid(pdg_reso).width)*1.e-3
                    elif gamma == "width_lower":
                        width = (Particle.from_pdgid(pdg_reso).width - Particle.from_pdgid(pdg_reso).width_lower)*1.e-3
                    else:
                        print(f"ERROR: gamma {gamma} not supported")
                        sys.exit()
                    for bpdf in bkg_pdf_reso: # loop over bkg pdfs
                        for mass_max in mass_maxs_reso: # loop over mass
                            for mass_min in mass_mins_reso:
                                if isparquet:
                                    data_hdl = DataHandler(df_sel_pt, var_name="delta_inv_mass_reso",
                                                           limits=[mass_min, mass_max],
                                                           nbins=100)
                                else:
                                    data_hdl = DataHandler(args.inputfile,
                                                           var_name="delta_inv_mass_reso",
                                                            limits=[mass_min, mass_max],
                                                           histoname=args.histoname[0])
                                #blockPrint() # TODO: add quiet option to F2MassFitter
                                fitter_reso = F2MassFitter(data_hdl,
                                                           name_signal_pdf=[spdf],
                                                           name_background_pdf=[bpdf],
                                                           name=f"{name_reso}_pt{pt_min}_{pt_max}_mass{mass_min}_{mass_max}_{spdf}_{bpdf}_{width}")
                                fitter_reso.set_signal_initpar(0, "gamma", width, fix=True)
                                fitter_reso.set_signal_initpar(0, "sigma", 0.0008, limits=[0.0004, 0.0020])
                                fitter_reso.set_particle_mass(0, mass=delta_mass, limits=[delta_mass-0.1, delta_mass+0.1])
                                fitter_reso.set_signal_initpar(0, "frac", 0.01, limits=[0., 1.])
                                fitter_reso.set_background_initpar(0, "mass", Particle.from_pdgid(pdg_v0).mass*1e-3)
                                fitter_reso.set_background_initpar(0, "c0", 100)
                                fitter_reso.set_background_initpar(0, "c1", 1)
                                fitter_reso.set_background_initpar(0, "c2", 1)
                                fitter_reso.mass_zfit()

                                rawy, rawy_unc = fitter_reso.get_raw_yield(0)
                                sign, sign_unc = fitter_reso.get_significance(0, nhwhm=3.)
                                mean, mean_unc = fitter_reso.get_signal_parameter(0, "mu")
                                sigma, sigma_unc = fitter_reso.get_signal_parameter(0, "sigma")
                                if fitter_reso.get_chi2_ndf() > 1.e+50: # bad fits may have chi2 = inf
                                    chi2 = 1.e+50
                                else:
                                    chi2 = fitter_reso.get_chi2_ndf()
                                if np.isnan(sign): # bad fits may have significance = nan
                                    sign = 0.
                                raw_yields.append(rawy)
                                raw_yields_unc.append(rawy_unc)
                                signif.append(sign)
                                signif_unc.append(sign_unc)
                                means.append(mean)
                                means_unc.append(mean_unc)
                                sigmas.append(sigma)
                                sigmas_unc.append(sigma_unc)
                                chi2s.append(chi2)
                                gamma_fits.append(width)

                                #enablePrint() # TODO: add quiet option to F2MassFitter
                                bar()
    chi2s_unc = [0. for i in range(len(chi2s))]
    gamma_fits_unc = [0. for i in range(len(gamma_fits))]

    #__________________________________________________________________________
    # Save results
    lims = [i for i in range(len(signif)+1)]
    max_ry = int(max(raw_yields)*1.2)
    file_root["h_rawyields"] = create_hist(lims, raw_yields, raw_yields_unc, label_pt="trial #")
    file_root["h_gamma"] = create_hist(lims, gamma_fits, gamma_fits_unc, label_pt="trial #")
    file_root["h_signif"] = create_hist(lims, signif, signif_unc, label_pt="trial #")
    file_root["h_means"] = create_hist(lims, means, means_unc, label_pt="trial #")
    file_root["h_sigmas"] = create_hist(lims, sigmas, sigmas_unc, label_pt="trial #")
    file_root["h_chi2"] = create_hist(lims, chi2s, chi2s_unc, label_pt="trial #")
    hsyst = Hist(hist.axis.Regular(max_ry, 0, max_ry, name="raw yield", label="raw yield"))
    for i, (rway, signif, chi2) in enumerate(zip(raw_yields, signif, chi2s)):
        if signif >= quality_criteria['significance'] and chi2 <= quality_criteria['chi2']:
            hsyst.fill(rway)
            print(f"trial {i} accepted")
            print(f"signif = {signif}, chi2 = {chi2}")
        else:
            print(f"WARNING: trial {i} discarded")
            print(f"signif = {signif}, chi2 = {chi2}")
    file_root["h_syst"] = hsyst
    file_root.close()
    input("Press enter to exit.")


main()