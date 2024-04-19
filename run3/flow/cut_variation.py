import ROOT
import yaml
import argparse
import numpy as np
import sys
from alive_progress import alive_bar
sys.path.append('..')
from flow_analysis_utils import get_vn_versus_mass, get_centrality_bins, compute_r2

def cut_var(config, an_res_file, centrality, outputdir, suffix, do_ep):
    with open(config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)
    pt_mins = config['ptmins']
    pt_maxs = config['ptmaxs']
    _, cent_bins = get_centrality_bins(centrality)
    det_A = config['detA']
    det_B = config['detB']
    det_C = config['detC']
    axis_cent = config['axes']['cent']
    axis_pt = config['axes']['pt']
    axis_mass = config['axes']['mass']
    axis_sp = config['axes']['sp']
    axis_deltaphi = config['axes']['deltaphi']
    inv_mass_bins = config['inv_mass_bins']
    axis_bdt_bkg = config['axes']['bdt_bkg']
    axis_bdt_sig = config['axes']['bdt_sig']
    bkg_cut_mins = config['cut_variation']['bdt_cut']['bkg']['min']
    bkg_cut_maxs = config['cut_variation']['bdt_cut']['bkg']['max']
    bkg_cut_steps = config['cut_variation']['bdt_cut']['bkg']['step']
    sig_cut_mins = config['cut_variation']['bdt_cut']['sig']['min']
    sig_cut_maxs = config['cut_variation']['bdt_cut']['sig']['max']
    sig_cut_steps = config['cut_variation']['bdt_cut']['sig']['step']

    infile = ROOT.TFile(an_res_file, 'READ')
    thnsparse = infile.Get('hf-task-flow-charm-hadrons_id10397/hSparseFlowCharm')

    outfile = ROOT.TFile(f'{outputdir}/cut_var{suffix}.root', 'RECREATE')
    hist_reso = ROOT.TH1F('hist_reso', 'hist_reso', len(cent_bins) - 1, cent_bins[0], cent_bins[-1])
    cent_min = cent_bins[0]
    cent_max = cent_bins[1]
    thnsparse_selcent = thnsparse.Clone(f'thnsparse_selcent{cent_min}_{cent_max}')
    thnsparse_selcent.GetAxis(axis_cent).SetRangeUser(cent_min, cent_max)
    reso = compute_r2(infile, cent_min, cent_max, det_A, det_B, det_C, do_ep)
    hist_reso.SetBinContent(hist_reso.FindBin(cent_min), reso)

    for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
        print(f'Processing pt bin {pt_min} - {pt_max}')
        inv_mass_bin = inv_mass_bins[ipt]
        outfile.mkdir(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
        thnsparse_selcent.GetAxis(axis_pt).SetRangeUser(pt_min, pt_max)
        bkg_cut_pt = [cut for cut in np.arange(bkg_cut_mins[ipt], bkg_cut_maxs[ipt], bkg_cut_steps[ipt])]
        sig_cut_pt = [cut for cut in np.arange(sig_cut_mins[ipt], sig_cut_maxs[ipt], sig_cut_steps[ipt])]
        outfile.cd(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
        with alive_bar(len(bkg_cut_pt) * len(sig_cut_pt), title='Processing BDT cuts') as bar:
            for _, (bkg_cut) in enumerate(bkg_cut_pt):
                for _, (sgn_cut) in enumerate(sig_cut_pt):
                    hist_out_label = f'cent{cent_min}_{cent_max}_pt{pt_min}_{pt_max}_bkg{bkg_cut}_sig{sgn_cut}'
                    thnsparse_selcent.GetAxis(axis_bdt_bkg).SetRangeUser(0, bkg_cut)
                    thnsparse_selcent.GetAxis(axis_bdt_sig).SetRangeUser(sgn_cut, 1.05)
                    hist_mass = thnsparse_selcent.Projection(axis_mass)
                    hist_mass.SetName(f'hist_mass_{hist_out_label}')
                    hist_mass.Write()
                    if do_ep:
                        hist_vn_ep = get_vn_versus_mass(thnsparse_selcent, inv_mass_bin,
                                                        axis_mass, axis_deltaphi, False)
                        hist_vn_ep.SetName(f'hist_vn_ep_{hist_out_label}')
                        hist_vn_ep.SetDirectory(0)
                        if reso > 0:
                            hist_vn_ep.Scale(1./reso)
                        hist_vn_ep.Write()
                    else:
                        hist_vn_sp = get_vn_versus_mass(thnsparse_selcent, inv_mass_bin, axis_mass, axis_sp)
                        hist_vn_sp.SetDirectory(0)
                        hist_vn_sp.SetName(f'hist_vn_sp_{hist_out_label}')
                        if reso > 0:
                            hist_vn_sp.Scale(1./reso)
                        hist_vn_sp.Write()
                bar()
        outfile.cd('..')
    outfile.cd()
    hist_reso.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("an_res_file", metavar="text",
                        default="an_res.root", help="input ROOT file with anres")
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k3050", help="centrality class")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--doEP",  action="store_true", default=False,
                        help="do EP resolution")
    args = parser.parse_args()

    cut_var(
        config=args.config,
        an_res_file=args.an_res_file,
        centrality=args.centrality,
        outputdir=args.outputdir,
        suffix=args.suffix,
        do_ep=args.doEP
    )

