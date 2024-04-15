import ROOT
import yaml
import argparse
import numpy as np
from alive_progress import alive_bar
from flow_analysis_utils import get_vn_versus_mass, get_centrality_bins, compute_r2

def check_anres(config, an_res_file, centrality, outputdir, suffix, do_ep):
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
    apply_btd_cuts = config['apply_btd_cuts']
    if apply_btd_cuts:
        bkg_ml_cuts = config['bkg_ml_cuts']
        sig_ml_cuts = config['sig_ml_cuts']

    infile = ROOT.TFile(an_res_file, 'READ')
    thnsparse = infile.Get('hf-task-flow-charm-hadrons/hSparseFlowCharm')

    outfile = ROOT.TFile(f'{outputdir}/proj{suffix}.root', 'RECREATE')
    hist_reso = ROOT.TH1F('hist_reso', 'hist_reso', len(cent_bins) - 1, cent_bins[0], cent_bins[-1])
    cent_min = cent_bins[0]
    cent_max = cent_bins[1]
    thnsparse_selcent = thnsparse.Clone(f'thnsparse_selcent{cent_min}_{cent_max}')
    thnsparse_selcent.GetAxis(axis_cent).SetRangeUser(cent_min, cent_max)
    reso = 1#compute_r2(infile, cent_min, cent_max, det_A, det_B, det_C, do_ep)
    hist_reso.SetBinContent(hist_reso.FindBin(cent_min), reso)

    with alive_bar(len(pt_mins)) as bar:
        for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
            print(f'Processing pt bin {pt_min} - {pt_max}')
            inv_mass_bin = inv_mass_bins[ipt]
            outfile.mkdir(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
            thnsparse_selcent.GetAxis(axis_pt).SetRangeUser(pt_min, pt_max)
            if apply_btd_cuts:
                print('Applying BDT cuts')
                thnsparse_selcent.GetAxis(axis_bdt_bkg).SetRangeUser(0, bkg_ml_cuts[0])
                thnsparse_selcent.GetAxis(axis_bdt_sig).SetRangeUser(sig_ml_cuts[0], 1)
            
            if do_ep:
                hist_vn_ep = get_vn_versus_mass(thnsparse_selcent, inv_mass_bin,
                                                axis_mass, axis_deltaphi, False)
                hist_vn_ep.SetName(f'hist_vn_ep_pt{pt_min}_{pt_max}')
                hist_vn_ep.SetDirectory(0)
                if reso > 0:
                    hist_vn_ep.Scale(1./reso)
                outfile.cd(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
                hist_vn_ep.Write()
            else:
                hist_vn_sp = get_vn_versus_mass(thnsparse_selcent, inv_mass_bin, axis_mass, axis_sp)
                hist_vn_sp.SetDirectory(0)
                hist_vn_sp.SetName(f'hist_vn_sp_pt{pt_min}_{pt_max}')
                if reso > 0:
                    hist_vn_sp.Scale(1./reso)
                outfile.cd(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
                hist_vn_sp.Write()
            hist_mass = thnsparse_selcent.Projection(axis_mass)
            hist_mass.SetName(f'hist_mass_cent{cent_min}_{cent_max}_pt{pt_min}_{pt_max}')
            hist_mass.Write()
            bar()
            outfile.cd('..')
    outfile.cd()
    hist_reso.Write()
    outfile.Close()


if __name__ == "__main__":
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

    check_anres(
        config=args.config,
        an_res_file=args.an_res_file,
        centrality=args.centrality,
        outputdir=args.outputdir,
        suffix=args.suffix,
        do_ep=args.doEP
    )
