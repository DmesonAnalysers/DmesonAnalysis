'''
Script for the variation with the BDT cuts.
python cut_variation.py config.yaml an_res.root -c k3050 -r resolution.root -o path/to/output -s text
'''
import ROOT
import yaml
import argparse
import numpy as np
import os
import sys
from alive_progress import alive_bar
sys.path.append('..')
from flow_analysis_utils import get_vn_versus_mass, get_centrality_bins, compute_r2

def cut_var(config, an_res_file, centrality, resolution, outputdir, suffix):
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
    inv_mass_bins = config['inv_mass_bins']
    axis_bdt_bkg = config['axes']['bdt_bkg']
    axis_bdt_sig = config['axes']['bdt_sig']
    bkg_cut_mins = config['cut_variation']['bdt_cut']['bkg']['min']
    bkg_cut_maxs = config['cut_variation']['bdt_cut']['bkg']['max']
    bkg_cut_steps = config['cut_variation']['bdt_cut']['bkg']['step']
    sig_cut_mins = config['cut_variation']['bdt_cut']['sig']['min']
    sig_cut_maxs = config['cut_variation']['bdt_cut']['sig']['max']
    sig_cut_steps = config['cut_variation']['bdt_cut']['sig']['step']

    # get resolution
    resoFile = ROOT.TFile(resolution, 'READ')
    histo_reso = resoFile.Get(f'{det_A}_{det_B}_{det_C}/histo_reso')
    histo_reso.SetName('hist_reso')
    histo_reso.SetDirectory(0)
    reso = histo_reso.GetBinContent(1)

    infile = ROOT.TFile(an_res_file, 'READ')
    thnsparse = infile.Get('hf-task-flow-charm-hadrons/hSparseFlowCharm')
    print(infile.GetName())

    cent_min = cent_bins[0]
    cent_max = cent_bins[1]
    thnsparse_selcent = thnsparse.Clone(f'thnsparse_selcent{cent_min}_{cent_max}')
    thnsparse_selcent.GetAxis(axis_cent).SetRangeUser(cent_min, cent_max)

    os.makedirs(f'{outputdir}/proj', exist_ok=True)

    for ipt in range(len(pt_mins)):
        if bkg_cut_steps[0] != 0:
            bkg_cut_pt = [cut for cut in np.arange(bkg_cut_mins[ipt], bkg_cut_maxs[ipt], bkg_cut_steps[ipt])]
        else:
            bkg_cut_pt = [bkg_cut_maxs[ipt]]
        sig_cut_pt = [cut for cut in np.arange(sig_cut_mins[ipt], sig_cut_maxs[ipt], sig_cut_steps[ipt])]

    with alive_bar(len(bkg_cut_pt) * len(sig_cut_pt), title='Processing BDT cuts') as bar:
        for _, (bkg_cut) in enumerate(bkg_cut_pt):
            for iCut, (sgn_cut) in enumerate(sig_cut_pt):

                outfile = ROOT.TFile(f'{outputdir}/proj/proj_{suffix}_{iCut:02d}.root', 'RECREATE')
                outfile.mkdir(f'cent_bins{cent_min}_{cent_max}')
                outfile.cd(f'cent_bins{cent_min}_{cent_max}')

                thnsparse_selcent.GetAxis(axis_bdt_bkg).SetRangeUser(0, bkg_cut)
                thnsparse_selcent.GetAxis(axis_bdt_sig).SetRangeUser(sgn_cut, 1)

                for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):

                    outfile.mkdir(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
                    outfile.cd(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')

                    inv_mass_bin = inv_mass_bins[ipt]
                    thnsparse_selcent.GetAxis(axis_pt).SetRangeUser(pt_min, pt_max)

                    hist_mass = thnsparse_selcent.Projection(axis_mass)
                    hist_mass.SetName(f'hist_mass_cent{cent_min}_{cent_max}_pt{pt_min}_{pt_max}')
                    hist_mass.Write()

                    hist_vn_sp = get_vn_versus_mass(thnsparse_selcent, inv_mass_bin, axis_mass, axis_sp)
                    hist_vn_sp.SetDirectory(0)
                    hist_vn_sp.SetName(f'hist_vn_sp_pt{pt_min}_{pt_max}')
                    if reso > 0:
                        hist_vn_sp.Scale(1./reso)
                    hist_vn_sp.Write()

                outfile.cd()
                histo_reso.Write()
                outfile.Close()
                bar()
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("an_res_file", metavar="text",
                        default="an_res.root", help="input ROOT file with anres")
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k3050", help="centrality class")
    parser.add_argument("--resolution", "-r", metavar="text",
                        default="reso.root", help="resolution file")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    cut_var(
        config=args.config,
        an_res_file=args.an_res_file,
        centrality=args.centrality,
        resolution=args.resolution,
        outputdir=args.outputdir,
        suffix=args.suffix
    )

