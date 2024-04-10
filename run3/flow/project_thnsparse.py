import ROOT
import yaml
import argparse
import numpy as np
import sys
from alive_progress import alive_bar
from flow_analysis_utils import get_vn_versus_mass

def compute_r2(reso_file, cent_min, cent_max, detA, detB, detC, do_ep):

    if do_ep:
        hist_name = 'hf-task-flow-charm-hadrons/epReso/hEpReso'
    else:
        hist_name = 'hf-task-flow-charm-hadrons/spReso/hSpReso'
    
    print(f'Computing resolution for {detA}, {detB}, {detC}')
    print(f'Centrality bin: {cent_min} - {cent_max}')
    print(f'Hist name: {hist_name}{detA}{detB}')
    print(f'Hist name: {hist_name}{detA}{detC}')
    print(f'Hist name: {hist_name}{detB}{detC}')
    detA_detB = reso_file.Get(f'{hist_name}{detA}{detB}')
    detA_detC = reso_file.Get(f'{hist_name}{detA}{detC}')
    detB_detC = reso_file.Get(f'{hist_name}{detB}{detC}')

    cent_bin_min = detA_detB.GetXaxis().FindBin(cent_min)
    cent_bin_max = detA_detB.GetXaxis().FindBin(cent_max)

    proj_detA_detB = detA_detB.ProjectionY(f'{hist_name}{detA}{detB}_proj{cent_min}_{cent_max}', cent_bin_min, cent_bin_max)
    proj_detA_detC = detA_detC.ProjectionY(f'{hist_name}{detA}{detC}_proj{cent_min}_{cent_max}', cent_bin_min, cent_bin_max)
    proj_detB_detC = detB_detC.ProjectionY(f'{hist_name}{detB}{detC}_proj{cent_min}_{cent_max}', cent_bin_min, cent_bin_max)

    average_detA_detB = proj_detA_detB.GetMean()
    average_detA_detC = proj_detA_detC.GetMean()
    average_detB_detC = proj_detB_detC.GetMean()

    reso = np.sqrt((average_detA_detB * average_detA_detC) / average_detB_detC) if average_detB_detC > 0 else 1.e+100000
    return reso

def check_anres(config, an_res_file, outputdir, suffix, do_ep):
    with open(config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)

    pt_mins = config['ptmins']
    pt_maxs = config['ptmaxs']
    cent_bins = config['centrality_bins']
    det_A = config['detA']
    det_B = config['detB']
    det_C = config['detC']
    axis_cent = config['axes']['cent']
    axis_pt = config['axes']['pt']
    axis_mass = config['axes']['mass']
    axis_sp = config['axes']['sp']
    axis_deltaphi = config['axes']['deltaphi']
    meson = config['Dmeson']
    #axis_bdt_bkg = config['axes']['bdt_bkg']
    #axis_bdt_sig = config['axes']['bdt_sig']

    # TODO: move this to the config file
    if meson == 'Dplus':
        inv_mass_bins = [1.74, 1.78, 1.82, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.90, 1.92, 1.96, 2.00]
        #[1.74, 1.75, 1.76, 1.77, 1.78, 1.79, 1.80, 1.81, 1.82, 1.83, 1.835, 1.84, 1.845, 1.85, 1.855, 1.86, 1.865, 1.87, 1.875, 1.88, 1.885, 1.89, 1.895, 1.90, 1.905, 1.91, 1.915, 1.92, 1.93, 1.94, 1.95, 1.96, 1.97, 1.98, 1.99, 2.00]
    elif meson == 'Ds':
        inv_mass_bins = [1.8, 1.82, 1.84, 1.86, 1.88, 1.90, 1.92, 1.94, 1.96, 1.98, 2.00, 2.02, 2.04, 2.06, 2.08, 2.10]
        
    infile = ROOT.TFile(an_res_file, 'READ')
    thnsparse = infile.Get('hf-task-flow-charm-hadrons/hSparseFlowCharm')

    outfile = ROOT.TFile(f'{outputdir}/proj{suffix}.root', 'RECREATE')
    hist_reso = ROOT.TH1F('hist_reso', 'hist_reso', len(cent_bins) - 1, cent_bins[0], cent_bins[-1])
    for _, (cent_min, cent_max) in enumerate(zip(cent_bins[:-1], cent_bins[1:])):
        thnsparse_selcent = thnsparse.Clone(f'thnsparse_selcent{cent_min}_{cent_max}')
        thnsparse_selcent.GetAxis(axis_cent).SetRangeUser(cent_min, cent_max)
        reso = compute_r2(infile, cent_min, cent_max, det_A, det_B, det_C, do_ep)
        hist_reso.SetBinContent(hist_reso.FindBin(cent_min), reso)

        with alive_bar(len(pt_mins)) as bar:
            for _, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
                print(f'Processing pt bin {pt_min} - {pt_max}')
                #if meson == 'Dplus':
                    #if pt_min >= 4:
                inv_mass_bins = [1.74, 1.78, 1.82, 1.84, 1.85, 1.86, 1.87, 1.88, 1.89, 1.90, 1.92, 1.96, 2.00]
                    #else:
                    #inv_mass_bins = [1.74, 1.77, 1.80, 1.81, 1.82, 1.83, 1.835, 1.84, 1.845, 1.85, 1.855, 1.86, 1.865, 1.87, 1.875, 1.88, 1.885, 1.89, 1.895, 1.90, 1.905, 1.91, 1.915, 1.92, 1.93, 1.96, 1.99, 2.00]
                outfile.mkdir(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
                thnsparse_selcent.GetAxis(axis_pt).SetRangeUser(pt_min, pt_max)
                
                if do_ep:
                    hist_vn_ep = get_vn_versus_mass(thnsparse_selcent, inv_mass_bins,
                                                    axis_mass, axis_deltaphi, False)
                    hist_vn_ep.SetName(f'hist_vn_ep_pt{pt_min}_{pt_max}')
                    hist_vn_ep.SetDirectory(0)
                    hist_vn_ep.Scale(1./reso)
                    outfile.cd(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
                    hist_vn_ep.Write()
                else:
                    hist_vn_sp = get_vn_versus_mass(thnsparse_selcent, inv_mass_bins, axis_mass, axis_sp)
                    hist_vn_sp.SetDirectory(0)
                    hist_vn_sp.SetName(f'hist_vn_sp_pt{pt_min}_{pt_max}')
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
        outputdir=args.outputdir,
        suffix=args.suffix,
        do_ep=args.doEP
    )
