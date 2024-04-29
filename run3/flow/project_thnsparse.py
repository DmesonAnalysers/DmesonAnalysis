import ROOT
import yaml
import argparse
import sys
from alive_progress import alive_bar
from flow_analysis_utils import get_vn_versus_mass, get_centrality_bins, compute_r2, get_invmass_vs_deltaphi

def check_anres(config, an_res_file, centrality, wagon_id, outputdir, suffix, vn_method):
    with open(config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)
    # read configuration
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

    # sanity check
    if len(pt_mins) != len(pt_maxs) or len(pt_mins) != len(inv_mass_bins):
        sys.exit('\033[91m FATAL: pt_mins, pt_maxs, and inv_mass_bins must have the same length\033[0m')
    for inv_mass_bin in inv_mass_bins:
        for ibin, (bin_low, bin_high) in enumerate(zip(inv_mass_bin[:-1], inv_mass_bin[1:])):
            if bin_low > bin_high:
                sys.exit(f'\033[91m FATAL: bin_low > bin_high for bin {ibin} in inv_mass_bins\033[0m')

    # BDT cuts
    axis_bdt_bkg = config['axes']['bdt_bkg']
    axis_bdt_sig = config['axes']['bdt_sig']
    apply_btd_cuts = config['apply_btd_cuts']
    if apply_btd_cuts:
        bkg_ml_cuts = config['bkg_ml_cuts']
        sig_ml_cuts = config['sig_ml_cuts']

    # load input file
    infile = ROOT.TFile(an_res_file, 'READ')
    if wagon_id != '':
        wagon_id = f'_id{wagon_id}'
    thnsparse = infile.Get(f'hf-task-flow-charm-hadrons{wagon_id}/hSparseFlowCharm')

    # output file
    outfile = ROOT.TFile(f'{outputdir}/proj{suffix}.root', 'RECREATE')
    hist_reso = ROOT.TH1F('hist_reso', 'hist_reso', len(cent_bins) - 1, cent_bins[0], cent_bins[-1])
    cent_min = cent_bins[0]
    cent_max = cent_bins[1]
    thnsparse_selcent = thnsparse.Clone(f'thnsparse_selcent{cent_min}_{cent_max}')
    thnsparse_selcent.GetAxis(axis_cent).SetRangeUser(cent_min, cent_max)
    reso = 1#compute_r2(infile, wagon_id, cent_min, cent_max, det_A, det_B, det_C, vn_method)
    hist_reso.SetBinContent(hist_reso.FindBin(cent_min), reso)
    vn_axis = axis_sp if vn_method == 'sp' else axis_deltaphi

    with alive_bar(len(pt_mins), title='Processing pt bins') as bar:
        # loop over pt bins
        for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
            inv_mass_bin = inv_mass_bins[ipt]
            outfile.mkdir(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
            thnsparse_selcent.GetAxis(axis_pt).SetRangeUser(pt_min, pt_max)
            # apply BDT cuts
            if apply_btd_cuts:
                print('Applying BDT cuts')
                thnsparse_selcent.GetAxis(axis_bdt_bkg).SetRangeUser(0, bkg_ml_cuts[0])
                thnsparse_selcent.GetAxis(axis_bdt_sig).SetRangeUser(sig_ml_cuts[0], 1)

            # create output histograms
            outfile.cd(f'cent_bins{cent_min}_{cent_max}/pt_bins{pt_min}_{pt_max}')
            # compute vn versus mass
            if vn_method == 'deltaphi':
                hist_mass_inplane, \
                hist_mass_outplane = get_invmass_vs_deltaphi(thnsparse_selcent,
                                                             axis_deltaphi,
                                                             axis_mass)
                hist_mass_inplane.SetName(f'hist_mass_inplane_cent{cent_min}_{cent_max}_pt{pt_min}_{pt_max}')
                hist_mass_outplane.SetName(f'hist_mass_outplane_cent{cent_min}_{cent_max}_pt{pt_min}_{pt_max}')
                hist_mass_inplane.SetDirectory(0)
                hist_mass_outplane.SetDirectory(0)
                hist_mass_inplane.Write()
                hist_mass_outplane.Write()
            else:
                hist_vn = get_vn_versus_mass(thnsparse_selcent, inv_mass_bin,
                                             axis_mass, vn_axis, False)
                hist_vn.SetName(f'hist_vn_{vn_method}_pt{pt_min}_{pt_max}')
                hist_vn.SetDirectory(0)
                if reso > 0:
                    hist_vn.Scale(1./reso)
                    hist_vn.Write()
                hist_mass = thnsparse_selcent.Projection(axis_mass)
                hist_mass.SetName(f'hist_mass_cent{cent_min}_{cent_max}_pt{pt_min}_{pt_max}')
                hist_mass.Write()
            bar()
            outfile.cd('..')

    # save output
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
    parser.add_argument("--wagon_id", "-w", metavar="text",
                        default="", help="wagon ID", required=False)
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--vn_method", "-vn", metavar="text",
                        default="sp", help="vn technique (sp, ep, deltaphi)")
    args = parser.parse_args()

    check_anres(
        config=args.config,
        an_res_file=args.an_res_file,
        centrality=args.centrality,
        wagon_id=args.wagon_id,
        outputdir=args.outputdir,
        suffix=args.suffix,
        vn_method=args.vn_method
    )