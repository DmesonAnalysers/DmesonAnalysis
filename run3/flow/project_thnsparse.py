import ROOT
from ROOT import TFile
import yaml
import argparse
import sys
import numpy as np
from alive_progress import alive_bar
from flow_analysis_utils import get_vn_versus_mass, get_centrality_bins, compute_r2, get_invmass_vs_deltaphi, get_occupancy, get_evselbits

def check_anres(config, an_res_file, centrality, resolution, wagon_id, 
                outputdir, suffix, vn_method, use_preprocessed):
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
    use_inv_mass_bins = config['use_inv_mass_bins']
    if use_inv_mass_bins:
        print('\033[93m WARNING: Using inv_mass_bins from '+ an_res_file +'\033[0m')
        inv_mass_bins = []
    # sanity check
    else:
        if len(pt_mins) != len(pt_maxs) or len(pt_mins) != len(inv_mass_bins):
            sys.exit('\033[91m FATAL: pt_mins, pt_maxs, and inv_mass_bins must have the same length\033[0m')
    for inv_mass_bin in inv_mass_bins:
        for ibin, (bin_low, bin_high) in enumerate(zip(inv_mass_bin[:-1], inv_mass_bin[1:])):
            if bin_low > bin_high:
                sys.exit(f'\033[91m FATAL: bin_low > bin_high for bin {ibin} in inv_mass_bins\033[0m')

    # resolution (TODO: add the possibility to use a file with resolutions for deltacent bins)
    if '.root' in resolution:
        reso_file = ROOT.TFile(resolution, 'READ')
        hist_reso = reso_file.Get(f'{det_A}_{det_B}_{det_C}/histo_reso_delta_cent')
        if not hist_reso:
            sys.exit(f'\033[91m FATAL: Could not find hist_reso in {resolution}\033[0m')
        if ((vn_method in ['sp', 'ep']) and (vn_method not in resolution)) or (vn_method == 'deltaphi' and 'ep' not in resolution):
            sys.exit(f'\033[91m FATAL: mimsatch between vn_method and resolution file\033[0m')
        reso = 0.
        for ibin in range(1, hist_reso.GetNbinsX() + 1):
            reso += hist_reso.GetBinContent(ibin)
        reso /= hist_reso.GetNbinsX()
        reso_file.Close()
    else:
        reso = float(resolution)

    # BDT cuts
    axis_bdt_bkg = config['axes']['bdt_bkg']
    axis_bdt_sig = config['axes']['bdt_sig']
    apply_btd_cuts = config['apply_btd_cuts']
    if apply_btd_cuts:
        bkg_ml_cuts = config['bkg_ml_cuts']
        sig_ml_cuts = config['sig_ml_cuts']

    # output file
    outfile = ROOT.TFile(f'{outputdir}/proj{suffix}.root', 'RECREATE')
    hist_reso = ROOT.TH1F('hist_reso', 'hist_reso', len(cent_bins) - 1, cent_bins[0], cent_bins[-1])
    cent_min = cent_bins[0]
    cent_max = cent_bins[1]
    hist_reso.SetBinContent(hist_reso.FindBin(cent_min), reso)
    vn_axis = axis_sp if vn_method == 'sp' else axis_deltaphi
    if use_inv_mass_bins:
        inv_mass_bins = [[] for _ in range(len(pt_mins))]

    # load input file
    thnsparse_list, thnsparse_selcent_list, thnsparse_selcents = [], [], []
    if use_preprocessed:
        for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
            infile = TFile.Open(f"{config['skimDir']}/AnalysisResults_pt_{int(pt_min*10)}_{int(pt_max*10)}.root", 'r')
            thnsparse_selcent_list.append(infile.Get(f'hf-task-flow-charm-hadrons/hSparseFlowCharm'))
            infile.Close()
        axis_bdt_sig = config['axestokeep'].index('score_FD')
        axis_sp = config['axestokeep'].index('sp')
        axis_mass = config['axestokeep'].index('Mass')
        vn_axis = axis_sp if vn_method == 'sp' else vn_axis
    else:
        for file in an_res_file:
            infile = ROOT.TFile(file, 'READ')
            if wagon_id != '':
                wagon_id = f'_id{wagon_id}'
            thnsparse_list.append(infile.Get(f'hf-task-flow-charm-hadrons{wagon_id}/hSparseFlowCharm'))
            if not thnsparse_list[-1]:
                sys.exit(f'\033[91m FATAL: Could not find thnsparse in {an_res_file}\033[0m')

        # sanity check of the pt bins
        bin_edges = [thnsparse_list[-1].GetAxis(axis_pt).GetBinLowEdge(bin) for bin in range(1, thnsparse_list[-1].GetAxis(axis_pt).GetNbins() + 1)]
        bin_edges.append(thnsparse_list[-1].GetAxis(axis_pt).GetBinUpEdge(thnsparse_list[-1].GetAxis(axis_pt).GetNbins()))
        if any(pt_min not in bin_edges or pt_max not in bin_edges for pt_min, pt_max in zip(pt_mins, pt_maxs)):
            sys.exit('\033[91m FATAL: Too granular pt bins, pt_min or pt_max not in the bin edges\033[0m')

        for thnsparse in thnsparse_list:
            thnsparse_selcent_list.append(thnsparse.Clone(f'thnsparse_selcent{cent_min}_{cent_max}'))
            thnsparse_selcent_list[-1].GetAxis(axis_cent).SetRangeUser(cent_min, cent_max)

    with alive_bar(len(pt_mins), title='Processing pt bins') as bar:
        # loop over pt bins
        for ipt, (pt_min, pt_max) in enumerate(zip(pt_mins, pt_maxs)):
            outfile.mkdir(f'cent_bins{cent_min}_{cent_max}/pt_bins{int(pt_min*10)}_{int(pt_max*10)}')
            if use_preprocessed:
                if config['apply_btd_cuts']:
                    thnsparse_selcent_list[ipt].GetAxis(axis_bdt_sig).SetRangeUser(sig_ml_cuts[ipt], 1)
                thnsparse_selcents = [thnsparse_selcent_list[ipt]]
                if use_inv_mass_bins:
                    inv_mass_bins_ipt = []
                    rebin = config['Rebin'][ipt]
                    hmass_dummy = thnsparse_selcent.Projection(axis_mass)
                    hmass_dummy = hmass_dummy.Rebin(rebin, f'hmass_dummy_{ipt}')
                    for ibin in range(1, hmass_dummy.GetNbinsX() + 1):
                        inv_mass_bins_ipt.append(hmass_dummy.GetBinLowEdge(ibin))
                        if ibin == hmass_dummy.GetNbinsX():
                            inv_mass_bins_ipt.append(hmass_dummy.GetBinLowEdge(ibin) + hmass_dummy.GetBinWidth(ibin))
                    del hmass_dummy
                    inv_mass_bins[ipt] = inv_mass_bins_ipt
                inv_mass_bin = inv_mass_bins[ipt]
            else:
                for iThn, thnsparse_selcent in enumerate(thnsparse_selcent_list):
                    thnsparse_selcent.GetAxis(axis_pt).SetRangeUser(pt_min, pt_max)
                    if use_inv_mass_bins:
                        inv_mass_bins_ipt = []
                        rebin = config['Rebin'][ipt]
                        hmass_dummy = thnsparse_selcent.Projection(axis_mass)
                        hmass_dummy = hmass_dummy.Rebin(rebin, f'hmass_dummy_{ipt}')
                        for ibin in range(1, hmass_dummy.GetNbinsX() + 1):
                            inv_mass_bins_ipt.append(hmass_dummy.GetBinLowEdge(ibin))
                            if ibin == hmass_dummy.GetNbinsX():
                                inv_mass_bins_ipt.append(hmass_dummy.GetBinLowEdge(ibin) + hmass_dummy.GetBinWidth(ibin))
                        del hmass_dummy
                        inv_mass_bins[ipt] = inv_mass_bins_ipt
                    inv_mass_bin = inv_mass_bins[ipt]
                    # apply BDT cuts (TODO: save the bdt score distribution after the cuts in the output file)
                    if apply_btd_cuts:
                        print('\033[93m WARNING: Applying BDT cuts\033[0m')
                        thnsparse_selcent.GetAxis(axis_bdt_bkg).SetRangeUser(0, bkg_ml_cuts[ipt])
                        thnsparse_selcent.GetAxis(axis_bdt_sig).SetRangeUser(sig_ml_cuts[ipt], 1)
                    thnsparse_selcents.append(thnsparse_selcent)
                    
            # create output histograms
            outfile.cd(f'cent_bins{cent_min}_{cent_max}/pt_bins{int(pt_min*10)}_{int(pt_max*10)}')
            # compute vn versus mass
            if vn_method == 'deltaphi':
                hist_mass_inplane, \
                hist_mass_outplane = get_invmass_vs_deltaphi(thnsparse_selcents,
                                                             axis_deltaphi,
                                                             axis_mass)
                hist_mass_inplane.SetName(f'hist_mass_inplane_proj_pt{int(pt_min*10)}_{int(pt_max*10)}')
                hist_mass_outplane.SetName(f'hist_mass_outplane_proj_pt{int(pt_min*10)}_{int(pt_max*10)}')
                hist_mass_inplane.SetDirectory(0)
                hist_mass_outplane.SetDirectory(0)
                hist_mass_inplane.Write()
                hist_mass_outplane.Write()
            else:
                hist_vn = get_vn_versus_mass(thnsparse_selcents, inv_mass_bin,
                                             axis_mass, vn_axis, False)
                hist_vn.SetName(f'hist_vn_{vn_method}_proj_pt{int(pt_min*10)}_{int(pt_max*10)}')
                hist_vn.SetDirectory(0)
                if reso > 0:
                    hist_vn.Scale(1./reso)
                    hist_vn.Write()
                for iThn, thnsparse_selcent in enumerate(thnsparse_selcents):
                    hist_mass_temp = thnsparse_selcent.Projection(axis_mass)
                    hist_mass_temp.SetName(f'hist_mass_proj_pt{int(pt_min*10)}_{int(pt_max*10)}_{iThn}')
                    if iThn == 0:
                        hist_mass = hist_mass_temp.Clone(f'hist_mass_proj_pt{int(pt_min*10)}_{int(pt_max*10)}')
                        hist_mass.SetDirectory(0)
                        hist_mass.Reset()
                    hist_mass.Add(hist_mass_temp)
                hist_mass.Write()
            
            if config['axes'].get('occupancy'):
                hist_occ = get_occupancy(thnsparse_selcents, config['axes']['occupancy'], False)
                hist_occ.SetName(f'hist_occ_pt{int(pt_min*10)}_{int(pt_max*10)}')
                hist_occ.SetDirectory(0)
                hist_occ.Write()

            if config['axes'].get('evselbits'):
                hist_evselbits = get_evselbits(thnsparse_selcents, config['axes']['evselbits'], False)
                hist_evselbits.SetName(f'hist_evselbits_pt{int(pt_min*10)}_{int(pt_max*10)}')
                hist_evselbits.SetDirectory(0)
                hist_evselbits.Write()
            
            bar()
            outfile.cd('..')

    # save output
    outfile.cd()
    hist_reso.Write('histo_reso_delta_cent')
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument('an_res_file', metavar='file.root', 
                        nargs='+', help='input ROOT files with anres')
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k3050", help="centrality class")
    parser.add_argument("--resolution", "-r",  default=1.,
                        help="resolution file/value", required=False)
    parser.add_argument("--wagon_id", "-w", metavar="text",
                        default="", help="wagon ID", required=False)
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--vn_method", "-vn", metavar="text",
                        default="sp", help="vn technique (sp, ep, deltaphi)")
    parser.add_argument('--preprocessed', action='store_true', 
                        help='Determines whether the sparses are pre-processed')
    args = parser.parse_args()

    check_anres(
        config=args.config,
        an_res_file=args.an_res_file,
        centrality=args.centrality,
        resolution=args.resolution,
        wagon_id=args.wagon_id,
        outputdir=args.outputdir,
        suffix=args.suffix,
        vn_method=args.vn_method,
        use_preprocessed=args.preprocessed
    )