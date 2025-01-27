'''
This sricpt is used to pre-process a/multi large AnRes.root for the BDT training:
    - split the input by pT
    - obtain the sigma from prompt enhance sample
python3 pre_process.py config_pre.yml AnRes_1.root AnRes_2.root --pre --sigma  
'''
import os
import sys
import yaml
import numpy as np
import array
import ROOT
from ROOT import TFile
import argparse
import itertools
import concurrent.futures
sys.path.append("/Users/mcosti/Analysis/DmesonAnalysis/run3/flow/BDT/")
sys.path.append("/Users/mcosti/Analysis/DmesonAnalysis/run3/flow/")
from flow_analysis_utils import get_centrality_bins
from sparse_dicts import get_sparses

def cook_thnsparse(thnsparse_list, ptmins, ptmaxs, axestokeep):
    '''
    Split the input THnSparse by pT bins and project onto axes to keep.

    Input:
        - thnsparse_list (list): List of input THnSparse objects.
        - ptmins (list): List of minimum pT values for each bin.
        - ptmaxs (list): List of maximum pT values for each bin.
        - axestokeep (list): List of axes to keep in the projection.

    Returns:
        - dict: Dictionary of projected THnSparse objects for each pT bin.
    '''
    sparses = {}
    for iThn, (sparse_key, sparse) in enumerate(thnsparse_list.items()):
        for iPt in range(0, len(ptmins)):
            binMin = sparse.GetAxis(1).FindBin(ptmins[iPt]*1.00001)
            binMax = sparse.GetAxis(1).FindBin(ptmaxs[iPt]*0.99999)
            sparse.GetAxis(1).SetRange(binMin, binMax)
            thn_proj = sparse.Projection(len(axestokeep), array.array('i', axestokeep), 'O')
            
            if iThn == 0:
                sparses[iPt] = thn_proj
            else:
                sparses[iPt].Add(thn_proj)
    return sparses

# def rebin_sparse(sparse, ipt, axvars, config):
    # if 'Mass' in axvars:
    #     bin_diffs = np.round(np.diff(config['inv_mass_bins'][ipt]), 2)
    #     print(f"bin_diffs: {bin_diffs}")
    #     min_diff = min(bin_diffs)
    #     print(f"min_diff: {min_diff}")
    #     divided_diffs = bin_diffs / min_diff
    #     print(f"divided_diffs: {divided_diffs}")
    #     are_all_integers = np.all(divided_diffs == np.floor(divided_diffs))
    #     print(f"are_all_integers: {are_all_integers}")
    #     rebin_factor = min_diff / sparse.GetAxis(0).GetBinWidth(1)
    #     iteration = 0
    #     while rebin_factor != np.floor(rebin_factor) or iteration<=5:
    #         iteration = iteration+1
    #         rebin_factor = min_diff / (2*iteration*sparse.GetAxis(0).GetBinWidth(1))
    #         print(f"sparse.GetAxis(0).GetBinWidth(1): {sparse.GetAxis(0).GetBinWidth(1)}")
    #         print(f"rebin_factor: {rebin_factor}")
            # floored_rebin = np.floor(rebin_factor)
            # print(f"floored_rebin: {floored_rebin}")
    # if 'score_FD' in axvars:
    #     if config['minimisation']['correlated']:
    #         bin_width = config['cut_variation']['corr_bdt_cut']['sig']['step'][ipt]
    #         rebin_fd = bin_width / (sparse.GetAxis(2).GetBinWidth(1))
    #         if rebin_fd == np.floor(rebin_fd):
    #             rebinned_sparse = sparse.Rebin(array.array('i', [4,1,int(rebin_fd)]))
    #         print(f"sparse.GetNbins(): {sparse.GetNbins()}")
    #         print(f"rebinned_sparse.GetNbins(): {rebinned_sparse.GetNbins()}")
    # return rebinned_sparse

def pre_process(config, ptmins, ptmaxs, centmin, centmax, axestokeep, outputDir):
    
    # Load the ThnSparse
    thnsparse_list, _, _, sparse_axes = get_sparses(config, True, False)

    out_file = TFile('Projections.root', 'recreate')
    for isparse, (key, sparse) in enumerate(thnsparse_list.items()):
        if 'Flow' in key:
            out_file.mkdir(f'Flow_{isparse}')
            out_file.cd(f'Flow_{isparse}')
            for idim in range(sparse.GetNdimensions()):
                histo = sparse.Projection(idim)
                histo.SetName(sparse.GetAxis(idim).GetName())
                histo.SetTitle(sparse.GetAxis(idim).GetTitle())
                histo.Write()
        
    def process_pt_bin(iPt, ptmin, ptmax, centmin, centmax, bkg_max_cut, thnsparse_list, axestokeep, outputDir):
        
        # add possibility to apply cuts for different variables
        for iThn, (sparse_key, sparse) in enumerate(thnsparse_list.items()):
            print(sparse)
            sparse.GetAxis(sparse_axes['Flow']['Pt']).SetRangeUser(ptmin, ptmax)
            sparse.GetAxis(sparse_axes['Flow']['cent']).SetRangeUser(centmin, centmax)
            sparse.GetAxis(sparse_axes['Flow']['score_bkg']).SetRangeUser(0, bkg_max_cut)
            
            print(axestokeep)
            thn_proj = sparse.Projection(len(axestokeep), array.array('i', [sparse_axes['Flow'][axtokeep] for axtokeep in axestokeep]), 'O')
            thn_proj.SetName(sparse.GetName())
            
            if iThn == 0:
                processed_sparse = thn_proj.Clone()
            else:
                processed_sparse.Add(thn_proj)
        
        if config.get('RebinSparse'):
            rebin_factors = [config['RebinSparse'][axtokeep] for axtokeep in axestokeep]
            sparse = sparse.Rebin(array.array('i', rebin_factors))
            
        os.makedirs(f'{outputDir}/pre/AnRes', exist_ok=True)
        outFile = ROOT.TFile(f'{outputDir}/pre/AnRes/AnalysisResults_pt_{int(ptmin*10)}_{int(ptmax*10)}.root', 'RECREATE')
        outFile.mkdir('hf-task-flow-charm-hadrons')
        outFile.cd('hf-task-flow-charm-hadrons')
        processed_sparse.Write()
        outFile.Close()
        
        out_file.mkdir(f'Flow_{iPt}_ipt_{ptmin}_{ptmax}')
        out_file.cd(f'Flow_{iPt}_ipt_{ptmin}_{ptmax}')
        for idim in range(processed_sparse.GetNdimensions()):
            histo = processed_sparse.Projection(idim)
            histo.SetName(processed_sparse.GetAxis(idim).GetName())
            histo.SetTitle(processed_sparse.GetAxis(idim).GetTitle())
            histo.Write()
        
        del processed_sparse
        
        print(f'Finished processing pT bin {ptmin} - {ptmax}')

    bkg_maxs = config['cut_variation']['corr_bdt_cut']['bkg_max']
    # Loop over each pt bin in parallel
    max_workers = 6
    with concurrent.futures.ThreadPoolExecutor(max_workers) as executor:
        tasks = [executor.submit(process_pt_bin, iPt, ptmin, ptmax, centmin, centmax, bkg_maxs[iPt], thnsparse_list, axestokeep, outputDir) for iPt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs))]
        for task in concurrent.futures.as_completed(tasks):
            task.result()
        
def get_sigma(preFiles, config_pre, centrality, resolution, outputDir, skip_projection=False):
    with open(config_pre, 'r') as cfgPre:
        config = yaml.safe_load(cfgPre)

    apply_btd_cuts = config['apply_btd_cuts']    
    
    if not apply_btd_cuts:
        print("Warning: 'apply_btd_cuts' is not enabled in the config_pre.")
    
    if not skip_projection:
        skip_proj = ''
    else:
        skip_proj = '--skip_projection'
    
    os.makedirs(f'{outputDir}/pre/sigma', exist_ok=True)

    str_preFiles = ' '.join(preFiles)
    command = (f"python3 ../flow/run_full_flow_analysis.py {config_pre} {str_preFiles} \
                -c {centrality} -o {outputDir}/pre/sigma \
                -s sigma -v sp --r {resolution} \
                --skip_efficiency --skip_resolution \
                {skip_proj}")
    os.system(f'{command}')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument('config_pre', metavar='text', 
                        default='config_pre.yml', help='configuration file')
    parser.add_argument('an_res_file', metavar='text', 
                        nargs='+', help='input ROOT files with anres')
    parser.add_argument('--pre', action='store_true', help='pre-process the AnRes.root')
    parser.add_argument('--sigma', action='store_true', help='get the sigma')
    parser.add_argument('--skip_projection', '-sp', action='store_true', help='skip the projection')
    args = parser.parse_args()

    if not args.pre and not args.sigma:
        print('Please specify the action to perform.')
        sys.exit(1)

    with open(args.config_pre, 'r') as cfgPre:
        config = yaml.safe_load(cfgPre)
  
    # Load the configuration
    ptmins = config['ptmins']
    ptmaxs = config['ptmaxs']
    axestokeep = config['axestokeep']
    outputDir = config['outputDir']
    
    centMin, centMax = get_centrality_bins(config['centrality'])[1]
    
    if args.pre:
        pre_process(config, ptmins, ptmaxs, centMin, centMax, axestokeep, outputDir)
    
    if args.sigma:
        
        centrality = config['centrality']
        resolution = config['resolution']
        
        if os.path.exists(f'{outputDir}/pre'):
            preFiles = [f'{outputDir}/pre/AnRes/AnalysisResults_pt{iFile}.root' for iFile in range(len(ptmins))]
        else:
            raise ValueError(f'No eff folder found in {outputDir}')
        preFiles.sort()
        
        # you have to know the sigma from the differet prompt enhance samples is stable first
        get_sigma(preFiles, args.config_pre, centrality, resolution, outputDir, skip_projection=args.skip_projection)
        
        
        
        
        
        
        
        
        
def rebin_sparses(config, ptmins, ptmaxs, varstokeep, outputDir):
    
    # Load the ThnSparse
    thnsparse_list, _, _, axes_dict = get_sparses(config, True, False)

    def process_pt_bin(iPt, ptmin, ptmax, centmax, centmin, bkg_max_cut, thnsparse_list, axes, axestokeep, axesbins, outputDir):
        pre_thnsparse = None
        # add posibility to apply cuts for different variables
        print(f"CIAOOO: {thnsparse_list}")
        print(f"CIAOOO: {axes}")
        bin_edges_arrays = [axbin for axbin in axesbins]
        print(f"bin_edges_arrays: {bin_edges_arrays}")
        
        
        list_ranges = [list(range(1,len(bin))) for bin in axesbins]
        result = list(itertools.product(*list_ranges))
        n_dims = len(bin_edges_arrays)  # Number of dimensions
        n_bins = [len(edges) - 1 for edges in bin_edges_arrays]  # Number of bins per dimension
        xmin = [min(edges) for edges in bin_edges_arrays]
        xmax = [max(edges) for edges in bin_edges_arrays]
        print(f"list_ranges: {list_ranges}")
        # print(result)
        print(f"n_bins: {n_bins}")
        print(f"n_dims: {n_dims}")
        print(f"array.array('d', xmin): {array.array('d', xmin)}")
        print(f"array.array('d', xmax): {array.array('d', xmax)}")
        sparse_refilled = ROOT.THnSparseD(
            "sparse_refilled", "Sparse Refilled",
            n_dims,
            array.array('i', n_bins),  # Number of bins per dimension
            array.array('d', xmin),   # Minimum values
            array.array('d', xmax)    # Maximum values
        )

        # Step 4: Set custom bin edges for each dimension
        for dim in range(n_dims):
            axis = sparse_refilled.GetAxis(dim)
            axis.Set(len(bin_edges_arrays[dim]) - 1, array.array('d', bin_edges_arrays[dim]))




        for iThn, (sparse_key, sparse) in enumerate(thnsparse_list.items()):
            print(f"thnsparse: {sparse}")


            sparse.GetAxis(axes['Flow']['Pt']).SetRangeUser(ptmin, ptmax)
            sparse.GetAxis(axes['Flow']['cent']).SetRangeUser(centmin, centmax)
            sparse.GetAxis(axes['Flow']['score_bkg']).SetRangeUser(0, bkg_max_cut)

            for ibin in result:
                print("\n")
                print(f"ibin: {ibin}")
                for iax, ax in enumerate(axestokeep):
                    # print(f"bin_edges_arrays[idim]: {bin_edges_arrays[idim]}")
                    # print(f"bin_edges_arrays[idim][ibin[idim]-1]: {bin_edges_arrays[idim][ibin[idim]-1]}")
                    # print(f"bin_edges_arrays[idim][ibin[idim]]: {bin_edges_arrays[idim][ibin[idim]]}")
                    sparse.GetAxis(axes_dict['Flow'][ax]).SetRangeUser(bin_edges_arrays[iax][ibin[iax]-1], bin_edges_arrays[iax][ibin[iax]])
                
                # print(f"sparse.Projection(0).Integral(): {sparse.Projection(0).Integral()}")
                sparse_refilled.SetBinContent(np.asarray(ibin, 'i'), sparse.Projection(0).Integral())


            thn_proj = sparse.Projection(len(axestokeep), array.array('i', axes['Flow'][varstokeep]), 'O')
            thn_proj.SetName(sparse.GetName())
            
            if iThn == 0:
                pre_sparse = thn_proj.Clone()
            else:
                pre_sparse.Add(thn_proj)
        
        os.makedirs(f'{outputDir}/pre/AnRes', exist_ok=True)
        outFile = ROOT.TFile(f'{outputDir}/pre/AnRes/AnalysisResults_pt{iPt}.root', 'RECREATE')
        outFile.mkdir('hf-task-flow-charm-hadrons')
        outFile.cd('hf-task-flow-charm-hadrons')
        pre_sparse.Write()
        outFile.Close()
        
        del pre_sparse
        
        print(f'Finished processing pT bin {ptmin} - {ptmax}')
   
    mass_pt_bins = config['inv_mass_bins']
    bkg_maxs = config['cut_variation']['corr_bdt_cut']['bkg_max']
    sig_cuts = config['cut_variation']['corr_bdt_cut']['sig']
    fd_bins = [np.arange(sig_cuts['min'][iPt], sig_cuts['max'][iPt], sig_cuts['step'][iPt]).tolist() for iPt in range(len(ptmins))] + [1]
    sp_bins = np.arange(config['sp_bins']['min'], config['sp_bins']['max'], config['sp_bins']['step']).tolist()
    new_bins = [[mass_bin, fd_bin, sp_bins] for mass_bin, fd_bin in zip(mass_pt_bins, fd_bins)]
    
    max_workers = 6
    # Loop over each pt bin in parallel
    with concurrent.futures.ThreadPoolExecutor(max_workers) as executor:
        tasks = [executor.submit(process_pt_bin, iPt, ptmin, ptmax, 30, 40, bkg_maxs[iPt], thnsparse_list, axes_dict, varstokeep, new_bins[iPt], outputDir) for iPt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs))]
        for task in concurrent.futures.as_completed(tasks):
            task.result()