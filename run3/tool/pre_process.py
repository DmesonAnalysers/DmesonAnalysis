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
script_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.append(f"{script_dir}/../flow/")
sys.path.append(f"{script_dir}/../flow/BDT")
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
    thnsparses = {}
    for iThn, thnsparse in enumerate(thnsparse_list):
        #TODO: add possibility to apply cuts for different variables
        for iPt in range(0, len(ptmins)):
            binMin = thnsparse.GetAxis(1).FindBin(ptmins[iPt]*1.00001)
            binMax = thnsparse.GetAxis(1).FindBin(ptmaxs[iPt]*0.99999)
            thnsparse.GetAxis(1).SetRange(binMin, binMax)
            thn_proj = thnsparse.Projection(len(axestokeep), array.array('i', axestokeep), 'O')
            
            if iThn == 0:
                thnsparses[iPt] = thn_proj
            else:
                thnsparses[iPt].Add(thn_proj)
    return thnsparses

def pre_process(config, ptmins, ptmaxs, centmin, centmax, axestokeep, outputDir):
    
    # Load the ThnSparse
    thnsparse_list, _, _, sparse_axes = get_sparses(config, True, False, False, config['flow_files'])

    os.makedirs(f'{outputDir}/pre/AnRes', exist_ok=True)
    out_file = TFile(f'{outputDir}/pre/AnRes/Projections_{centmin}_{centmax}_{ptmins}_{ptmaxs}.root', 'recreate')
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
        print(f'Processing pT bin {ptmin} - {ptmax}, cent {centmin}-{centmax}')
        # add possibility to apply cuts for different variables
        for iThn, (sparse_key, sparse) in enumerate(thnsparse_list.items()):
            cloned_sparse = sparse.Clone()
            cloned_sparse.GetAxis(sparse_axes['Flow']['Pt']).SetRangeUser(ptmin, ptmax)
            cloned_sparse.GetAxis(sparse_axes['Flow']['cent']).SetRangeUser(centmin, centmax)
            cloned_sparse.GetAxis(sparse_axes['Flow']['score_bkg']).SetRangeUser(0, bkg_max_cut)
            thn_proj = cloned_sparse.Projection(len(axestokeep), array.array('i', [sparse_axes['Flow'][axtokeep] for axtokeep in axestokeep]), 'O')
            thn_proj.SetName(cloned_sparse.GetName())
            
            if iThn == 0:
                processed_sparse = thn_proj.Clone()
            else:
                processed_sparse.Add(thn_proj)
        
        if config.get('RebinSparse'):
            rebin_factors = [config['RebinSparse'][axtokeep] for axtokeep in axestokeep]
            processed_sparse = processed_sparse.Rebin(array.array('i', rebin_factors))
        
        outFile = ROOT.TFile(f'{outputDir}/pre/AnRes/AnalysisResults_pt_{int(ptmin*10)}_{int(ptmax*10)}.root', 'recreate')
        outFile.mkdir('hf-task-flow-charm-hadrons')
        outFile.cd('hf-task-flow-charm-hadrons')
        processed_sparse.Write('hSparseFlowCharm')
        outFile.Close()
        
        out_file.mkdir(f'Flow_pt_{ptmin}_{ptmax}')
        out_file.cd(f'Flow_pt_{ptmin}_{ptmax}')
        for idim in range(processed_sparse.GetNdimensions()):
            histo = processed_sparse.Projection(idim)
            histo.SetName(processed_sparse.GetAxis(idim).GetName())
            histo.SetTitle(processed_sparse.GetAxis(idim).GetTitle())
            histo.Write()
        
        del processed_sparse
        
        print(f'Finished processing pT bin {ptmin} - {ptmax}')

    bkg_maxs = config['bkg_cuts']
    # Loop over each pt bin in parallel
    max_workers = 12 # hyperparameter
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

def process_pt_bin_Singlecut(iPt, ptmin, ptmax, centMin, centMax, config, bkg_max_cut, sig_mins, sig_maxs, thnsparse_list, sparse_axes, axestokeep, outputDir):

    print(f'Processing pT bin {ptmin} - {ptmax}, cent {centMin}-{centMax}')

    # add possibility to apply cuts for different variables
    processed_sparses = []
    for iThn, (sparse_key, sparse) in enumerate(thnsparse_list.items()):
        cloned_sparse = sparse.Clone()
        cloned_sparse.GetAxis(sparse_axes['Flow']['Pt']).SetRangeUser(ptmin, ptmax)
        cloned_sparse.GetAxis(sparse_axes['Flow']['cent']).SetRangeUser(centMin, centMax)
        cloned_sparse.GetAxis(sparse_axes['Flow']['score_bkg']).SetRangeUser(0, bkg_max_cut)
        
        temp_thn_projs = []
        for iSig, (sig_min, sig_max) in enumerate(zip(sig_mins, sig_maxs)):
            temp_cloned_sparse = cloned_sparse.Clone()
            temp_cloned_sparse.GetAxis(sparse_axes['Flow']['score_FD']).SetRangeUser(sig_min, sig_max)
            temp_thn_projs.append(temp_cloned_sparse.Projection(len(axestokeep), array.array('i', [sparse_axes['Flow'][axtokeep] for axtokeep in axestokeep]), 'O'))
            temp_thn_projs[-1].SetName(cloned_sparse.GetName() + f'_sig_{iSig}')
            temp_cloned_sparse.Delete()
            del temp_cloned_sparse
            
        # delete the cloned sparse
        cloned_sparse.Delete()
        del cloned_sparse
        
        if iThn == 0:
            for iSig, thn_proj in enumerate(temp_thn_projs):
                processed_sparse = thn_proj.Clone()
                processed_sparses.append(processed_sparse)
                temp_thn_projs[iSig].Delete()
        else:
            for iSig, thn_proj in enumerate(temp_thn_projs):
                processed_sparses[iSig].Add(thn_proj)
                temp_thn_projs[iSig].Delete()

        del temp_thn_projs
    
        if config.get('RebinSparse'):
            rebin_factors = array.array('i', [config['RebinSparse'][axtokeep] for axtokeep in axestokeep])
            if -1 not in rebin_factors:
                processed_sparse = processed_sparse.Rebin(len(rebin_factors), rebin_factors)

    for iSig, processed_sparse in enumerate(processed_sparses):
        outFile = ROOT.TFile(f'{outputDir}/pre_sys/AnRes/{iSig:02d}/AnalysisResults_pt_{int(ptmin*10)}_{int(ptmax*10)}.root', 'recreate')
        outFile.mkdir('hf-task-flow-charm-hadrons')
        outFile.cd('hf-task-flow-charm-hadrons')
        processed_sparse.Write('hSparseFlowCharm')
        outFile.Close()
        processed_sparses[iSig].Delete()
    del processed_sparse

def pre_sys_process(config, ptmins, ptmaxs, centmin, centmax, axestokeep, outputDir):
    
    os.makedirs(f'{outputDir}/pre_sys/AnRes', exist_ok=True)
    
    # Load the ThnSparse
    thnsparse_list, _, _, sparse_axes = get_sparses(config, True, False, False, config['flow_files'])

    bkg_cuts = config['bdt_cut']['bkg_cuts']
    sig_mins = config['bdt_cut']['sig_mins']
    sig_maxs = config['bdt_cut']['sig_maxs']

    mCutset = max(len(sig_min) for sig_min in sig_mins)
    for iCut in range(mCutset):
        os.makedirs(f'{outputDir}/pre_sys/AnRes/{iCut:02d}', exist_ok=True)

    # Loop over each pt bin in parallel
    max_workers = 12 # hyperparameter
    args = [(iPt, ptmin, ptmax, centmin, centmax, config, bkg_cuts[iPt], sig_mins[iPt], sig_maxs[iPt], thnsparse_list, sparse_axes, axestokeep, outputDir) for iPt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs))]
    with concurrent.futures.ProcessPoolExecutor(max_workers=4) as executor:
        tasks = executor.map(process_pt_bin_Singlecut, *zip(*args))
        for result in tasks:
            result

    for iPt in range(len(ptmins)):
        if len(sig_mins[iPt]) < mCutset:
            available_file_index = len(sig_mins[iPt])
            for iCut in range(available_file_index, mCutset):
                os.system(f'cp -r {outputDir}/pre_sys/AnRes/{(available_file_index-1):02d}/AnalysisResults_pt_{int(ptmins[iPt]*10)}_{int(ptmaxs[iPt]*10)}.root {outputDir}/pre_sys/AnRes/{iCut:02d}/')
    # with concurrent.futures.ThreadPoolExecutor(max_workers) as executor:
    #     tasks = [executor.submit(process_pt_bin_Singlecut, iPt, ptmin, ptmax, centmin, centmax, bkg_maxs[iPt], sig_mins[iPt], sig_maxs[iPt], thnsparse_list, sparse_axes, axestokeep, outputDir) for iPt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs))]
    #     for task in concurrent.futures.as_completed(tasks):
    #         task.result()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument('config_pre', metavar='text', 
                        default='config_pre.yml', help='configuration file')
    parser.add_argument('--out_dir', metavar='text', default="", 
                        help='output directory for projected .root files')
    parser.add_argument('--pre', action='store_true', help='pre-process the AnRes.root')
    parser.add_argument('--sigma', action='store_true', help='get the sigma')
    parser.add_argument('--pre_sys', action='store_true', help='pre-process the AnRes.root for systematic')
    parser.add_argument('--skip_projection', '-sp', action='store_true', help='skip the projection')
    parser.add_argument("--suffix", "-s", metavar="text", default="", help="suffix for output files")
    args = parser.parse_args()

    if not args.pre and not args.sigma and not args.pre_sys:
        print('Please specify the action to perform.')
        sys.exit(1)

    with open(args.config_pre, 'r') as cfgPre:
        config = yaml.safe_load(cfgPre)
  
    # Load the configuration
    ptmins = config['ptmins']
    ptmaxs = config['ptmaxs']
    axestokeep = config['axestokeep']
    outputDir = args.out_dir if args.out_dir != "" else config['skim_out_dir'] 
    
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
        
    if args.pre_sys:
        pre_sys_process(config, ptmins, ptmaxs, centMin, centMax, axestokeep, outputDir)
        os.system(f'cp {args.config_pre} {outputDir}/pre_sys/AnRes')
