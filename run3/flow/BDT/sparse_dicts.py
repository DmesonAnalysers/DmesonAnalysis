import ROOT
import yaml
from ROOT import TFile # pyright: ignore # type: ignore

def get_sparses(config, get_data, get_mc_reco, get_mc_gen, anres_files=[], preprocessed=False, preprocess_dir='', debug=False):
    """Load the sparses and axes infos

    Args:
        config (dict): the flow config dictionary
        get_data (bool): whether to get the data
        get_mc_reco (bool): whether to get the mc reco level
        get_mc_gen (bool): whether to get the mc gen level
        anres_files (list, optional): a list of AnRes.root. Defaults to [].
        preprocessed (bool, optional): whether to use the pre-selected AnRes.root. Defaults to False.
        preprocess_dir (str, optional): path of the config_pre.yaml. Defaults to ''.
        debug (bool, optional): print debug info. Defaults to False.

    Outputs:
        sparsesFlow: thnSparse in the flow task
        sparsesReco: thnSparse of reco level from the D meson task
        sparsesGen: thnSparse of gen level from the D meson task
        axes_dict (dict): dictionary of the axes for each sparse
    """
    
    
    sparsesFlow, sparsesReco, sparsesGen, axes_dict = {}, {}, {}, {}    
    
    if get_data:
        if preprocessed:

            with open(f"{config['skim_out_dir']}/config_pre.yml", 'r') as CfgPre:
                config_pre = yaml.safe_load(CfgPre)

            axes_dict['Flow'] = {ax: iax for iax, ax in enumerate(config_pre['axestokeep'])}
            for ptmin, ptmax in zip(config['ptmins'], config['ptmaxs']):
                print(f"Loading flow sparse from file: {preprocess_dir}/AnalysisResults_pt_{int(ptmin*10)}_{int(ptmax*10)}.root")
                infileflow = TFile(f"{preprocess_dir}/AnalysisResults_pt_{int(ptmin*10)}_{int(ptmax*10)}.root")
                sparsesFlow[f'Flow_{ptmin*10}_{ptmax*10}'] = infileflow.Get('hf-task-flow-charm-hadrons/hSparseFlowCharm')
                infileflow.Close()
        else:
            axes_dict['Flow'] = {
                'Mass': 0,
                'Pt': 1,
                'cent': 2,
                'sp': 3,
                'score_bkg': 4,
                'score_FD': 5,
                'occ': 6
            }
            # REVIEW: I would suggest to separete the config_flow and config_pre
            # and load the flow files from the arguments
            for ifile, file in enumerate(anres_files):
                print(f"Loading flow sparse from file: {file}")
                infileflow = TFile(file)
                sparsesFlow[f'Flow_{ifile}'] = infileflow.Get('hf-task-flow-charm-hadrons/hSparseFlowCharm')
                infileflow.Close()
       
    if get_mc_gen or get_mc_reco:
        print(f"Loading mc sparse from: {config['eff_filename']}")
        infiletask = TFile(config['eff_filename'])
    
    if get_mc_reco: 
        if config['Dmeson'] == 'Dzero':
            axes_reco = {
                'score_bkg': 0,
                'score_prompt': 1,
                'score_FD': 2,
                'Mass': 3,
                'Pt': 4,
                'y': 5,
                'cand_type': 6,
                'pt_bmoth': 7,
                'origin': 8,
                'npvcontr': 9,
                'cent': 10,
                'occ': 11,
            }
            print(f"infiletask: {infiletask}")
            sparsesReco['RecoPrompt'] = infiletask.Get('hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type')
            print(f"sparsesReco['RecoPrompt']: {sparsesReco['RecoPrompt']}")
            print('\n')
            sparsesReco['RecoPrompt'].GetAxis(axes_reco['origin']).SetRange(2, 2)  # make sure it is prompt
            sparsesReco['RecoPrompt'].GetAxis(axes_reco['cand_type']).SetRange(1, 2) # make sure it is signal
            axes_dict['RecoPrompt'] = axes_reco
            sparsesReco['RecoFD'] = infiletask.Get('hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type')
            sparsesReco['RecoFD'].GetAxis(axes_reco['origin']).SetRange(3, 3)  # make sure it is non-prompt
            sparsesReco['RecoFD'].GetAxis(axes_reco['cand_type']).SetRange(1, 2)  # make sure it is signal
            axes_dict['RecoFD'] = axes_reco
            sparsesReco['RecoRefl'] = infiletask.Get('hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type')
            sparsesReco['RecoRefl'].GetAxis(axes_reco['cand_type']).SetRange(3, 4)  # make sure it is reflection
            axes_dict['RecoRefl'] = axes_reco
            sparsesReco['RecoReflPrompt'] = infiletask.Get('hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type')
            sparsesReco['RecoReflPrompt'].GetAxis(axes_reco['cand_type']).SetRange(3, 4)  # make sure it is reflection
            sparsesReco['RecoReflPrompt'].GetAxis(axes_reco['origin']).SetRange(2, 2) # make sure it is prompt
            axes_dict['RecoReflPrompt'] = axes_reco
            sparsesReco['RecoReflFD'] = infiletask.Get('hf-task-d0/hBdtScoreVsMassVsPtVsPtBVsYVsOriginVsD0Type')
            sparsesReco['RecoReflFD'].GetAxis(axes_reco['cand_type']).SetRange(3, 4)  # make sure it is reflection
            sparsesReco['RecoReflFD'].GetAxis(axes_reco['origin']).SetRange(3, 3) # make sure it is FD
            axes_dict['RecoReflFD'] = axes_reco
            #TODO: safety checks for Dmeson reflecton and secondary peak
        elif config['Dmeson'] == 'Dplus':
            sparsesReco['RecoPrompt'] = infiletask.Get('hf-task-dplus/hSparseMassPrompt')
            axes_dict['RecoPrompt'] = {
                'Mass': 0,
                'Pt': 1,
                'score_bkg': 2,
                'score_prompt': 3,
                'score_FD': 4,
                'cent': 5,
                'occ': 6,
            }
            sparsesReco['RecoFD'] = infiletask.Get('hf-task-dplus/hSparseMassFD')
            axes_dict['RecoFD'] = {
                'Mass': 0,
                'Pt': 1,
                'pt_bmoth': 2,
                'flag_bhad': 3,
                'score_bkg': 4,
                'score_prompt': 5,
                'score_FD': 6,
                'cent': 7,
                'occ': 8
            }
        elif config['Dmeson'] == 'Ds':
            sparsesReco['RecoPrompt'] = infiletask.Get('hf-task-ds/MC/Ds/Prompt/hSparseMass')
            axes_dict['RecoPrompt'] = {
                'Mass': 0,
                'Pt': 1,
                'cent': 3,
                'npvcontr': 4,
                'score_bkg': 5,
                'score_prompt': 6,
                'score_FD': 7,
                'occ': 8,
            }
            sparsesReco['RecoFD'] = infiletask.Get('hf-task-ds/MC/Ds/NonPrompt/hSparseMass')
            axes_dict['RecoFD'] = {
                'Mass': 0,
                'Pt': 1,
                'cent': 2,
                'score_bkg': 3,
                'score_prompt': 4,
                'score_FD': 5,
                'pt_bmoth': 6,
                'flag_bhad': 7,
                'occ': 8
            }

    if get_mc_gen: 
        print(f"Loading mc gen sparse from: {config['eff_filename']}")
        infiletask = TFile(config['eff_filename'])
        if config['Dmeson'] == 'Dzero':
            axes_gen = {
                'Pt': 0,
                'pt_bmoth': 1,
                'y': 2,
                'origin': 3,
                'npvcontr': 4,
                'cent': 5,
                'occ': 6
            }
            sparsesGen['GenPrompt'] = infiletask.Get('hf-task-d0/hSparseAcc')
            sparsesGen['GenPrompt'].GetAxis(axes_gen['origin']).SetRange(2, 2)  # make sure it is prompt
            axes_dict['GenPrompt'] = axes_gen 
            print(f"sparseGenPrompt: {sparsesGen['GenPrompt']}")
            sparsesGen['GenFD'] = infiletask.Get('hf-task-d0/hSparseAcc')
            sparsesGen['GenFD'].GetAxis(axes_gen['origin']).SetRange(3, 3)  # make sure it is non-prompt
            axes_dict['GenFD'] = axes_gen
            print(f"sparseGenFD: {sparsesGen['GenFD']}")
            #TODO: safety checks for Dmeson reflecton and secondary peak
        elif config['Dmeson'] == 'Dplus':
            sparsesGen['GenPrompt'] = infiletask.Get('hf-task-dplus/hSparseMassGenPrompt')
            axes_dict['GenPrompt'] = {
                'Pt': 0,
                'y': 1,
                'cent': 2,
                'occ': 3
            }   
            print(f"sparseGenPrompt: {sparsesGen['GenPrompt']}")
            sparsesGen['GenFD'] = infiletask.Get('hf-task-dplus/hSparseMassGenFD')
            axes_dict['GenFD'] = {
                'Pt': 0,
                'y': 1,
                'pt_bmoth': 2,
                'flag_bhad': 3,
                'cent': 4,
                'occ': 5
            }
            print(f"sparseGenFD: {sparsesGen['GenFD']}")
        elif config['Dmeson'] == 'Ds':
            sparsesGen['GenPrompt'] = infiletask.Get('hf-task-ds/MC/Ds/Prompt/hSparseGen')
            axes_dict['GenPrompt'] = {
                'Pt': 0,
                'y': 1,
                'npvcontr': 2,
                'cent': 3,
                'occ': 4
            }   
            print(f"sparseGenPrompt: {sparsesGen['GenPrompt']}")
            sparsesGen['GenFD'] = infiletask.Get('hf-task-ds/MC/Ds/NonPrompt/hSparseGen')
            axes_dict['GenFD'] = {
                'Pt': 0,
                'y': 1,
                'cent': 2,
                'pt_bmoth': 3,
                'flag_bhad': 4,
                'occ': 5
            }
            print(f"sparseGenFD: {sparsesGen['GenFD']}")

    if get_mc_gen or get_mc_reco:
        infiletask.Close()

    print(f"Loaded sparses!")
    if debug:
        print('\n')
        print('###############################################################')
        for key, value in axes_dict.items():
            print(f"{key}:")
            for sub_key, sub_value in value.items():
                print(f"    {sub_key}: {sub_value}")
        print('###############################################################')
        print('\n')

    return sparsesFlow, sparsesReco, sparsesGen, axes_dict
