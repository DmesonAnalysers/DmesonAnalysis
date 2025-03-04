import yaml
import ROOT
import copy
import os
import numpy as np
import argparse
import itertools
from alive_progress import alive_bar # type: ignore
class DebugConfig:
    debug = False  # DebugConfig.debug = True to enable debug

def modify_yaml(input_config, output_dir, modifications_config):
    with open(input_config, 'r') as f:
        data = yaml.safe_load(f)

    with open(modifications_config, 'r') as f:
        modif = yaml.safe_load(f)

    infile_name = os.path.basename(input_config)
    infile_name = infile_name.replace(".yml", "")
    output_dir_folder = os.path.join(output_dir, "config_syst")
    if not os.path.exists(output_dir_folder):
        os.makedirs(output_dir_folder)

    #____________________________________________________________________________________
    # Collect all possible combinations of the modifications
    bkg_vn, invm_bins, mass_min, mass_max, bkg_func, \
    rebins, use_invm_bins, ptmins_maxs = ([] for i in range(8))
    if modif['BkgFuncVn']:
        bkg_vn = modif['BkgFuncVn']
    if modif['BkgFunc']:
        bkg_func = modif['BkgFunc']
    if modif['inv_mass_bins']:
        for step in modif['inv_mass_bins']['steps']:
            invm_bins.append(np.arange(modif['inv_mass_bins']['min'], modif['inv_mass_bins']['max'], step).tolist())
    if modif['MassMin']:
        mass_min = modif['MassMin']
    if modif['MassMax']:
        mass_max = modif['MassMax']
    if modif['Rebin']:
        rebins = modif['Rebin']
    if modif['use_inv_mass_bins']:
        use_invm_bins = modif['use_inv_mass_bins']
    if modif['ptmins'] and modif['ptmaxs']:
        for ptmin, ptmax in zip(modif['ptmins'], modif['ptmaxs']):
            ptmins_maxs.append([ptmin, ptmax])

    #____________________________________________________________________________________
    # Create a new yaml file for each combination of the modifications
    # Multitrial systematics
    is_multitrial = False
    for bkg in bkg_vn:
        for bkgf in bkg_func:
            for reb in rebins:
                for use_invm_bin in use_invm_bins:
                    for _, (invm_bin, step) in enumerate(zip(invm_bins, modif['inv_mass_bins']['steps'])):
                        for mmin in mass_min:
                            for mmax in mass_max:
                                new_data = data.copy()
                                for i in range(len(new_data['BkgFuncVn'])):
                                    new_data['BkgFuncVn'][i] = bkg
                                for i in range(len(new_data['BkgFunc'])):
                                    new_data['BkgFunc'][i] = bkgf
                                for i in range(len(new_data['Rebin'])):
                                    new_data['Rebin'][i] = reb
                                new_data['use_inv_mass_bins'] = use_invm_bin
                                for j in range(len(new_data['inv_mass_bins'])):
                                    new_data['inv_mass_bins'][j] = [l for l in invm_bin]
                                for i in range(len(new_data['MassMin'])):
                                    new_data['MassMin'][i] = mmin
                                for i in range(len(new_data['MassMax'])):
                                    new_data['MassMax'][i] = mmax

                                        
                                outfile_name = os.path.join(output_dir_folder, f'{infile_name}_bkgvn{bkg}_bkg{bkgf}_reb{reb}_invm{step}_useinvm{use_invm_bin}_min{mmin}_max{mmax}.yml')

                                with open(outfile_name, 'w') as f:
                                    yaml.dump(new_data, f, default_flow_style=False)
                                if not is_multitrial:
                                    is_multitrial = True
    if is_multitrial: # If there are multitrial systematics, return, since we don't need to do the rest
        return
   
    # Reso systematics
    for ptmin_max in ptmins_maxs:
        new_data = data.copy()
        new_data['ptmins'] = [ptmin_max[0]]
        new_data['ptmaxs'] = [ptmin_max[1]]
        new_data['use_inv_mass_bins'] = True

        outfile_name = os.path.join(output_dir_folder,
                                    f'{infile_name}_ptmin{ptmin_max[0]}_ptmax{ptmin_max[1]}.yml')

        with open(outfile_name, 'w') as f:
            yaml.dump(new_data, f, default_flow_style=False)


def cook_inv_mass_bins(nPtBins, terms_inv_mass_bins_lower, terms_inv_mass_bins_upper, steps):
    '''
    Cook the inv_mass_bins for each pt bin: range, min, max
    
    args:
    nPtBins: number of pt bins
    terms_inv_mass_bins_lower: list of first mass bins for each pt bin
    terms_inv_mass_bins_upper: list of last mass bins for each pt bin
    steps: list of steps of mass bins for each pt bin
    
    return:
    terms_inv_mass_bins: list of mass bins for each pt bin
    terms_MassMin: list of min mass bins for each pt bin
    terms_MassMax: list of max mass bins for each pt bin
    
    '''
    terms_inv_mass_bins, terms_MassMin, terms_MassMax = [], [], []
    for step in steps:
        for term_inv_mass_bins_lower in terms_inv_mass_bins_lower:
            for term_inv_mass_bins_upper in terms_inv_mass_bins_upper:
                term_inv_mass_bins, term_MassMin, term_MassMax = [], [], []
                for iPt in range(nPtBins):
                    term_inv_mass_bins.append(np.arange(term_inv_mass_bins_lower[iPt], term_inv_mass_bins_upper[iPt], step).tolist())
                    # [-1] means the current pt bin: [0] means the first mass bin for this term, [-1] means the last mass bin for this term
                    term_MassMin.append(term_inv_mass_bins[-1][0])
                    term_MassMax.append(term_inv_mass_bins[-1][-1])
                terms_inv_mass_bins.append(term_inv_mass_bins)
                terms_MassMin.append(term_MassMin)
                terms_MassMax.append(term_MassMax)
    return terms_inv_mass_bins, terms_MassMin, terms_MassMax

def find_threshold(nPtBins, default_values, threshold_values):
    lower_thresholds = []
    upper_thresholds = []
    center_values = []
    
    threshold_values = sorted(threshold_values)
    
    for iPt in range(nPtBins):
        default_value = default_values[iPt]
        center_values.append(default_value)
        lower_threshold, upper_threshold = default_value, default_value
        # if fit option can't be dependent on pt, then forcelly set the threshold
        for threshold_value in threshold_values:
            if threshold_value < default_value:
                lower_threshold = threshold_value
            elif threshold_value > default_value:
                upper_threshold = threshold_value
                break
        lower_thresholds.append(lower_threshold)
        upper_thresholds.append(upper_threshold)
    return lower_thresholds, upper_thresholds, center_values

def find_2threshold(nPtBins, default_values, threshold_values):
    lower_thresholds = []
    lower_thresholds_2 = []
    upper_thresholds = []
    upper_thresholds_2 = []
    center_values = []
    
    threshold_values = sorted(threshold_values)
    
    for iPt in range(nPtBins):
        default_value = default_values[iPt]
        lower_threshold, upper_threshold = default_value, default_value
        lower_threshold_2, upper_threshold_2 = default_value, default_value
        center_value = default_value
        # if fit option can't be dependent on pt, then forcelly set the threshold
        for threshold_value in threshold_values:
            if threshold_value < default_value:
                lower_threshold_2 = lower_threshold
                lower_threshold = threshold_value
            elif threshold_value > default_value:
                upper_threshold = threshold_value
                break
        for threshold_value in reversed(threshold_values):
            if threshold_value > default_value:
                upper_threshold_2 = upper_threshold
                upper_threshold = threshold_value
            elif threshold_value < default_value:
                lower_threshold = threshold_value
                break
        lower_thresholds_2.append(lower_threshold_2)
        lower_thresholds.append(lower_threshold)
        upper_thresholds_2.append(upper_threshold_2)
        upper_thresholds.append(upper_threshold)
        center_values.append(center_value)
    return lower_thresholds, upper_thresholds, lower_thresholds_2, upper_thresholds_2, center_values

def clean_flow_configs(flow_configs, output_dir):
    for config_name, config in flow_configs.items():
        config['out_dir'] = f'{output_dir}/trails/all_pt'
        config['suffix'] = config_name
        config['skim_out_dir'] = f'{output_dir}'
        config['minimisation']['correlated'] = False
        config['minimisation']['combined'] = True
        config['nworkers'] = 1
        config['FixSigma'] = 1
        config['FixSigmaFromFile'] = ''
        if config['minimisation'].get('skip_cuts', []):
            config['minimisation'].pop('skip_cuts')
        if config['minimisation'].get('systematics', []):
            config['minimisation'].pop('systematics')
    return flow_configs

def generate_flow_config_variations(flow_configs, multi_terms, multi_terms_name, debug=False):

    # The number of terms in multi_terms should be the same

    debug = DebugConfig.debug
    # define a new dict to store the new flow configs
    new_flow_configs = {}

    # key, value = next(iter(flow_configs.items()))
    # new_flow_configs[key] = value

    # define two lists to store the new flow configs and their names
    flow_configs_list = []
    flow_configs_name = []
    
    # loop over the flow configs
    for flow_config_name, flow_config in flow_configs.items():
        temp_flow_configs_name = []
        temp_flow_configs = []
        # loop for each fit option term
        for iTerm, (terms, terms_name) in enumerate(zip(multi_terms, multi_terms_name)):
            # loop for the variation range of the term
            for term_index, term in enumerate(terms):
                # if it is the first term, create a new flow config
                if iTerm == 0:
                    temp_flow_configs.append(copy.deepcopy(flow_config))
                    temp_flow_configs_name.append(flow_config_name + f'_{terms_name}-{term_index}')
                else:
                    temp_flow_configs_name[term_index] = temp_flow_configs_name[term_index] + f'_{terms_name}-{term_index}'
                # update the fit option term with the variation
                temp_flow_configs[term_index][terms_name] = term
        flow_configs_list.extend(temp_flow_configs)
        flow_configs_name.extend(temp_flow_configs_name)
    # update the flow configs dict with the new flow configs
    new_flow_configs.update(dict(zip(flow_configs_name, flow_configs_list)))

    if debug:
        print(f'After adding {multi_terms_name} variations ({len(terms)}): {len(new_flow_configs)}')
    return new_flow_configs

def generate_flow_config_variations_add(flow_configs, multi_terms, multi_terms_name, debug=False):
    debug = DebugConfig.debug
    
    # define two lists to store the new flow configs and their names
    flow_configs_list = []
    flow_configs_name = []
    
    # loop over the flow configs
    for flow_config_name, flow_config in flow_configs.items():
        temp_flow_configs_name = []
        temp_flow_configs = []
        # loop for each fit option term
        for iTerm, (terms, terms_name) in enumerate(zip(multi_terms, multi_terms_name)):
            # loop for the variation range of the term
            for term_index, term in enumerate(terms):
                # if it is the first term, create a new flow config
                if iTerm == 0:
                    temp_flow_configs.append(copy.deepcopy(flow_config))
                    temp_flow_configs_name.append(flow_config_name + f'_{terms_name}-{term_index}')
                else:
                    temp_flow_configs_name[term_index] = temp_flow_configs_name[term_index] + f'_{terms_name}-{term_index}'
                # update the fit option term with the variation
                temp_flow_configs[term_index][terms_name] = term
        flow_configs_list.extend(temp_flow_configs)
        flow_configs_name.extend(temp_flow_configs_name)
    # update the flow configs dict with the new flow configs
    flow_configs.update(dict(zip(flow_configs_name, flow_configs_list)))
    
    if debug:
        print(f'After cumulatively adding {multi_terms_name} variations ({len(terms)}): {len(flow_configs)}')
    return flow_configs

def combination_fit_option(config_flow_name, cfg_flow, nPtBins, cfg_mod, output_dir):
    '''
    TEMPLETE:

    terms = [option1, option2, ...]
    terms_name = string of the option name
    multi_terms = [terms1, terms2, ...]
    multi_terms_name = [terms_name1, terms_name2, ...]
    
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Example:
    new_flow_configs = generate_flow_config_variations(old_flow_configs, multi_terms=[multi_terms], multi_terms_name=[multi_terms_name])
    
    len(new_flow_configs) == len(terms1) * len(old_flow_configs)
    
    ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    Example:
    new_flow_configs = generate_flow_config_variations_add(old_flow_configs, multi_terms=[multi_terms], multi_terms_name=[multi_terms_name])
    
    len(new_flow_configs) == len(terms1) * len(old_flow_configs) + len(old_flow_configs)
    
    '''
    DebugConfig.debug = True
    
    fit_option_dict = {}
    fit_opts_dependent_pt = []

    for opt_name, opt_values in cfg_mod['fit_options'].items():
        fit_option_dict[opt_name] = opt_values

    flow_configs = {} # key: config name, value: config dict
    flow_configs_default_mass_bins = {} # key: config name, value: config dict
    flow_configs[config_flow_name] = copy.deepcopy(cfg_flow)
    flow_configs_default_mass_bins[config_flow_name + '-default'] = copy.deepcopy(cfg_flow)

    # sigma upper, median, lower
    if fit_option_dict['Sigma']['FixSigma'] == 0:
        terms_FixSigma = [[0 for _ in range(nPtBins)] for _ in range(3)]
        flow_configs = generate_flow_config_variations(flow_configs, multi_terms=[terms_FixSigma], multi_terms_name=['FixSigma'])
    else:
        # sigma_file = ROOT.TFile(fit_option_dict['Sigma']['FixSigmaFromFile'])
        # hSigma = sigma_file.Get('hSigmaSimFit')
        # hSigma.SetDirectory(0)
        # sigma_file.Close()
        terms_FixSigma = [1 for _ in range(3)]
        terms_FixSigmaFromFile = ['' for  _ in range(3)]
        if fit_option_dict['Sigma']['FixSigma'] == -1:
            # Sigma_uppers = [hSigma.GetBinContent(iPt+1) + hSigma.GetBinError(iPt+1) for iPt in range(nPtBins)]
            # for iSig, Sigma_upper in enumerate(Sigma_uppers):
            #     if Sigma_upper > 0.055:
            #         Sigma_uppers[iSig] = 0.055
            # Sigma_median = [hSigma.GetBinContent(iPt+1) for iPt in range(nPtBins)]
            # Sigma_lowers = [hSigma.GetBinContent(iPt+1) - hSigma.GetBinError(iPt+1) for iPt in range(nPtBins)]
            Sigma_uppers = fit_option_dict['Sigma']['upper']
            Sigma_median = fit_option_dict['Sigma']['med']
            Sigma_lowers = fit_option_dict['Sigma']['lower']
            terms_sigma = [Sigma_lowers, Sigma_median, Sigma_uppers]
            flow_configs = generate_flow_config_variations(flow_configs, multi_terms=[terms_sigma], multi_terms_name=['Sigma'])
            flow_configs_default_mass_bins = generate_flow_config_variations_add(flow_configs_default_mass_bins, multi_terms=[terms_sigma, terms_FixSigmaFromFile], multi_terms_name=['Sigma'])
        # else:
        #     Sigma_uppers = [hSigma.GetBinContent(iPt+1) * (1 + fit_option_dict['Sigma']['FixSigma']) for iPt in range(nPtBins)]
        #     for iSig, Sigma_upper in enumerate(Sigma_uppers):
        #         if Sigma_upper > 0.055:
        #             Sigma_uppers[iSig] = 0.055
        #     Sigma_median = [hSigma.GetBinContent(iPt+1) for iPt in range(nPtBins)]
        #     Sigma_lowers = [hSigma.GetBinContent(iPt+1) * (1 - fit_option_dict['Sigma']['FixSigma']) for iPt in range(nPtBins)]
        #     terms_sigma = [Sigma_lowers, Sigma_median, Sigma_uppers]
            # flow_configs = generate_flow_config_variations(flow_configs, multi_terms=[terms_sigma], multi_terms_name=['Sigma'])
            # flow_configs_default_mass_bins = generate_flow_config_variations_add(flow_configs_default_mass_bins, multi_terms=[terms_sigma], multi_terms_name=['Sigma'])
    fit_opts_dependent_pt.append('Sigma')

    # delete the original config
    flow_configs_default_mass_bins.pop(config_flow_name + '-default', None)

    # bkg function for vn
    terms_bkg_func_vn = [[bkg_func_vn for _ in range(nPtBins)] for bkg_func_vn in fit_option_dict['BkgFuncVn']]
    flow_configs = generate_flow_config_variations(flow_configs, multi_terms=[terms_bkg_func_vn], multi_terms_name=['BkgFuncVn'])
    flow_configs_default_mass_bins = generate_flow_config_variations_add(flow_configs_default_mass_bins, multi_terms=[terms_bkg_func_vn], multi_terms_name=['BkgFuncVn'])
    fit_opts_dependent_pt.append('BkgFuncVn')
    
    # bkg function for mass
    terms_bkg_func = [[bkg_func for _ in range(nPtBins)] for bkg_func in fit_option_dict['BkgFunc']]
    flow_configs = generate_flow_config_variations(flow_configs, multi_terms=[terms_bkg_func], multi_terms_name=['BkgFunc'])
    flow_configs_default_mass_bins = generate_flow_config_variations_add(flow_configs_default_mass_bins, multi_terms=[terms_bkg_func], multi_terms_name=['BkgFunc'])
    fit_opts_dependent_pt.append('BkgFunc')
    
    # rebin
    terms_rebin = find_threshold(nPtBins, cfg_flow['Rebin'], fit_option_dict['Rebin'])
    flow_configs = generate_flow_config_variations(flow_configs, multi_terms=[terms_rebin], multi_terms_name=['Rebin'])
    flow_configs_default_mass_bins = generate_flow_config_variations_add(flow_configs_default_mass_bins, multi_terms=[terms_rebin], multi_terms_name=['Rebin'])
    fit_opts_dependent_pt.append('Rebin')
    
    # inv mass bins
    default_inv_mass_bins_lower = [cfg_flow['inv_mass_bins'][iPt][0] for iPt in range(nPtBins)]
    default_inv_mass_bins_upper = [cfg_flow['inv_mass_bins'][iPt][-1] for iPt in range(nPtBins)]     
    terms_inv_mass_bins_lower = find_threshold(nPtBins, default_inv_mass_bins_lower, fit_option_dict['inv_mass_bins']['MassMin'])
    terms_inv_mass_bins_upper = find_threshold(nPtBins, default_inv_mass_bins_upper, fit_option_dict['inv_mass_bins']['MassMax'])
    terms_inv_mass_bins, terms_MassMin, terms_MassMax = cook_inv_mass_bins(nPtBins, terms_inv_mass_bins_lower, terms_inv_mass_bins_upper, fit_option_dict['inv_mass_bins']['steps'])
    flow_configs = generate_flow_config_variations(flow_configs=flow_configs, multi_terms=[terms_inv_mass_bins, terms_MassMin, terms_MassMax], multi_terms_name=['inv_mass_bins', 'MassMin', 'MassMax'])
    fit_opts_dependent_pt.append('MassMin')
    fit_opts_dependent_pt.append('MassMax')
    fit_opts_dependent_pt.append('inv_mass_bins')

    # use inv mass bins #! For D0, use_inv_mass_bins is same with rebin = 1 ==> commented
    # terms_use_inv_mass_bins = fit_option_dict['use_inv_mass_bins']
    # flow_configs = generate_flow_config_variations(flow_configs, multi_terms=[terms_use_inv_mass_bins], multi_terms_name=['use_inv_mass_bins'])
    # flow_configs_default_mass_bins = generate_flow_config_variations_add(flow_configs_default_mass_bins, multi_terms=[terms_use_inv_mass_bins], multi_terms_name=['use_inv_mass_bins'])

    flow_configs.update(flow_configs_default_mass_bins)
    return clean_flow_configs(flow_configs, output_dir), fit_opts_dependent_pt

def slice_single_pt(flow_configs, nPtBins, fit_opts_dependent_pt, output_dir):
    
    flow_configs_pt = {i: [] for i in range(nPtBins)}
    
    # loop over the flow configs
    for config_name, config in flow_configs.items():
        # get the header of the config name
        config_name_header = config_name.split('_')[0]
        # loop over the pt bins
        for iPt in range(nPtBins):
            # create a new config for each pt bin
            new_config = copy.deepcopy(config)
            # create a new config name for each pt bin
            new_config_name = config_name_header + f'_pt{config["ptmins"][iPt]}_{config["ptmaxs"][iPt]}'
            
            new_config['ptmins'] = [config['ptmins'][iPt]]
            new_config['ptmaxs'] = [config['ptmaxs'][iPt]]
            new_config['out_dir'] = f'{output_dir}/trails/pt_{int(config["ptmins"][iPt]*10)}_{int(config["ptmaxs"][iPt]*10)}'
            
            # update the options, that are dependent on pt, with the values for the current pt bin
            for fit_opt in fit_opts_dependent_pt:
                new_config[fit_opt] = [config[fit_opt][iPt]]
                if fit_opt != 'inv_mass_bins':
                    new_config_name = new_config_name + f'_{fit_opt}-{config[fit_opt][iPt]}'
            
            # split the bdt cut
            new_config['cut_variation']["uncorr_bdt_cut"]['bkg_max'] = [config['cut_variation']["uncorr_bdt_cut"]['bkg_max'][iPt]]
            uncorr_sig_cut = config['cut_variation']["uncorr_bdt_cut"]['sig'][iPt]
            new_config['cut_variation']["uncorr_bdt_cut"]['sig'] = [{'min': uncorr_sig_cut['min'], 'max': uncorr_sig_cut['max']}]
            new_config['cut_variation']["corr_bdt_cut"]['bkg_max'] = [config['cut_variation']["corr_bdt_cut"]['bkg_max'][iPt]]
            # in principle, does not need the corr_bdt_cut
            new_config['cut_variation']["corr_bdt_cut"]['sig']['max'] = [config['cut_variation']["corr_bdt_cut"]['sig']['max'][iPt]]
            new_config['cut_variation']["corr_bdt_cut"]['sig']['min'] = [config['cut_variation']["corr_bdt_cut"]['sig']['min'][iPt]]
            new_config['cut_variation']["corr_bdt_cut"]['sig']['step'] = [config['cut_variation']["corr_bdt_cut"]['sig']['step'][iPt]]
            
            flow_configs_pt[iPt].append({new_config_name: new_config})

    return flow_configs_pt

def produce_pre_config(cfg_flow, cfg_mod, output_dir):
    
    pre_config_dict = {}
    pre_config_dict['flow_files'] = cfg_flow['anresdir']
    pre_config_dict['ptmins'] = cfg_flow['ptmins']
    pre_config_dict['ptmaxs'] = cfg_flow['ptmaxs']
    pre_config_dict['centrality'] = cfg_flow['centrality']
    pre_config_dict['skim_out_dir'] = output_dir
    pre_config_dict['bdt_cut'] = {}
    pre_config_dict['bdt_cut']['bkg_cuts'] = [max(cfg_flow['cut_variation']['uncorr_bdt_cut']['bkg_max'][iPt]) for iPt in range(len(cfg_flow['ptmins']))]
    pre_config_dict['bdt_cut']['sig_mins'] = [cfg_flow['cut_variation']['uncorr_bdt_cut']['sig'][iPt]['min'] for iPt in range(len(cfg_flow['ptmins']))]
    pre_config_dict['bdt_cut']['sig_maxs'] = [cfg_flow['cut_variation']['uncorr_bdt_cut']['sig'][iPt]['max'] for iPt in range(len(cfg_flow['ptmins']))]
    # different bkg cut for a dedicated pt bin
    pre_config_dict['axestokeep'] = ['Mass', 'sp']
    pre_config_dict['RebinSparse'] = {}
    pre_config_dict['RebinSparse']['Mass'] = -1
    pre_config_dict['RebinSparse']['sp'] = -1
    pre_config_dict['RebinSparse']['score_bkg'] = -1
    pre_config_dict['RebinSparse']['score_FD'] = -1
    pre_config_dict['RebinSparse']['Pt'] = -1
    pre_config_dict['RebinSparse']['cent'] = -1

    return pre_config_dict

def modify_yaml_bdt(config_flow, config_mod, output_dir):
    
    with open(config_flow, 'r') as CfgFlow:
        cfg_flow = yaml.safe_load(CfgFlow)

    cfg_flow['skim_out_dir'] = f'{output_dir}'
    
    nPtBins = len(cfg_flow['ptmins'])

    with open(config_mod, 'r') as CfgMod:
        cfg_mod = yaml.safe_load(CfgMod)
        
    config_flow = os.path.basename(config_flow)
    config_flow = config_flow.replace(".yml", "")

    # all possible combinations of the modifications: {config_name: config}
    # # pt dependent fit options
    flow_configs, fit_opts_dependent_pt = combination_fit_option(config_flow, cfg_flow, nPtBins, cfg_mod, output_dir)

    with alive_bar(len(flow_configs), title='Writing yaml files, which contain all pt bins') as bar:
        os.makedirs(f'{output_dir}/config_sys/all_pt', exist_ok=True)
        for config_name, config in flow_configs.items():
            outfile_name = os.path.join(f'{output_dir}/config_sys/all_pt', f'{config_name}.yml')
            with open(outfile_name, 'w') as f:
                yaml.dump(config, f, default_flow_style=False)
            bar()

    # # slice the flow configs into single pt bins
    # flow_configs_pt = slice_single_pt(flow_configs, nPtBins, fit_opts_dependent_pt, output_dir)

    # with alive_bar(nPtBins * len(flow_configs_pt[0]), title='Writing yaml files, which contain single pt bins') as bar:
    #     for iPt, (ptmin, ptmax) in enumerate(zip(cfg_flow['ptmins'], cfg_flow['ptmaxs'])):
    #         print(f'Writing yaml files for pT bin {iPt}')
    #         os.makedirs(f'{output_dir}/config_sys/pt_{int(ptmin*10)}_{int(ptmax*10)}', exist_ok=True)
    #         for flow_configs_single_pt in flow_configs_pt[iPt]:
    #             for config_name, config in flow_configs_single_pt.items():
    #                 outfile_name = os.path.join(f'{output_dir}/config_sys/pt_{int(ptmin*10)}_{int(ptmax*10)}', f'{config_name}.yml')
    #                 with open(outfile_name, 'w') as f:
    #                     yaml.dump(config, f, default_flow_style=False)
    #                 bar()

    # #____________________________________________________________________________________________________________________________________________________
    # # config_pre for sysmatical
    pre_config = produce_pre_config(cfg_flow, cfg_mod, output_dir)
    
    with open(f'{output_dir}/config_sys/config_pre.yml', 'w') as f:
        yaml.dump(pre_config, f, default_flow_style=False)
    
    # #____________________________________________________________________________________________________________________________________________________
    # config_pre for uncorrelated and correlated
    os.makedirs(f'{output_dir}/config_sys/reference', exist_ok=True)

    # # uncorrelated config
    with open(f'{output_dir}/config_sys/reference/config_uncorrelated.yml', 'w') as f:
        cfg_uncorr = copy.deepcopy(cfg_flow)
        cfg_uncorr['out_dir'] = f'{output_dir}/pre_sys'
        cfg_uncorr['suffix'] = 'uncorr'
        cfg_uncorr['skim_out_dir'] = f'{output_dir}'
        cfg_uncorr['minimisation']['correlated'] = False
        cfg_uncorr['minimisation']['combined'] = False
        cfg_uncorr['nworkers'] = 25
        if cfg_uncorr['minimisation'].get('skip_cuts', []):
            cfg_uncorr['minimisation'].pop('skip_cuts')
        if cfg_uncorr['minimisation'].get('systematics', []):
            cfg_uncorr['minimisation'].pop('systematics')
        yaml.dump(cfg_uncorr, f, default_flow_style=False)

    # # correlated config correlated
    with open(f'{output_dir}/config_sys/reference/config_correlated.yml', 'w') as f:
        cfg_corr = copy.deepcopy(cfg_flow)
        cfg_corr['out_dir'] = f'{output_dir}/pre_sys'
        cfg_corr['suffix'] = 'corr'
        cfg_corr['skim_out_dir'] = f'{output_dir}'
        cfg_corr['minimisation']['correlated'] = True
        cfg_corr['minimisation']['combined'] = False
        cfg_corr['nworkers'] = 12
        yaml.dump(cfg_corr, f, default_flow_style=False)

    # # combined config
    with open(f'{output_dir}/config_sys/reference/config_combined.yml', 'w') as f:
        cfg_comb = copy.deepcopy(cfg_flow)
        cfg_comb['out_dir'] = f'{output_dir}/pre_sys'
        cfg_comb['suffix'] = 'combined'
        cfg_comb['minimisation']['correlated'] = False
        cfg_comb['minimisation']['combined'] = True
        cfg_comb['minimisation']['correlatedPath'] = f'{output_dir}/pre_sys/cutvar_corr'
        yaml.dump(cfg_comb, f, default_flow_style=False)
    pass

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('input_config', metavar='text', default='config_Ds_Fit.yml')
    parser.add_argument('--modifications_config', "-m", metavar='text', default='')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    parser.add_argument("--multitrial_bdt", "-mb", action="store_true", default=False,
                        help="multitrial systematics for BDT")
    args = parser.parse_args()

    if not args.multitrial_bdt:
        modify_yaml(args.input_config,
                    args.outputdir,
                    args.modifications_config)
    else:
        modify_yaml_bdt(args.input_config,
                         args.modifications_config,
                         args.outputdir)
