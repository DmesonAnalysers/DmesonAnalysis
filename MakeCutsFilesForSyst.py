'''
python script to create yaml files with set of cuts for cut-variation studies
'''

import os
import argparse
import copy
from itertools import product
import yaml

def get_variation_mult(edge, kind): # pylint: disable=too-many-return-statements
    if kind == 'loose_1':
        if edge == 'min':
            return -1.
        if edge == 'max':
            return 1.
    if kind == 'loose_2':
        if edge == 'min':
            return -2.
        if edge == 'max':
            return 2.
    if kind == 'tight_1':
        if edge == 'min':
            return 1.
        if edge == 'max':
            return -1.
    if kind == 'tight_2':
        if edge == 'min':
            return 2.
        if edge == 'max':
            return -2.
    return 0.

def check_value(new_value, cut_lim, histo_lim, edge):
    if edge == 'min':
        return cut_lim <= new_value < histo_lim
    if edge == 'max':
        return histo_lim < new_value <= cut_lim
    return False

def make_cuts():
    var_name = ['cospkphi3', 'cp', 'cpxy', 'd0', 'deltamKK', 'dl', 'dlxy', 'ndlxy', 'sigvtx', 'topo']
    var_key = ['CosPiKPhi3', 'CosP', 'CosPXY', 'd0', 'DeltaMassKK', 'DecL', 'DecLXY', 'NormDecLXY',
               'SigmaVtx', 'd0d0Exp']
    edge_to_vary = ['min', 'min', 'min', 'max', 'max', 'min', 'min', 'min', 'max', 'max']
    step_variation = [0.5, 0.5, 0.5, 10., 1., 10., 10., 1., 5., 0.5]
    histo_lims = [3.5, 100., 100., 0., 0., 105., 105., 10.5, 0., 0.]
    variation_kind = ['loose_1', 'loose_2', 'tight_1', 'tight_2']
    # [-1, -2, +1, +2]*step_variation added to central value

    in_dir = 'configfiles/'
    cut_file_central = 'cutset_3050_central_2018.yml'
    cut_file_loose = 'cutset_3050_loose_2018.yml'
    out_dir = 'configfiles/syst_cuts_Ds_3050/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(in_dir + cut_file_central, 'r') as cut_file_yml:
        cutset = yaml.load(cut_file_yml, yaml.FullLoader)
    with open(in_dir + cut_file_loose, 'r') as cut_file_loose_yml:
        cutset_loose = yaml.load(cut_file_loose_yml, yaml.FullLoader)

    cutset_central = copy.deepcopy(cutset)
    with open(out_dir + cut_file_central, 'w') as outfile:
        yaml.dump(cutset_central, outfile, default_flow_style=False)

    for name, key, edge, step, histo_lim in zip(var_name, var_key, edge_to_vary, step_variation, histo_lims):
        loose_values = cutset_loose['cutvars'][key][edge]
        for kind in variation_kind:
            cutset_mod = copy.deepcopy(cutset)
            mult_value = get_variation_mult(edge, kind)
            modified_list = []
            for value, cut_lim in zip(cutset_mod['cutvars'][key][edge], loose_values):
                new_value = value + step * mult_value
                if check_value(new_value, cut_lim, histo_lim, edge):
                    modified_list.append(new_value)
                else:
                    modified_list.append(value)
            cutset_mod['cutvars'][key][edge] = modified_list
            cut_file_mod = cut_file_central.replace('central_2018', name + '_' + kind)
            with open(out_dir + cut_file_mod, 'w') as outfile_mod:
                yaml.dump(cutset_mod, outfile_mod, default_flow_style=False)

def make_cuts_ml():
    var_key = ['ML_output_FD', 'ML_output_Bkg']
    var_tag = ['outFD', 'outBkg'] # used in file names to reduce length
    step_variation = [{"2": 0.01, "4": 0.01, "6": 0.01, "8": 0.01, "12": 0.01},
                      {"2": 0.002, "4": 0.002, "6": 0.002, "8": 0.002, "12": 0.002}]
    num_step_pos = 3
    num_step_neg = 3
    edge_to_vary = ['min', 'max']

    in_dir = 'configfiles/cutsets/Ds/pp/'
    cut_file_central = 'cutset_pp5TeV_FDen_conservative.yml'
    out_dir = 'configfiles/cutsets/Ds/pp/syst_cuts_TEST_GRID/'
    out_file_tag = 'cutset_pp5TeV_FDen'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
    with open(in_dir + cut_file_central, 'r') as cut_file_yml:
        cutset = yaml.load(cut_file_yml, yaml.FullLoader)

    #same number of steps for all variables
    neg_steps = [-i for i in range(1, num_step_neg + 1)]
    neg_steps = neg_steps[::-1]
    pos_steps = list(range(0, num_step_pos + 1))
    steps = neg_steps + pos_steps

    n_combinations = len(var_key)

    for prod_steps in product(steps, repeat=n_combinations):
        cutset_mod = copy.deepcopy(cutset)
        file_tag = str()
        for i, step in enumerate(prod_steps):
            modified_list = []
            cuts = cutset_mod['cutvars'][var_key[i]]
            for min_val, max_val, pt_min in zip(cuts['min'], cuts['max'], cutset_mod['cutvars']['Pt']['min']):
                if edge_to_vary[i] == 'min':
                    new_value = min_val + step * step_variation[i][f'{pt_min:.0f}']
                    if(new_value < 0. or new_value > max_val):
                        new_value = min_val
                else:
                    new_value = max_val + step * step_variation[i][f'{pt_min:.0f}']
                    if(new_value > 1. or new_value < min_val):
                        new_value = max_val
                modified_list.append(new_value)
            cuts[edge_to_vary[i]] = modified_list
            step_name = 'pos'
            if step < 0.:
                step_name = 'neg'
                step += num_step_neg + 1
            # more than 100 files unlikely
            name = f'_{var_tag[i]}_{step_name}{str(int(step)).zfill(2)}'
            file_tag += name
        cut_file_mod = f'{out_file_tag}_{file_tag}.yml'
        with open(out_dir + cut_file_mod, 'w') as outfile_mod:
            yaml.dump(cutset_mod, outfile_mod, default_flow_style=False)

def main():
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument("--ml", help="make cuts for ml", action="store_true")
    args = parser.parse_args()
    if args.ml:
        make_cuts_ml()
    else:
        make_cuts()

main()
