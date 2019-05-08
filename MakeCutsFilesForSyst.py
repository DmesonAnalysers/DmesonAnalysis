import os
import copy
import yaml

def get_variation_mult(edge, kind):
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
    var_name = ['cospkphi3', 'cp', 'cpxy', 'd0', 'deltamKK', 'dl', 'dlxy', 'ndlxy',
                'sigvtx', 'topo']
    var_key = ['CosPiKPhi3', 'CosP', 'CosPXY', 'd0', 'DeltaMassKK', 'DecL', 'DecLXY', 'NormDecLXY',
               'SigmaVtx', 'd0d0Exp']
    edge_to_vary = ['min', 'min', 'min', 'max', 'max', 'min', 'min', 'min', 'max', 'max']
    step_variation = [0.5, 0.5, 0.5, 10., 1., 10., 10., 1., 5., 0.5]
    histo_lims = [3.5, 100., 100., 0., 0., 105., 105., 10.5, 0., 0.]
    variation_kind = ['loose_1', 'loose_2', 'tight_1', 'tight_2'] # [-1, -2, +1, +2]*step_variation
                                                                  # added to central value
    in_dir = 'configfiles/'
    cut_file_central = 'cutset_3050_central_2018.yml'
    cut_file_loose = 'cutset_3050_loose_2018.yml'
    out_dir = 'configfiles/syst_cuts_Ds_3050/'

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    with open(in_dir + cut_file_central, 'r') as cut_file_yml:
        cutset = yaml.load(cut_file_yml)
    with open(in_dir + cut_file_loose, 'r') as cut_file_loose_yml:
        cutset_loose = yaml.load(cut_file_loose_yml)

    cutset_central = copy.deepcopy(cutset)
    with open(out_dir + cut_file_central, 'w') as outfile:
        yaml.dump(cutset_central, outfile, default_flow_style=False)

    for name, key, edge, step, histo_lim in zip(var_name, var_key, edge_to_vary,
                                                step_variation, histo_lims):
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

make_cuts()
