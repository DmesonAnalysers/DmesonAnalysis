'''
Script to create a yaml file with a set of cuts for ML
python3 make_yaml_for_ml.py config_flow.yml -o path/to/output -s text

'''
import yaml
import argparse
import os
import numpy as np

def make_combination(mass_axis, pt_axis, bkg_axis, sig_axis, ptmins, ptmaxs, sig_cut_file, bkg_cut_file):
    #TODO: add comments
    combinations = {}
    for iFile in range((len(sig_cut_file))):
        combinations[iFile] = {
            'icutset': iFile,
            'cutvars': {
                'InvMass': {
                    'axisnum': mass_axis,
                    'name': 'inv_mass'
                },
                'Pt': {
                    'axisnum': pt_axis,
                    'min': [i for i in ptmins],
                    'max': [j for j in ptmaxs],
                    'name': 'pt_cand'
                },
                'ML_output_Bkg': {
                    'axisnum': bkg_axis,
                    'min': [-1.0 for _ in range(len(bkg_cut_file[iFile]))],
                    'max': [float(j) for j in bkg_cut_file[iFile]],
                    'name': 'ML_output_Bkg'
                },
                'ML_output_FD': {
                    'axisnum': sig_axis,
                    'min': [float(i) for i in sig_cut_file[iFile]],
                    'max': [1.1 for _ in range(len(sig_cut_file[iFile]))],
                    'name': 'ML_output_FD'
                }
            }
        }
    return combinations

def make_yaml(flow_config, outputdir, suffix):
    with open(flow_config, 'r') as f:
        input = yaml.safe_load(f)

    os.makedirs(outputdir, exist_ok=True)

    # load the variable from the input config
    # axis
    axis = input['axes_mc']
    mass_axis = axis['mass']
    pt_axis = axis['pt']
    bkg_axis = axis['bdt_bkg']
    sig_axis = axis['bdt_sig']

    # pt
    ptmins = input['ptmins']
    ptmaxs = input['ptmaxs']

    # cut variation
    bkg_cut_mins = input['cut_variation']['bdt_cut']['bkg']['min']
    bkg_cut_maxs = input['cut_variation']['bdt_cut']['bkg']['max']
    bkg_cut_steps = input['cut_variation']['bdt_cut']['bkg']['step']
    sig_cut_mins = input['cut_variation']['bdt_cut']['sig']['min']
    sig_cut_maxs = input['cut_variation']['bdt_cut']['sig']['max']
    sig_cut_steps = input['cut_variation']['bdt_cut']['sig']['step']

    ## safely check
    if len(ptmins) != len(ptmaxs):
        raise ValueError(f'''The number of pt bins({len(ptmins)}, {len(ptmaxs)}),
                         bkg cuts({len(bkg_cut_mins)}, {len(bkg_cut_maxs)}),
                         and sig cuts({len(sig_cut_mins)}, {len(sig_cut_maxs)}) are not the same''')

    # create a new yaml file for each combination
    nCut = len(np.arange(sig_cut_mins[0], sig_cut_maxs[0], sig_cut_steps[0]))
    bkg_cut_pt, sig_cut_pt = {}, {}

    for ipt in range(len(ptmins)):
        sig_cut_pt[ipt] = list(np.arange(sig_cut_mins[ipt], sig_cut_maxs[ipt], sig_cut_steps[ipt]))
        bkg_cut_pt[ipt] = list(np.arange(bkg_cut_mins[ipt], bkg_cut_maxs[ipt], bkg_cut_steps[ipt])) if bkg_cut_steps[ipt] != 0 else [bkg_cut_maxs[ipt]] * nCut

    sig_cut_file = {iFile: [sig_cut_pt[ipt][iFile] for ipt in range(len(ptmins))] for iFile in range(nCut)}
    bkg_cut_file = {iFile: [bkg_cut_pt[ipt][iFile] for ipt in range(len(ptmins))] for iFile in range(nCut)}

    combinations = make_combination(mass_axis, pt_axis, bkg_axis, sig_axis, ptmins, ptmaxs, sig_cut_file, bkg_cut_file)

    os.makedirs(f'{outputdir}/config', exist_ok=True)
    for iFile in range(nCut):
        with open(f'{outputdir}/config/cutset_{suffix}_{iFile:02}.yml', 'w') as file:
            yaml.dump(combinations[iFile], file, default_flow_style=False)
    with open(f'{outputdir}/config/Output.yml', 'w') as file:
        yaml.dump("#output results", file, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('flow_config', metavar='text', default='config_flow.yml')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text", default="", help="suffix for output files")
    args = parser.parse_args()

    make_yaml(args.flow_config,
                args.outputdir,
                args.suffix)