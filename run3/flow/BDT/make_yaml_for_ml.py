'''
Script to create a yaml file with a set of cuts for ML
python3 make_yaml_for_ml.py config_flow.yml -o path/to/output -s text

'''
import yaml
import argparse
import os
import numpy as np
import sys
sys.path.append('..')
from flow_analysis_utils import get_cut_sets_config

def make_combination(mass_axis, pt_axis, bkg_axis, sig_axis, ptmins, ptmaxs, nCutSets,
                     sig_cut_lower_file, sig_cut_upper_file, bkg_cut_lower_file, bkg_cut_upper_file):
    '''
    Create a dictionary with the combination of cuts for each cutset

    Parameters:
    mass_axis: int
        axis number for the invariant mass
    pt_axis: int
        axis number for the pt
    bkg_axis: int
        axis number for the bkg output
    sig_axis: int
        axis number for the signal output
    ptmins: list
        list of lower pt edges
    ptmaxs: list
        list of upper pt edges
    sig_cut_file: dict
        dictionary with the signal cut for each cutset
    bkg_cut_file: dict
        dictionary with the background cut for each cutset

    Returns:
    combinations: dict
        dictionary with the combination of cuts for each cutset
    '''
    combinations = {}
    for iFile in range(nCutSets):
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
                    'min': [float(i) for i in bkg_cut_lower_file[iFile]],
                    'max': [float(j) for j in bkg_cut_upper_file[iFile]],
                    'name': 'ML_output_Bkg'
                },
                'ML_output_FD': {
                    'axisnum': sig_axis,
                    'min': [float(i) for i in sig_cut_lower_file[iFile]],
                    'max': [float(j) for j in sig_cut_upper_file[iFile]],
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

    ## safely check
    if len(ptmins) != len(ptmaxs):
        raise ValueError(f'''The number of pt bins({len(ptmins)}, {len(ptmaxs)} are not the same''')

    nCutSets, sig_cut_lower, sig_cut_upper, bkg_cut_lower, bkg_cut_upper = get_cut_sets_config(flow_config)
    
    maxCutSets = max(nCutSets)
    sig_cut_lower_file, sig_cut_upper_file, bkg_cut_lower_file, bkg_cut_upper_file = {}, {}, {}, {}
    for iCut in range(maxCutSets):
        sig_cut_lower_file[iCut], sig_cut_upper_file[iCut], bkg_cut_lower_file[iCut], bkg_cut_upper_file[iCut] = [], [], [], []
        for iPt in range(len(ptmins)):
            # consider the different number of cutsets for each pt bin
            if iCut < nCutSets[iPt]:
                sig_cut_lower_file[iCut].append(sig_cut_lower[iPt][iCut])
                sig_cut_upper_file[iCut].append(sig_cut_upper[iPt][iCut])
                bkg_cut_lower_file[iCut].append(bkg_cut_lower[iPt][iCut])
                bkg_cut_upper_file[iCut].append(bkg_cut_upper[iPt][iCut])
            else:
                sig_cut_lower_file[iCut].append(sig_cut_lower[iPt][nCutSets[iPt]-1])
                sig_cut_upper_file[iCut].append(sig_cut_upper[iPt][nCutSets[iPt]-1])
                bkg_cut_lower_file[iCut].append(bkg_cut_lower[iPt][nCutSets[iPt]-1])
                bkg_cut_upper_file[iCut].append(bkg_cut_upper[iPt][nCutSets[iPt]-1])

    combinations = make_combination(mass_axis, pt_axis, bkg_axis, sig_axis, ptmins, ptmaxs, maxCutSets,
                                    sig_cut_lower_file, sig_cut_upper_file, bkg_cut_lower_file, bkg_cut_upper_file)

    print(f'''
The axis number for the invariant mass is {mass_axis}
The axis number for the pt is {pt_axis}
The axis number for the bkg output is {bkg_axis}
The axis number for the signal output is {sig_axis}
''')
    for iFile in range(maxCutSets):
        print(f'''For cutset {iFile}:
        ptmin: {ptmins}
        ptmax: {ptmaxs}
        sig cut: {["{:.3f}".format(x) for x in sig_cut_lower_file[iFile]]}
                 {["{:.3f}".format(x) for x in sig_cut_upper_file[iFile]]}
        bkg cut: {["{:.3f}".format(x) for x in bkg_cut_lower_file[iFile]]}
                 {["{:.3f}".format(x) for x in bkg_cut_upper_file[iFile]]}
''')

    os.makedirs(f'{outputdir}/config', exist_ok=True)
    for iFile in range(maxCutSets):
        with open(f'{outputdir}/config/cutset_{suffix}_{iFile:02}.yml', 'w') as file:
            yaml.dump(combinations[iFile], file, default_flow_style=False)
    print(f'Yaml files are saved in {outputdir}/config')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('flow_config', metavar='text', default='config_flow.yml')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text", default="", help="suffix for output files")
    args = parser.parse_args()

    make_yaml(args.flow_config,
                args.outputdir,
                args.suffix)