import yaml
import os
import numpy as np
import argparse
import itertools

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
    bkg_vn, invm_bins, mass_min, mass_max, bkg_func, rebins = ([] for i in range(6))
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
    

    #____________________________________________________________________________________
    # Create a new yaml file for each combination of the modifications
    for bkg in bkg_vn:
        for bkgf in bkg_func:
            for reb in rebins:
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
                            for j in range(len(new_data['inv_mass_bins'])):
                                new_data['inv_mass_bins'][j] = [l for l in invm_bin]
                            for i in range(len(new_data['MassMin'])):
                                new_data['MassMin'][i] = mmin
                            for i in range(len(new_data['MassMax'])):
                                new_data['MassMax'][i] = mmax
                            outfile_name = os.path.join(output_dir_folder, f'{infile_name}_bkgvn{bkg}_bkg{bkgf}_{reb}_invm{step}_min{mmin}_max{mmax}.yml')

                            with open(outfile_name, 'w') as f:
                                yaml.dump(new_data, f, default_flow_style=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('input_config', metavar='text', default='config_Ds_Fit.yml')
    parser.add_argument('--modifications_config', "-m", metavar='text', default='')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    args = parser.parse_args()

    modify_yaml(args.input_config,
                args.outputdir,
                args.modifications_config)
