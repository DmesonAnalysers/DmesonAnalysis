import argparse
import yaml
import pandas as pd
import numpy as np
import uproot
import ROOT
import ctypes
from ROOT import TFile, TKDE, TCanvas, TH1D

def kde_producer(tree_file, var, pt_min, pt_max, flag, outfile='', tree_name='O2hfcanddplite'):
    
    print(f"Producing KDE from {tree_file} for var {var}, {pt_min} <= pt < {pt_max}, flag {flag}")
    # convert the tree_file to a pandas dataframe
    dfsData = []
    with uproot.open(f'{tree_file}') as f:
        for key in f.keys():
            if tree_name in key:
                dfData = f[key].arrays(library='pd')
                dfsData.append(dfData)      
    full_dataset = pd.concat([df for df in dfsData], ignore_index=True)
    pt_filtered_df = full_dataset.query(f"{pt_min} <= fPt < {pt_max}")
    filtered_df = pt_filtered_df.query(f"fFlagMcDecayChanRec == {flag} or fFlagMcDecayChanRec == {-flag}")
    var_values = filtered_df[f'{var}'].tolist()  # Or use `.tolist()` to get a list

    kde = TKDE(len(var_values), np.asarray(var_values, 'd'), min(var_values), max(var_values))
    kde_func = kde.GetFunction(1000)
    
    binned_var_values = TH1D(f'hBinned', f'hBinned', 100, min(var_values), max(var_values))
    for var_value in var_values:
        binned_var_values.Fill(var_value)
    
    max_content = 0
    for bin_idx in range(1, binned_var_values.GetNbinsX() + 1):
        bin_content = binned_var_values.GetBinContent(bin_idx)
        if bin_content > max_content:
            max_content = bin_content
            max_bin = bin_idx
    binned_var_values.Scale(kde_func.GetMaximum() / binned_var_values.GetBinContent(max_bin))
    
    if outfile != '':
        outfile.mkdir(f'KDE_pT_{pt_min}_{pt_max}_flag{flag}')
        outfile.cd(f'KDE_pT_{pt_min}_{pt_max}_flag{flag}')
        kde_func.Write()
        binned_var_values.Write()
    
    return kde, kde_func, binned_var_values

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--config", "-cfg", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("--var", "-v", metavar="text",
                        default="fM", help="variable of interest")
    parser.add_argument("--ptmin", "-pmin", metavar="text",
                        default="2.", help="min pt")
    parser.add_argument("--ptmax", "-pmax", metavar="text",
                        default="4.", help="max pt")
    parser.add_argument("--flag", "-f", metavar="chn flag",
                        default="2", help="channel flag")
    parser.add_argument("--input", "-in", metavar="path/input.root",
                        default="AnalysisResults.root", help="path to file containing tree")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    
    args = parser.parse_args()
    
    KDEs = []
    binned_histos = []
    if args.config != parser.get_default("config"):
        with open(args.config, 'r') as ymlCfgFile:
            config = yaml.load(ymlCfgFile, yaml.FullLoader)
        
    output_dir = config["outputdir"] if args.outputdir == parser.get_default("outputdir") else args.outputdir
    suffix = config["suffix"] if args.suffix == parser.get_default("suffix") else args.suffix
    outfile = ROOT.TFile(f'{output_dir}/kde_{suffix}.root', 'RECREATE')
        
    if args.config != parser.get_default("config"):
        for pt_low, pt_max in zip(config['pt_mins'], config['pt_maxs']):
            KDE, histo = kde_producer(config['input'], config['variable'], pt_low, 
                                           pt_max, config['chn_flag'], outfile)
    else:
        KDE, histo = kde_producer(args.input, args.var, args.ptmin,
                                       args.ptmax, args.flag, outfile)
    outfile.Close()    
