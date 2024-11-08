import argparse
import yaml
import pandas as pd
import numpy as np
import uproot
import ROOT
import ctypes
from ROOT import TFile, TKDE, TCanvas, TH1D

def kde_producer(tree, var, cent, pt_min, pt_max, inv_mass_bins, flag, outfile):
    
    # convert the tree to a pandas dataframe
    dfsData = []
    # print(tree)
    with uproot.open(f'{tree}') as f:
        # print(f.keys())
        for iKey, key in enumerate(f.keys()):
            if 'O2hfcanddplite' in key:
                # print(key)
                dfData = f[key].arrays(library='pd')
                dfsData.append(dfData)        

    full_dataset = pd.concat([df for df in dfsData], ignore_index=True)
    print(full_dataset['fFlagMcMatchRec'].tolist())
    print(full_dataset['fOriginMcRec'].tolist())
    print(full_dataset['fFlagMcDecayChanRec'].tolist())
    pt_filtered_df = full_dataset.query(f"{pt_min} <= fPt < {pt_max}")
    filtered_df = pt_filtered_df.query(f"fFlagMcMatchRec == {flag} or fFlagMcMatchRec == {-flag}")
    # filtered_df = pt_filtered_df.query(f"fFlagMcMatchRec == 1")
    var_values = filtered_df[f'{var}'].tolist()  # Or use `.tolist()` to get a list
    # print(var_values)
    
    kde = TKDE(len(var_values), np.asarray(var_values, 'd'), 1.7, 2.0)
    kde_func = kde.GetFunction(1000)
    
    binned_var_values = TH1D(f'hBinned_{var}_{flag}', f'hBinned_{var}_{flag}', 5000, 1, 3)
    for var_value in var_values:
        binned_var_values.Fill(var_value)
    
    cOverlap = TCanvas('cOverlap', 'cOverlap', 600, 600)
    cOverlap.cd()
    binned_var_values.Draw()
    kde_func.Draw('SAME')
    
    outfile.mkdir(f'pT_{pt_low}_{pt_max}')
    outfile.cd(f'pT_{pt_low}_{pt_max}')
    kde_func.Write()
    binned_var_values.Write()
    cOverlap.Write()
    
    return kde_func, binned_var_values 


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("--var", "-v", metavar="text",
                        default="fM", help="variable of interest")
    parser.add_argument("--config", "-cfg", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k2050", help="centrality class")
    parser.add_argument("--ptmin", "-pmin", metavar="text",
                        default="2.", help="min pt")
    parser.add_argument("--ptmax", "-pmax", metavar="text",
                        default="4.", help="max pt")
    parser.add_argument("--invmassbins", "-imbins", metavar="text",
                        default="[1, 1.5, 2]", help="inv mass bins")
    parser.add_argument("--flag", "-f", metavar="chn flag",
                        default="2", help="channel flag")
    parser.add_argument("--input", "-in", metavar="path/input.root",
                        default="AnalysisResults.root", help="path to tree")
    parser.add_argument("--wagon_id", "-w", metavar="text",
                        default="", help="wagon ID", required=False)
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    
    args = parser.parse_args()
    print(args)
    
    KDEs = []
    binned_histos = []
    if args.config != parser.get_default("config"):
        with open(args.config, 'r') as ymlCfgFile:
            config = yaml.load(ymlCfgFile, yaml.FullLoader)
        
    output_dir = config["outputdir"] if args.outputdir == parser.get_default("outputdir") else args.outputdir
    suffix = config["suffix"] if args.suffix == parser.get_default("suffix") else args.suffix
    outfile = ROOT.TFile(f'{output_dir}/kde_{suffix}.root', 'RECREATE')
        
    if args.config != parser.get_default("config"):
        for pt_low, pt_max, inv_mass_bins in zip(config['pt_mins'], config['pt_maxs'], config['inv_mass_bins']):
            KDE, histo = kde_producer(config['input'], config['variable'], config['centrality'], \
                                      pt_low, pt_max, inv_mass_bins, config['chn_flag'],
                                      outfile)
    else:
        KDE, histo = kde_producer(args.input, args.var, args.centrality, \
                                  args.ptmin, args.ptmax, args.invmassbins, args.flag,
                                  outfile)
    outfile.Close()    
    