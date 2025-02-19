import argparse
import yaml
import pandas as pd
import numpy as np
import uproot
import ROOT
import ctypes
from ROOT import TFile, TKDE, TCanvas, TH1D, TF1

def templ_producer_kde(tree_file, var, pt_min, pt_max, query, name, outfile='', tree_name='O2hfcanddplite'):
    
    print(f"Producing KDE from {tree_file} for var {var}, {pt_min} <= pt < {pt_max}, name {name}")
    print(f"-----> Query: {pt_min} <= fPt < {pt_max} and {query}")
    # convert the tree_file to a pandas dataframe
    dfsData = []
    with uproot.open(f'{tree_file}') as f:
        for key in f.keys():
            if tree_name in key:
                dfData = f[key].arrays(library='pd')
                dfsData.append(dfData)      
    full_dataset = pd.concat([df for df in dfsData], ignore_index=True)
    full_dataset.query(f"{pt_min} <= fPt < {pt_max} and {query}", inplace=True)
    var_values = full_dataset[f'{var}'].tolist()

    # kde = TKDE(len(var_values), np.asarray(var_values, 'd'), 0, 3)
    # kde_func = TF1("kde_name",kde,-10,10,0)
    # print(type(kde_func))
    # Assuming TKDE class is defined as such, and var_values is your data
    kde = TKDE(len(var_values), np.asarray(var_values, 'd'), 0, 3)

    # Use a lambda function to wrap the TKDE's Eval method (or similar)
    kde_func = kde.GetFunction(5000)
    # kde_func = TF1("kde_name", lambda x, par: kde.Eval(x[0]), -10, 10)

    # Print the type of the TF1 object
    # print(type(kde_func))
    # print(f'kde_func.Eval(2): {kde_func.Eval(2.)}')
    # file = TFile(f'/home/mdicosta/FlowDplus/FinalResults/templs_from_histo/cutvar_withtempls/ry/debug_{pt_min}_{pt_max}_{name}.root', 'recreate')
    # print('ciao1')
    # kde_func.Write()
    # print('ciao2')
    # file.Close()
    # kde_func = kde.GetFunction(1000)
    
    binned_var_values = TH1D(f'hBinned', f'hBinned', 3000, 0, 3)
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
        cOverlap = TCanvas('cOverlap', 'cOverlap', 600, 600)
        cOverlap.cd()
        binned_var_values.Draw()
        kde_func.Draw('same')
        outfile.mkdir(f'KDE_pT_{pt_min}_{pt_max}_{name}')
        outfile.cd(f'KDE_pT_{pt_min}_{pt_max}_{name}')
        kde.Write('kde')
        binned_var_values.Write()
        kde_func.Write()
        cOverlap.Write()
    
    return kde, kde_func, binned_var_values

def templ_producer_histo(tree_file, var, pt_min, pt_max, queries, names, relweights=[], outfile='', tree_name='O2hfcanddplite'):

    print(f"Producing KDE from {tree_file} for var {var}, {pt_min} <= pt < {pt_max}, names {names}")
    # convert the tree_file to a pandas dataframe
    dfsData = []
    print(f"tree_file: {tree_file}")
    with uproot.open(f'{tree_file}') as f:
        for key in f.keys():
            if tree_name in key:
                dfData = f[key].arrays(library='pd')
                dfsData.append(dfData)      
    df = pd.concat([df for df in dfsData], ignore_index=True)
    histos_templ = []
    for query, name in zip(queries, names):
        print(f"query: {query}")
        print(f"{pt_min} < fPt < {pt_max} and {query}")
        templ_df = df.query(f"{pt_min} < fPt < {pt_max} and {query}")[var].to_numpy()
        histos_templ.append(ROOT.TH1D(
            f"hist_templ_{name}_pt{pt_min:.1f}_{pt_max:.1f}",
            "#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27))
        for var_value in templ_df:
            histos_templ[-1].Fill(var_value)

    histo_comb = ROOT.TH1D(
        f"hist_templ_combined_pt{pt_min:.1f}_{pt_max:.1f}",
        "#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27)

    if relweights != []:
        for irelweight, histo_templ in zip(relweights, histos_templ):
            histo_comb.Add(histo_templ, irelweight)
    else:
        for irelweight, histo_templ in zip(relweights, histos_templ):
            histo_comb.Add(histo_templ, 1)

    histo_comb_smoothened = histo_comb.Clone(f"{histo_comb.GetName()}_smooth")
    histo_comb_smoothened.Smooth(100)

    if outfile != '':
        outfile.mkdir(f'hTempl_pT_{pt_min}_{pt_max}')
        outfile.cd(f'hTempl_pT_{pt_min}_{pt_max}')
        for hist in histos_templ:
            hist.Write()
        histo_comb.Write()
        histo_comb_smoothened.Write()

    return histo_comb

def get_templates_weights(tree_file, pt_min, pt_max, sgn_weight, templ_weights, names, outfile = ''):
    """
    """

    dfsData = []
    print(f"tree_file: {tree_file}")
    with uproot.open(f'{tree_file}') as f:
        for key in f.keys():
            if tree_name in key:
                dfData = f[key].arrays(library='pd')
                dfsData.append(dfData)      
    df = pd.concat([df for df in dfsData], ignore_index=True)

    df_bkg = df.query("abs(fFlagMcMatchRec) == 4")
    df_signal = df.query("abs(fFlagMcMatchRec) == 1")

    hist_frac_bkg_to_signal = ROOT.TH1D("hist_frac_bkg_to_signal",
                                        ";#it{p}_{T} (GeV/#it{c});bkg corr / signal",
                                        len(pt_bins)-1, np.array(pt_bins, dtype=np.float64))
    
    hist_templ_comb = ROOT.TH1D(f"hist_templ_comb_pt{pt_min:.1f}_{pt_max:.1f}",
                                ";#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27)
    
    hist_templs, df_templs = [], []
    for name, weight in zip(names, templ_weights):
        hist_templs.append(ROOT.TH1D(f"hist_templ_{name}_pt{pt_min:.1f}_{pt_max:.1f}",
                                     ";#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27))
        df_templs.append(df_bkg.query(f"{pt_min} < fPt < {pt_max} and abs(fFlagMcDecayChanRec) > 2"))
        for mass in df_templs[-1]["fM"].to_numpy():
            hist_templs[-1].Fill(mass)
        hist_templ_comb.Add(hist_templs[-1], weight)
            
    hist_signal = ROOT.TH1D(f"hist_signal_pt{pt_min:.1f}_{pt_max:.1f}",
                            ";#it{M}(K#pi#pi) (GeV/#it{c})", 600, 1.67, 2.27)
    df_pt_signal = df_signal.query(f"{pt_min} < fPt < {pt_max}")
    for mass in df_pt_signal["fM"].to_numpy():
        hist_signal[ipt].Fill(mass)
    hist_signal.Scale(sgn_weight)

    # Ds/D+ is underestimated in pythia CRMode2
    hist_frac_bkg_to_signal.SetBinContent(ipt+1,
                                          hist_templ_comb.Integral() / hist_signal.Integral())

    if outfile != ''
        hist_frac_bkg_to_signal.Write()
        hist_signal.Write()
        for hist in hist_templs:
            hist.Write()
            hist_smooth = hist.Clone(f"{hist.GetName()}_smooth")
            hist_smooth.Smooth(100)
            hist_smooth.Write()

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
    
    
        # for iPt, (bkgStr, sgnStr, bkgVnStr) in enumerate(zip(BkgFuncStr, SgnFuncStr, BkgFuncVnStr)):
    #     if 'IncludeTempls' in fitConfig and fitConfig['IncludeTempls'][iPt]:
    #         for iFlag, flag in enumerate(fitConfig['TemplsFlags']):
    #             if fitConfig['TemplInputType'][iFlag] == 'Grid':
    #                 KDEtemplates[iPt][iFlag], kde_func, hRebinnedHistos[iPt][iFlag] = kde_producer(fitConfig['TemplsInputs'][iFlag], 'fM', ptMins[iPt], ptMaxs[iPt], flag,
    #                                                                     templatesFile, fitConfig['TemplsTreeNames'][iFlag])
    #             elif fitConfig['TemplInputType'][iFlag] == 'File':
    #                 print('GETTING TEMPLATE FROM FILE')
    #                 KDEtemplates[iPt][iFlag], kde_func, hRebinnedHistos[iPt][iFlag] = kde_producer(fitConfig['TemplsInputs'][iFlag], 'massPiKaPi', 
    #                                                                                                ptMins[iPt], ptMaxs[iPt], 633,
    #                                                                                                templatesFile, 'treeD', "sumPdgDaus == 211+211+211")
    #             else:
    #                 print(f'ERROR: incorrect setting for including KDEs in fit! Exit!')
    #                 sys.exit()

    #             templatesFile.cd()
    #             hRebinnedHistos[iPt][iFlag].Write(f"hBinned_pt_{iPt}")
    #             cTemplOverlap[iPt][iFlag] = TCanvas(f'cOverlap_{iPt}_{flag}', f'cOverlap_{iPt}_{flag}', 600, 600)
    #             cTemplOverlap[iPt][iFlag].cd()
    #             hRebinnedHistos[iPt][iFlag].Draw()
    #             kde_func.Draw('same')
    #             cTemplOverlap[iPt][iFlag].Write()
    #     else:
    #         KDEtemplates[iPt], kde_func, hRebinnedHistos[iPt], cTemplOverlap[iPt] = None, None, None, None
    #         continue
        
    # templatesFile.Close()
    # KDEtemplatesFuncts = []
    # print(f"KDEtemplates: {KDEtemplates}")
    # if KDEtemplates is not None:
    #     for KDEtemplatesPt in KDEtemplates:
    #         if KDEtemplatesPt is not None:
    #             KDEtemplatesFuncts.append([KDE.GetFunction() for KDE in KDEtemplatesPt])
    #         else:
    #             KDEtemplatesFuncts.append(None)