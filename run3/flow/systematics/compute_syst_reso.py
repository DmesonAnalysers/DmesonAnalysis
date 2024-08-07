import sys
import os
import numpy as np
import argparse
import yaml
from ROOT import TFile, TCanvas, TH1F, kBlack, kSpring, kOrange, kAzure, kGray, kRed, TLegend, TGraphAsymmErrors
sys.path.append('../')
from flow_analysis_utils import get_centrality_bins


def compute_syst_reso(config_modify, config_def, rypathsyst,  cent, reso_file, outputdir):

    #______________________________________________________________________________________
    # Read configuration
    with open(config_modify, 'r') as ymlCfgFile:
        config_mod = yaml.load(ymlCfgFile, yaml.FullLoader)
    ptmins = config_mod['ptmins']
    ptmaxs = config_mod['ptmaxs']

    with open(config_def, 'r') as ymlCfgFile:
        config_def = yaml.load(ymlCfgFile, yaml.FullLoader)
    det_A = config_def['detA']
    det_B = config_def['detB']
    det_C = config_def['detC']

    #______________________________________________________________________________________
    # Collect all .root files from rypathsyst
    hraw_yields = {}
    hresolutions = {} # only needed to collect cenrtalities
    for _, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs)):
        pt_label = f"ptmin{ptmin}_ptmax{ptmax}"
        hraw_yields[pt_label] = []
        hresolutions[pt_label] = []
        for file in os.listdir(rypathsyst):
            if file.endswith('.root') and f'ptmin{ptmin}_ptmax{ptmax}' in file:
                print(f'INFO: Found {file}')
                tfile = TFile().Open(os.path.join(rypathsyst, file))
                # check if the file has the right keys
                if not tfile.GetListOfKeys().Contains('hRawYieldsSimFit'):
                    print(f'ERROR: {file} does not contain hRawYieldsSimFit')
                    continue
                # get last 5 characters of the file name
                file_name = file.replace('.root', '')
                file_suffix = file_name[-5:]
                hraw_yields[pt_label].append(tfile.Get('hRawYieldsSimFit'))
                hresolutions[pt_label].append(tfile.Get('hist_reso'))
                hraw_yields[pt_label][-1].SetDirectory(0)
                hresolutions[pt_label][-1].SetDirectory(0)
                hraw_yields[pt_label][-1].SetName(f'hRawYieldsSimFit_{ptmin}_{ptmax}_{file_suffix}')
                tfile.Close()
    if not hraw_yields:
        sys.exit(f'ERROR: No .root files found in {rypathsyst}')

    #______________________________________________________________________________________
    # Get raw yields vs centrality
    hist_ry_vs_cent = {}
    _, centMinMax = get_centrality_bins(cent)
    for ip, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs)):
        pt_label = f"ptmin{ptmin}_ptmax{ptmax}"
        hraw_yields_pt = hraw_yields[pt_label]
        hreso_pt = hresolutions[pt_label]
        hist_ry_vs_cent[pt_label] = TH1F(f'hRawYieldsSimFit_vs_cent_{ptmin}_{ptmax}', '',
                                         centMinMax[1] - centMinMax[0],
                                         centMinMax[0], centMinMax[1])
        for _, (hry, hreso) in enumerate(zip(hraw_yields_pt, hreso_pt)):
            # get centrality bin center
            bincentreso = hreso.GetBinCenter(1)
            # get the bin index in the raw yield histogram
            binindexry = hist_ry_vs_cent[pt_label].FindBin(bincentreso)
            hist_ry_vs_cent[pt_label].SetBinContent(binindexry, hry.GetBinContent(1))
            hist_ry_vs_cent[pt_label].SetBinError(binindexry, hry.GetBinError(1))
        hist_ry_vs_cent[pt_label].SetDirectory(0)

    #______________________________________________________________________________________
    # Plot the raw yields vs centrality
    cols = [kOrange+1, kAzure+4, kRed+1, kGray+1]
    canv = TCanvas('canv', 'canv', 800, 800)
    canv.cd()
    canv.SetLeftMargin(0.15)
    legw = TLegend(0.6, 0.6, 0.9, 0.9)
    legw.SetBorderSize(0)
    legw.SetFillStyle(0)
    legw.SetTextSize(0.03)
    legw.SetTextFont(42)
    for ip, hist in enumerate(hist_ry_vs_cent.values()):
        hist.SetLineColor(cols[ip])
        hist.SetMarkerColor(cols[ip])
        hist.SetMarkerStyle(20)
        hist.SetMarkerSize(1)
        hist.SetLineWidth(2)
        hist.SetStats(0)
        if ip != len(hist_ry_vs_cent) - 1:
            legw.AddEntry(hist, f'{ptmins[ip]} < #it{{p}}_{{T}} < {ptmaxs[ip]} GeV/#it{{c}}', 'lp')
        if ip == 0:
            hist.GetYaxis().SetRangeUser(hist.GetMinimum()*0.3, hist.GetMaximum()*1.1)
            hist.Draw()
        hist.GetXaxis().SetTitle('centrality (%)')
        hist.GetYaxis().SetTitle('raw yield')
        hist.Draw('same')
    legw.Draw()
    canv.Update()

    #______________________________________________________________________________________
    # Load resolution file
    reso_file = TFile().Open(reso_file)
    # check if the file has the right keys
    if not reso_file.GetListOfKeys().Contains(f'{det_A}_{det_B}_{det_C}'):
        sys.exit(f'ERROR: {reso_file} does not contain {det_A}_{det_B}_{det_C}')
    hist_reso_def_deltacent = reso_file.Get(f'{det_A}_{det_B}_{det_C}/histo_reso_delta_cent') # only needed to collect reference resolution
    hist_reso_def_deltacent.SetDirectory(0)
    hist_reso_def = reso_file.Get(f'{det_A}_{det_B}_{det_C}/histo_reso')
    hist_reso_def.SetDirectory(0)
    reso_file.Close()

    #______________________________________________________________________________________
    # Compute systematics
    hreso_weighted_mean = {}
    hreso_weighted_mean_err = {}
    sum_weight = {}
    hreso_arithm_mean = 0
    counter = 0
    for icent in range(1, hist_reso_def.GetNbinsX() + 1):
        reso = hist_reso_def.GetBinContent(icent)
        if reso == 0:
            continue
        counter += 1
        bincent = hist_reso_def.GetBinCenter(icent)

        # Arithmetic mean
        hreso_arithm_mean += reso

        # Loop over pt bins
        for _, (ptmin, ptmax) in enumerate(zip(ptmins[:-1], ptmaxs[:-1])):
            pt_label = f"ptmin{ptmin}_ptmax{ptmax}"

            if icent == 2:
                sum_weight[pt_label] = 0
                hreso_weighted_mean[pt_label] = 0
                hreso_weighted_mean_err[pt_label] = 0

            binindexry = hist_ry_vs_cent[pt_label].FindBin(bincent)
            weigth = hist_ry_vs_cent[pt_label].GetBinContent(binindexry)
            sum_weight[pt_label] += weigth

            hreso_weighted_mean[pt_label] += reso * weigth
            hreso_weighted_mean_err[pt_label] += (reso * weigth)**2

    for _, (ptmin, ptmax) in enumerate(zip(ptmins[:-1], ptmaxs[:-1])):
        pt_label = f"ptmin{ptmin}_ptmax{ptmax}"
        hreso_weighted_mean[pt_label] /= sum_weight[pt_label]
    hreso_arithm_mean /= counter

    canvsyst = TCanvas('canvsyst', 'canvsyst', 800, 800)
    canvsyst.cd()
    canvsyst.SetLeftMargin(0.2)
    canvsyst.SetRightMargin(0.05)
    leg = TLegend(0.2, 0.6, 0.9, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetTextFont(42)
    greso_def_deltacent = TGraphAsymmErrors()
    greso_def_deltacent.SetPoint(0, (centMinMax[1]+centMinMax[0]) / 2,
                                 1)
    greso_def_deltacent.SetPointError(0, (centMinMax[1]-centMinMax[0]) / 2,
                                         (centMinMax[1]-centMinMax[0]) / 2,
                                         0,0)
    greso_def_deltacent.SetLineColor(kBlack)
    greso_def_deltacent.SetMarkerColor(kBlack)
    greso_def_deltacent.SetMarkerStyle(20)
    greso_def_deltacent.SetMarkerSize(1)
    greso_def_deltacent.SetLineWidth(2)
    greso_def_deltacent.SetTitle('')
    greso_def_deltacent.GetXaxis().SetTitle('centrality (%)')
    greso_def_deltacent.GetYaxis().SetTitle('#it{R}_{2}{SP} / #it{R}_{2}{SP}^{ref.}')
    greso_def_deltacent.GetYaxis().SetTitleOffset(1.9)
    greso_def_deltacent.GetYaxis().SetDecimals()
    greso_def_deltacent.GetYaxis().SetRangeUser(0.95, 1.08)
    greso_def_deltacent.Draw('')
    leg.AddEntry(greso_def_deltacent, 'Reference resolution', 'lp')
    print(f'INFO: Reference resolution: {hist_reso_def_deltacent.GetBinContent(1)}')

    greso_ari = greso_def_deltacent.Clone('hist_reso_ari')
    greso_ari.SetPoint(0, (centMinMax[1]+centMinMax[0]) / 2,
                       hreso_arithm_mean / hist_reso_def_deltacent.GetBinContent(1))
    greso_ari.SetLineColor(kSpring+2)
    greso_ari.SetMarkerColor(kSpring+2)
    greso_ari.SetMarkerStyle(20)
    greso_ari.SetMarkerSize(1)
    greso_ari.SetLineWidth(2)
    greso_ari.Draw('same p')
    leg.AddEntry(greso_ari, 'Arithmetic mean resolution', 'lp')
    gist_reso_weight = {}
    print(f'INFO: Arithmetic mean resolution: {hreso_arithm_mean}')

    for ipt, (ptmin, ptmax) in enumerate(zip(ptmins[:-1], ptmaxs[:-1])):
        pt_label = f"ptmin{ptmin}_ptmax{ptmax}"
        gist_reso_weight[pt_label] = greso_def_deltacent.Clone('hist_reso_weight_' + pt_label)
        gist_reso_weight[pt_label].SetPoint(0, (centMinMax[1]+centMinMax[0]) / 2,
                                            hreso_weighted_mean[pt_label] / hist_reso_def_deltacent.GetBinContent(1))
        gist_reso_weight[pt_label].SetLineColor(cols[ipt])
        gist_reso_weight[pt_label].SetMarkerColor(cols[ipt])
        gist_reso_weight[pt_label].SetMarkerStyle(20)
        gist_reso_weight[pt_label].SetMarkerSize(1)
        gist_reso_weight[pt_label].SetLineWidth(2)
        gist_reso_weight[pt_label].Draw('same p')
        leg.AddEntry(gist_reso_weight[pt_label], f'Weighted mean resolution {ptmin} < #it{{p}}_{{T}} < {ptmax} GeV/#it{{c}}', 'lp')
        print(f'INFO: Weighted mean resolution for {pt_label}: {hreso_weighted_mean[pt_label]}')
    leg.Draw()
    canvsyst.Update()
    input()

    #______________________________________________________________________________________
    # Save output
    outputfile = os.path.join(outputdir, f'syst_reso_{cent}.root')
    outfile = TFile(outputfile, 'RECREATE')
    hist_reso_def_deltacent.Write()
    greso_ari.Write()
    for hist in hist_ry_vs_cent.values():
        hist.Write()
    for hist in gist_reso_weight.values():
        hist.Write()
    canv.Write()
    canvsyst.Write()
    outfile.Close()
    print(f'INFO: Output saved in {outputfile}')

    canv.SaveAs(f'{outputdir}/SystReso_centweights.pdf')
    canvsyst.SaveAs(f'{outputdir}/SystReso_vs_cent.pdf')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('config_modifications', metavar='text', default='config_modifications.yml')
    parser.add_argument('config_default', metavar='text', default='config_default.yml')
    parser.add_argument('rypathsyst', metavar='text',
                        default='path to the directory containing the .root files from multitrial')
    parser.add_argument('--centClass', "-c", metavar="text", default=".", help="centrality class")
    parser.add_argument('reso_file', metavar='text', default='resolution .root file')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    args = parser.parse_args()

    compute_syst_reso(args.config_modifications,
                      args.config_default,
                      args.rypathsyst,
                      args.centClass,
                      args.reso_file,
                      args.outputdir)
