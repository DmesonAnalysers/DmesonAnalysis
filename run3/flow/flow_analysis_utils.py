'''
Analysis utilities for flow analysis
'''
import ROOT
import sys
import numpy as np

def get_vn_versus_mass(thnSparse, inv_mass_bins, mass_axis, vn_axis, debug=False):
    '''
    Project vn versus mass

    Input:
        - thnSparse:
            THnSparse, input THnSparse obeject (already projected in centrality and pt)
        - inv_mass_bins:
            list of floats, bin edges for the mass axis
        - mass_axis:
            int, axis number for mass
        - vn_axis:
            int, axis number for vn
        - debug:
            bool, if True, create a debug file with the projections (default: False)

    Output:
        - hist_mass_proj:
            TH1D, histogram with vn as a function of mass
    '''
    hist_vn_proj = thnSparse.Projection(vn_axis, mass_axis)
    hist_mass_proj = thnSparse.Projection(mass_axis)
    hist_mass_proj.Reset()
    invmass_bins = np.array(inv_mass_bins)
    hist_mass_proj = ROOT.TH1D('hist_mass_proj', 'hist_mass_proj', len(invmass_bins)-1, invmass_bins)
    
    if debug:
        outfile = ROOT.TFile('debug.root', 'RECREATE')
    
    for i in range(hist_mass_proj.GetNbinsX()):
        bin_low = hist_vn_proj.GetXaxis().FindBin(invmass_bins[i])
        bin_high = hist_vn_proj.GetXaxis().FindBin(invmass_bins[i+1])
        profile = hist_vn_proj.ProfileY(f'profile_{bin_low}_{bin_high}', bin_low, bin_high)
        mean_sp = profile.GetMean()
        mean_sp_err = profile.GetMeanError()
        hist_mass_proj.SetBinContent(i+1, mean_sp)
        hist_mass_proj.SetBinError(i+1, mean_sp_err)
    
    if debug:
        hist_vn_proj.Write()
        hist_mass_proj.Write()
        outfile.Close()
    
    return hist_mass_proj

def get_resolution(resolution_file_name, dets, centMin, centMax, doEP=False):
    '''
    Compute resolution for SP or EP method

    Input:
        - resolution_file_name: 
            str, path to the resolution file
        - dets: 
            list of strings, subsystems to compute the resolution
            if 2 subsystems are given, the resolution is computed as the product of the two
            if 3 subsystems are given, the resolution is computed as the product of the first two divided by the third
        - centMin:
            int, minimum centrality bin
        - centMax:
            int, maximum centrality bin
        - doEP:
            bool, if True, compute EP resolution
            if False, compute SP resolution (default)

    Output:
        - histo_reso:
            TH1D, histogram with the resolution value as a function of centrality
        - histo_dets:
            list of TH1D, list of projections used to compute the resolution
        - histo_means:
            list of TH1D, list of histograms with the mean value of the projections as a function of centrality
    '''
    infile = ROOT.TFile(resolution_file_name, 'READ')
    histo_projs, histo_dets, histo_means = [], [], []
    detA = dets[0]
    detB = dets[1]
    if len(dets) == 3:
        detC = dets[2]
    dets = [f'{detA}{detB}', f'{detA}{detC}', f'{detB}{detC}'] if len(dets) == 3 else [f'{detA}{detB}']
    
    # set path and prefix
    if doEP:
        path = 'hf-task-flow-charm-hadrons/epReso/'
        prefix = 'EpReso'
    else:
        path = 'hf-task-flow-charm-hadrons/spReso/'
        prefix = 'SpReso'
   
    # collect the qvecs and the prepare histo for mean and resolution
    for det in dets:
        histo_dets.append(infile.Get(f'{path}h{prefix}{det}'))
        histo_dets[-1].SetDirectory(0)
        histo_dets[-1].SetName(f'h{prefix}{det}')
        histo_means.append(histo_dets[-1].ProjectionX(f'proj_{histo_dets[-1].GetName()}_mean'))
        histo_means[-1].SetDirectory(0)
        histo_means[-1].Reset()
        histo_projs.append([])

        # collect projections
        for cent in range(centMin, centMax):
            bin_cent_low = histo_dets[-1].GetXaxis().FindBin(cent) # common binning
            bin_cent_high = histo_dets[-1].GetXaxis().FindBin(cent)
            histo_projs[-1].append(histo_dets[-1].ProjectionY(f'proj_{histo_dets[-1].GetName()}_{cent}_{cent}', bin_cent_low, bin_cent_high))
            histo_projs[-1][-1].SetDirectory(0)

        # Appllying absolute value to the projections
        for ihist, _ in enumerate(histo_projs[-1]): histo_means[-1].SetBinContent(ihist+1, histo_projs[-1][ihist].GetMean())
    infile.Close()

    histo_reso = ROOT.TH1F('', '', 100, 0, 100)
    histo_reso.SetDirectory(0)
    for icent in range(centMin, centMax):
        histo_reso.SetBinContent(icent+1, compute_resolution([histo_means[i].GetBinContent(icent+1) for i in range(len(dets))]))

    return histo_reso, histo_dets, histo_means


def compute_resolution(subMean):
    '''
    Compute resolution for SP or EP method

    Input:
        - subMean:
            list of floats, list of mean values of the projections

    Output:
        - resolution:
            float, resolution value
    '''
    if len(subMean) == 1:
        resolution =  subMean[0]
        if resolution <= 0:
            return 0
        else:
            return np.sqrt(resolution)
    elif len(subMean) == 3:
        resolution = (subMean[2] * subMean[1]) / subMean[0] if subMean[0] != 0 else 0
        if resolution <= 0:
            return 0
        else:
            return np.sqrt(resolution)
    else:
        print('ERROR: dets must be a list of 2 or 3 subsystems')
        sys.exit(1)