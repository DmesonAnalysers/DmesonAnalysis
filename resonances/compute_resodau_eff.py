'''
python script to calculate the efficiency of the resonances
run: python ComputeResoDauEff.py config_proj_reso.yml
'''
import sys
import os
import argparse
import yaml
import numpy as np
from alive_progress import alive_bar
from particle import Particle
from ROOT import TFile, TH3F, TH2F, TCanvas, TLatex, gROOT
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle

def sparse_loader(infile, dirname, listname, sparsetype, verbose=False):
    '''
    Load the THnSparse from the task file

    Parameters:
    infile (str): name of the input file
    dirname (str): name of the directory
    listname (str): name of the list
    sparsetype (str): name of the THnSparse

    Returns:
    sparse (THnSparse): THnSparse
    '''
    if verbose:
        print(f'Loading {listname} from {dirname} in {infile}')
    infile = TFile(infile)
    indir = infile.Get(dirname)
    if not indir:
        print(f'Cannot find directory {dirname} in {infile}')
        return None
    inlist = indir.Get(listname)
    if not inlist:
        print(f'Cannot find list {listname} in {dirname} in {infile}')
        return None
    sparse = inlist.FindObject(sparsetype)
    if not sparse:
        print(f'Cannot find {sparsetype} in {listname} in {dirname} in {infile}')
        return None

    return sparse

def decode_cuts_to_axis(selections, verbose=False):
    '''
    Decode the cuts to the axis of the THnSparse

    Parameters:
    selections (dict): dictionary with the cuts
    verbose (bool): print the cuts (default: False)

    Returns:
    axisnum (int): number of the axis to apply the cut
    cutvar_min (float or list): minimum value of the cut variable or list of minimum values
    cutvar_max (float or list): maximum value of the cut variable or list of maximum values
    pt_mins (list): list of minimum values of the pt, if empty the cut is not pt dependent (default: [])
    pt_maxs (list): list of maximum values of the pt, if empty the cut is not pt dependent (default: [])
    '''
    selsToaxis = {}
    for key in selections:
        selsToaxis[key] = {}
        selsToaxis[key]['axisnum'] = None
        selsToaxis[key]['cutvar_min'] = None
        selsToaxis[key]['cutvar_max'] = None
        selsToaxis[key]['pt_mins'] = []
        selsToaxis[key]['pt_maxs'] = []
        if key == 'pt_min':
            selsToaxis[key]['axisnum'] = 0
            selsToaxis[key]['cutvar_min'] = selections[key]
            selsToaxis[key]['cutvar_max'] = None
        elif 'delta_mass' in key:
            # TODO: implement the delta mass cut with nsigma
            continue
        elif key == 'mass_min':
            selsToaxis[key]['axisnum'] = 3
            selsToaxis[key]['cutvar_min'] = selections['mass_min']
        elif key == 'mass_max':
            selsToaxis[key]['axisnum'] = 3
            selsToaxis[key]['cutvar_max'] = selections['mass_max']
        elif key == 'BDT':
            if selections[key]['BDT_bkg']:
                selsToaxis['BDT_bkg'] = {}
                selsToaxis['BDT_bkg']['axisnum'] = 4
                selsToaxis['BDT_bkg']['cutvar_max'] = selections[key]['BDT_bkg']
                selsToaxis['BDT_bkg']['cutvar_min'] = [0. for i in range(len(selections[key]['pt_mins']))]
                selsToaxis['BDT_bkg']['pt_mins'] = selections[key]['pt_mins']
                selsToaxis['BDT_bkg']['pt_maxs'] = selections[key]['pt_maxs']
            if selections[key]['BDT_prompt']:
                selsToaxis['BDT_prompt'] = {}
                selsToaxis['BDT_prompt']['axisnum'] = 5
                selsToaxis['BDT_prompt']['cutvar_min'] = selections[key]['BDT_prompt']
                selsToaxis['BDT_prompt']['cutvar_max'] = [1. for i in range(len(selections[key]['pt_mins']))]
                selsToaxis['BDT_prompt']['pt_mins'] = selections[key]['pt_mins']
                selsToaxis['BDT_prompt']['pt_maxs'] = selections[key]['pt_maxs']
        elif key == 'cosp_min':
            selsToaxis[key]['axisnum'] = 4
            selsToaxis[key]['cutvar_min'] = selections[key]
            selsToaxis[key]['cutvar_max'] = None
        elif key == 'declen_xy_min':
            selsToaxis[key]['axisnum'] = 5
            selsToaxis[key]['cutvar_min'] = selections[key]
            selsToaxis[key]['cutvar_max'] = None
        elif key == 'dca_dau_min':
            selsToaxis[key]['axisnum'] = 6
            selsToaxis[key]['cutvar_min'] = selections[key]
            selsToaxis[key]['cutvar_max'] = None
        else:
            print(f'Unrecognized cut on {key}')

    key_to_remove = []
    for key in selsToaxis:
        if selsToaxis[key]['axisnum'] is None:
            key_to_remove.append(key)
    for key in key_to_remove:
        selsToaxis.pop(key)

    if verbose:
        print(selsToaxis)
    return selsToaxis

def apply_analysis_cuts_to_sparse(sparse, sels_decoded, pt_min=None, pt_max=None, verbose=False):
    '''
    Apply the analysis cuts to the THnSparse

    Parameters:
    sparse (THnSparse): THnSparse
    sels_decoded (dict): dictionary with the cuts
    pt_min (float): minimum pt of the pt-dependent cuts
    pt_max (float): maximum pt of the pt-dependent cuts
    verbose (bool): print the cuts applied (default: False)

    Returns:
    sparse (THnSparse): THnSparse with the analysis cuts applied
    '''
    for key in sels_decoded.keys():
        axisnum = sels_decoded[key]['axisnum']
        cut_var_min = sels_decoded[key]['cutvar_min']
        cut_var_max = sels_decoded[key]['cutvar_max']
        pt_mins = sels_decoded[key]['pt_mins']
        pt_maxs = sels_decoded[key]['pt_maxs']

        if axisnum is None:
            print(f'No analysis cuts to apply to {sparse.GetName()}')
            return sparse

        if (pt_mins == [] or pt_maxs == []):
            if cut_var_min:
                bin_min = sparse.GetAxis(axisnum).FindBin(cut_var_min * 1.0001)
            else:
                bin_min = sparse.GetAxis(axisnum).FindBin(sparse.GetAxis(axisnum).GetXmin())
            if cut_var_max:
                bin_max = sparse.GetAxis(axisnum).FindBin(cut_var_max * 0.9999)
            else:
                bin_max = sparse.GetAxis(axisnum).FindBin(sparse.GetAxis(axisnum).GetXmax())
            if verbose:
                print(f'Applying cuts to {sparse.GetName()} on axis {axisnum} = [{cut_var_min}, {cut_var_max}] --> [{bin_min}, {bin_max}]')
            sparse.GetAxis(axisnum).SetRange(bin_min, bin_max)
        else:
            if verbose:
                print(f'Applying pt dependent analysis cuts to {sparse.GetName()} on axis {axisnum}')
            # TODO: extend to broeader pt bins
            for ipt, (PtMin, PtMax) in enumerate(zip(pt_mins, pt_maxs)):
                if (pt_min >= PtMin and pt_max <= PtMax):
                    lower_sel = cut_var_min[ipt]
                    higher_sel = cut_var_max[ipt]
                    break
            bin_min = sparse.GetAxis(axisnum).FindBin(lower_sel * 1.0001)
            bin_max = sparse.GetAxis(axisnum).FindBin(higher_sel * 0.9999)
            bin_pt_min = sparse.GetAxis(0).FindBin(pt_min * 1.0001)
            bin_pt_max = sparse.GetAxis(0).FindBin(pt_max * 0.9999)
            sparse.GetAxis(0).SetRange(bin_pt_min, bin_pt_max)
            sparse.GetAxis(axisnum).SetRange(bin_min, bin_max)
            sparse.GetAxis(axisnum).SetRange(bin_min, bin_max)

    return sparse

def compute_efficiencies(config, trigger, pdg_d, pdg_v0, output_dir, suffix, mult_weights):
    '''
    function for efficiency computation

    Parameters
    ----------

    - config (str): path of config file
    - trigger (str): trigger class (MB or HM)
    - pdg_d (int): PDG code of the D meson
    - pdg_v0 (int): PDG code of the V0
    - output_dir (str): output directory
    - suffix (str): suffix to append to output files
    - mult_weights (str): multitplicity weights (all, cand, candinmass)
    '''

    SetGlobalStyle(padtopmargin=0.05, padrightmargin=0.15, padleftmargin=0.15,
                   padbottommargin=0.15, palette=55, labelsize=0.03, titlesize=0.03,
                   labeloffset=0.01, titleoffsety=1.8, titleoffsetx=1.8, titleoffsetz=1.,
                   opttitle=0, optstat=0)
    gROOT.SetBatch(True)

    with open(config, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    input_files = inputCfg['inputs']['mc'][trigger]
    origins = ['Prompt', 'NonPrompt']

    if trigger not in ['MB', 'HM']:
        raise ValueError(f'Wrong trigger: {trigger}. Allowed values at the moment: MB, HM')

    # define titles and labels
    if pdg_d == 411:
        name_d = "Dplus"
        name_decay = "DplustoKpipi"
        label_d = "D^{+}"
    elif pdg_d == 413:
        name_d = "Dstar"
        name_decay = "DstartoD0pi"
        label_d = "D*^{+}"
    else:
        print(f"ERROR: D pdg code {pdg_d} not supported")

    if pdg_v0 == 310:
        name_v0 = "K0S"
        label_v0 = "K_{S}^{0}"
    elif pdg_v0 == 3122:
        name_v0 = "Lambda"
        label_v0 = "#Lambda"
    else:
        print(f"ERROR: V0 pdg code {pdg_v0} not supported")

    pdg_reso = -1
    if pdg_d == 411:
        if pdg_v0 == 310:
            pdg_reso = 435
        else:
            print(f"ERROR: combination of D and V0 pdg codes {pdg_d}-{pdg_v0} not supported")
    elif pdg_d == 413:
        if pdg_v0 == 310:
            pdg_reso = 10433
        else:
            print(f"ERROR: combination of D and V0 pdg codes {pdg_d}-{pdg_v0} not supported")
    name_reso = Particle.from_pdgid(pdg_reso).name

    print('Configuration file loaded successfully.')
    print('Resonance to be analysed: ', name_reso)
    print('Input files: ', input_files)
    print('Trigger: ', trigger)
    print('_______________________________________________________')

    print(f'\033[1mComputing efficiency for:  {name_reso} --> {name_d} + {name_v0}\033[0m')

    #__________________________________________________________________________
    # Loading the THnSparse
    sparseDmeson, sparseV0 = {}, {}
    print('\033[1mLoading the THnSparse...\033[0m\r', end='')
    for file in input_files:
        mult_weights_suffix = ""
        if mult_weights != "":
            mult_weights_suffix = f"_multweights_{mult_weights}"
        indir = f"PWGHF_D2H_HFResoBuilder{name_decay}{name_d}{name_v0}_{trigger}{mult_weights_suffix}"
        inlist = f"coutput{name_decay}{name_d}{name_v0}_{trigger}{mult_weights_suffix}"
        sparseDmeson, sparseV0 = {}, {}

        for status in ['Reco', 'Gen']:
            sparseDmeson[status] = {}
            sparseV0[status] = {}
            for origin in origins:
                if status == 'Reco':
                    sparseDmeson[status][origin] = sparse_loader(file, indir, inlist, f'h{origin}DmesonMCReco_{pdg_d}')
                else:
                    sparseDmeson[status][origin] = sparse_loader(file, indir, inlist, f'h{origin}DmesonMCGen_{pdg_d}')

            sparseV0[status] = {}
            if status == 'Reco':
                sparseV0[status] = sparse_loader(file, indir, inlist, f'hV0MCReco_{pdg_v0}')
            else:
                sparseV0[status] = sparse_loader(file, indir, inlist, f'hV0MCGen_{pdg_v0}')

    print('\033[1mLoading the THnSparse: complete.\033[0m')

    #__________________________________________________________________________
    # Computing the efficiency
    outfile_name = os.path.join(
        output_dir, f'effmap_{name_d}_{name_v0}_{trigger}{mult_weights_suffix}{suffix}.root')
    outfile = TFile(outfile_name, 'recreate')
    outfile.cd()
    pt_bins_d = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24] # FIXME: hardcoded for now
    npt_bins_d = len(pt_bins_d) - 1
    pt_bins_d = np.array(pt_bins_d, dtype=np.float32)
    pt_bins_v0 = [0, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24] # FIXME: hardcoded for now
    npt_bins_v0 = len(pt_bins_v0) - 1
    pt_bins_v0 = np.array(pt_bins_v0, dtype=np.float32)
    y_bins = [-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0] # FIXME: hardcoded for now
    ny_bins = len(y_bins) - 1
    y_bins = np.array(y_bins, dtype=np.float32)
    phi_bins = [np.linspace(0., 2*np.pi, 10)] # hardcoded for now
    phi_bins = [round(phi, 2) for phi in np.concatenate(phi_bins)]
    nphi_bins = len(phi_bins) - 1
    phi_bins = np.array(phi_bins, dtype=np.float32)
    selections = inputCfg['selections']

    # mesons
    sels = selections[pdg_d][trigger]
    sels_decoded = decode_cuts_to_axis(sels, verbose=False)

    h3Eff_mes, h3Num, h3Den, hEff_pt_phi_mes, hEff_pt_y_mes, hEff_phi_y_mes =  ([] for i in range(6))
    canvas_mes = TCanvas('cResoEffMeson', '', 2000, 900)
    canvas_mes.Divide(2)

    for iorigin, (origin) in enumerate(origins):
        canvas_mes.cd(iorigin+1)
        canvas_mes.cd(iorigin+1).SetTitle(f'{name_d}_{origin}')
        h3Eff_mes.append(TH3F(f'h3Eff_mes{name_d}{origin}',
                          ';#it{p}_{T} (GeV/#it{c});#it{y}; #phi',
                          npt_bins_d, pt_bins_d, ny_bins, y_bins,
                          nphi_bins, phi_bins))
        h3Num.append(TH3F(f'h3Num_mes{name_d}{origin}',
                          ';#it{p}_{T} (GeV/#it{c});#it{y}; #phi',
                          npt_bins_d, pt_bins_d, ny_bins, y_bins,
                          nphi_bins, phi_bins))
        h3Den.append(TH3F(f'h3Den_mes{name_d}{origin}',
                          ';#it{p}_{T} (GeV/#it{c});#it{y}; #phi',
                          npt_bins_d, pt_bins_d, ny_bins, y_bins,
                          nphi_bins, phi_bins))
        h3Eff_mes[-1].GetXaxis().CenterTitle()
        h3Eff_mes[-1].GetYaxis().CenterTitle()
        h3Eff_mes[-1].GetZaxis().CenterTitle()
        hEff_pt_phi_mes.append(TH2F(f'hEff{name_d}{origin}_pt_phi',
                               f'{name_d} {origin} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #phi',
                               npt_bins_d, pt_bins_d, nphi_bins, phi_bins))
        hEff_pt_y_mes.append(TH2F(f'hEff{name_d}{origin}_pt_y',
                             f'{name_d} {origin} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #it{{y}}',
                             npt_bins_d, pt_bins_d, ny_bins, y_bins))
        hEff_phi_y_mes.append(TH2F(f'hEff{name_d}{origin}_phi_y',
                              f'{name_d} {origin} efficiency; #phi; #it{{y}}',
                              nphi_bins, phi_bins, ny_bins, y_bins))
        h3Eff_mes[-1].SetDirectory(0)
        h3Eff_mes[-1].SetTitle(f'{name_d}_{origin}_{status}')

        with alive_bar(len(pt_bins_d)-1, bar='smooth', spinner='classic',
                       title=f'Filling {name_d} {origin} efficiency') as bar:

            for ipt, (pt_min, pt_max) in enumerate(zip(pt_bins_d[:-1], pt_bins_d[1:])):
                lower_pt_thr = sels_decoded['pt_min']['cutvar_min']
                if pt_min < lower_pt_thr:
                    continue
                sparseDmeson['Reco'][origin] = apply_analysis_cuts_to_sparse(
                    sparseDmeson['Reco'][origin],
                    sels_decoded,
                    pt_min=pt_min,
                    pt_max=pt_max,
                    verbose=False
                )
                pt_bin_min = sparseDmeson['Gen'][origin].GetAxis(0).FindBin(pt_min*1.0001)
                pt_bin_max = sparseDmeson['Gen'][origin].GetAxis(0).FindBin(pt_max*0.9999)
                sparseDmeson['Gen'][origin].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                sparseDmeson['Reco'][origin].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                for iphi, (phi_min, phi_max) in enumerate(zip(phi_bins[:-1], phi_bins[1:])):
                    phi_bin_min = sparseDmeson['Gen'][origin].GetAxis(2).FindBin(phi_min*1.0001)
                    phi_bin_max = sparseDmeson['Gen'][origin].GetAxis(2).FindBin(phi_max*0.9999)
                    sparseDmeson['Gen'][origin].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                    sparseDmeson['Reco'][origin].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                    for iy, (y_min, y_max) in enumerate(zip(y_bins[:-1], y_bins[1:])):
                        y_bin_min = sparseDmeson['Gen'][origin].GetAxis(1).FindBin(y_min*1.0001)
                        y_bin_max = sparseDmeson['Gen'][origin].GetAxis(1).FindBin(y_max*0.9999)
                        sparseDmeson['Gen'][origin].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                        sparseDmeson['Reco'][origin].GetAxis(1).SetRange(y_bin_min, y_bin_max)

                        hNum = sparseDmeson['Reco'][origin].Projection(0, 1, 2)
                        hDen = sparseDmeson['Gen'][origin].Projection(0, 1, 2)
                        hNum.SetDirectory(0)
                        hDen.SetDirectory(0)
                        num = hNum.Integral()
                        den = hDen.Integral()
                        eff = num/den if den > 0 else 0
                        h3Num[-1].SetBinContent(ipt+1, iy+1, iphi+1, num)
                        h3Den[-1].SetBinContent(ipt+1, iy+1, iphi+1, den)
                        hEff_pt_phi_mes[-1].SetBinContent(ipt+1, iphi+1, eff)
                        hEff_pt_y_mes[-1].SetBinContent(ipt+1, iy+1, eff)
                        hEff_phi_y_mes[-1].SetBinContent(iphi+1, iy+1, eff)
                        del hNum, hDen
                bar()
        h3Eff_mes[-1] = h3Num[-1].Clone(f'h3Eff{name_d}{origin}')
        h3Eff_mes[-1].Divide(h3Eff_mes[-1], h3Den[-1], 1, 1, 'B')
        hy_pt_y_num = h3Num[-1].Project3D('yxe')
        hy_pt_y_den = h3Den[-1].Project3D('yxe')
        hy_pt_y_eff = hy_pt_y_num.Clone(f'hy_pt_y_eff{name_d}{origin}')
        hy_pt_y_eff.Divide(hy_pt_y_eff, hy_pt_y_den, 1, 1, 'B')
        h_pt_num = hy_pt_y_num.ProjectionX()
        h_pt_den = hy_pt_y_den.ProjectionX()
        h_pt_eff = h_pt_num.Clone(f'h_pt_eff{name_d}{origin}')
        h_pt_eff.Divide(h_pt_eff, h_pt_den, 1, 1, 'B')
        h3Eff_mes[-1].Draw('BOX2 colz')
        latLabel = TLatex()
        latLabel.SetNDC()
        latLabel.SetTextSize(0.05)
        latLabel.DrawLatex(0.2, 0.95, f'{origin} {label_d}')
        outfile.mkdir(f'{name_d}_{origin}')
        outfile.cd(f'{name_d}_{origin}')
        h3Den[-1].Write()
        h3Num[-1].Write()
        h3Eff_mes[-1].Write()
        hy_pt_y_eff.Write(f'h3Dy_pt_y_eff_{name_d}_{origin}')
        h_pt_eff.Write(f'h_pt_eff_{name_d}_{origin}')
        hEff_pt_phi_mes[-1].Write()
        hEff_pt_y_mes[-1].Write()
        hEff_phi_y_mes[-1].Write()
        canvas_mes.Update()

    outfile.cd()
    canvas_mes.Write()

    # v0s
    sels = selections[pdg_v0][trigger]
    sels_decoded = decode_cuts_to_axis(sels, verbose=False)

    h3Eff_v0, h3Num_v0, h3Den_v0, hEff_pt_phi_v0, hEff_pt_y_v0, hEff_phi_y_v0 = ([] for _ in range(6))
    canvas_v0 = TCanvas('cResoEffv0', '', 1600, 900)
    outfile.mkdir(f'{name_v0}_{name_d}')
    h3Eff_v0.append(TH3F(f'h3Eff_{name_d}_{name_v0}',
                    f'{name_v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                    npt_bins_v0, pt_bins_v0, ny_bins, y_bins,
                    nphi_bins, phi_bins))
    h3Num_v0.append(TH3F(f'h3Num_{name_d}_{name_v0}',
                         f'{name_v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                         npt_bins_v0, pt_bins_v0, ny_bins, y_bins,
                         nphi_bins, phi_bins))
    h3Den_v0.append(TH3F(f'h3Dum_{name_d}_{name_v0}',
                         f'{name_v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                         npt_bins_v0, pt_bins_v0, ny_bins, y_bins,
                         nphi_bins, phi_bins))
    hEff_pt_phi_v0.append(TH2F(f'hEff_pt_phi_{name_d}_{name_v0}',
                          f'{name_v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #phi',
                          npt_bins_v0, pt_bins_v0, nphi_bins, phi_bins))
    hEff_pt_y_v0.append(TH2F(f'hEff_pt_y_{name_d}_{name_v0}',
                             f'{name_v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #it{{y}}',
                             npt_bins_v0, pt_bins_v0, ny_bins, y_bins))
    hEff_phi_y_v0.append(TH2F(f'hEff_phi_y_{name_d}_{name_v0}',
                              f'{name_v0} efficiency;#phi; #it{{y}}',
                              nphi_bins, phi_bins, ny_bins, y_bins))
    outfile.cd(f'{name_v0}_{name_d}')
    with alive_bar(len(pt_bins_v0)-1, bar='smooth',
                   spinner='classic', title=f'Filling {name_v0} efficiency') as bar:
        for ipt, (pt_min, pt_max) in enumerate(zip(pt_bins_v0[:-1], pt_bins_v0[1:])):
            lower_pt_thr = sels_decoded['pt_min']['cutvar_min']
            if pt_min < lower_pt_thr:
                print(f'WARNING: {name_v0} efficiency for pt < {lower_pt_thr} GeV/c is not computed')
                continue
            sparseV0['Reco'] = apply_analysis_cuts_to_sparse(
                sparseV0['Reco'],
                sels_decoded,
                pt_min=pt_min,
                pt_max=pt_max,
                verbose=False
            )
            pt_bin_min = sparseV0['Reco'].GetAxis(0).FindBin(pt_min*1.0001)
            pt_bin_max = sparseV0['Reco'].GetAxis(0).FindBin(pt_max*0.9999)
            sparseV0['Reco'].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
            sparseV0['Gen'].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
            for iphi, (phi_min, phi_max) in enumerate(zip(phi_bins[:-1], phi_bins[1:])):
                phi_bin_min = sparseV0['Reco'].GetAxis(2).FindBin(phi_min*1.0001)
                phi_bin_max = sparseV0['Reco'].GetAxis(2).FindBin(phi_max*0.9999)
                sparseV0['Reco'].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                sparseV0['Gen'].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                for iy, (y_min, y_max) in enumerate(zip(y_bins[:-1], y_bins[1:])):
                    y_bin_min = sparseV0['Reco'].GetAxis(1).FindBin(y_min*1.0001)
                    y_bin_max = sparseV0['Reco'].GetAxis(1).FindBin(y_max*0.9999)
                    sparseV0['Reco'].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                    sparseV0['Gen'].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                    hNum = sparseV0['Reco'].Projection(0, 1, 2)
                    hNum.SetDirectory(0)
                    hDen = sparseV0['Gen'].Projection(0, 1, 2)
                    hDen.SetDirectory(0)
                    num = hNum.Integral()
                    den = hDen.Integral()
                    eff = num/den if den > 0 else 0
                    h3Den_v0[-1].SetBinContent(ipt+1, iy+1, iphi+1, den)
                    h3Num_v0[-1].SetBinContent(ipt+1, iy+1, iphi+1, num)
                    hEff_pt_phi_v0[-1].SetBinContent(ipt+1, iphi+1, eff)
                    hEff_pt_y_v0[-1].SetBinContent(ipt+1, iy+1, eff)
                    hEff_phi_y_v0[-1].SetBinContent(iphi+1, iy+1, eff)
                    del hNum, hDen
            bar()
        h3Den_v0[-1].Write(f'hDen_{name_d}_{name_v0}')
        h3Num_v0[-1].Write(f'hNum_{name_d}_{name_v0}')
        h3Eff_v0[-1] = h3Num_v0[-1].Clone(f'hEff_{name_d}_{name_v0}')
        h3Eff_v0[-1].Divide(h3Num_v0[-1], h3Den_v0[-1], 1, 1, 'B')
        h_pt_y_num = h3Num_v0[-1].Project3D('yxe')
        h_pt_y_den = h3Den_v0[-1].Project3D('yxe')
        h_pt_y_eff = h_pt_y_num.Clone(f'hEff_pt_y_{name_d}_{name_v0}')
        h_pt_y_eff.Divide(h_pt_y_eff, h_pt_y_den, 1, 1, 'B')
        h_pt_y_eff.Write(f'hEff_pt_y_{name_d}_{name_v0}')
        h_pt_num = h3Num_v0[-1].ProjectionX()
        h_pt_den = h3Den_v0[-1].ProjectionX()
        h_pt_eff = h_pt_num.Clone(f'hEff_pt_{name_d}_{name_v0}')
        h_pt_eff.Divide(h_pt_eff, h_pt_den, 1, 1, 'B')
        h_pt_eff.Write(f'hEff_pt_{name_d}_{name_v0}')
        h3Eff_v0[-1].Write()
        hEff_pt_phi_v0[-1].Write()
        hEff_pt_y_v0[-1].Write()
        hEff_phi_y_v0[-1].Write()
        outfile.cd('..')
    canvas_v0.cd()
    h3Eff_mes[-1].Draw('BOX2 colz same')
    latlabel = TLatex()
    latlabel.SetNDC()
    latlabel.SetTextSize(0.04)
    latlabel.DrawLatex(0.15, 0.85, f'{label_v0}')
    canvas_v0.Update()

    outfile.cd()
    canvas_v0.Write()

    outfile.Close()
    print(f'Output written to {outfile_name}')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config_proj.yml", help="input config file")
    parser.add_argument("--trigger", "-t", metavar="text",
                        default="HM", required=True, help="trigger class (options: [HM, MB]")
    parser.add_argument("--pdg_D", "-d", type=int, required=True, default=413,
                        help="pdg code of the D meson (options: [411, 413])")
    parser.add_argument("--pdg_V0", "-v0", type=int, required=True, default=310,
                        help="pdg code of the V0 (options: [310, 3122])")
    parser.add_argument("--output_dir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--mult_weights", "-m", metavar="text",
                        default="", help="multiplicity weights (all, cand, candinmass)")
    args = parser.parse_args()

    compute_efficiencies(
        args.config,
        args.trigger,
        args.pdg_D,
        args.pdg_V0,
        args.output_dir,
        args.suffix,
        args.mult_weights
    )
