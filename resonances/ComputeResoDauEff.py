# python script to calculate the efficiency of the resonances
# run: python ComputeResoDauEff.py config_proj_reso.yml
import sys
import argparse
import yaml
import numpy as np
from alive_progress import alive_bar
from particle import Particle
from ROOT import TFile, TH3F, TH2F, TCanvas, kRed
sys.path.insert(0, '..')
#from utils.TaskFileLoader import LoadSingleSparseFromTask, LoadNormObjFromTask, LoadSingleSparseFromTask
from utils.StyleFormatter import SetGlobalStyle, LatLabel

SetGlobalStyle(padtopmargin=0.05, padrightmargin=0.05, padleftmargin=0.15,
               padbottommargin=0.15, palette=55)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='config_proj_reso.yml',
                    help='config file name with root input files')
args = parser.parse_args()

def main():
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    input_files = inputCfg['eff']['input_files']
    input_dirs = inputCfg['eff']['input_dirs']
    input_list = [dir[dir.find('Builder'):].replace('Builder', '') for dir in input_dirs]
    mult = inputCfg['eff']['mult']
    mesons = inputCfg['eff']['mesons']
    v0s = inputCfg['eff']['v0s']
    mesons_code = {'Dplus': 411, 'Dstar': 413}
    v0s_code = {'K0S': 310}
    for mes in mesons:
        if mes not in mesons_code.keys():
            raise ValueError(f'Wrong meson name: {mes}.\
                             Allowed values at the moment: {mesons_code.keys()}')
    for v0 in v0s:
        if v0 not in ['K0S']:
            raise ValueError(f'Wrong V0 name: {v0}.\
                             Allowed values at the moment: {v0s_code.keys()}')
    if mult not in ['MB', 'HM']:
        raise ValueError(f'Wrong multiplicity: {mult}.\
                         Allowed values at the moment: MB, HM')
    mes_label = {'Dplus': 'D^{+}', 'Dstar': 'D^{*+}'}
    v0s_label = {'K0S': 'K^{0}_{S}'}
    origins = ['Prompt', 'NonPrompt']
    for i, (dir, list) in enumerate(zip(input_dirs, input_list)):
        input_dirs[i] = f'{dir}_{mult}'
        input_list[i] = f'coutput{list}_{mult}'

    #__________________________________________________________________________
    # Loading the THnSparse
    sparseDmeson, sparseV0 = {}, {}
    print('\033[1mLoading the THnSparse...\033[0m\r', end='')
    for file in input_files: # loop over the input files
        for i, (dir, list) in enumerate(zip(input_dirs, input_list)):
            for meson in mesons:
                if meson not in dir:
                    continue
                meson_code = mesons_code[meson]
                sparseDmeson[meson] = {}
                sparseV0[meson] = {}
                for status in ['Reco', 'Gen']:
                    sparseDmeson[meson][status] = {}
                    sparseV0[meson][status] = {}
                    for origin in origins:
                        if status == 'Reco':
                            sparseDmeson[meson][status][origin] = SparseLoader(file, dir, list,
                                                                               f'h{origin}DmesonMCReco_{meson_code}')
                        else:
                            sparseDmeson[meson][status][origin] = SparseLoader(file, dir, list,
                                                                               f'h{origin}DmesonMCGen_{meson_code}')
                    for v0 in v0s_code.keys():
                        v0_code = v0s_code[v0]
                        sparseV0[meson][status][v0] = {}
                        if status == 'Reco':
                            sparseV0[meson][status][v0] = SparseLoader(file, dir, list, f'hV0MCReco_{v0_code}')
                        else:
                            sparseV0[meson][status][v0] = SparseLoader(file, dir, list, f'hV0MCGen_{v0_code}')
    print('\033[1mLoading the THnSparse: complete.\033[0m')

    #__________________________________________________________________________
    # Applying the analysis cuts to the THnSparse
    selections = inputCfg['selections']

    # meson
    print('\033[1mApplying the analysis selections to the THnSparse...\033[0m\r', end='')
    for meson in ['Dplus', 'Dstar']:
        meson_code = mesons_code[meson]
        sels = selections[meson_code][mult]
        sels_decoded = DecodeCutsToAxis(sels)
        for origin in origins:
            for key in sels_decoded.keys():
                axisnum = sels_decoded[key]['axisnum']
                cut_var_min = sels_decoded[key]['cutvar_min']
                cut_var_max = sels_decoded[key]['cutvar_max']
                pt_mins = sels_decoded[key]['pt_mins']
                pt_maxs = sels_decoded[key]['pt_maxs']
                sparseDmeson[meson]['Reco'][origin] = ApplyAnalysesCutsToSparse(sparseDmeson[meson]['Reco'][origin],
                                                                                axisnum,
                                                                                cut_var_min,
                                                                                cut_var_max,
                                                                                pt_mins,
                                                                                pt_maxs)
    # v0s
    for v0 in v0s:
        v0_code = v0s_code[v0]
        sels = selections[v0_code][mult]
        sels_decoded = DecodeCutsToAxis(sels)
        for key in sels_decoded.keys():
            axisnum = sels_decoded[key]['axisnum']
            cut_var_min = sels_decoded[key]['cutvar_min']
            cut_var_max = sels_decoded[key]['cutvar_max']
            pt_mins = sels_decoded[key]['pt_mins']
            pt_maxs = sels_decoded[key]['pt_maxs']
            for meson in ['Dplus', 'Dstar']:
                for v0 in v0s_code.keys():
                    sparseV0[meson]['Reco'][v0] = ApplyAnalysesCutsToSparse(sparseV0[meson]['Reco'][v0],
                                                                            axisnum,
                                                                            cut_var_min,
                                                                            cut_var_max,
                                                                            pt_mins,
                                                                            pt_maxs)
    print('\033[1mApplying the analysis selections to the THnSparse: complete.\033[0m')

    #__________________________________________________________________________
    # Computing the efficiency
    outfile = TFile(inputCfg['eff']['outfile'], 'recreate')
    outfile.cd()
    pt_bins = [0, 1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 24, 36] # hardcoded for now
    npt_bins = len(pt_bins) - 1
    pt_bins = np.array(pt_bins, dtype=np.float32)
    y_bins = [-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0] # hardcoded for now
    ny_bins = len(y_bins) - 1
    y_bins = np.array(y_bins, dtype=np.float32)
    phi_bins = [np.linspace(0., 2*np.pi, 10)] # hardcoded for now
    phi_bins = [round(phi, 2) for phi in np.concatenate(phi_bins)]
    nphi_bins = len(phi_bins) - 1
    phi_bins = np.array(phi_bins, dtype=np.float32)

    # mesons
    h3Eff_mes, hEff_pt_phi_mes, hEff_pt_y_mes, hEff_phi_y_mes = [], [], [], []
    canvas_mes = TCanvas('cResoEffMeson', '', 1600, 900)
    canvas_mes.Divide(len(mesons), 2)
    counter = 0
    for _, (meson) in enumerate(mesons):
        for origin in origins:
            counter += 1
            canvas_mes.cd(counter)
            canvas_mes.cd(counter).SetTitle(f'{meson}_{origin}')
            h3Eff_mes.append(TH3F(f'h3Eff_mes{meson}{origin}',
                              f'{meson} {origin} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                              npt_bins, pt_bins, ny_bins, y_bins,
                              nphi_bins, phi_bins))
            hEff_pt_phi_mes.append(TH2F(f'hEff{meson}{origin}_pt_phi',
                                    f'{meson} {origin} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #phi',
                                    npt_bins, pt_bins, nphi_bins, phi_bins))
            hEff_pt_y_mes.append(TH2F(f'hEff{meson}{origin}_pt_y',
                                  f'{meson} {origin} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #it{{y}}',
                                  npt_bins, pt_bins, ny_bins, y_bins))
            hEff_phi_y_mes.append(TH2F(f'hEff{meson}{origin}_phi_y',
                                   f'{meson} {origin} efficiency; #phi; #it{{y}}',
                                   nphi_bins, phi_bins, ny_bins, y_bins))
            h3Eff_mes[-1].SetDirectory(0)
            h3Eff_mes[-1].SetTitle(f'{meson}_{origin}_{status}')
            with alive_bar(len(pt_bins)-1, bar='smooth', spinner='classic',
                           title=f'Filling {meson} {origin} efficiency') as bar:
                for ipt, (pt_min, pt_max) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
                    pt_bin_min = sparseDmeson[meson]['Gen'][origin].GetAxis(0).FindBin(pt_min*
                                                                                       1.0001)
                    pt_bin_max = sparseDmeson[meson]['Gen'][origin].GetAxis(0).FindBin(pt_max*
                                                                                       0.9999)
                    for iphi, (phi_min, phi_max) in enumerate(zip(phi_bins[:-1], phi_bins[1:])):
                        phi_bin_min = sparseDmeson[meson]['Gen'][origin].GetAxis(2).FindBin(phi_min*1.0001)
                        phi_bin_max = sparseDmeson[meson]['Gen'][origin].GetAxis(2).FindBin(phi_max*0.9999)
                        for iy, (y_min, y_max) in enumerate(zip(y_bins[:-1], y_bins[1:])):
                            y_bin_min = sparseDmeson[meson]['Gen'][origin].GetAxis(1).FindBin(y_min*1.0001)
                            y_bin_max = sparseDmeson[meson]['Gen'][origin].GetAxis(1).FindBin(y_max*0.9999)
                            sparseDmeson[meson]['Gen'][origin].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                            sparseDmeson[meson]['Gen'][origin].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                            sparseDmeson[meson]['Gen'][origin].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                            sparseDmeson[meson]['Reco'][origin].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                            sparseDmeson[meson]['Reco'][origin].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                            sparseDmeson[meson]['Reco'][origin].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                            hNum = sparseDmeson[meson]['Reco'][origin].Projection(0, 1, 2)
                            hDen = sparseDmeson[meson]['Gen'][origin].Projection(0, 1, 2)
                            hNum.SetDirectory(0)
                            hDen.SetDirectory(0)
                            num = hNum.Integral()
                            den = hDen.Integral()
                            eff = num/den if den > 0 else 0
                            h3Eff_mes[-1].SetBinContent(ipt+1, iy+1, iphi+1, eff)
                            hEff_pt_phi_mes[-1].SetBinContent(ipt+1, iphi+1, eff)
                            hEff_pt_y_mes[-1].SetBinContent(ipt+1, iy+1, eff)
                            hEff_phi_y_mes[-1].SetBinContent(iphi+1, iy+1, eff)
                            del hNum, hDen
                    bar()
            h3Eff_mes[-1].Draw('BOX2 colz same')
            LatLabel(f'{origin} {mes_label[meson]}', 0.2, 0.9, textsize=0.05)
            outfile.mkdir(f'{meson}_{origin}')
            outfile.cd(f'{meson}_{origin}')
            h3Eff_mes[-1].Write()
            hEff_pt_phi_mes[-1].Write()
            hEff_pt_y_mes[-1].Write()
            hEff_phi_y_mes[-1].Write()
            outfile.cd('..')
            canvas_mes.Update()
    outfile.cd('..')
    canvas_mes.Write()
    del sparseDmeson

    # v0s
    h3Eff_v0, hEff_pt_phi_v0, hEff_pt_y_v0, hEff_phi_y_v0 = [], [], [], []
    canvas_v0 = TCanvas('cResoEffv0', '', 1600, 900)
    canvas_v0.Divide(len(mesons))
    for v0 in v0s:
        for imes, (mes) in enumerate(mesons):
            canvas_v0.cd(imes+1)
            outfile.mkdir(f'{v0}_{mes}')
            h3Eff_v0.append(TH3F(f'h3Eff_{mes}_{v0}',
                                 f'{v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                                 npt_bins, pt_bins, ny_bins, y_bins,
                                 nphi_bins, phi_bins))
            hEff_pt_phi_v0.append(TH2F(f'hEff_pt_phi_{mes}_{v0}',
                                  f'{v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #phi',
                                  npt_bins, pt_bins, nphi_bins, phi_bins))
            hEff_pt_y_v0.append(TH2F(f'hEff_pt_y_{mes}_{v0}',
                                     f'{v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #it{{y}}',
                                     npt_bins, pt_bins, ny_bins, y_bins))
            hEff_phi_y_v0.append(TH2F(f'hEff_phi_y_{mes}_{v0}',
                                      f'{v0} efficiency;#phi; #it{{y}}',
                                      nphi_bins, phi_bins, ny_bins, y_bins))
            outfile.cd(f'{v0}_{mes}')
            with alive_bar(len(pt_bins)-1, bar='smooth',
                           spinner='classic', title=f'Filling {mes} {v0} efficiency') as bar:
                for ipt, (pt_min, pt_max) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
                    pt_bin_min = sparseV0[mes]['Reco'][v0].GetAxis(0).FindBin(pt_min*1.0001)
                    pt_bin_max = sparseV0[mes]['Reco'][v0].GetAxis(0).FindBin(pt_max*0.9999)
                    sparseV0[mes]['Reco'][v0].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                    sparseV0[mes]['Gen'][v0].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                    for iphi, (phi_min, phi_max) in enumerate(zip(phi_bins[:-1], phi_bins[1:])):
                        phi_bin_min = sparseV0[mes]['Reco'][v0].GetAxis(2).FindBin(phi_min*1.0001)
                        phi_bin_max = sparseV0[mes]['Reco'][v0].GetAxis(2).FindBin(phi_max*0.9999)
                        sparseV0[mes]['Reco'][v0].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                        sparseV0[mes]['Gen'][v0].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                        for iy, (y_min, y_max) in enumerate(zip(y_bins[:-1], y_bins[1:])):
                            y_bin_min = sparseV0[mes]['Reco'][v0].GetAxis(1).FindBin(y_min*1.0001)
                            y_bin_max = sparseV0[mes]['Reco'][v0].GetAxis(1).FindBin(y_max*0.9999)
                            sparseV0[mes]['Reco'][v0].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                            sparseV0[mes]['Gen'][v0].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                            hNum = sparseV0[mes]['Reco'][v0].Projection(0, 1, 2)
                            hNum.SetDirectory(0)
                            hDen = sparseV0[mes]['Gen'][v0].Projection(0, 1, 2)
                            hDen.SetDirectory(0)
                            num = hNum.Integral()
                            den = hDen.Integral()
                            eff = num/den if den > 0 else 0
                            h3Eff_v0[-1].SetBinContent(ipt+1, iy+1, iphi+1, eff)
                            hEff_pt_phi_v0[-1].SetBinContent(ipt+1, iphi+1, eff)
                            hEff_pt_y_v0[-1].SetBinContent(ipt+1, iy+1, eff)
                            hEff_phi_y_v0[-1].SetBinContent(iphi+1, iy+1, eff)
                            del hNum, hDen
                    bar()
                h3Eff_v0[-1].Write()
                hEff_pt_phi_v0[-1].Write()
                hEff_pt_y_v0[-1].Write()
                hEff_phi_y_v0[-1].Write()
                outfile.cd('..')
        h3Eff_mes[-1].Draw('BOX2 colz same')
        LatLabel(f'{v0s_label[v0]}', 0.2, 0.8, textsize=0.04)
    del sparseV0

    canvas_mes.SaveAs(f'{inputCfg["eff"]["outfile"].replace(".root", "")}_Dmeson_efficiency.pdf')
    canvas_v0.SaveAs(f'{inputCfg["eff"]["outfile"].replace(".root", "")}_V0s_efficiency.pdf')

    input('Press enter to exit')


def SparseLoader(infile, dirname, listname, sparsetype, verbose=False):
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

def DecodeCutsToAxis(selections):
    '''
    Decode the cuts to the axis of the THnSparse

    Parameters:
    selections (dict): dictionary with the cuts

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
    return selsToaxis

def ApplyAnalysesCutsToSparse(sparse, axisnum, cutvar_min, cutvar_max, pt_mins=[], pt_maxs=[], verbose=False):
    '''
    Apply the analysis cuts to the THnSparse

    Parameters:
    sparse (THnSparse): THnSparse
    axisnum (int): number of the axis to apply the cut
    cutvar_min (float or list): minimum value of the cut variable or list of minimum values
    cutvar_max (float or list): maximum value of the cut variable or list of maximum values
    pt_mins (list): list of minimum values of the pt, if empty the cut is not pt dependent (default: [])
    pt_maxs (list): list of maximum values of the pt, if empty the cut is not pt dependent (default: [])

    Returns:
    sparse (THnSparse): THnSparse with the analysis cuts applied
    '''
    if axisnum is None:
        print(f'No analysis cuts to apply to {sparse.GetName()}')
        return sparse
        
    if (pt_mins == [] or pt_maxs == []):
        if cutvar_min:
            bin_min = sparse.GetAxis(axisnum).FindBin(cutvar_min * 1.0001)
        else:
            bin_min = sparse.GetAxis(axisnum).FindBin(sparse.GetAxis(axisnum).GetXmin())
        if cutvar_max:
            bin_max = sparse.GetAxis(axisnum).FindBin(cutvar_max * 0.9999)
        else:
            bin_max = sparse.GetAxis(axisnum).FindBin(sparse.GetAxis(axisnum).GetXmax())
        if verbose: print(f'Applying cuts to {sparse.GetName()} on axis {axisnum} = [{cutvar_min}, {cutvar_max}]')
        sparse.GetAxis(axisnum).SetRange(bin_min, bin_max)
    else:
        if verbose: print(f'Applying pt dependent analysis cuts to {sparse.GetName()} on axis {axisnum}')
        for ipt, (ptmin, ptmax) in enumerate(zip(pt_mins, pt_maxs)):
            if verbose: print(f' pt = [{ptmin}, {ptmax}], cut = [{cutvar_min[ipt]}, {cutvar_max[ipt]}]')
            bin_min = sparse.GetAxis(axisnum).FindBin(cutvar_min[ipt] * 1.0001)
            bin_max = sparse.GetAxis(axisnum).FindBin(cutvar_max[ipt] * 0.9999)
            sparse.GetAxis(axisnum).SetRange(bin_min, bin_max)
            sparse.GetAxis(axisnum).SetRange(bin_min, bin_max)
    return sparse

main()
