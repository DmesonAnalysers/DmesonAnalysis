'''
python script to calculate the efficiency of the resonances
run: python ComputeResoDauEff.py config_proj_reso.yml
'''
import sys
import argparse
import yaml
import numpy as np
from alive_progress import alive_bar
from particle import Particle
from ROOT import TFile, TH3F, TH2F, TCanvas, TLatex
sys.path.insert(0, '..')
#from utils.TaskFileLoader import LoadSingleSparseFromTask, LoadNormObjFromTask, LoadSingleSparseFromTask
from utils.StyleFormatter import SetGlobalStyle

SetGlobalStyle(padtopmargin=0.05, padrightmargin=0.15, padleftmargin=0.15,
               padbottommargin=0.15, palette=55, labelsize=0.03, titlesize=0.03,
               labeloffset=0.01, titleoffsety=1.8, titleoffsetx=1.8, titleoffsetz=1.,
               opttitle=0, optstat=0)
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
    weights = inputCfg['eff']['weights']
    resonances = inputCfg['eff']['resonances']
    origins = ['Prompt', 'NonPrompt']
    ptshapes = ['Dsptshape', 'Lcptshape'] # 'Dsptshape' or 'Lcptshape'

    if mult not in ['MB', 'HM']:
        raise ValueError(f'Wrong multiplicity: {mult}.\
                         Allowed values at the moment: MB, HM')
    for i, (dir, list) in enumerate(zip(input_dirs, input_list)):
        input_dirs[i] = f'{dir}_{mult}{weights}'
        input_list[i] = f'coutput{list}_{mult}{weights}'

    mesons_code = {'Dplus': 411, 'Dstar': 413}
    v0s_code = {'K0S': 310}
    reso_label = {10433: Particle.from_pdgid(10433).name, 435: Particle.from_pdgid(435).name}
    mes_label = {'Dplus': 'D^{+}', 'Dstar': 'D^{*+}'}
    v0s_label = {'K0S': 'K^{0}_{S}'}

    print('Configuration file loaded successfully.')
    print('Resonances to be analysed: ', resonances)
    print('Input files: ', input_files)
    print('Input dirs: ', input_dirs)
    print('Input list: ', input_list)
    print('Multiplicity: ', mult)
    print('Weights: ', weights)
    print('_______________________________________________________')
        
    for ireso, (reso) in enumerate(resonances):
        meson, v0 = None, None

        if reso == 10433:
            meson = 'Dstar'
            v0 = 'K0S'
        elif reso == 435:
            meson = 'Dplus'
            v0 = 'K0S'
        else:
            raise ValueError(f'Wrong resonance: {reso}.\
                             Allowed values at the moment: 10433, 435')
        print(f'\033[1mComputing efficiency for:  {reso_label[reso]} --> {mes_label[meson]} + {v0s_label[v0]}\033[0m')

        #__________________________________________________________________________
        # Loading the THnSparse
        sparseDmeson, sparseV0 = {}, {}
        print('\033[1mLoading the THnSparse...\033[0m\r', end='')
        for file in input_files: # extend to multiple files
            dir = input_dirs[ireso]
            list = input_list[ireso]
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
                v0_code = v0s_code[v0]
                sparseV0[meson][status][v0] = {}
                if status == 'Reco':
                    sparseV0[meson][status][v0] = SparseLoader(file, dir, list, f'hV0MCReco_{v0_code}')
                else:
                    sparseV0[meson][status][v0] = SparseLoader(file, dir, list, f'hV0MCGen_{v0_code}')
        print('\033[1mLoading the THnSparse: complete.\033[0m')

        #__________________________________________________________________________
        # Computing the efficiency
        outfile_name = f'{inputCfg["eff"]["outfile"]}'.replace('.root', f'_{reso}{weights}.root')
        outfile = TFile(outfile_name, 'recreate')
        outfile.cd()
        pt_bins = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 16, 24] # hardcoded for now
        npt_bins = len(pt_bins) - 1
        pt_bins = np.array(pt_bins, dtype=np.float32)
        y_bins = [-1.0, -0.8, -0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, 0.8, 1.0] # hardcoded for now
        ny_bins = len(y_bins) - 1
        y_bins = np.array(y_bins, dtype=np.float32)
        phi_bins = [np.linspace(0., 2*np.pi, 10)] # hardcoded for now
        phi_bins = [round(phi, 2) for phi in np.concatenate(phi_bins)]
        nphi_bins = len(phi_bins) - 1
        phi_bins = np.array(phi_bins, dtype=np.float32)
        selections = inputCfg['selections']

        # mesons
        meson_code = mesons_code[meson]
        sels = selections[meson_code][mult]
        sels_decoded = DecodeCutsToAxis(sels, verbose=False)

        h3Eff_mes, h3Num, h3Den,\
        hEff_pt_phi_mes, hEff_pt_y_mes,\
        hEff_phi_y_mes = ([] for i in range(6))
        canvas_mes = TCanvas('cResoEffMeson', '', 2000, 900)
        canvas_mes.Divide(2)
    
        for iorigin, (origin) in enumerate(origins):
            canvas_mes.cd(iorigin+1)
            canvas_mes.cd(iorigin+1).SetTitle(f'{meson}_{origin}')
            h3Eff_mes.append(TH3F(f'h3Eff_mes{meson}{origin}',
                              f';#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                              npt_bins, pt_bins, ny_bins, y_bins,
                              nphi_bins, phi_bins))
            h3Num.append(TH3F(f'h3Num_mes{meson}{origin}',
                              f';#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                              npt_bins, pt_bins, ny_bins, y_bins,
                              nphi_bins, phi_bins))
            h3Den.append(TH3F(f'h3Den_mes{meson}{origin}',
                              f';#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                              npt_bins, pt_bins, ny_bins, y_bins,
                              nphi_bins, phi_bins))
            h3Eff_mes[-1].GetXaxis().CenterTitle()
            h3Eff_mes[-1].GetYaxis().CenterTitle()
            h3Eff_mes[-1].GetZaxis().CenterTitle()
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
                    lower_pt_thr = sels_decoded['pt_min']['cutvar_min']
                    if pt_min < lower_pt_thr:
                        print(f'WARNING: {meson} {origin} efficiency for pt < {lower_pt_thr} GeV/c is not computed')
                        continue
                    sparseDmeson[meson]['Reco'][origin] = ApplyAnalysesCutsToSparse(sparseDmeson[meson]['Reco'][origin],
                                                                                    sels_decoded,
                                                                                    pt_min=pt_min,
                                                                                    pt_max=pt_max,
                                                                                    verbose=False)
                    pt_bin_min = sparseDmeson[meson]['Gen'][origin].GetAxis(0).FindBin(pt_min*
                                                                                       1.0001)
                    pt_bin_max = sparseDmeson[meson]['Gen'][origin].GetAxis(0).FindBin(pt_max*
                                                                                       0.9999)
                    sparseDmeson[meson]['Gen'][origin].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                    sparseDmeson[meson]['Reco'][origin].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                    for iphi, (phi_min, phi_max) in enumerate(zip(phi_bins[:-1], phi_bins[1:])):
                        phi_bin_min = sparseDmeson[meson]['Gen'][origin].GetAxis(2).FindBin(phi_min*1.0001)
                        phi_bin_max = sparseDmeson[meson]['Gen'][origin].GetAxis(2).FindBin(phi_max*0.9999)
                        sparseDmeson[meson]['Gen'][origin].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                        sparseDmeson[meson]['Reco'][origin].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                        for iy, (y_min, y_max) in enumerate(zip(y_bins[:-1], y_bins[1:])):
                            y_bin_min = sparseDmeson[meson]['Gen'][origin].GetAxis(1).FindBin(y_min*1.0001)
                            y_bin_max = sparseDmeson[meson]['Gen'][origin].GetAxis(1).FindBin(y_max*0.9999)
                            sparseDmeson[meson]['Gen'][origin].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                            sparseDmeson[meson]['Reco'][origin].GetAxis(1).SetRange(y_bin_min, y_bin_max)

                            hNum = sparseDmeson[meson]['Reco'][origin].Projection(0, 1, 2)
                            hDen = sparseDmeson[meson]['Gen'][origin].Projection(0, 1, 2)
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
            h3Eff_mes[-1] = h3Num[-1].Clone(f'h3Eff{meson}{origin}')
            h3Eff_mes[-1].Divide(h3Eff_mes[-1], h3Den[-1], 1, 1, 'B')
            hy_pt_y_num = h3Num[-1].Project3D('yxe')
            hy_pt_y_den = h3Den[-1].Project3D('yxe')
            hy_pt_y_eff = hy_pt_y_num.Clone(f'hy_pt_y_eff{meson}{origin}')
            hy_pt_y_eff.Divide(hy_pt_y_eff, hy_pt_y_den, 1, 1, 'B')
            h_pt_num = hy_pt_y_num.ProjectionX()
            h_pt_den = hy_pt_y_den.ProjectionX()
            h_pt_eff = h_pt_num.Clone(f'h_pt_eff{meson}{origin}')
            h_pt_eff.Divide(h_pt_eff, h_pt_den, 1, 1, 'B')
            h3Eff_mes[-1].Draw('BOX2 colz')
            latLabel = TLatex()
            latLabel.SetNDC()
            latLabel.SetTextSize(0.05)
            latLabel.DrawLatex(0.2, 0.95, f'{origin} {mes_label[meson]}')
            outfile.mkdir(f'{meson}_{origin}')
            outfile.cd(f'{meson}_{origin}')
            h3Den[-1].Write()
            h3Num[-1].Write()
            h3Eff_mes[-1].Write()
            hy_pt_y_eff.Write(f'h3Dy_pt_y_eff_{meson}_{origin}')
            h_pt_eff.Write(f'h_pt_eff_{meson}_{origin}')
            hEff_pt_phi_mes[-1].Write()
            hEff_pt_y_mes[-1].Write()
            hEff_phi_y_mes[-1].Write()
            outfile.cd('..')
            canvas_mes.Update()
            #input('Press enter to continue')
        outfile.cd('..')
        canvas_mes.Write()
        del sparseDmeson

        # v0s
        v0_code = v0s_code[v0]
        v0_label = v0s_label[v0]
        sels = selections[v0_code][mult]
        sels_decoded = DecodeCutsToAxis(sels, verbose=False)

        h3Eff_v0, h3Num_v0, h3Den_v0,\
        hEff_pt_phi_v0, hEff_pt_y_v0,\
        hEff_phi_y_v0 = ([] for i in range(6))
        canvas_v0 = TCanvas('cResoEffv0', '', 1600, 900)
        canvas_v0.cd()
        outfile.mkdir(f'{v0}_{meson}')
        h3Eff_v0.append(TH3F(f'h3Eff_{meson}_{v0}',
                             f'{v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                             npt_bins, pt_bins, ny_bins, y_bins,
                             nphi_bins, phi_bins))
        h3Num_v0.append(TH3F(f'h3Num_{meson}_{v0}',
                             f'{v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                             npt_bins, pt_bins, ny_bins, y_bins,
                             nphi_bins, phi_bins))
        h3Den_v0.append(TH3F(f'h3Dum_{meson}_{v0}',
                             f'{v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}});#it{{y}}; #phi',
                             npt_bins, pt_bins, ny_bins, y_bins,
                             nphi_bins, phi_bins))
        hEff_pt_phi_v0.append(TH2F(f'hEff_pt_phi_{meson}_{v0}',
                              f'{v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #phi',
                              npt_bins, pt_bins, nphi_bins, phi_bins))
        hEff_pt_y_v0.append(TH2F(f'hEff_pt_y_{meson}_{v0}',
                                 f'{v0} efficiency;#it{{p}}_{{T}} (GeV/#it{{c}}); #it{{y}}',
                                 npt_bins, pt_bins, ny_bins, y_bins))
        hEff_phi_y_v0.append(TH2F(f'hEff_phi_y_{meson}_{v0}',
                                  f'{v0} efficiency;#phi; #it{{y}}',
                                  nphi_bins, phi_bins, ny_bins, y_bins))
        outfile.cd(f'{v0}_{meson}')
        with alive_bar(len(pt_bins)-1, bar='smooth',
                       spinner='classic', title=f'Filling {meson} {v0} efficiency') as bar:
            for ipt, (pt_min, pt_max) in enumerate(zip(pt_bins[:-1], pt_bins[1:])):
                lower_pt_thr = sels_decoded['pt_min']['cutvar_min']
                if pt_min < lower_pt_thr:
                    print(f'WARNING: {v0} efficiency for pt < {lower_pt_thr} GeV/c is not computed')
                    continue
                sparseV0[meson]['Reco'][v0] = ApplyAnalysesCutsToSparse(sparseV0[meson]['Reco'][v0],
                                                                        sels_decoded,
                                                                        pt_min=pt_min,
                                                                        pt_max=pt_max,
                                                                        verbose=False)
                pt_bin_min = sparseV0[meson]['Reco'][v0].GetAxis(0).FindBin(pt_min*1.0001)
                pt_bin_max = sparseV0[meson]['Reco'][v0].GetAxis(0).FindBin(pt_max*0.9999)
                sparseV0[meson]['Reco'][v0].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                sparseV0[meson]['Gen'][v0].GetAxis(0).SetRange(pt_bin_min, pt_bin_max)
                for iphi, (phi_min, phi_max) in enumerate(zip(phi_bins[:-1], phi_bins[1:])):
                    phi_bin_min = sparseV0[meson]['Reco'][v0].GetAxis(2).FindBin(phi_min*1.0001)
                    phi_bin_max = sparseV0[meson]['Reco'][v0].GetAxis(2).FindBin(phi_max*0.9999)
                    sparseV0[meson]['Reco'][v0].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                    sparseV0[meson]['Gen'][v0].GetAxis(2).SetRange(phi_bin_min, phi_bin_max)
                    for iy, (y_min, y_max) in enumerate(zip(y_bins[:-1], y_bins[1:])):
                        y_bin_min = sparseV0[meson]['Reco'][v0].GetAxis(1).FindBin(y_min*1.0001)
                        y_bin_max = sparseV0[meson]['Reco'][v0].GetAxis(1).FindBin(y_max*0.9999)
                        sparseV0[meson]['Reco'][v0].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                        sparseV0[meson]['Gen'][v0].GetAxis(1).SetRange(y_bin_min, y_bin_max)
                        hNum = sparseV0[meson]['Reco'][v0].Projection(0, 1, 2)
                        hNum.SetDirectory(0)
                        hDen = sparseV0[meson]['Gen'][v0].Projection(0, 1, 2)
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
            h3Den_v0[-1].Write(f'hDen_{meson}_{v0}')
            h3Num_v0[-1].Write(f'hNum_{meson}_{v0}')
            h3Eff_v0[-1] = h3Num_v0[-1].Clone(f'hEff_{meson}_{v0}')
            h3Eff_v0[-1].Divide(h3Num_v0[-1], h3Den_v0[-1], 1, 1, 'B')
            h_pt_y_num = h3Num_v0[-1].Project3D('yxe')
            h_pt_y_den = h3Den_v0[-1].Project3D('yxe')
            h_pt_y_eff = h_pt_y_num.Clone(f'hEff_pt_y_{meson}_{v0}')
            h_pt_y_eff.Divide(h_pt_y_eff, h_pt_y_den, 1, 1, 'B')
            h_pt_y_eff.Write(f'hEff_pt_y_{meson}_{v0}')
            h_pt_num = h3Num_v0[-1].ProjectionX()
            h_pt_den = h3Den_v0[-1].ProjectionX()
            h_pt_eff = h_pt_num.Clone(f'hEff_pt_{meson}_{v0}')
            h_pt_eff.Divide(h_pt_eff, h_pt_den, 1, 1, 'B')
            h_pt_eff.Write(f'hEff_pt_{meson}_{v0}')
            h3Eff_v0[-1].Write()
            hEff_pt_phi_v0[-1].Write()
            #hEff_pt_y_v0[-1].Write()
            hEff_phi_y_v0[-1].Write()
            outfile.cd('..')
        h3Eff_mes[-1].Draw('BOX2 colz same')
        latlabel = TLatex()
        latlabel.SetNDC()
        latlabel.SetTextSize(0.04)
        latlabel.DrawLatex(0.15, 0.85, f'{v0_label}')
        del sparseV0

        outfile.Close()
        print(f'Output written to {outfile_name}')

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

def DecodeCutsToAxis(selections, verbose=False):
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

def ApplyAnalysesCutsToSparse(sparse, sels_decoded, pt_min=None, pt_max=None, verbose=False):
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

main()
