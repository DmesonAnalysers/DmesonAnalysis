'''
Script to project the MC distributions and apply the pt weights from the AnRes.root of Dtask
python3 proj_thn_mc.py config_flow.yml config_cutset.yml -o path/to/output -s text
                                                        --ptWeights path/to/file histName 
                                                        --ptWeightsB path/to/file histName
'''
import ROOT
import uproot
import yaml
import argparse
import sys
import os
from ROOT import gROOT, TFile
from alive_progress import alive_bar
from scipy.interpolate import InterpolatedUnivariateSpline
from sparse_dicts import get_sparses
sys.path.append('..')
from flow_analysis_utils import get_vn_versus_mass, get_centrality_bins

### please fill your path of DmesonAnalysis
sys.path.append('../../..')

def proj_data(cent_min, cent_max, pt_min, pt_max, sparseFlow, axes, inv_mass_bins, reso):
    hist_mass = sparseFlow.Projection(axes['Flow']['mass'])
    hist_mass.SetName(f'hist_mass_cent{cent_min}_{cent_max}_pt{pt_min}_{pt_max}')
    hist_mass.Write()

    # project the vn
    hist_vn_sp = get_vn_versus_mass(sparseFlow, inv_mass_bins, axes['Flow']['mass'], axes['Flow']['sp'])
    hist_vn_sp.SetDirectory(0)
    hist_vn_sp.SetName(f'hist_vn_sp_pt{pt_min}_{pt_max}')
    if reso > 0:
        hist_vn_sp.Scale(1./reso)
    print(f"hist_vn_sp.Integral(): {hist_vn_sp.Integral()}")
    hist_vn_sp.Write()

def proj_MC(config, cutsetConfig, ptweights, ptweightsB, outputdir, suffix):
def proj_mc_reco(ptLowLabel, ptHighLabel, config, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB):
    # read configuration
    particleName = config['Dmeson']
    enableRef = config['enableRef']
    enableSecPeak = config['enableSecPeak']
    Bspeciesweights = config['Bspeciesweights'] if 'Bspeciesweights' in config else None

    #TODO: safety checks for Dmeson reflecton and secondary peak

    # check the arguments
    if ptweights == []:
        print('\033[91m WARNING: pt weights will not be provided! \033[0m')
        ptweights = None
    if ptweightsB == []:
        print('\033[91m WARNING: B weights will not not be provided! \033[0m')
        ptweightsB = None
    if Bspeciesweights == None:
        print('\033[91m WARNING: B species weights will not be provided! \033[0m')
    enableRef = config.get('enableRef')
    enableSecPeak = config.get('enableSecPeak')

    for iVar in ('mass', 'pt'):
        hVarRefl, hVarReflPrompt, hVarReflFD, hVarPrompt, hVarFD = None, None, None, None, None
        if 'RecoAll' in sparsesReco:
            hVar = sparsesReco['RecoAll'].Projection(axes['RecoAll'][iVar])
            hVar.SetName(f'h{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            outfile.cd()
            hVar.Write()

        if iVar != 'pt':
            hVarPrompt = sparsesReco['RecoPrompt'].Projection(axes['RecoPrompt'][iVar])
            hVarPrompt.SetName(f'hPrompt{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarFD = sparsesReco['RecoFD'].Projection(axes['RecoFD'][iVar])
            hVarFD.SetName(f'hFD{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        else:
            ### no pt weights
            if not ptWeights and not ptWeightsB:
                hVarPrompt = sparsesReco['RecoPrompt'].Projection(axes['RecoPrompt'][iVar])
                hVarFD = sparsesReco['RecoFD'].Projection(axes['RecoFD'][iVar])

        ### pt weights for prompt
        if ptWeights:
            hVarPrompt = sparsesReco['RecoPrompt'].Projection(axes['RecoPrompt'][iVar])
            for iBin in range(1, hVarPrompt.GetNbinsX()+1):
                if hVarPrompt.GetBinContent(iBin) > 0.:
                    relStatUnc = hVarPrompt.GetBinError(iBin) / hVarPrompt.GetBinContent(iBin)
                    ptCent = hVarPrompt.GetBinCenter(iBin)
                    hVarPrompt.SetBinContent(iBin, hVarPrompt.GetBinContent(iBin) * sPtWeights(ptCent))
                    hVarPrompt.SetBinError(iBin, hVarPrompt.GetBinContent(iBin) * relStatUnc)

        ### pt weights for non-prompt, but no B species weights
        if ptWeightsB and not Bspeciesweights:
            hPtBvsPtD = sparsesReco['RecoFD'].Projection(2, axes['RecoFD'][iVar])
            for iPtD in range(1, hPtBvsPtD.GetXaxis().GetNbins()+1):
                for iPtB in range(1, hPtBvsPtD.GetYaxis().GetNbins()+1):
                    ptCentB = hPtBvsPtD.GetYaxis().GetBinCenter(iPtB)
                    origContent = hPtBvsPtD.GetBinContent(iPtD, iPtB)
                    origError = hPtBvsPtD.GetBinError(iPtD, iPtB)
                    weight = 0
                    if sPtWeightsB(ptCentB) > 0:
                        weight = sPtWeightsB(ptCentB)
                    content = origContent * weight
                    error = 0
                    if origContent > 0:
                        error = origError / origContent * content
                    hPtBvsPtD.SetBinContent(iPtD, iPtB, content)
                    hPtBvsPtD.SetBinError(iPtD, iPtB, error)
            hVarFD = hPtBvsPtD.ProjectionX(f'hFD{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}',
                                            0, hPtBvsPtD.GetYaxis().GetNbins()+1, 'e')
        ### pt weights from B for non-prompt and B species weights
        elif ptWeightsB and Bspeciesweights:
            hPtBvsBspecievsPtD = sparsesReco['RecoFD'].Projection(axes['RecoFD'][iVar], 3, 2)
            for iPtD in range(1, hPtBvsBspecievsPtD.GetXaxis().GetNbins()+1):
                for iBspecie in range(1, hPtBvsBspecievsPtD.GetYaxis().GetNbins()+1):
                    for iPtB in range(1, hPtBvsBspecievsPtD.GetZaxis().GetNbins()+1):
                        ptCentB = hPtBvsBspecievsPtD.GetZaxis().GetBinCenter(iPtB)
                        origContent = hPtBvsBspecievsPtD.GetBinContent(iPtD, iBspecie, iPtB)
                        origError = hPtBvsBspecievsPtD.GetBinError(iPtD, iBspecie, iPtB)
                        weight = Bspeciesweights[iBspecie-1]
                        if sPtWeightsB(ptCentB) > 0:
                            weight *= sPtWeightsB(ptCentB)
                        content = origContent * weight
                        error = 0
                        if origContent > 0:
                            error = origError / origContent * content
                        hPtBvsBspecievsPtD.SetBinContent(iPtD, iBspecie, iPtB, content)
                        hPtBvsBspecievsPtD.SetBinError(iPtD, iBspecie, iPtB, error)
            hVarFD = hPtBvsBspecievsPtD.ProjectionX(f'hFD{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}',
                                                    0, hPtBvsBspecievsPtD.GetYaxis().GetNbins()+1,
                                                    0, hPtBvsBspecievsPtD.GetZaxis().GetNbins()+1, 'e')
        ### only B species weights
        elif Bspeciesweights:
            hBspecievsPtD = sparsesReco['RecoFD'].Projection(3, axes['RecoFD'][iVar])
            for iPtD in range(1, hBspecievsPtD.GetXaxis().GetNbins()+1):
                for iBspecie in range(1, hBspecievsPtD.GetYaxis().GetNbins()+1):
                    origContent = hBspecievsPtD.GetBinContent(iPtD, iBspecie)
                    origError = hBspecievsPtD.GetBinError(iPtD, iBspecie)
                    weight = Bspeciesweights[iBspecie-1]
                    content = origContent * weight
                    error = 0
                    if origContent > 0:
                        error = origError / origContent * content
                    hBspecievsPtD.SetBinContent(iPtD, iBspecie, content)
                    hBspecievsPtD.SetBinError(iPtD, iBspecie, error)
            hVarFD = hBspecievsPtD.ProjectionX(f'hFD{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}',
                                            0, hBspecievsPtD.GetYaxis().GetNbins()+1, 'e')
        ### use the pt weights from prompt for non-prompt
        elif ptWeights: # if pt weights for prompt are present apply them
            for iBin in range(1, hVarFD.GetNbinsX()+1):
                    if hVarFD.GetBinContent(iBin) > 0.:
                        relStatUnc = hVarFD.GetBinError(iBin) / hVarFD.GetBinContent(iBin)
                        ptCent = hVarFD.GetBinCenter(iBin)
                        hVarFD.SetBinContent(iBin, hVarFD.GetBinContent(iBin) * sPtWeights(ptCent))
                        hVarFD.SetBinError(iBin, hVarFD.GetBinContent(iBin) * relStatUnc)

        ## write the output           
        hVarPrompt.SetName(f'hPrompt{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarPrompt.Write()
        hVarFD.SetName(f'hFD{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarFD.Write()
        if config.get('enableRef'):
            hVarRefl = sparsesReco['RecoRefl'].Projection(axes['RecoRefl'][iVar])
            hVarRefl.SetName(f'hVarRefl{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarRefl.Write()
            hVarReflPrompt = sparsesReco['RecoReflPrompt'].Projection(axes['RecoReflPrompt'][iVar])
            hVarReflPrompt.SetName(f'hVarReflPrompt{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarReflPrompt.Write()
            hVarReflFD = sparsesReco['RecoReflFD'].Projection(axes['RecoReflFD'][iVar])
            hVarReflFD.SetName(f'hVarReflFD{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarReflFD.Write()
        if config.get('enableSecPeak'):
            hVarPromptSecPeak = sparsesReco['RecoSecPeakPrompt'].Projection(axes['RecoSecPeakPrompt'][iVar])
            hVarPromptSecPeak.SetName(f'hPromptSecPeak{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarPromptSecPeak.Write()
            hVarFDSecPeak = sparsesReco['RecoSecPeakFD'].Projection(axes['RecoSecPeakFD'][iVar])
            hVarFDSecPeak.SetName(f'hFDSecPeak{iVar}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarFDSecPeak.Write()

def proj_mc_gen(ptLowLabel, ptHighLabel, config, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB):

    ## no pt weights
    if not ptWeights and not ptWeightsB:
        hGenPtPrompt = sparsesGen['GenPrompt'].Projection(axes['GenPrompt']['pt'])
        hGenPtFD = sparsesGen['GenFD'].Projection(axes['GenFD']['pt'])

    ## apply pt weights for prompt
    if ptWeights:
        hGenPtPrompt = sparsesGen['GenPrompt'].Projection(axes['GenPrompt']['pt'])
        for iBin in range(1, hGenPtPrompt.GetNbinsX()+1):
            if hGenPtPrompt.GetBinContent(iBin) > 0:
                relStatUnc = hGenPtPrompt.GetBinError(iBin) / hGenPtPrompt.GetBinContent(iBin)
                ptCent = hGenPtPrompt.GetBinCenter(iBin)
                hGenPtPrompt.SetBinContent(iBin, hGenPtPrompt.GetBinContent(iBin) * sPtWeights(ptCent))
                hGenPtPrompt.SetBinError(iBin, hGenPtPrompt.GetBinContent(iBin) * relStatUnc)

    ## apply pt weights for non-prompt
    ### pt weights from B for non-prompt, but no B species weights
    if ptWeightsB and not Bspeciesweights:
        hPtBvsPtGenD = sparsesGen['GenFD'].Projection(2, 0)
        for iPtD in range(1, hPtBvsPtGenD.GetXaxis().GetNbins()+1):
            for iPtB in range(1, hPtBvsPtGenD.GetYaxis().GetNbins()+1):
                ptCentB = hPtBvsPtGenD.GetYaxis().GetBinCenter(iPtB)
                origContent = hPtBvsPtGenD.GetBinContent(iPtD, iPtB)
                origError = hPtBvsPtGenD.GetBinError(iPtD, iPtB)
                weight = 0
                if sPtWeightsB(ptCent) > 0:
                    weight = sPtWeightsB(ptCentB)
                content = hPtBvsPtGenD.GetBinContent(iPtD, iPtB) * weight
                error = 0
                if origContent > 0:
                    error = origError / origContent * content
                hPtBvsPtGenD.SetBinContent(iPtD, iPtB, content)
                hPtBvsPtGenD.SetBinError(iPtD, iPtB, error)
        hGenPtFD = hPtBvsPtGenD.ProjectionX(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}',
                                            0, hPtBvsPtGenD.GetYaxis().GetNbins()+1, 'e')
    ### pt weights from B for non-prompt and B species weights
    elif ptWeightsB and Bspeciesweights:
        hPtBvsBspecievsPtGenD = sparsesGen['GenFD'].Projection(0, 3, 2)
        for iPtD in range(1, hPtBvsBspecievsPtGenD.GetXaxis().GetNbins()+1):
            for iBspecie in range(1, hPtBvsBspecievsPtGenD.GetYaxis().GetNbins()+1):
                for iPtB in range(1, hPtBvsBspecievsPtGenD.GetZaxis().GetNbins()+1):
                    ptCentB = hPtBvsBspecievsPtGenD.GetZaxis().GetBinCenter(iPtB)
                    origContent = hPtBvsBspecievsPtGenD.GetBinContent(iPtD, iBspecie, iPtB)
                    origError = hPtBvsBspecievsPtGenD.GetBinError(iPtD, iBspecie, iPtB)
                    weight = Bspeciesweights[iBspecie-1]
                    if sPtWeightsB(ptCentB) > 0:
                        weight *= sPtWeightsB(ptCentB)
                    content = origContent * weight
                    error = 0
                    if origContent > 0:
                        error = origError / origContent * content
                    hPtBvsBspecievsPtGenD.SetBinContent(iPtD, iBspecie, iPtB, content)
                    hPtBvsBspecievsPtGenD.SetBinError(iPtD, iBspecie, iPtB, error)
        hGenPtFD = hPtBvsBspecievsPtGenD.ProjectionX(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}',
                                                    0, hPtBvsBspecievsPtGenD.GetYaxis().GetNbins()+1,
                                                    0, hPtBvsBspecievsPtGenD.GetZaxis().GetNbins()+1, 'e')
    ### only B species weights
    elif Bspeciesweights:
        hBspecievsPtGenD = sparsesGen['GenFD'].Projection(3, 0)
        for iPtD in range(1, hBspecievsPtGenD.GetXaxis().GetNbins()+1):
            for iBspecie in range(1, hBspecievsPtGenD.GetYaxis().GetNbins()+1):
                origContent = hBspecievsPtGenD.GetBinContent(iPtD, iBspecie)
                origError = hBspecievsPtGenD.GetBinError(iPtD, iBspecie)
                weight = Bspeciesweights[iBspecie-1]
                content = origContent * weight
                error = 0
                if origContent > 0:
                    error = origError / origContent * content
                hBspecievsPtGenD.SetBinContent(iPtD, iBspecie, content)
                hBspecievsPtGenD.SetBinError(iPtD, iBspecie, error)
        hGenPtFD = hBspecievsPtGenD.ProjectionX(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}',
                                                0, hBspecievsPtGenD.GetYaxis().GetNbins()+1, 'e')
    ### use the pt weights from prompt for non-prompt
    elif ptWeights: # if pt weights for prompt are present apply them
        hGenPtFD = sparsesGen['GenFD'].Projection(axes['GenFD']['pt'])
                            
    ## print the output
    hGenPtPrompt.SetName(f'hPromptGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
    hGenPtPrompt.Write()
    hGenPtFD.SetName(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
    hGenPtFD.Write()
    if config.get('enableSecPeak'):
        hGenPtPromptSecPeak = sparsesGen['GenSecPeakPrompt'].Projection(axes['GenSecPeakPrompt']['pt'])
        hGenPtPromptSecPeak.SetName(f'hPromptSecPeakGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hGenPtPromptSecPeak.Write()
        hGenPtFDSecPeak = sparsesGen['GenSecPeakFD'].Projection(axes['GenSecPeakFD']['pt'])
        hGenPtFDSecPeak.SetName(f'hFDSecPeakGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hGenPtFDSecPeak.Write()

def pt_weights_info(config):
    # compute info for pt weights
    if config.get('ptWeights'):
        ptWeights = uproot.open(config['ptweightPath'])[config['ptweightName']]
        bins = ptWeights.axis(0).edges()
        ptCentW = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
        sPtWeights = InterpolatedUnivariateSpline(ptCentW, ptWeights.values())
    else:
        print('\033[91m WARNING: pt weights will not be provided! \033[0m')
        ptWeights = None
        sPtWeights = None

    if config.get('ptWeightsB'):
        ptWeightsB = uproot.open(config['ptweightBPath'])[config['ptweightBName']]
        bins = ptWeightsB.axis(0).edges()
        ptCentWB = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
        sPtWeightsB = InterpolatedUnivariateSpline(ptCentWB, ptWeightsB.values())
    else:
        print('\033[91m WARNING: B weights will not not be provided! \033[0m')
        ptWeightsB = None
        sPtWeightsB = None

    if config.get('Bspeciesweights'):
        Bspeciesweights = config['Bspeciesweights']
    else:
        print('\033[91m WARNING: B species weights will not be provided! \033[0m')
        Bspeciesweights = None
    
    return ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="flow configuration file")
    parser.add_argument('cutsetConfig', metavar='text',
                        default='cutsetConfig.yaml', help='cutset configuration file')
    parser.add_argument("--ptweights", "-w", metavar="text", nargs=2, required=False,
                        default=[], help="path to pt weights file and histogram name")
    parser.add_argument("--ptweightsB", "-wb", metavar="text", nargs=2, required=False,
                        default=[], help="path to pt weightsB file and histogram name")
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k3050", help="centrality class")
    parser.add_argument("--resolution", "-r", metavar="text",
                        default="reso.root", help="resolution file")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    with open(args.config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)


    os.makedirs(f'{args.outputdir}/proj', exist_ok=True)
    outfile = ROOT.TFile(f'{args.outputdir}/proj/proj_{args.suffix}.root', 'RECREATE')
    
    outfile_dir = 'hf-candidate-creator-2prong' if config['Dmeson'] == 'Dzero' else 'hf-candidate-creator-3prong'
    infilemc = TFile.Open(config['MC_filename'], 'r')
    histo_cent = infilemc.Get(f'{outfile_dir}/hSelCollisionsCent')
    resofile = TFile.Open(args.resolution, 'r')
    try:
        det_A = config['detA']
        det_B = config['detB']
        det_C = config['detC']
        histo_reso = resofile.Get(f'{det_A}_{det_B}_{det_C}/histo_reso')
        histo_reso.SetName('hist_reso')
        histo_reso.SetDirectory(0)
        reso = histo_reso.GetBinContent(1)
    except:
        histo_reso = resofile.Get(f'hf-task-flow-charm-hadrons/spReso/hSpReso{det_B}{det_C}')
        histo_reso.SetName('hist_reso')
        histo_reso.SetDirectory(0)
        reso = histo_reso.GetBinContent(1)
    
    outfile.mkdir(outfile_dir)
    outfile.cd(outfile_dir)
    histo_reso.Write()
    histo_cent.Write()
    resofile.Close()
    infilemc.Close()
    
    # quit()

    # load thnsparse
    sparseFlow, sparsesReco, sparsesGen, axes = get_sparses(config)
    cent_min, cent_max = get_centrality_bins(args.centrality)[1]
    sparseFlow.GetAxis(axes['Flow']['cent']).SetRangeUser(cent_min, cent_max)
    for key, iSparse in sparsesGen.items():
        iSparse.GetAxis(axes[key]['cent']).SetRangeUser(cent_min, cent_max)

    with open(args.cutsetConfig, 'r') as ymlCutSetFile:
        cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
    cutVars = cutSetCfg['cutvars']
    print(f"cutVars: {cutVars}")

    # compute info for pt weights
    ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB = pt_weights_info(config)

    with alive_bar(len(cutVars['pt']['min']), title='Processing pT bins') as bar:
        for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['pt']['min'], cutVars['pt']['max'])):
            print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
            ptLowLabel = ptMin * 10
            ptHighLabel = ptMax * 10
            
            ## apply cuts
            for iVar in cutVars:
                sparseFlow.GetAxis(axes['Flow'][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                for key, iSparse in sparsesReco.items():
                    iSparse.GetAxis(axes[key][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                if iVar == 'pt':
                    for key, iSparse in sparsesGen.items():
                        iSparse.GetAxis(axes[key][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                if iVar == 'score_FD' or iVar == 'score_bkg':
                    print(f'{iVar}: {cutVars[iVar]["min"][iPt]} < {iVar} < {cutVars[iVar]["max"][iPt]}')

            outfile.mkdir(f'cent_bins{cent_min}_{cent_max}/pt_bins{ptMin}_{ptMax}')
            outfile.cd(f'cent_bins{cent_min}_{cent_max}/pt_bins{ptMin}_{ptMax}')
            proj_data(cent_min, cent_max, ptMin, ptMax, sparseFlow, axes, config['inv_mass_bins'][iPt], reso)
            print(f"Projected data!")
            # outfile.cd()
            proj_mc_reco(ptLowLabel, ptHighLabel, config, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB)
            print(f"Projected mc reco!")
            proj_mc_gen(ptLowLabel, ptHighLabel, config, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB)
            print(f"Projected mc gen!")
            
            bar()
    
    outfile.Close()
    outfile.Close()
