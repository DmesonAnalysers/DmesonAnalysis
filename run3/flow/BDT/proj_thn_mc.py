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

def proj_data(sparse_flow, ptMin, ptMax, axes, inv_mass_bins, reso):

    if isinstance(sparse_flow, dict):
        for isparse, (key, sparse) in enumerate(sparse_flow.items()):
            hist_mass_temp = sparse.Projection(axes['Flow']['Mass'])
            hist_mass_temp.SetName(f'hist_mass_proj_{isparse}')
            hist_mass_temp.SetDirectory(0)

            if isparse == 0:
                hist_mass = hist_mass_temp.Clone('hist_mass_proj')
                hist_mass.SetDirectory(0)
                hist_mass.Reset()

            hist_mass.Add(hist_mass_temp)

        hist_vn_sp = get_vn_versus_mass(list(sparse_flow.values()), inv_mass_bins, axes['Flow']['Mass'], axes['Flow']['sp'])
        hist_vn_sp.SetDirectory(0)
        if reso > 0:
            hist_vn_sp.Scale(1./reso)
    else:
        hist_mass = sparse_flow.Projection(axes['Flow']['Mass'])
        hist_vn_sp = get_vn_versus_mass(sparse_flow, inv_mass_bins, axes['Flow']['Mass'], axes['Flow']['sp'])
        hist_vn_sp.SetDirectory(0)
        if reso > 0:
            hist_vn_sp.Scale(1./reso)

    hist_mass.Write(f'hist_mass_proj_pt{int(ptMin*10)}_{int(ptMax*10)}')
    hist_vn_sp.Write(f'hist_vn_sp_proj_pt{int(ptMin*10)}_{int(ptMax*10)}')

def proj_mc_reco(config, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB):
    
    ### project mass
    if 'RecoAll' in sparsesReco:
        hMass = sparsesReco['RecoAll'].Projection(axes['RecoAll']['Mass'])
        hPt = sparsesReco['RecoAll'].Projection(axes['RecoAll']['Pt'])

    hMassPrompt = sparsesReco['RecoPrompt'].Projection(axes['RecoPrompt']['Mass'])
    hMassFD = sparsesReco['RecoFD'].Projection(axes['RecoFD']['Mass'])    

    ### project pt
    hPtRefl, hPtReflPrompt, hPtReflFD, hPtPrompt, hPtFD = None, None, None, None, None

    ### no pt weights
    if not ptWeights and not ptWeightsB:
        hPtPrompt = sparsesReco['RecoPrompt'].Projection(axes['RecoPrompt']['Pt'])
        hPtFD = sparsesReco['RecoFD'].Projection(axes['RecoFD']['Pt'])

    ### pt weights for prompt
    if ptWeights:
        hPtPrompt = sparsesReco['RecoPrompt'].Projection(axes['RecoPrompt']['Pt'])
        for iBin in range(1, hPtPrompt.GetNbinsX()+1):
            if hPtPrompt.GetBinContent(iBin) > 0.:
                relStatUnc = hPtPrompt.GetBinError(iBin) / hPtPrompt.GetBinContent(iBin)
                ptCent = hPtPrompt.GetBinCenter(iBin)
                hPtPrompt.SetBinContent(iBin, hPtPrompt.GetBinContent(iBin) * sPtWeights(ptCent))
                hPtPrompt.SetBinError(iBin, hPtPrompt.GetBinContent(iBin) * relStatUnc)
        ### initially use prompt weights for FD, eventually overwrite
        hPtFD = sparsesReco['RecoFD'].Projection(axes['RecoFD']['Pt'])
        for iBin in range(1, hPtFD.GetNbinsX()+1):
            if hPtFD.GetBinContent(iBin) > 0.:
                relStatUnc = hPtFD.GetBinError(iBin) / hPtFD.GetBinContent(iBin)
                ptCent = hPtFD.GetBinCenter(iBin)
                hPtFD.SetBinContent(iBin, hPtFD.GetBinContent(iBin) * sPtWeights(ptCent))
                hPtFD.SetBinError(iBin, hPtFD.GetBinContent(iBin) * relStatUnc)
    ### pt weights for non-prompt, but no B species weights
    if ptWeightsB and not Bspeciesweights:
        hPtBvsPtD = sparsesReco['RecoFD'].Projection(axes['RecoFD']['Pt'], axes['RecoFD']['pt_bmoth'])
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
        hPtFD = hPtBvsPtD.ProjectionX(f'hFDPt', 0, hPtBvsPtD.GetYaxis().GetNbins()+1, 'e')
    ### pt weights from B for non-prompt and B species weights
    elif ptWeightsB and Bspeciesweights:
        hPtBvsBspecievsPtD = sparsesReco['RecoFD'].Projection(axes['RecoFD']['Pt'], axes['RecoFD']['flag_bhad'], axes['RecoFD']['pt_bmoth'])
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
        hPtFD = hPtBvsBspecievsPtD.ProjectionX(f'hFDPt', 0, hPtBvsBspecievsPtD.GetYaxis().GetNbins()+1,
                                                0, hPtBvsBspecievsPtD.GetZaxis().GetNbins()+1, 'e')
    ### only B species weights
    elif Bspeciesweights:
        hBspecievsPtD = sparsesReco['RecoFD'].Projection(axes['RecoFD']['Pt'], axes['RecoFD']['flag_bhad'])
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
        hPtFD = hBspecievsPtD.ProjectionX(f'hFDPt', 0, hBspecievsPtD.GetYaxis().GetNbins()+1, 'e')

    ## write the output      
    hMassPrompt.Write(f'hPromptMass')
    hMassFD.Write(f'hFDMass')
    hPtPrompt.Write(f'hPromptPt')
    hPtFD.Write(f'hFDPt')

    if 'RecoAll' in sparsesReco:
        hMass.Write(f'hRecoAllMass')
        hPt.Write(f'RecoAllPt')
    if config.get('enableRef'):
        hMassRefl = sparsesReco['RecoRefl'].Projection(axes['RecoRefl']['Mass'])
        hMassRefl.Write(f'hReflMass')
        hMassReflPrompt = sparsesReco['RecoReflPrompt'].Projection(axes['RecoReflPrompt']['Mass'])
        hMassReflPrompt.Write(f'hReflPromptMass')
        hMassReflFD = sparsesReco['RecoReflFD'].Projection(axes['RecoReflFD']['Mass'])
        hMassReflFD.Write(f'hReflFDMass')
        hPtRefl = sparsesReco['RecoRefl'].Projection(axes['RecoRefl']['Pt'])
        hPtRefl.Write(f'hReflPt')
        hPtReflPrompt = sparsesReco['RecoReflPrompt'].Projection(axes['RecoReflPrompt']['Pt'])
        hPtReflPrompt.Write(f'hReflPromptPt')
        hPtReflFD = sparsesReco['RecoReflFD'].Projection(axes['RecoReflFD']['Pt'])
        hPtReflFD.Write(f'hReflFDPt')
    if config.get('enableSecPeak'):
        hMassPromptSecPeak = sparsesReco['RecoSecPeakPrompt'].Projection(axes['RecoSecPeakPrompt']['Mass'])
        hMassPromptSecPeak.Write(f'hPromptSecPeakMass')
        hMassFDSecPeak = sparsesReco['RecoSecPeakFD'].Projection(axes['RecoSecPeakFD']['Mass'])
        hMassFDSecPeak.Write(f'hFDSecPeakMass')
        hPtPromptSecPeak = sparsesReco['RecoSecPeakPrompt'].Projection(axes['RecoSecPeakPrompt']['Pt'])
        hPtPromptSecPeak.Write(f'hPromptSecPeakPt')
        hPtFDSecPeak = sparsesReco['RecoSecPeakFD'].Projection(axes['RecoSecPeakFD']['Pt'])
        hPtFDSecPeak.Write(f'hFDSecPeakPt')

def proj_mc_gen(config, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB):

    if config.get('enableSecPeak'):
        hGenPtPromptSecPeak = sparsesGen['GenSecPeakPrompt'].Projection(axes['GenSecPeakPrompt']['Pt'])
        hGenPtPromptSecPeak.Write(f'hPromptSecPeakGenPt')
        hGenPtFDSecPeak = sparsesGen['GenSecPeakFD'].Projection(axes['GenSecPeakFD']['Pt'])
        hGenPtFDSecPeak.Write(f'hFDSecPeakGenPt')

    ## no pt weights
    # if not ptWeights and not ptWeightsB:
    #     print('NO PT WEIGHTS FOR GENERATED')
    #     hGenPtPrompt = sparsesGen['GenPrompt'].Projection(axes['GenPrompt']['Pt'])
    #     # hGenPtFD = sparsesGen['GenFD'].Projection(axes['GenFD']['Pt'])

    ## apply pt weights for prompt
    if ptWeights:
        hGenPtPrompt = sparsesGen['GenPrompt'].Projection(axes['GenPrompt']['Pt'])
        for iBin in range(1, hGenPtPrompt.GetNbinsX()+1):
            if hGenPtPrompt.GetBinContent(iBin) > 0:
                relStatUnc = hGenPtPrompt.GetBinError(iBin) / hGenPtPrompt.GetBinContent(iBin)
                ptCent = hGenPtPrompt.GetBinCenter(iBin)
                hGenPtPrompt.SetBinContent(iBin, hGenPtPrompt.GetBinContent(iBin) * sPtWeights(ptCent))
                hGenPtPrompt.SetBinError(iBin, hGenPtPrompt.GetBinContent(iBin) * relStatUnc)
    else:
        hGenPtPrompt = sparsesGen['GenPrompt'].Projection(axes['GenPrompt']['Pt'])

    ## apply pt weights for non-prompt
    ### pt weights from B for non-prompt, but no B species weights
    if ptWeightsB and not Bspeciesweights:
        hPtBvsPtGenD = sparsesGen['GenFD'].Projection(axes['GenFD']['Pt'], axes['GenFD']['pt_bmoth'])
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
        hGenPtFD = hPtBvsPtGenD.ProjectionX(f'hFDGenPt', 0, hPtBvsPtGenD.GetYaxis().GetNbins()+1, 'e')
    ### pt weights from B for non-prompt and B species weights
    elif ptWeightsB and Bspeciesweights:
        hPtBvsBspecievsPtGenD = sparsesGen['GenFD'].Projection(axes['GenFD']['Pt'], axes['GenFD']['flag_bhad'], axes['GenFD']['pt_bmoth'])
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
        hGenPtFD = hPtBvsBspecievsPtGenD.ProjectionX(f'hFDGenPt', 0, hPtBvsBspecievsPtGenD.GetYaxis().GetNbins()+1,
                                                     0, hPtBvsBspecievsPtGenD.GetZaxis().GetNbins()+1, 'e')
    ### only B species weights
    elif Bspeciesweights:
        hBspecievsPtGenD = sparsesGen['GenFD'].Projection(axes['GenFD']['Pt'], axes['GenFD']['flag_bhad'])
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
        hGenPtFD = hBspecievsPtGenD.ProjectionX(f'hFDGenPt', 0, hBspecievsPtGenD.GetYaxis().GetNbins()+1, 'e')
    ### no pt weights for generated FD
    else:
        hGenPtFD = sparsesGen['GenFD'].Projection(axes['GenFD']['Pt'])
                            
    ## write the output
    hGenPtPrompt.Write(f'hPromptGenPt')
    hGenPtFD.Write(f'hFDGenPt')

def pt_weights_info(ptweights, ptweightsB):
    # compute info for pt weights
    if ptweights != []:
        ptWeights = uproot.open(ptweights[0])[ptweights[1]]
        bins = ptWeights.axis(0).edges()
        ptCentW = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
        sPtWeights = InterpolatedUnivariateSpline(ptCentW, ptWeights.values())
    else:
        print('\033[91m WARNING: pt weights will not be provided! \033[0m')
        ptWeights = None
        sPtWeights = None

    if ptweightsB != []:
        ptWeightsB = uproot.open(ptweightsB[0])[ptweightsB[1]]
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
    parser.add_argument('--preprocessed', action='store_true', 
                        help='Determines whether the sparses are pre-processed')
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
    
    cent, (cent_min, cent_max) = get_centrality_bins(args.centrality)
    outfile_dir = 'hf-candidate-creator-2prong' if config['Dmeson'] == 'Dzero' else 'hf-candidate-creator-3prong'
    infilemc = TFile.Open(config['MC_filename'], 'r')
    histo_cent = infilemc.Get(f'{outfile_dir}/hSelCollisionsCent')
    histo_cent.GetXaxis().SetRangeUser(cent_min, cent_max)
    resofile = TFile.Open(args.resolution, 'r')
    try:
        det_A = config['detA']
        det_B = config['detB']
        det_C = config['detC']
        histo_reso = resofile.Get(f'{det_A}_{det_B}_{det_C}/histo_reso_delta_cent')
        histo_reso.SetName('histo_reso_delta_cent')
        histo_reso.SetDirectory(0)
        reso = histo_reso.GetBinContent(1)
    except:
        histo_reso = resofile.Get(f'hf-task-flow-charm-hadrons/spReso/hSpReso{det_B}{det_C}')
        histo_reso.SetName('histo_reso_delta_cent')
        histo_reso.SetDirectory(0)
        reso = histo_reso.GetBinContent(1)
    
    outfile.cd()
    histo_reso.Write()
    outfile.mkdir(outfile_dir)
    outfile.cd(outfile_dir)
    histo_cent.Write()
    resofile.Close()
    infilemc.Close()

    # load thnsparse
    sparsesFlow, sparsesReco, sparsesGen, axes = get_sparses(config, True, True, True, args.preprocessed, f'{config["skimDir"]}')
    if not args.preprocessed:
        for key, iSparse in sparsesFlow.items():
            iSparse.GetAxis(axes['Flow']['cent']).SetRangeUser(cent_min, cent_max)
    for key, iSparse in sparsesGen.items():
        iSparse.GetAxis(axes[key]['cent']).SetRangeUser(cent_min, cent_max)
    for key, iSparse in sparsesReco.items():
        iSparse.GetAxis(axes[key]['cent']).SetRangeUser(cent_min, cent_max)

    with open(args.cutsetConfig, 'r') as ymlCutSetFile:
        cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
    cutVars = cutSetCfg['cutvars']

    # compute info for pt weights
    ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB = pt_weights_info(args.ptweights, args.ptweightsB)

    with alive_bar(len(cutVars['Pt']['min']), title='Processing pT bins') as bar:
        for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
            print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
            ptLowLabel = ptMin * 10
            ptHighLabel = ptMax * 10
            outfile.mkdir(f'cent_bins{cent}/pt_bins{int(ptLowLabel)}_{int(ptHighLabel)}')
            outfile.cd(f'cent_bins{cent}/pt_bins{int(ptLowLabel)}_{int(ptHighLabel)}')
    
            print(f"sparsesFlow: {sparsesFlow}")
            if args.preprocessed:
                print('PREPROCESSED')
                sparsesFlow[f"Flow_{ptLowLabel}_{ptHighLabel}"].GetAxis(axes['Flow']['score_FD']).SetRangeUser(cutVars['score_FD']['min'][iPt], cutVars['score_FD']['max'][iPt])
                proj_data(sparsesFlow[f"Flow_{ptLowLabel}_{ptHighLabel}"], ptMin, ptMax, axes, config['inv_mass_bins'][iPt], reso)
                file_proj = TFile(f'{args.outputdir}/proj/ProjFlow_{args.suffix}.root', 'recreate')
                for idim in range(sparsesFlow[f"Flow_{ptLowLabel}_{ptHighLabel}"].GetNdimensions()):
                    histo = sparsesFlow[f"Flow_{ptLowLabel}_{ptHighLabel}"].Projection(idim)
                    histo.Write()
                file_proj.Close()
                outfile.cd(f'cent_bins{cent}/pt_bins{int(ptLowLabel)}_{int(ptHighLabel)}')
                print(f"Projected data!")
            
            if not args.preprocessed:
                print('NOT PREPROCESSED')
                for iSparse, (key, sparse) in enumerate(sparsesFlow.items()):
                    for iVar in cutVars:
                        sparse.GetAxis(axes['Flow'][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                proj_data(sparsesFlow, ptMin, ptMax, axes, config['inv_mass_bins'][iPt], reso)
                print(f"Projected data!")
            
            for iVar in cutVars:
                for key, iSparse in sparsesReco.items():
                    iSparse.GetAxis(axes[key][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                if iVar == 'Pt':
                    for key, iSparse in sparsesGen.items():
                        iSparse.GetAxis(axes[key][iVar]).SetRangeUser(cutVars[iVar]['min'][iPt], cutVars[iVar]['max'][iPt])
                if iVar == 'score_FD' or iVar == 'score_bkg':
                    print(f'{iVar}: {cutVars[iVar]["min"][iPt]} < {iVar} < {cutVars[iVar]["max"][iPt]}')

            proj_mc_reco(config, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB)
            print(f"Projected mc reco!")
            proj_mc_gen(config, ptWeights, ptWeightsB, Bspeciesweights, sPtWeights, sPtWeightsB)
            print(f"Projected mc gen!")
            
            bar()
    
    outfile.Close()
    outfile.Close()
