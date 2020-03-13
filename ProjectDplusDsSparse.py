'''
python script for the projection of D+ and Ds+ mesons THnSparses
run: python ProjectDplusDsSparse.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import argparse
import yaml
from ROOT import TFile, TH1F  # pylint: disable=import-error,no-name-in-module
from utils.TaskFileLoader import LoadSparseFromTask, LoadNormObjFromTask
from utils.AnalysisUtils import MergeHists

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                    help='input file with cut set')
parser.add_argument('outFileName', metavar='text', default='outFileName.root',
                    help='output root file name')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

infilenames = inputCfg['filename']
if not isinstance(infilenames, list):
    infilenames = [infilenames]
isMC = inputCfg['isMC']
enableSecPeak = inputCfg['enableSecPeak']

for iFile, infilename in enumerate(infilenames):
    if iFile == 0:
        sparseReco, sparseGen = LoadSparseFromTask(infilename, inputCfg)
        hEv, normCounter = LoadNormObjFromTask(infilename, inputCfg)
    else:
        sparseRecoPart, sparseGenPart = LoadSparseFromTask(infilename, inputCfg)
        hEvPart, normCounterPart = LoadNormObjFromTask(infilename, inputCfg)
        for sparsetype in sparseRecoPart:
            sparseReco[sparsetype].Add(sparseRecoPart[sparsetype])
        for sparsetype in sparseGenPart:
            sparseGen[sparsetype].Add(sparseGenPart[sparsetype])
        hEv.Add(hEvPart)
        normCounter.Add(normCounterPart)

with open(args.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']

# dicts of TH1
all_dict = {'InvMass': [], 'Pt': []}
prompt_dict = {'InvMass': [], 'Pt': []}
fd_dict = {'InvMass': [], 'Pt': []}
prompt_gen_list = []
fd_gen_list = []
prompt_dict_secpeak = {'InvMass': [], 'Pt': []}
fd_dict_secpeak = {'InvMass': [], 'Pt': []}
prompt_gen_list_secpeak = []
fd_gen_list_secpeak = []

outfile = TFile(args.outFileName, 'recreate')

for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
    print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
    ptLowLabel = ptMin * 10
    ptHighLabel = ptMax * 10
    for iVar in cutVars:
        if iVar == 'InvMass':
            continue
        binMin = sparseReco['RecoAll'].GetAxis(
            cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min'][iPt] * 1.0001)
        binMax = sparseReco['RecoAll'].GetAxis(
            cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max'][iPt] * 0.9999)
        sparseReco['RecoAll'].GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
        if isMC:
            sparseReco['RecoPrompt'].GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
            sparseReco['RecoFD'].GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
            if enableSecPeak:
                sparseReco['RecoSecPeakPrompt'].GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
                sparseReco['RecoSecPeakFD'].GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)

    for iVar in ('InvMass', 'Pt'):
        varName = cutVars[iVar]['name']
        hVar = sparseReco['RecoAll'].Projection(cutVars[iVar]['axisnum'])
        hVar.SetName(f'h{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        outfile.cd()
        all_dict[iVar].append(hVar)
        hVar.Write()
        if isMC:
            hVarPrompt = sparseReco['RecoPrompt'].Projection(cutVars[iVar]['axisnum'])
            hVarPrompt.SetName(f'hPrompt%{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            prompt_dict[iVar].append(hVarPrompt)
            hVarPrompt.Write()
            hVarFD = sparseReco['RecoFD'].Projection(cutVars[iVar]['axisnum'])
            hVarFD.SetName(f'hFD{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            fd_dict[iVar].append(hVarFD)
            hVarFD.Write()
            if enableSecPeak:
                hVarPromptSecPeak = sparseReco['RecoSecPeakPrompt'].Projection(cutVars[iVar]['axisnum'])
                hVarPromptSecPeak.SetName(f'hPromptSecPeak{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
                prompt_dict_secpeak[iVar].append(hVarPromptSecPeak)
                hVarPromptSecPeak.Write()
                hVarFDSecPeak = sparseReco['RecoSecPeakFD'].Projection(cutVars[iVar]['axisnum'])
                hVarFDSecPeak.SetName(f'hFDSecPeak{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
                fd_dict_secpeak[iVar].append(hVarFDSecPeak)
                hVarFDSecPeak.Write()
    if isMC:
        binGenMin = sparseGen['GenPrompt'].GetAxis(0).FindBin(ptMin*1.0001)
        binGenMax = sparseGen['GenPrompt'].GetAxis(0).FindBin(ptMax*0.9999)
        sparseGen['GenPrompt'].GetAxis(0).SetRange(binGenMin, binGenMax)
        sparseGen['GenFD'].GetAxis(0).SetRange(binGenMin, binGenMax)
        hGenPtPrompt = sparseGen['GenPrompt'].Projection(0)
        hGenPtPrompt.SetName(f'hPromptGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        prompt_gen_list.append(hGenPtPrompt)
        hGenPtPrompt.Write()
        hGenPtFD = sparseGen['GenFD'].Projection(0)
        hGenPtFD.SetName(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        fd_gen_list.append(hGenPtFD)
        hGenPtFD.Write()
        if enableSecPeak:
            sparseGen['GenSecPeakPrompt'].GetAxis(0).SetRange(binGenMin, binGenMax)
            sparseGen['GenSecPeakFD'].GetAxis(0).SetRange(binGenMin, binGenMax)
            hGenPtPromptSecPeak = sparseGen['GenSecPeakPrompt'].Projection(0)
            hGenPtPromptSecPeak.SetName(f'hPromptSecPeakGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            prompt_gen_list_secpeak.append(hGenPtPromptSecPeak)
            hGenPtPromptSecPeak.Write()
            hGenPtFDSecPeak = sparseGen['GenSecPeakFD'].Projection(0)
            hGenPtFDSecPeak.SetName(f'hFDSecPeakGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            fd_gen_list_secpeak.append(hGenPtFDSecPeak)
            hGenPtFDSecPeak.Write()

    for iVar in cutVars:
        sparseReco['RecoAll'].GetAxis(
            cutVars[iVar]['axisnum']).SetRange(-1, -1)
        if isMC:
            sparseReco['RecoPrompt'].GetAxis(
                cutVars[iVar]['axisnum']).SetRange(-1, -1)
            sparseReco['RecoFD'].GetAxis(
                cutVars[iVar]['axisnum']).SetRange(-1, -1)
            if enableSecPeak:
                sparseReco['RecoSecPeakPrompt'].GetAxis(
                    cutVars[iVar]['axisnum']).SetRange(-1, -1)
                sparseReco['RecoSecPeakFD'].GetAxis(
                    cutVars[iVar]['axisnum']).SetRange(-1, -1)

for iPt in range(0, len(cutVars['Pt']['min']) - 1):
    ptLowLabel = cutVars['Pt']['min'][iPt] * 10
    ptHighLabel = cutVars['Pt']['max'][iPt+1] * 10
    for iVar in ('InvMass', 'Pt'):
        varName = cutVars[iVar]['name']
        hVar_merged = MergeHists([all_dict[iVar][iPt], all_dict[iVar][iPt+1]])
        hVar_merged.SetName(f'h{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVar_merged.Write()
        if isMC:
            hVarPrompt_merged = MergeHists([prompt_dict[iVar][iPt], prompt_dict[iVar][iPt+1]])
            hVarPrompt_merged.SetName(f'hPrompt{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarPrompt_merged.Write()
            hVarFD_merged = MergeHists([fd_dict[iVar][iPt], fd_dict[iVar][iPt+1]])
            hVarFD_merged.SetName(f'hFD{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarFD_merged.Write()
            if enableSecPeak:
                hVarPrompt_secpeak_merged = MergeHists([prompt_dict_secpeak[iVar][iPt],
                                                        prompt_dict_secpeak[iVar][iPt+1]])
                hVarPrompt_secpeak_merged.SetName(f'hPromptSecPeak{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
                hVarPrompt_secpeak_merged.Write()
                hVarFD_secpeak_merged = MergeHists([fd_dict_secpeak[iVar][iPt], fd_dict_secpeak[iVar][iPt+1]])
                hVarFD_secpeak_merged.SetName(f'hFDSecPeak{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
                hVarFD_secpeak_merged.Write()

    if isMC:
        hVarPromptGen_merged = MergeHists([prompt_gen_list[iPt], prompt_gen_list[iPt+1]])
        hVarPromptGen_merged.SetName(f'hPromptGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarPromptGen_merged.Write()
        hVarPromptFD_merged = MergeHists([fd_gen_list[iPt], fd_gen_list[iPt+1]])
        hVarPromptFD_merged.SetName(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarPromptFD_merged.Write()
        if enableSecPeak:
            hVarPromptGen_secpeak_merged = MergeHists([prompt_gen_list_secpeak[iPt], prompt_gen_list_secpeak[iPt+1]])
            hVarPromptGen_secpeak_merged.SetName(f'hPromptGenSecPeakPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarPromptGen_secpeak_merged.Write()
            hVarFDGen_secpeak_merged = MergeHists([fd_gen_list_secpeak[iPt], fd_gen_list_secpeak[iPt+1]])
            hVarFDGen_secpeak_merged.SetName(f'hFDGenSecPeakPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarFDGen_secpeak_merged.Write()

hEvForNorm = TH1F("hEvForNorm", ";;Number of events", 2, 0., 2.)
hEvForNorm.GetXaxis().SetBinLabel(1, "norm counter")
hEvForNorm.GetXaxis().SetBinLabel(2, "accepted events")
hEvForNorm.SetBinContent(1, normCounter.GetNEventsForNorm())
for iBin in range(1, hEv.GetNbinsX() + 1):
    binLabel = hEv.GetXaxis().GetBinLabel(iBin)
    if 'isEvSelected' in binLabel or 'accepted' in binLabel:
        hEvForNorm.SetBinContent(2, hEv.GetBinContent(iBin))
        break

hEvForNorm.Write()
outfile.Close()
