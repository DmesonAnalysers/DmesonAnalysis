'''
python script for the projection of D+ and Ds+ mesons THnSparses
run: python ProjectDplusDsSparse.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import sys
import argparse
import yaml
from ROOT import TFile, TH1F, TList  # pylint: disable=import-error,no-name-in-module
from TaskFileLoader import LoadSparseFromTask, LoadNormObjFromTask


def merge_hist(first_hist, second_hist):
    h_merged = first_hist.Clone()
    merge_list = TList()
    merge_list.Add(second_hist)
    h_merged.Merge(merge_list)
    return h_merged


def main():  # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    PARSER = argparse.ArgumentParser(description='Arguments to pass')
    PARSER.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                        help='config file name with root input files')
    PARSER.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                        help='input file with cut set')
    PARSER.add_argument('outFileName', metavar='text', default='outFileName.root',
                        help='output root file name')
    ARGS = PARSER.parse_args()

    with open(ARGS.cfgFileName, 'r') as ymlCfgFile:
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
            sparseRecoPart, sparseGenPart = LoadSparseFromTask(
                infilename, inputCfg)
            hEvPart, normCounterPart = LoadNormObjFromTask(
                infilename, inputCfg)
            for sparsetype in sparseRecoPart:
                sparseReco[sparsetype].Add(sparseRecoPart[sparsetype])
            for sparsetype in sparseGenPart:
                sparseGen[sparsetype].Add(sparseGenPart[sparsetype])
            hEv.Add(hEvPart)
            normCounterPart.Add(normCounterPart)

    with open(ARGS.cutSetFileName, 'r') as ymlCutSetFile:
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

    outfile = TFile(ARGS.outFileName, 'RECREATE')

    for iPt in range(0, len(cutVars['Pt']['min'])):
        print("Projecting distributions for %0.1f < pT < %0.1f GeV/c" % (cutVars['Pt']['min'][iPt],
                                                                         cutVars['Pt']['max'][iPt]))
        for iVar in cutVars:
            if iVar == 'InvMass':
                continue
            binMin = sparseReco['RecoAll'].GetAxis(
                cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min'][iPt]*1.0001)
            binMax = sparseReco['RecoAll'].GetAxis(
                cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max'][iPt]*0.9999)
            sparseReco['RecoAll'].GetAxis(
                cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
            if isMC:
                sparseReco['RecoPrompt'].GetAxis(
                    cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
                sparseReco['RecoFD'].GetAxis(
                    cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
                if enableSecPeak:
                    sparseReco['RecoSecPeakPrompt'].GetAxis(
                        cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
                    sparseReco['RecoSecPeakFD'].GetAxis(
                        cutVars[iVar]['axisnum']).SetRange(binMin, binMax)

        for iVar in ('InvMass', 'Pt'):
            hVar = sparseReco['RecoAll'].Projection(cutVars[iVar]['axisnum'])
            hVar.SetName('h%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                            cutVars['Pt']['max'][iPt]*10))
            outfile.cd()
            all_dict[iVar].append(hVar)
            hVar.Write()
            if isMC:
                hVarPrompt = sparseReco['RecoPrompt'].Projection(
                    cutVars[iVar]['axisnum'])
                hVarPrompt.SetName('hPrompt%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                            cutVars['Pt']['max'][iPt]*10))
                prompt_dict[iVar].append(hVarPrompt)
                hVarPrompt.Write()
                hVarFD = sparseReco['RecoFD'].Projection(
                    cutVars[iVar]['axisnum'])
                hVarFD.SetName('hFD%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                    cutVars['Pt']['max'][iPt]*10))
                fd_dict[iVar].append(hVarFD)
                hVarFD.Write()
                if enableSecPeak:
                    hVarPromptSecPeak = sparseReco['RecoSecPeakPrompt'].Projection(
                        cutVars[iVar]['axisnum'])
                    hVarPromptSecPeak.SetName('hPromptSecPeak%s_%0.f_%0.f' % (cutVars[iVar]['name'],
                                                                              cutVars['Pt']['min'][iPt]*10,
                                                                              cutVars['Pt']['max'][iPt]*10))
                    prompt_dict_secpeak[iVar].append(hVarPromptSecPeak)
                    hVarPromptSecPeak.Write()
                    hVarFDSecPeak = sparseReco['RecoSecPeakFD'].Projection(
                        cutVars[iVar]['axisnum'])
                    hVarFDSecPeak.SetName('hFDSecPeak%s_%0.f_%0.f' % (cutVars[iVar]['name'],
                                                                      cutVars['Pt']['min'][iPt]*10,
                                                                      cutVars['Pt']['max'][iPt]*10))
                    fd_dict_secpeak[iVar].append(hVarFDSecPeak)
                    hVarFDSecPeak.Write()
        if isMC:
            binGenMin = sparseGen['GenPrompt'].GetAxis(
                0).FindBin(cutVars['Pt']['min'][iPt]*1.0001)
            binGenMax = sparseGen['GenPrompt'].GetAxis(
                0).FindBin(cutVars['Pt']['max'][iPt]*0.9999)
            sparseGen['GenPrompt'].GetAxis(0).SetRange(binGenMin, binGenMax)
            sparseGen['GenFD'].GetAxis(0).SetRange(binGenMin, binGenMax)
            hGenPtPrompt = sparseGen['GenPrompt'].Projection(0)
            hGenPtPrompt.SetName('hPromptGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                             cutVars['Pt']['max'][iPt]*10))
            prompt_gen_list.append(hGenPtPrompt)
            hGenPtPrompt.Write()
            hGenPtFD = sparseGen['GenFD'].Projection(0)
            hGenPtFD.SetName('hFDGenPt_%0.f_%0.f' % (
                cutVars['Pt']['min'][iPt]*10, cutVars['Pt']['max'][iPt]*10))
            fd_gen_list.append(hGenPtFD)
            hGenPtFD.Write()
            if enableSecPeak:
                sparseGen['GenSecPeakPrompt'].GetAxis(
                    0).SetRange(binGenMin, binGenMax)
                sparseGen['GenSecPeakFD'].GetAxis(
                    0).SetRange(binGenMin, binGenMax)
                hGenPtPromptSecPeak = sparseGen['GenSecPeakPrompt'].Projection(
                    0)
                hGenPtPromptSecPeak.SetName('hPromptSecPeakGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                               cutVars['Pt']['max'][iPt]*10))
                prompt_gen_list_secpeak.append(hGenPtPromptSecPeak)
                hGenPtPromptSecPeak.Write()
                hGenPtFDSecPeak = sparseGen['GenSecPeakFD'].Projection(0)
                hGenPtFDSecPeak.SetName('hFDSecPeakGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                       cutVars['Pt']['max'][iPt]*10))
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
        for iVar in ('InvMass', 'Pt'):
            hVar_merged = merge_hist(
                all_dict[iVar][iPt], all_dict[iVar][iPt+1])
            hVar_merged.SetName('h%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                   cutVars['Pt']['max'][iPt+1]*10))
            hVar_merged.Write()
            if isMC:
                hVarPrompt_merged = merge_hist(
                    prompt_dict[iVar][iPt], prompt_dict[iVar][iPt+1])
                hVarPrompt_merged.SetName('hPrompt%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                                   cutVars['Pt']['max'][iPt+1]*10))
                hVarPrompt_merged.Write()
                hVarFD_merged = merge_hist(
                    fd_dict[iVar][iPt], fd_dict[iVar][iPt+1])
                hVarFD_merged.SetName('hFD%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                           cutVars['Pt']['max'][iPt+1]*10))
                hVarFD_merged.Write()
                if enableSecPeak:
                    hVarPrompt_secpeak_merged = merge_hist(prompt_dict_secpeak[iVar][iPt],
                                                           prompt_dict_secpeak[iVar][iPt+1])
                    hVarPrompt_secpeak_merged.SetName('hPromptSecPeak%s_%0.f_%0.f' % (cutVars[iVar]['name'],
                                                                                      cutVars['Pt']['min'][iPt]*10,
                                                                                      cutVars['Pt']['max'][iPt+1]*10))
                    hVarPrompt_secpeak_merged.Write()
                    hVarFD_secpeak_merged = merge_hist(
                        fd_dict_secpeak[iVar][iPt], fd_dict_secpeak[iVar][iPt+1])
                    hVarFD_secpeak_merged.SetName('hFDSecPeak%s_%0.f_%0.f' % (cutVars[iVar]['name'],
                                                                              cutVars['Pt']['min'][iPt]*10,
                                                                              cutVars['Pt']['max'][iPt+1]*10))
                    hVarFD_secpeak_merged.Write()

        if isMC:
            hVarPromptGen_merged = merge_hist(
                prompt_gen_list[iPt], prompt_gen_list[iPt+1])
            hVarPromptGen_merged.SetName('hPromptGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                     cutVars['Pt']['max'][iPt+1]*10))
            hVarPromptGen_merged.Write()
            hVarPromptFD_merged = merge_hist(
                fd_gen_list[iPt], fd_gen_list[iPt+1])
            hVarPromptFD_merged.SetName('hFDGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                cutVars['Pt']['max'][iPt+1]*10))
            hVarPromptFD_merged.Write()
            if enableSecPeak:
                hVarPromptGen_secpeak_merged = merge_hist(
                    prompt_gen_list_secpeak[iPt], prompt_gen_list_secpeak[iPt+1])
                hVarPromptGen_secpeak_merged.SetName('hPromptGenSecPeakPt_%0.f_%0.f' %
                                                     (cutVars['Pt']['min'][iPt]*10, cutVars['Pt']['max'][iPt+1]*10))
                hVarPromptGen_secpeak_merged.Write()
                hVarFDGen_secpeak_merged = merge_hist(
                    fd_gen_list_secpeak[iPt], fd_gen_list_secpeak[iPt+1])
                hVarFDGen_secpeak_merged.SetName('hFDGenSecPeakPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                                cutVars['Pt']['max'][iPt+1]*10))
                hVarFDGen_secpeak_merged.Write()

    hEvForNorm = TH1F("hEvForNorm", ";;Number of events", 2, 0., 2.)
    hEvForNorm.GetXaxis().SetBinLabel(1, "norm counter")
    hEvForNorm.GetXaxis().SetBinLabel(2, "accepted events")
    hEvForNorm.SetBinContent(1, normCounter.GetNEventsForNorm())
    hEvForNorm.SetBinContent(2, hEv.GetBinContent(5))

    hEvForNorm.Write()
    outfile.Close()


main()
