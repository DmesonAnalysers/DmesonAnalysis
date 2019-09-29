#*******************************************************************************************#
# python script for the projection of D+ and Ds+ mesons THnSparses                          #
# run: python ProjectDplusDsSparse.py cfgFileName.yml cutSetFileName.yml outFileName.root   #
# author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino                            #
#*******************************************************************************************#
import sys
import yaml
from ROOT import TFile, TH1F, TList # pylint: disable=import-error,no-name-in-module

def merge_hist(first_hist, second_hist):
    h_merged = first_hist.Clone()
    merge_list = TList()
    merge_list.Add(second_hist)
    h_merged.Merge(merge_list)
    return h_merged

def main(): # pylint: disable=too-many-locals,too-many-branches,too-many-statements
    cfgFileName = sys.argv[1]
    cutSetFileName = sys.argv[2]
    outFileName = sys.argv[3]

    with open(cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    isMC = inputCfg['isMC']
    infilenames = inputCfg['filename']
    if not isinstance(infilenames, list):
        infilenames = [infilenames]

    for iFile, infilename in enumerate(infilenames):
        infileData = TFile(infilename)
        indirData = infileData.Get(inputCfg['dirname'])
        inlistData = indirData.Get(inputCfg['listname'])
        enableSecPeak = inputCfg['enableSecPeak']
        if iFile == 0:
            sMassPtCutVars = inlistData.FindObject(inputCfg['sparsenameAll'])
            normCounter = indirData.Get(inputCfg['normname'])
            hEv = inlistData.FindObject(inputCfg['histoevname'])
        else:
            sMassPtCutVars.Add(inlistData.FindObject(inputCfg['sparsenameAll']))
            normCounter.Add(indirData.Get(inputCfg['normname']))
            hEv.Add(inlistData.FindObject(inputCfg['histoevname']))
        if isMC:
            if iFile == 0:
                sMassPtCutVarsPrompt = inlistData.FindObject(inputCfg['sparsenamePrompt'])
                sMassPtCutVarsFD = inlistData.FindObject(inputCfg['sparsenameFD'])
                sGenPrompt = inlistData.FindObject(inputCfg['sparsenameGenPrompt'])
                sGenFD = inlistData.FindObject(inputCfg['sparsenameGenFD'])
                if enableSecPeak:
                    sMassPtCutVarsPromptSecPeak = inlistData.FindObject(inputCfg['sparsenamePromptSecPeak'])
                    sMassPtCutVarsFDSecPeak = inlistData.FindObject(inputCfg['sparsenameFDSecPeak'])
                    sGenPromptSecPeak = inlistData.FindObject(inputCfg['sparsenameGenPromptSecPeak'])
                    sGenFDSecPeak = inlistData.FindObject(inputCfg['sparsenameGenFDSecPeak'])
            else:
                sMassPtCutVarsPrompt.Add(inlistData.FindObject(inputCfg['sparsenamePrompt']))
                sMassPtCutVarsFD.Add(inlistData.FindObject(inputCfg['sparsenameFD']))
                sGenPrompt.Add(inlistData.FindObject(inputCfg['sparsenameGenPrompt']))
                sGenFD.Add(inlistData.FindObject(inputCfg['sparsenameGenFD']))
                if enableSecPeak:
                    sMassPtCutVarsPromptSecPeak.Add(inlistData.FindObject(inputCfg['sparsenamePromptSecPeak']))
                    sMassPtCutVarsFDSecPeak.Add(inlistData.FindObject(inputCfg['sparsenameFDSecPeak']))
                    sGenPromptSecPeak.Add(inlistData.FindObject(inputCfg['sparsenameGenPromptSecPeak']))
                    sGenFDSecPeak.Add(inlistData.FindObject(inputCfg['sparsenameGenFDSecPeak']))

    with open(cutSetFileName, 'r') as ymlCutSetFile:
        cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)

    cutVars = cutSetCfg['cutvars']

    #dicts of TH1
    all_dict = {'InvMass': [], 'Pt': []}
    prompt_dict = {'InvMass': [], 'Pt': []}
    fd_dict = {'InvMass': [], 'Pt': []}
    prompt_gen_list = []
    fd_gen_list = []
    prompt_dict_secpeak = {'InvMass': [], 'Pt': []}
    fd_dict_secpeak = {'InvMass': [], 'Pt': []}
    prompt_gen_list_secpeak = []
    fd_gen_list_secpeak = []

    outfile = TFile(outFileName, 'RECREATE')

    for iPt in range(0, len(cutVars['Pt']['min'])):
        print("Projecting distributions for %0.1f < pT < %0.1f GeV/c" % (cutVars['Pt']['min'][iPt],
                                                                         cutVars['Pt']['max'][iPt]))
        for iVar in cutVars:
            if iVar == 'InvMass':
                continue
            binMin = sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['min'][iPt]*1.0001)
            binMax = sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).FindBin(cutVars[iVar]['max'][iPt]*0.9999)
            sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
            if isMC:
                sMassPtCutVarsPrompt.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
                sMassPtCutVarsFD.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
                if enableSecPeak:
                    sMassPtCutVarsPromptSecPeak.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)
                    sMassPtCutVarsFDSecPeak.GetAxis(cutVars[iVar]['axisnum']).SetRange(binMin, binMax)

        for iVar in ('InvMass', 'Pt'):
            hVar = sMassPtCutVars.Projection(cutVars[iVar]['axisnum'])
            hVar.SetName('h%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                            cutVars['Pt']['max'][iPt]*10))
            outfile.cd()
            all_dict[iVar].append(hVar)
            hVar.Write()
            if isMC:
                hVarPrompt = sMassPtCutVarsPrompt.Projection(cutVars[iVar]['axisnum'])
                hVarPrompt.SetName('hPrompt%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                            cutVars['Pt']['max'][iPt]*10))
                prompt_dict[iVar].append(hVarPrompt)
                hVarPrompt.Write()
                hVarFD = sMassPtCutVarsFD.Projection(cutVars[iVar]['axisnum'])
                hVarFD.SetName('hFD%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                    cutVars['Pt']['max'][iPt]*10))
                fd_dict[iVar].append(hVarFD)
                hVarFD.Write()
                if enableSecPeak:
                    hVarPromptSecPeak = sMassPtCutVarsPromptSecPeak.Projection(cutVars[iVar]['axisnum'])
                    hVarPromptSecPeak.SetName('hPromptSecPeak%s_%0.f_%0.f' % (cutVars[iVar]['name'],
                                                                              cutVars['Pt']['min'][iPt]*10,
                                                                              cutVars['Pt']['max'][iPt]*10))
                    prompt_dict_secpeak[iVar].append(hVarPromptSecPeak)
                    hVarPromptSecPeak.Write()
                    hVarFDSecPeak = sMassPtCutVarsFDSecPeak.Projection(cutVars[iVar]['axisnum'])
                    hVarFDSecPeak.SetName('hFDSecPeak%s_%0.f_%0.f' % (cutVars[iVar]['name'],
                                                                      cutVars['Pt']['min'][iPt]*10,
                                                                      cutVars['Pt']['max'][iPt]*10))
                    fd_dict_secpeak[iVar].append(hVarFDSecPeak)
                    hVarFDSecPeak.Write()
        if isMC:
            binGenMin = sGenPrompt.GetAxis(0).FindBin(cutVars['Pt']['min'][iPt]*1.0001)
            binGenMax = sGenPrompt.GetAxis(0).FindBin(cutVars['Pt']['max'][iPt]*0.9999)
            sGenPrompt.GetAxis(0).SetRange(binGenMin, binGenMax)
            sGenFD.GetAxis(0).SetRange(binGenMin, binGenMax)
            hGenPtPrompt = sGenPrompt.Projection(0)
            hGenPtPrompt.SetName('hPromptGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                             cutVars['Pt']['max'][iPt]*10))
            prompt_gen_list.append(hGenPtPrompt)
            hGenPtPrompt.Write()
            hGenPtFD = sGenFD.Projection(0)
            hGenPtFD.SetName('hFDGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10, cutVars['Pt']['max'][iPt]*10))
            fd_gen_list.append(hGenPtFD)
            hGenPtFD.Write()
            if enableSecPeak:
                sGenPromptSecPeak.GetAxis(0).SetRange(binGenMin, binGenMax)
                sGenFDSecPeak.GetAxis(0).SetRange(binGenMin, binGenMax)
                hGenPtPromptSecPeak = sGenPromptSecPeak.Projection(0)
                hGenPtPromptSecPeak.SetName('hPromptSecPeakGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                               cutVars['Pt']['max'][iPt]*10))
                prompt_gen_list_secpeak.append(hGenPtPromptSecPeak)
                hGenPtPromptSecPeak.Write()
                hGenPtFDSecPeak = sGenFDSecPeak.Projection(0)
                hGenPtFDSecPeak.SetName('hFDSecPeakGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                       cutVars['Pt']['max'][iPt]*10))
                fd_gen_list_secpeak.append(hGenPtFDSecPeak)
                hGenPtFDSecPeak.Write()

        for iVar in cutVars:
            sMassPtCutVars.GetAxis(cutVars[iVar]['axisnum']).SetRange(-1, -1)
            if isMC:
                sMassPtCutVarsPrompt.GetAxis(cutVars[iVar]['axisnum']).SetRange(-1, -1)
                sMassPtCutVarsFD.GetAxis(cutVars[iVar]['axisnum']).SetRange(-1, -1)
                if enableSecPeak:
                    sMassPtCutVarsPromptSecPeak.GetAxis(cutVars[iVar]['axisnum']).SetRange(-1, -1)
                    sMassPtCutVarsFDSecPeak.GetAxis(cutVars[iVar]['axisnum']).SetRange(-1, -1)

    for iPt in range(0, len(cutVars['Pt']['min']) - 1):
        for iVar in ('InvMass', 'Pt'):
            hVar_merged = merge_hist(all_dict[iVar][iPt], all_dict[iVar][iPt+1])
            hVar_merged.SetName('h%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                   cutVars['Pt']['max'][iPt+1]*10))
            hVar_merged.Write()
            if isMC:
                hVarPrompt_merged = merge_hist(prompt_dict[iVar][iPt], prompt_dict[iVar][iPt+1])
                hVarPrompt_merged.SetName('hPrompt%s_%0.f_%0.f' % (cutVars[iVar]['name'], cutVars['Pt']['min'][iPt]*10,
                                                                   cutVars['Pt']['max'][iPt+1]*10))
                hVarPrompt_merged.Write()
                hVarFD_merged = merge_hist(fd_dict[iVar][iPt], fd_dict[iVar][iPt+1])
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
                    hVarFD_secpeak_merged = merge_hist(fd_dict_secpeak[iVar][iPt], fd_dict_secpeak[iVar][iPt+1])
                    hVarFD_secpeak_merged.SetName('hFDSecPeak%s_%0.f_%0.f' % (cutVars[iVar]['name'],
                                                                              cutVars['Pt']['min'][iPt]*10,
                                                                              cutVars['Pt']['max'][iPt+1]*10))
                    hVarFD_secpeak_merged.Write()

        if isMC:
            hVarPromptGen_merged = merge_hist(prompt_gen_list[iPt], prompt_gen_list[iPt+1])
            hVarPromptGen_merged.SetName('hPromptGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                     cutVars['Pt']['max'][iPt+1]*10))
            hVarPromptGen_merged.Write()
            hVarPromptFD_merged = merge_hist(fd_gen_list[iPt], fd_gen_list[iPt+1])
            hVarPromptFD_merged.SetName('hFDGenPt_%0.f_%0.f' % (cutVars['Pt']['min'][iPt]*10,
                                                                cutVars['Pt']['max'][iPt+1]*10))
            hVarPromptFD_merged.Write()
            if enableSecPeak:
                hVarPromptGen_secpeak_merged = merge_hist(prompt_gen_list_secpeak[iPt], prompt_gen_list_secpeak[iPt+1])
                hVarPromptGen_secpeak_merged.SetName('hPromptGenSecPeakPt_%0.f_%0.f' %
                                                     (cutVars['Pt']['min'][iPt]*10, cutVars['Pt']['max'][iPt+1]*10))
                hVarPromptGen_secpeak_merged.Write()
                hVarFDGen_secpeak_merged = merge_hist(fd_gen_list_secpeak[iPt], fd_gen_list_secpeak[iPt+1])
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
