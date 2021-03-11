'''
python script for the projection of D+ and Ds+ mesons THnSparses
run: python ProjectDplusDsSparse.py cfgFileName.yml cutSetFileName.yml outFileName.root
                                    [--ptweights PtWeightsFileName.root histoName]
                                    [--ptweightsB PtWeightsFileName.root histoName]

if the --ptweights argument is provided, pT weights will be applied to prompt and FD pT distributions
if the --ptweightsB argument is provided, pT weights will be applied to FD pT distributions instead of
those for the prompt
'''

import argparse
import yaml
import uproot
from scipy.interpolate import InterpolatedUnivariateSpline
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
parser.add_argument('--ptweights', metavar=('text', 'text'), nargs=2, required=False,
                    help='First path of the pT weights file, second name of the pT weights histogram')
parser.add_argument('--ptweightsB', metavar=('text', 'text'), nargs=2, required=False,
                    help='First path of the pT weights file, second name of the pT weights histogram')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

infilenames = inputCfg['filename']
if not isinstance(infilenames, list):
    infilenames = [infilenames]
enableSecPeak = inputCfg['enableSecPeak']
isMC = inputCfg['isMC']
isRedVar = inputCfg['isReducedVariables']
isWithBinfo = inputCfg['isWithBinfo']
shiftForRedVar = inputCfg['shiftForRedVariables']
if not isMC:
    if args.ptweights:
        print('WARNING: pt weights will not be applied since it is not MC')
        args.ptweights = None
    if args.ptweightsB:
        print('WARNING: ptB weights will not be applied since it is not MC')
        args.ptweightsB = None
if isRedVar:
    print(('WARNING: option for reduced number of variables in THnSparse set to true.'
           'If this is not the case, the code will work producing wrong results'))

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

refSparse = 'RecoAll'
if isMC:
    refSparse = 'RecoPrompt'

# compute pt weights
if args.ptweights:
    ptWeights = uproot.open(args.ptweights[0])[args.ptweights[1]]
    bins = ptWeights.axis(0).edges()
    ptCentW = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
    sPtWeights = InterpolatedUnivariateSpline(ptCentW, ptWeights.values())
    if not args.ptweightsB:
        sPtWeightsGenDfromB = sPtWeights
        sPtWeightsRecoDfromB = sPtWeights

if args.ptweightsB:
    ptWeightsB = uproot.open(args.ptweightsB[0])[args.ptweightsB[1]]
    bins = ptWeightsB.axis(0).edges()
    ptCentWB = [(bins[iBin]+bins[iBin+1])/2 for iBin in range(len(bins)-1)]
    sPtWeightsB = InterpolatedUnivariateSpline(ptCentWB, ptWeightsB.values())
    hPtBvsPtGenD = sparseGen['GenFD'].Projection(2, 0).ProfileX()
    if isWithBinfo:
        hPtBvsPtRecoD = sparseReco['RecoFD'].Projection(2, 0).ProfileX()
    averagePtBvsPtGen, averagePtBvsPtReco = [], []
    for iPt in range(1, hPtBvsPtGenD.GetNbinsX()+1):
        averagePtBvsPtGen.append(hPtBvsPtGenD.GetBinContent(iPt))
        if isWithBinfo:
            averagePtBvsPtReco.append(hPtBvsPtRecoD.GetBinContent(iPt))
        else:
            averagePtBvsPtReco.append(hPtBvsPtGenD.GetBinContent(iPt))
    aPtGenWeightsB = [sPtWeightsB(pt) for pt in averagePtBvsPtGen]
    aPtRecoWeightsB = [sPtWeightsB(pt) for pt in averagePtBvsPtReco]
    sPtWeightsGenDfromB = InterpolatedUnivariateSpline(ptCentW, aPtGenWeightsB)
    sPtWeightsRecoDfromB = InterpolatedUnivariateSpline(ptCentW, aPtRecoWeightsB)

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
        axisNum = cutVars[iVar]['axisnum']
        if axisNum >= 2: # check if axis is a cut variable (not inv. mass or pt)
            if isRedVar:
                axisNum -= shiftForRedVar
            if isWithBinfo:
                axisNum += 2
        binMin = sparseReco[refSparse].GetAxis(axisNum).FindBin(cutVars[iVar]['min'][iPt] * 1.0001)
        binMax = sparseReco[refSparse].GetAxis(axisNum).FindBin(cutVars[iVar]['max'][iPt] * 0.9999)
        if 'RecoAll' in sparseReco:
            sparseReco['RecoAll'].GetAxis(axisNum).SetRange(binMin, binMax)
        if isMC:
            sparseReco['RecoPrompt'].GetAxis(axisNum).SetRange(binMin, binMax)
            sparseReco['RecoFD'].GetAxis(axisNum).SetRange(binMin, binMax)
            if enableSecPeak:
                sparseReco['RecoSecPeakPrompt'].GetAxis(axisNum).SetRange(binMin, binMax)
                sparseReco['RecoSecPeakFD'].GetAxis(axisNum).SetRange(binMin, binMax)

    for iVar in ('InvMass', 'Pt'):
        varName = cutVars[iVar]['name']
        axisNum = cutVars[iVar]['axisnum']
        if 'RecoAll' in sparseReco:
            hVar = sparseReco['RecoAll'].Projection(axisNum)
            hVar.SetName(f'h{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            outfile.cd()
            all_dict[iVar].append(hVar)
            hVar.Write()
        if isMC:
            hVarPrompt = sparseReco['RecoPrompt'].Projection(axisNum)
            # apply pt weights
            if iVar == 'Pt' and args.ptweights:
                for iBin in range(1, hVarPrompt.GetNbinsX()+1):
                    if hVarPrompt.GetBinContent(iBin) > 0.:
                        relStatUnc = hVarPrompt.GetBinError(iBin) / hVarPrompt.GetBinContent(iBin)
                        ptCent = hVarPrompt.GetBinCenter(iBin)
                        hVarPrompt.SetBinContent(iBin, hVarPrompt.GetBinContent(iBin) * sPtWeights(ptCent))
                        hVarPrompt.SetBinError(iBin, hVarPrompt.GetBinContent(iBin) * relStatUnc)
            hVarPrompt.SetName(f'hPrompt{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            prompt_dict[iVar].append(hVarPrompt)
            hVarPrompt.Write()
            hVarFD = sparseReco['RecoFD'].Projection(axisNum)
            # apply pt weights
            if iVar == 'Pt' and (args.ptweightsB or args.ptweights):
                for iBin in range(1, hVarFD.GetNbinsX()+1):
                    if hVarFD.GetBinContent(iBin) > 0.:
                        relStatUnc = hVarFD.GetBinError(iBin) / hVarFD.GetBinContent(iBin)
                        ptCent = hVarFD.GetBinCenter(iBin)
                        hVarFD.SetBinContent(iBin, hVarFD.GetBinContent(iBin) * sPtWeightsRecoDfromB(ptCent))
                        hVarFD.SetBinError(iBin, hVarFD.GetBinContent(iBin) * relStatUnc)
            hVarFD.SetName(f'hFD{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            fd_dict[iVar].append(hVarFD)
            hVarFD.Write()
            if enableSecPeak:
                hVarPromptSecPeak = sparseReco['RecoSecPeakPrompt'].Projection(axisNum)
                hVarPromptSecPeak.SetName(f'hPromptSecPeak{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
                prompt_dict_secpeak[iVar].append(hVarPromptSecPeak)
                hVarPromptSecPeak.Write()
                hVarFDSecPeak = sparseReco['RecoSecPeakFD'].Projection(axisNum)
                hVarFDSecPeak.SetName(f'hFDSecPeak{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
                fd_dict_secpeak[iVar].append(hVarFDSecPeak)
                hVarFDSecPeak.Write()
    if isMC:
        binGenMin = sparseGen['GenPrompt'].GetAxis(0).FindBin(ptMin*1.0001)
        binGenMax = sparseGen['GenPrompt'].GetAxis(0).FindBin(ptMax*0.9999)
        sparseGen['GenPrompt'].GetAxis(0).SetRange(binGenMin, binGenMax)
        sparseGen['GenFD'].GetAxis(0).SetRange(binGenMin, binGenMax)
        hGenPtPrompt = sparseGen['GenPrompt'].Projection(0)
        # apply pt weights
        if args.ptweights:
            for iBin in range(1, hGenPtPrompt.GetNbinsX()+1):
                if hGenPtPrompt.GetBinContent(iBin) > 0:
                    relStatUnc = hGenPtPrompt.GetBinError(iBin) / hGenPtPrompt.GetBinContent(iBin)
                    ptCent = hGenPtPrompt.GetBinCenter(iBin)
                    hGenPtPrompt.SetBinContent(iBin, hGenPtPrompt.GetBinContent(iBin) * sPtWeights(ptCent))
                    hGenPtPrompt.SetBinError(iBin, hGenPtPrompt.GetBinContent(iBin) * relStatUnc)
        hGenPtPrompt.SetName(f'hPromptGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        prompt_gen_list.append(hGenPtPrompt)
        hGenPtPrompt.Write()
        hGenPtFD = sparseGen['GenFD'].Projection(0)
        # apply pt weights
        if args.ptweightsB or args.ptweights:
            for iBin in range(1, hGenPtFD.GetNbinsX()+1):
                if hGenPtFD.GetBinContent(iBin) > 0:
                    relStatUnc = hGenPtFD.GetBinError(iBin) / hGenPtFD.GetBinContent(iBin)
                    ptCent = hGenPtFD.GetBinCenter(iBin)
                    hGenPtFD.SetBinContent(iBin, hGenPtFD.GetBinContent(iBin) * sPtWeightsGenDfromB(ptCent))
                    hGenPtFD.SetBinError(iBin, hGenPtFD.GetBinContent(iBin) * relStatUnc)
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
        axisNum = cutVars[iVar]['axisnum']
        if axisNum >=2: # check if axis is a cut variable (not inv. mass or pt)
            if isRedVar:
                axisNum -= shiftForRedVar
            if isWithBinfo:
                axisNum += 2
        if 'RecoAll' in sparseReco:
            sparseReco['RecoAll'].GetAxis(axisNum).SetRange(-1, -1)
        if isMC:
            sparseReco['RecoPrompt'].GetAxis(axisNum).SetRange(-1, -1)
            sparseReco['RecoFD'].GetAxis(axisNum).SetRange(-1, -1)
            if enableSecPeak:
                sparseReco['RecoSecPeakPrompt'].GetAxis(axisNum).SetRange(-1, -1)
                sparseReco['RecoSecPeakFD'].GetAxis(axisNum).SetRange(-1, -1)

for iPt in range(0, len(cutVars['Pt']['min']) - 1):
    ptLowLabel = cutVars['Pt']['min'][iPt] * 10
    ptHighLabel = cutVars['Pt']['max'][iPt+1] * 10
    for iVar in ('InvMass', 'Pt'):
        varName = cutVars[iVar]['name']
        if 'RecoAll' in sparseReco:
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
        hVarFDGen_merged = MergeHists([fd_gen_list[iPt], fd_gen_list[iPt+1]])
        hVarFDGen_merged.SetName(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarFDGen_merged.Write()
        if enableSecPeak:
            hVarPromptGen_secpeak_merged = MergeHists([prompt_gen_list_secpeak[iPt], prompt_gen_list_secpeak[iPt+1]])
            hVarPromptGen_secpeak_merged.SetName(f'hPromptGenSecPeakPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarPromptGen_secpeak_merged.Write()
            hVarFDGen_secpeak_merged = MergeHists([fd_gen_list_secpeak[iPt], fd_gen_list_secpeak[iPt+1]])
            hVarFDGen_secpeak_merged.SetName(f'hFDGenSecPeakPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarFDGen_secpeak_merged.Write()

# merge all pT bins
ptLowLabel = cutVars['Pt']['min'][0] * 10
ptHighLabel = cutVars['Pt']['max'][-1] * 10
for iVar in ('InvMass', 'Pt'):
    varName = cutVars[iVar]['name']
    if 'RecoAll' in sparseReco:
        hVar_merged_allPt = MergeHists(all_dict[iVar])
        hVar_merged_allPt.SetName(f'h{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVar_merged_allPt.Write()
    if isMC:
        hVarPrompt_merged_allPt = MergeHists(prompt_dict[iVar])
        hVarPrompt_merged_allPt.SetName(f'hPrompt{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarPrompt_merged_allPt.Write()
        hVarFD_merged_allPt = MergeHists(fd_dict[iVar])
        hVarFD_merged_allPt.SetName(f'hFD{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarFD_merged_allPt.Write()
        if enableSecPeak:
            hVarPrompt_secpeak_merged_allPt = MergeHists(prompt_dict_secpeak[iVar])
            hVarPrompt_secpeak_merged_allPt.SetName(f'hPromptSecPeak{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarPrompt_secpeak_merged_allPt.Write()
            hVarFD_secpeak_merged_allPt = MergeHists(fd_dict_secpeak[iVar])
            hVarFD_secpeak_merged_allPt.SetName(f'hFDSecPeak{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hVarFD_secpeak_merged_allPt.Write()
if isMC:
    hVarPromptGen_merged_allPt = MergeHists(prompt_gen_list)
    hVarPromptGen_merged_allPt.SetName(f'hPromptGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
    hVarPromptGen_merged_allPt.Write()
    hVarFDGen_merged_allPt = MergeHists(fd_gen_list)
    hVarFDGen_merged_allPt.SetName(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
    hVarFDGen_merged_allPt.Write()
    if enableSecPeak:
        hVarPromptGen_secpeak_merged_allPt = MergeHists(prompt_gen_list_secpeak)
        hVarPromptGen_secpeak_merged_allPt.SetName(f'hPromptGenSecPeakPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarPromptGen_secpeak_merged_allPt.Write()
        hVarFDGen_secpeak_merged_allPt = MergeHists(fd_gen_list_secpeak)
        hVarFDGen_secpeak_merged_allPt.SetName(f'hFDGenSecPeakPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hVarFDGen_secpeak_merged_allPt.Write()

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
