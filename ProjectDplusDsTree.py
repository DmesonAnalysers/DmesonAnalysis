'''
python script for the projection of of D+ and Ds+ particles TTrees
run: python ProjectDplusDsTree.py cfgFileName.yml cutSetFileName.yml outFileName.root
                                  [--ptweights PtWeightsFileName.root histoName]
                                  [--ptweightsB PtWeightsFileName.root histoName]
                                  [--std]

if the --ptweights argument is provided, pT weights will be applied to prompt and FD pT distributions
if the --ptweightsB argument is provided, pT weights will be applied to FD pT distributions instead of
those for the prompt

--std, used to apply standard analysis cuts on tree (account for differences in conventions)
'''

import sys
import argparse
import yaml
import numpy as np
import uproot
from scipy.interpolate import InterpolatedUnivariateSpline
from root_numpy import fill_hist
from ROOT import TFile, TH1F, TDatabasePDG # pylint: disable=import-error,no-name-in-module
from utils.TaskFileLoader import LoadNormObjFromTask, LoadSparseFromTask
from utils.DfUtils import FilterBitDf, LoadDfFromRootOrParquet
from utils.AnalysisUtils import MergeHists, ApplySplineFuncToColumn

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
parser.add_argument('--std', help='adapt to std. analysis cuts', action='store_true')
args = parser.parse_args()

#config with input file details
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
inFileNames = inputCfg['filename']
if not isinstance(inFileNames, list):
    inFileNames = [inFileNames]
isMC = inputCfg['isMC']
if not isMC:
    if args.ptweights:
        print('WARNING: pt weights will not be applied since it is not MC')
        args.ptweights = None
    if args.ptweightsB:
        print('WARNING: ptB weights will not be applied since it is not MC')
        args.ptweightsB = None

#define filter bits
bitSignal = 0
bitPrompt = 2
bitFD = 3
bitRefl = 4
# define mass binning
particle = inputCfg['tree']['particle']
if particle == 'Ds':
    mD = TDatabasePDG.Instance().GetParticle(431).Mass()
elif particle == 'Dplus':
    mD = TDatabasePDG.Instance().GetParticle(411).Mass()
elif particle == 'Lc':
    mD = TDatabasePDG.Instance().GetParticle(4122).Mass()
else:
    print('Error: only Dplus, Ds particles and Lc supported. Exit!')
    sys.exit()
massBins = 500
massLimLow = mD - 0.25
massLimHigh = mD + 0.25

# selections to be applied
with open(args.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']
selToApply = []
for iPt, _ in enumerate(cutVars['Pt']['min']):
    selToApply.append('')
    for varName in cutVars:
        if varName == 'InvMass':
            continue
        if selToApply[iPt] != '':
            selToApply[iPt] += ' & '
        if args.std and varName == 'CosPiKPhi3':
            selToApply[iPt] += '~'
        selToApply[iPt] += f"({cutVars[varName]['min'][iPt]}<{cutVars[varName]['name']}<{cutVars[varName]['max'][iPt]})"

# dicts of TH1
allDict = {'InvMass': [], 'Pt': []}
promptDict = {'InvMass': [], 'Pt': []}
FDDict = {'InvMass': [], 'Pt': []}
promptGenList = []
FDGenList = []
# TODO: add second peak histograms for Ds

outFile = TFile(args.outFileName, 'recreate')

# load objects from task outputs
for iFile, inFileName in enumerate(inFileNames):
    if iFile == 0:
        hEv, normCounter = LoadNormObjFromTask(inFileName, inputCfg)
        if isMC:
            _, sparseGen = LoadSparseFromTask(inFileName, inputCfg) #only gen sparses used
    else:
        hEvPart, normCounterPart = LoadNormObjFromTask(inFileName, inputCfg)
        hEv.Add(hEvPart)
        normCounter.Add(normCounterPart)
        if isMC:
            _, sparseGenPart = LoadSparseFromTask(inFileName, inputCfg) #only gen sparses used
            for sparseType in sparseGenPart:
                sparseGen[sparseType].Add(sparseGenPart[sparseType])


# define pT binning (from gen sparses if MC)
if isMC:
    nPtBins = sparseGen['GenPrompt'].GetAxis(0).GetNbins()
    ptLimLow = sparseGen['GenPrompt'].GetAxis(0).GetBinLowEdge(1)
    ptLimHigh = sparseGen['GenPrompt'].GetAxis(0).GetBinLowEdge(nPtBins) + \
        sparseGen['GenPrompt'].GetAxis(0).GetBinWidth(nPtBins)
else:
    nPtBins = 500
    ptLimLow = 0.
    ptLimHigh = 50.
ptBinWidth = (ptLimHigh-ptLimLow) / nPtBins

# load trees
if isMC:
    dataFramePrompt = LoadDfFromRootOrParquet(inputCfg['tree']['filenamePrompt'], inputCfg['tree']['dirname'],
                                              inputCfg['tree']['treename'])
    if 'cand_type' in dataFramePrompt.columns: #if not filtered tree, select only prompt and not reflected
        dataFramePrompt = FilterBitDf(dataFramePrompt, 'cand_type', [bitSignal, bitPrompt], 'and')
        dataFramePrompt = FilterBitDf(dataFramePrompt, 'cand_type', [bitRefl], 'not')
    dataFramePrompt.reset_index(inplace=True)

    dataFrameFD = LoadDfFromRootOrParquet(inputCfg['tree']['filenameFD'], inputCfg['tree']['dirname'],
                                          inputCfg['tree']['treename'])
    if 'cand_type' in dataFrameFD.columns: #if not filtered tree, select only FD and not reflected
        dataFrameFD = FilterBitDf(dataFrameFD, 'cand_type', [bitSignal, bitFD], 'and')
        dataFrameFD = FilterBitDf(dataFrameFD, 'cand_type', [bitRefl], 'not')
    dataFrameFD.reset_index(inplace=True)

    # compute pt weights
    if args.ptweights:
        ptWeights = uproot.open(args.ptweights[0])[args.ptweights[1]]
        ptCentW = [(ptWeights.edges[iBin]+ptWeights.edges[iBin+1])/2 for iBin in range(len(ptWeights.edges)-1)]
        sPtWeights = InterpolatedUnivariateSpline(ptCentW, ptWeights.values)
        dataFramePrompt['pt_weights'] = ApplySplineFuncToColumn(dataFramePrompt, 'pt_cand', sPtWeights, 0, 50)
        if not args.ptweightsB:
            dataFrameFD['pt_weights'] = ApplySplineFuncToColumn(dataFrameFD, 'pt_cand', sPtWeights, 0, 50)
            sPtWeightsDfromB = sPtWeights

    if args.ptweightsB:
        ptWeightsB = uproot.open(args.ptweights[0])[args.ptweights[1]]
        ptCentWB = [(ptWeightsB.edges[iBin]+ptWeightsB.edges[iBin+1])/2 for iBin in range(len(ptWeights.edges)-1)]
        sPtWeightsB = InterpolatedUnivariateSpline(ptCentWB, ptWeightsB.values)
        dataFrameFD['pt_weights'] = ApplySplineFuncToColumn(dataFrameFD, 'pt_B', sPtWeightsB, 0, 50)
        # average correction for gen part since tree not available (--> good approximation)
        hPtBvsPtGenD = sparseGen['GenFD'].Projection(2, 0).ProfileX()
        ptCentGen, averagePtBvsPtGen = [], []
        for iPt in range(1, hPtBvsPtGenD.GetNbinsX()+1):
            ptCentGen.append(hPtBvsPtGenD.GetBinCenter(iPt))
            averagePtBvsPtGen.append(hPtBvsPtGenD.GetBinContent(iPt))
        aPtGenWeightsB = list(sPtWeightsB(averagePtBvsPtGen))
        sPtWeightsDfromB = InterpolatedUnivariateSpline(ptCentGen, aPtGenWeightsB)

    for (cuts, ptMin, ptMax) in zip(selToApply, cutVars['Pt']['min'], cutVars['Pt']['max']):
        print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
        ptLowLabel = ptMin * 10
        ptHighLabel = ptMax * 10

        # gen histos from sparses
        binGenMin = sparseGen['GenPrompt'].GetAxis(0).FindBin(ptMin * 1.0001)
        binGenMax = sparseGen['GenPrompt'].GetAxis(0).FindBin(ptMax * 0.9999)
        sparseGen['GenPrompt'].GetAxis(0).SetRange(binGenMin, binGenMax)
        sparseGen['GenFD'].GetAxis(0).SetRange(binGenMin, binGenMax)

        hGenPtPrompt = sparseGen['GenPrompt'].Projection(0)
        hGenPtPrompt.Sumw2()
        if args.ptweights:
            for iPt in range(1, hGenPtPrompt.GetNbinsX()+1):
                if hGenPtPrompt.GetBinContent(iPt) > 0:
                    relStatUnc = hGenPtPrompt.GetBinError(iPt) / hGenPtPrompt.GetBinContent(iPt)
                    ptCent = hGenPtPrompt.GetBinCenter(iPt)
                    hGenPtPrompt.SetBinContent(iPt, hGenPtPrompt.GetBinContent(iPt) * sPtWeights(ptCent))
                    hGenPtPrompt.SetBinError(iPt, hGenPtPrompt.GetBinContent(iPt) * relStatUnc)
        hGenPtPrompt.SetName(f'hPromptGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        promptGenList.append(hGenPtPrompt)

        hGenPtFD = sparseGen['GenFD'].Projection(0)
        hGenPtFD.Sumw2()
        if args.ptweights or args.ptweightsB:
            for iPt in range(1, hGenPtFD.GetNbinsX()+1):
                if hGenPtFD.GetBinContent(iPt) > 0:
                    relStatUnc = hGenPtFD.GetBinError(iPt) / hGenPtFD.GetBinContent(iPt)
                    ptCent = hGenPtFD.GetBinCenter(iPt)
                    hGenPtFD.SetBinContent(iPt, hGenPtFD.GetBinContent(iPt) * sPtWeightsDfromB(ptCent))
                    hGenPtFD.SetBinError(iPt, hGenPtFD.GetBinContent(iPt) * relStatUnc)
        hGenPtFD.SetName(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        FDGenList.append(hGenPtFD)

        # reco histos from trees
        dataFramePromptSel = dataFramePrompt.astype(float).query(cuts)
        dataFrameFDSel = dataFrameFD.astype(float).query(cuts)
        hPtPrompt = TH1F(f'hPromptPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}', '', nPtBins, ptLimLow, ptLimHigh)
        hInvMassPrompt = TH1F(f'hPromptMass_{ptLowLabel:.0f}_{ptHighLabel:.0f}', '', massBins, massLimLow, massLimHigh)
        hPtFD = TH1F(f'hFDPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}', '', nPtBins, ptLimLow, ptLimHigh)
        hInvMassFD = TH1F(f'hFDMass_{ptLowLabel:.0f}_{ptHighLabel:.0f}', '', massBins, massLimLow, massLimHigh)

        if args.ptweights:
            hTmp = hPtPrompt.Clone('hTmp')
            fill_hist(hTmp, dataFramePromptSel['pt_cand'].values) # for stat unc
            fill_hist(hPtPrompt, dataFramePromptSel['pt_cand'].values, weights=dataFramePromptSel['pt_weights'].values)
            for iPt in range(1, hTmp.GetNbinsX()+1):
                if hTmp.GetBinContent(iPt) == 0.:
                    hPtPrompt.SetBinError(iPt, 0.)
                else:
                    hPtPrompt.SetBinError(iPt, 1./np.sqrt(hTmp.GetBinContent(iPt))*hPtPrompt.GetBinContent(iPt))
        else:
            fill_hist(hPtPrompt, dataFramePromptSel['pt_cand'].values)
            hPtPrompt.Sumw2()
        fill_hist(hInvMassPrompt, dataFramePromptSel['inv_mass'].values)

        if args.ptweightsB or args.ptweights:
            hTmp = hPtFD.Clone('hTmp')
            fill_hist(hTmp, dataFrameFDSel['pt_cand'].values) # for stat unc
            fill_hist(hPtFD, dataFrameFDSel['pt_cand'].values, weights=dataFrameFDSel['pt_weights'].values)
            for iPt in range(1, hTmp.GetNbinsX()+1):
                if hTmp.GetBinContent(iPt) == 0.:
                    hPtFD.SetBinError(iPt, 0.)
                else:
                    hPtFD.SetBinError(iPt, 1./np.sqrt(hTmp.GetBinContent(iPt))*hPtFD.GetBinContent(iPt))
        else:
            fill_hist(hPtFD, dataFrameFDSel['pt_cand'].values)
            hPtFD.Sumw2()
        fill_hist(hInvMassFD, dataFrameFDSel['inv_mass'].values)

        promptDict['InvMass'].append(hInvMassPrompt)
        promptDict['Pt'].append(hPtPrompt)
        FDDict['InvMass'].append(hInvMassFD)
        FDDict['Pt'].append(hPtFD)
        outFile.cd()
        hGenPtPrompt.Write()
        hGenPtFD.Write()
        hPtPrompt.Write()
        hInvMassPrompt.Write()
        hPtFD.Write()
        hInvMassFD.Write()

    # merge adiacent pt bin histograms
    for iPt in range(0, len(cutVars['Pt']['min']) - 1):
        ptLowLabel = cutVars['Pt']['min'][iPt] * 10
        ptHighLabel = cutVars['Pt']['max'][iPt+1] * 10
        for iVar in ('InvMass', 'Pt'):
            if iVar == 'Pt':
                varName = iVar
            else:
                varName = 'Mass'
            hPromptMerged = MergeHists([promptDict[iVar][iPt], promptDict[iVar][iPt+1]])
            hPromptMerged.SetName(f'hPrompt{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hPromptMerged.Write()
            hFDMerged = MergeHists([FDDict['Pt'][iPt], FDDict['Pt'][iPt+1]])
            hFDMerged.SetName(f'hFD{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
            hFDMerged.Write()
        hPtPromptGenMerged = MergeHists([promptGenList[iPt], promptGenList[iPt+1]])
        hPtPromptGenMerged.SetName(f'hPromptGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hPtPromptGenMerged.Write()
        hPtFDGenMerged = MergeHists([FDGenList[iPt], FDGenList[iPt+1]])
        hPtFDGenMerged.SetName(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hPtFDGenMerged.Write()

else:
    dataFrame = LoadDfFromRootOrParquet(inputCfg['tree']['filenameAll'], inputCfg['tree']['dirname'],
                                        inputCfg['tree']['treename'])

    for (cuts, ptMin, ptMax) in zip(selToApply, cutVars['Pt']['min'], cutVars['Pt']['max']):
        print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
        ptLowLabel = ptMin * 10
        ptHighLabel = ptMax * 10
        dataFrameSel = dataFrame.astype(float).query(cuts)
        hPt = TH1F(f'hPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}', '', nPtBins, ptLimLow, ptLimHigh)
        hInvMass = TH1F(f'hMass_{ptLowLabel:.0f}_{ptHighLabel:.0f}', '', massBins, massLimLow, massLimHigh)
        fill_hist(hPt, dataFrameSel['pt_cand'].values)
        fill_hist(hInvMass, dataFrameSel['inv_mass'].values)
        allDict['InvMass'].append(hInvMass)
        allDict['Pt'].append(hPt)
        outFile.cd()
        hPt.Write()
        hInvMass.Write()

    # merge adiacent pt bin histograms
    for iPt in range(0, len(cutVars['Pt']['min']) - 1):
        ptLowLabel = cutVars['Pt']['min'][iPt] * 10
        ptHighLabel = cutVars['Pt']['max'][iPt+1] * 10
        hPtMerged = MergeHists([allDict['Pt'][iPt], allDict['Pt'][iPt+1]])
        hPtMerged.SetName(f'hPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hPtMerged.Write()
        hInvMassMerged = MergeHists([allDict['InvMass'][iPt], allDict['InvMass'][iPt+1]])
        hInvMassMerged.SetName(f'hMass_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hInvMassMerged.Write()

# merge all pT bins
ptLowLabel = cutVars['Pt']['min'][0] * 10
ptHighLabel = cutVars['Pt']['max'][-1] * 10
for iVar in ('InvMass', 'Pt'):
    if iVar == 'Pt':
        varName = iVar
    else:
        varName = 'Mass'
    if not isMC:
        hAllMergedAllPt = MergeHists(allDict[iVar])
        hAllMergedAllPt.SetName(f'h{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hAllMergedAllPt.Write()
    else:
        hPromptMergedAllPt = MergeHists(promptDict[iVar])
        hPromptMergedAllPt.SetName(f'hPrompt{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hPromptMergedAllPt.Write()
        hFDMergedAllPt = MergeHists(FDDict[iVar])
        hFDMergedAllPt.SetName(f'hFD{varName}_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
        hFDMergedAllPt.Write()
if isMC:
    hPromptGenMergedAllPt = MergeHists(promptGenList)
    hPromptGenMergedAllPt.SetName(f'hPromptGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
    hPromptGenMergedAllPt.Write()
    hFDGenMergedAllPt = MergeHists(FDGenList)
    hFDGenMergedAllPt.SetName(f'hFDGenPt_{ptLowLabel:.0f}_{ptHighLabel:.0f}')
    hFDGenMergedAllPt.Write()

# normalisation
hEvForNorm = TH1F("hEvForNorm", ";;Number of events", 2, 0., 2.)
hEvForNorm.GetXaxis().SetBinLabel(1, "norm counter")
hEvForNorm.GetXaxis().SetBinLabel(2, "accepted events")
hEvForNorm.SetBinContent(1, normCounter.GetNEventsForNorm())
for iBin in range(1, hEv.GetNbinsX() + 1):
    binLabel = hEv.GetXaxis().GetBinLabel(iBin)
    if 'isEvSelected' in binLabel or 'accepted' in binLabel:
        hEvForNorm.SetBinContent(2, hEv.GetBinContent(iBin))
        break
outFile.cd()
hEvForNorm.Write()
outFile.Close()
