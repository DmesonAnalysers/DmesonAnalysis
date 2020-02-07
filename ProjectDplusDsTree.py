'''
python script for the projection of of D+ and Ds+ mesons TTrees
run: python ProjectDplusDsTree.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import argparse
import yaml
from root_numpy import fill_hist
from ROOT import TFile, TH1F, TDatabasePDG # pylint: disable=import-error,no-name-in-module
from utils.TaskFileLoader import LoadNormObjFromTask, LoadSparseFromTask
from utils.DfUtils import FilterBitDf, LoadDfFromRootOrParquet


bitSignal = 0
bitPrompt = 2
bitFD = 3
bitRefl = 4

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                    help='input file with cut set')
parser.add_argument('outFileName', metavar='text', default='outFileName.root',
                    help='output root file name')
args = parser.parse_args()

#config with input file details
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
inFileNames = inputCfg['filename']
if not isinstance(inFileNames, list):
    inFileNames = [inFileNames]
isMC = inputCfg['isMC']

#define mass binning
meson = inputCfg['tree']['meson']
if meson == 'Ds':
    mD = TDatabasePDG.Instance().GetParticle(431).Mass()
elif meson == 'Dplus':
    mD = TDatabasePDG.Instance().GetParticle(411).Mass()
else:
    print('Error: only Dplus and Ds mesons supported. Exit!')
    exit()

#selections to be applied
with open(args.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']
selToApply = []
for iPt, _ in enumerate(cutVars['Pt']['min']):
    selToApply.append('')
    for iVar, varName in enumerate(cutVars):
        if varName == 'InvMass':
            continue
        if selToApply[iPt] != '':
            selToApply[iPt] += ' & '
        selToApply[iPt] += f"{cutVars[varName]['min'][iPt]}<{cutVars[varName]['name']}<{cutVars[varName]['max'][iPt]}"

# dicts of TH1
allDict = {'InvMass': [], 'Pt': []}
promptDict = {'InvMass': [], 'Pt': []}
FDDict = {'InvMass': [], 'Pt': []}
promptGenList = []
FDGenList = []

outFile = TFile(args.outFileName, 'recreate')

#load objects from task outputs
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

#load trees
if isMC:
    dataFramePrompt = LoadDfFromRootOrParquet(inputCfg['tree']['filenamePrompt'], inputCfg['tree']['dirname'],
                                              inputCfg['tree']['treename'])

    if 'cand_type' in dataFramePrompt.columns: #if not filtered tree, select only FD and not reflected
        dataFramePrompt = FilterBitDf(dataFramePrompt, 'cand_type', [bitSignal, bitPrompt], 'and')
        dataFramePrompt = FilterBitDf(dataFramePrompt, 'cand_type', [bitRefl], 'not')

    dataFrameFD = LoadDfFromRootOrParquet(inputCfg['tree']['filenameFD'], inputCfg['tree']['dirname'],
                                          inputCfg['tree']['treename'])

    if 'cand_type' in dataFrameFD.columns: #if not filtered tree, select only FD and not reflected
        dataFrameFD = FilterBitDf(dataFrameFD, 'cand_type', [bitSignal, bitFD], 'and')
        dataFrameFD = FilterBitDf(dataFrameFD, 'cand_type', [bitRefl], 'not')

    for iPt, (cuts, ptMin, ptMax) in enumerate(zip(selToApply, cutVars['Pt']['min'], cutVars['Pt']['max'])):
        print("Projecting distributions for %0.1f < pT < %0.1f GeV/c" % (ptMin, ptMax))
        #gen histos from sparses
        binGenMin = sparseGen['GenPrompt'].GetAxis(0).FindBin(ptMin*1.0001)
        binGenMax = sparseGen['GenPrompt'].GetAxis(0).FindBin(ptMax*0.9999)
        sparseGen['GenPrompt'].GetAxis(0).SetRange(binGenMin, binGenMax)
        sparseGen['GenFD'].GetAxis(0).SetRange(binGenMin, binGenMax)
        hGenPtPrompt = sparseGen['GenPrompt'].Projection(0)
        hGenPtPrompt.SetName('hPromptGenPt_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10))
        promptGenList.append(hGenPtPrompt)
        hGenPtFD = sparseGen['GenFD'].Projection(0)
        hGenPtFD.SetName('hFDGenPt_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10))
        FDGenList.append(hGenPtFD)
        #reco histos from trees
        dataFramePromptSel = dataFramePrompt.astype(float).query(cuts)
        dataFrameFDSel = dataFrameFD.astype(float).query(cuts)
        hPtPrompt = TH1F('hPromptPt_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10), '', 500, 0., 50.)
        hInvMassPrompt = TH1F('hPromptMass_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10), '', 400, mD-0.2, mD+0.2)
        hPtFD = TH1F('hFDPt_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10), '', 500, 0., 50.)
        hInvMassFD = TH1F('hFDMass_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10), '', 400, mD-0.2, mD+0.2)
        fill_hist(hPtPrompt, dataFramePromptSel['pt_cand'])
        fill_hist(hInvMassPrompt, dataFramePromptSel['inv_mass'])
        fill_hist(hPtFD, dataFrameFDSel['pt_cand'])
        fill_hist(hInvMassFD, dataFrameFDSel['inv_mass'])
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
else:

    dataFrame = LoadDfFromRootOrParquet(inputCfg['tree']['filenameAll'], inputCfg['tree']['dirname'],
                                        inputCfg['tree']['treename'])

    for iPt, (cuts, ptMin, ptMax) in enumerate(zip(selToApply, cutVars['Pt']['min'], cutVars['Pt']['max'])):
        print("Projecting distributions for %0.1f < pT < %0.1f GeV/c" % (ptMin, ptMax))
        dataFrameSel = dataFrame.astype(float).query(cuts)
        hPt = TH1F('hPt_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10), '', 500, 0., 50.)
        hInvMass = TH1F('hMass_{0:.0f}_{1:.0f}'.format(ptMin*10, ptMax*10), '', 400, mD-0.2, mD+0.2)
        fill_hist(hPt, dataFrameSel['pt_cand'])
        fill_hist(hInvMass, dataFrameSel['inv_mass'])
        allDict['InvMass'].append(hInvMass)
        allDict['Pt'].append(hPt)

        outFile.cd()
        hPt.Write()
        hInvMass.Write()

#normalisation
hEvForNorm = TH1F("hEvForNorm", ";;Number of events", 2, 0., 2.)
hEvForNorm.GetXaxis().SetBinLabel(1, "norm counter")
hEvForNorm.GetXaxis().SetBinLabel(2, "accepted events")
hEvForNorm.SetBinContent(1, normCounter.GetNEventsForNorm())
for iBin in range(1, hEv.GetNbinsX()+1):
    binLabel = hEv.GetXaxis().GetBinLabel(iBin)
    if 'isEvSelected' in binLabel or 'accepted' in binLabel:
        hEvForNorm.SetBinContent(2, hEv.GetBinContent(iBin))
        break

outFile.cd()
hEvForNorm.Write()
outFile.Close()
