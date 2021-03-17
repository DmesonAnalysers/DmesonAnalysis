'''
python script for the projection of the D-meson mass spectrum from femto task
'''

import sys
import argparse
import numpy as np
import yaml
from ROOT import TFile, TH1F # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.AnalysisUtils import MergeHists #pylint: disable=wrong-import-position,import-error

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFileName', metavar='text', default='inFileName.root',
                    help='root input file name')
parser.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                    help='input file with cut set')
parser.add_argument('outFileName', metavar='text', default='outFileName.root',
                    help='output root file name')
parser.add_argument('--prefix', metavar='text', default='MB',
                    help='prefix for directory inside task output file')
parser.add_argument('--suffix', metavar='text', default='0',
                    help='suffix for directory inside task output file')
parser.add_argument('--HFsuffix', metavar='text', default='',
                    help='HF suffix for directory inside task output file')
args = parser.parse_args()

with open(args.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']

dirName = f'{args.prefix}_CharmFemto_{args.HFsuffix}DChargedQA{args.suffix}'
listName = f'{dirName}/{args.prefix}_CharmFemto_{args.HFsuffix}DChargedQA{args.suffix}'
print(f'Read input file: {args.inFileName}')
print(f'           list: {listName}')

inFile = TFile.Open(args.inFileName)
inList = inFile.Get(listName)
hDminusMassVsPt = inList.FindObject('fHistDminusInvMassPt')
hDplusMassVsPt = inList.FindObject('fHistDplusInvMassPt')
hMassVsPt = hDminusMassVsPt.Clone()
if hDminusMassVsPt.GetYaxis().GetXmin() == hDplusMassVsPt.GetXaxis().GetXmin():
    hMassVsPt.Add(hDplusMassVsPt)
else:
    print('WARNING: using only D- because of different binning between D+ and D- histograms!')

hMass, hPt = [], []
for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
    ptBinMin = hMassVsPt.GetXaxis().FindBin(ptMin*1.0001)
    ptBinMax = hMassVsPt.GetXaxis().FindBin(ptMax*0.9999)
    hMassVsPt.GetXaxis().SetRange(ptBinMin, ptBinMax)
    hMass.append(hMassVsPt.ProjectionY(f'hMass_{ptMin*10:.0f}_{ptMax*10:.0f}'))
    hPt.append(hMassVsPt.ProjectionX(f'hPt_{ptMin*10:.0f}_{ptMax*10:.0f}'))
hMass.append(MergeHists(hMass))
hPt.append(MergeHists(hPt))
hMass[-1].SetName(f'hMass_{cutVars["Pt"]["min"][0]*10:.0f}_{cutVars["Pt"]["max"][-1]*10:.0f}')
hPt[-1].SetName(f'hPt{cutVars["Pt"]["min"][0]*10:.0f}_{cutVars["Pt"]["max"][-1]*10:.0f}')

dirName = f'{args.prefix}_CharmFemto_QA0'
listName = f'{dirName}/{args.prefix}_CharmFemto_QA0'

inListEvents = inFile.Get(listName)
inListAliEventCuts = inListEvents.FindObject('AliEventCuts')
hEvents = inListAliEventCuts.FindObject('fCutStats')
nEvents = hEvents.GetBinContent(hEvents.GetNbinsX())

hEvForNorm = TH1F("hEvForNorm", ";;Number of events", 2, 0., 2.)
hEvForNorm.GetXaxis().SetBinLabel(1, "norm counter")
hEvForNorm.GetXaxis().SetBinLabel(2, "accepted events")
hEvForNorm.SetBinContent(1, nEvents) # Norm counter not present in femto task for the time being
hEvForNorm.SetBinContent(2, nEvents)

outFile = TFile(args.outFileName, 'recreate')
for hM, hP in zip(hMass, hPt):
    hM.Write()
    hP.Write()
hEvForNorm.Write()
outFile.Close()

print(f'Projected histograms saved in output file: {args.outFileName}')
