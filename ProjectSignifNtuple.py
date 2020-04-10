'''
python script for the projection of the Tupla generated in ScanSelections.py script
run: python ScanSelectionTree.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import sys
import argparse
import itertools
import yaml
import numpy as np
from root_numpy import fill_hist
from ROOT import TFile, TH1F, TCanvas, TNtuple, TDirectoryFile  # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, kGreen, kOpenCrossX
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.DfUtils import LoadDfFromRootOrParquet

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                    help='input file with cut set')
parser.add_argument('outFileName', metavar='text', default='outFileName.root',
                    help='output root file name')
args = parser.parse_args()

#config input file and df definition
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
inFileNames = inputCfg['infiles']['name']
dfSignif = LoadDfFromRootOrParquet(inputCfg['infiles']['name'], inputCfg['infiles']['dirname'],
                                   inputCfg['infiles']['treename'])
dfSignif['Pt'] = dfSignif.apply(lambda row: (row.PtMin + row.PtMax) / 2, axis=1)
VarDrawList = inputCfg['VarDrawList']
if not isinstance(VarDrawList, list):
    VarDrawList = [VarDrawList]

#selections to be applied
with open(args.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']
if not 'ML_output_Bkg' or not 'ML_output_FD' in cutVars:
    print('\t\t---Warning: no ML Bkg or FD output cut was provided. Are you sure you want to continue?---\n')
selToApply = []
counter = 0
for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
    selToApply.append('')
    for iVar, varName in enumerate(cutVars):
        if selToApply[counter] != '':
            selToApply[counter] += ' & '
        selToApply[counter] += f"{cutVars[varName]['min'][counter]} < {cutVars[varName]['name']} <= {cutVars[varName]['max'][counter]}"
    counter += 1

#output file preparation
outFile = TFile(args.outFileName, 'RECREATE')
outFile.cd()
TProject = TCanvas('DsNtupleProjOverPt', '', 1920, 1080)
TProject.Divide(2, round(len(VarDrawList)/2))
hProject = []

#output histos
nbins = len(cutVars['Pt']['max'])
xbins = cutVars['Pt']['min'] + list(set(cutVars['Pt']['max']) - set(cutVars['Pt']['min']))
for iVar, VartoDraw in enumerate(VarDrawList):
    hProject.append(TH1F(f'hProject{VartoDraw}', f'{VartoDraw} over p''_{T}'f'; p''_{T}'' [GeV/c] ;'f'{VartoDraw}',
                    nbins, np.asarray(xbins, float)))
    SetObjectStyle(hProject[iVar], color=kGreen-iVar, markerstyle=kOpenCrossX, markersize=1.5, linewidh=2, linestyle=7)
    for iPt, _ in enumerate(cutVars['Pt']['min']):
        dfSignifSel = dfSignif.query(selToApply[iPt])
        hProject[iVar].SetBinContent(iPt+1, dfSignifSel[f'{VartoDraw}'])
        if VartoDraw == 'EffAccFD' or VartoDraw == 'EffAccPrompt' or VartoDraw == 'S':
            hProject[iVar].SetBinError(iPt+1, dfSignifSel[f'{VartoDraw}Error'])          
    TProject.cd(iVar+1)
    hProject[iVar].DrawCopy()
    TProject.Update()
    TProject.Modified()
TProject.Write()
if inputCfg['saveaspdf']:
    TProject.SaveAs(f'DsNtupleProjOverPt.pdf')

#saving output and closing
outFile.Write()
outFile.Close()
input('Press enter to exit')