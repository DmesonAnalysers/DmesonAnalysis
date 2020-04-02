'''
python script for the projection of the Tupla generated in ScanSelections.py script
run: python ScanSelectionTree.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import sys
import argparse
import itertools
import six
import numpy as np
import yaml
from root_numpy import fill_hist
from ROOT import TFile, TH1F, TH2F, TCanvas, TNtuple, TDirectoryFile  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                    help='input file with cut set')
parser.add_argument('outFileName', metavar='text', default='outFileName.root',
                    help='output root file name')
args = parser.parse_args()

#config input file and tuple extraction
with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
inFileNames = inputCfg['infilename']
infile = TFile(inFileNames)
tSignif = infile.Get('tSignif')
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
for iPt, _ in enumerate(cutVars['PtMin']['min']):
    for iBkg, (Bkg_score_min, Bkg_score_max)   in enumerate(zip(cutVars['ML_output_Bkg']['min'], cutVars['ML_output_Bkg']['max'])):
        for iFD, (FD_score_min, FD_score_max) in enumerate(zip(cutVars['ML_output_FD']['min'], cutVars['ML_output_FD']['max'])):
            selToApply.append('')
            for iVar, varName in enumerate(cutVars):
                if selToApply[counter] != '':
                    selToApply[counter] += ' & '
                if varName == 'PtMin':
                    selToApply[counter] += f"{cutVars[varName]['name']}>={cutVars[varName]['min'][iPt]}"
                elif varName == 'PtMax':
                    selToApply[counter] += f"{cutVars[varName]['name']}<={cutVars[varName]['max'][iPt]}"
                elif varName == 'ML_output_Bkg':
                    selToApply[counter] += f" {cutVars[varName]['name']}>={Bkg_score_min} & {cutVars[varName]['name']}<={Bkg_score_max}"
                elif varName == 'ML_output_FD':
                    selToApply[counter] += f"{cutVars[varName]['name']}>={FD_score_min} & {cutVars[varName]['name']}<={FD_score_max}"
            counter += 1

#output file
outFile = TFile(args.outFileName, 'RECREATE')
outFile.cd()
outDirPros = TDirectoryFile('DsNTupleProjection', 'DsNTupleProjection')
outDirPros.Write()
outDirProjection = []
cDist = []

#loop over cut and TNtuple projection
counter = 0 
for iPt, (PtMin, PtMax) in enumerate(zip(cutVars['PtMin']['min'], cutVars['PtMax']['max'])):
    print(f'Projecting tSignif for pT{PtMin}-{PtMax}')
    outDirPros.cd()
    outDirProjection.append(TDirectoryFile(f'pT{PtMin}-{PtMax}', f'pT{PtMin}-{PtMax}'))
    outDirProjection[iPt].Write()
    outDirProjection[iPt].cd()
    for iBkg, (Bkg_out_min, Bkg_out_max) in enumerate(zip(cutVars['ML_output_Bkg']['min'], cutVars['ML_output_Bkg']['max'])):
        for iFD, (FD_out_min, FD_out_max) in enumerate(zip(cutVars['ML_output_FD']['min'],cutVars['ML_output_FD']['max'])):
            outDirProjection[iPt].mkdir(f'ML_output_Bkg{Bkg_out_min}-{Bkg_out_max}-ML_output_FD{FD_out_min}-{FD_out_max}')
            outDirProjection[iPt].cd(f'ML_output_Bkg{Bkg_out_min}-{Bkg_out_max}-ML_output_FD{FD_out_min}-{FD_out_max}')
            cDist.append(TCanvas(f'pT{PtMin}-{PtMax}_ML_output_Bkg{Bkg_out_min}-{Bkg_out_max}-ML_output_FD{FD_out_min}-{FD_out_max}', '', 1920, 1080))
            cDist[counter].Divide(2, round(len(VarDrawList)/2))
            for iVar, VartoDraw in enumerate(VarDrawList):
                cDist[counter].cd(iVar+1)
                tSignif.Draw(f'{VartoDraw}>>h{VartoDraw}_pT{PtMin}-{PtMax}_ML_output_Bkg{Bkg_out_min}-{Bkg_out_max}-ML_output_FD{FD_out_min}-{FD_out_max}',
                             f'{selToApply[counter]}', 'histo')
                #TODO: need to check if there's a problem with cuts or is just too poor dataset
            cDist[counter].Write()
            if inputCfg['saveaspdf']:
                outDirPros.cd()
                cDist[counter].SaveAs(f'DsNtupleProj_pT{PtMin}-{PtMax}_ML_output_Bkg{Bkg_out_min}-{Bkg_out_max}-ML_output_FD{FD_out_min}-{FD_out_max}.pdf')
            counter += 1
outFile.Write()
outFile.Close()
if six.PY2:
    raw_input('Press enter to exit')
elif six.PY3:
    input('Press enter to exit')