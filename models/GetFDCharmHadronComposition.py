'''
python script to get the fraction of the non-prompt charmed hadrons coming from each B hadron
'''

import argparse
import numpy as np
from ROOT import TFile #pylint: disable=import-error,no-name-in-module

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFileName', metavar='text', default='predictions.root',
                    help='input root or txt file with non-prompt FONLL predictions')
args = parser.parse_args()

charmHadronNames = ['D0   ', 'Dplus', 'Ds   ', 'Lc   ', 'Dstar']
beautyHadronNames = ['B0', 'Bplus', 'Bs', 'Lb']
histNames = ['hBRHbtoD0', 'hBRHbtoDplus', 'hBRHbtoDs', 'hBRHbtoLc', 'hBRHbtoDstar']

brMatrix = np.zeros((len(charmHadronNames), len(beautyHadronNames)))
ffArray = np.zeros(len(beautyHadronNames))

#hHbFF
inFile = TFile.Open(args.inFileName)

ffHisto = inFile.Get('hHbFF')
for iBin in range(ffHisto.GetNbinsX()):
    ffArray[iBin] = ffHisto.GetBinContent(iBin+1)

for iHadron, histName in enumerate(histNames):
    histo = inFile.Get(histName)
    for iBin in range(histo.GetNbinsX()):
        brMatrix[iHadron][iBin] = histo.GetBinContent(iBin+1)

contrMatrix = np.multiply(brMatrix, ffArray)
normContrMatrix = (contrMatrix.T/contrMatrix.sum(axis=1)).T

print('     ', *beautyHadronNames, sep='       ')
for iHadron, hadronName in enumerate(charmHadronNames):
    print(hadronName, normContrMatrix[iHadron])
