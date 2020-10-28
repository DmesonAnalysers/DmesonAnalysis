'''
python script to modify the prompt Lc FONLL prediction as D0 FONLL x measured Lc/D0 ratio
'''

import sys
import os
import argparse
from ROOT import TObject, TCanvas, TFile, TMath #pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kFullSquare, kOpenCircle #pylint: disable=import-error,no-name-in-module
sys.path.append('..')
#pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

SetGlobalStyle(padbottommargin=0.14)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inFileName', metavar='text', default='predictions.root',
                    help='input root file with FONLL predictions')
parser.add_argument('--BRD0', type=float, default=0.0389,
                    help='BR for D0 -> Kpi channel. This must correspond to the one used in the original file!')
parser.add_argument('--BRpKpi', type=float, default=0.0623,
                    help='BR for Lc -> pKpi channel. It is convenient to have the same as for FD')
parser.add_argument('--BRpK0s', type=float, default=0.0158,
                    help='BR for Lc -> pK0s channel. It is convenient to have the same as for FD')
args = parser.parse_args()

FONLLVar = ['central', 'min', 'max']
hPromptD0, hPromptLcpK0s, hPromptLcpKpi, hPromptLcpiL = ({} for _ in range(4))
inFile = TFile.Open(args.inFileName)
for var in FONLLVar:
    hPromptD0[var] = inFile.Get(f'hD0Kpipred_{var}')
    hPromptD0[var].Scale(1./args.BRD0)
    hPromptLcpK0s[var] = inFile.Get(f'hLcK0sppred_{var}')
    hPromptLcpKpi[var] = inFile.Get(f'hLcpkpipred_{var}')

for iPt in range(1, hPromptD0['central'].GetNbinsX()+1):
    for var in FONLLVar:
        pT = hPromptD0[var].GetBinCenter(iPt)
        crossSec = hPromptD0[var].GetBinContent(iPt)
        factor = 0.11 + 4.68735e-01 * TMath.Gaus(pT, 1, 4.90037) # parametrised from 5.02 TeV measurement
        hPromptLcpK0s[var].SetBinContent(iPt, factor * crossSec * args.BRpK0s)
        hPromptLcpKpi[var].SetBinContent(iPt, factor * crossSec * args.BRpKpi)

outFileName = args.inFileName.replace('.root', '_PromptLcMod.root')
os.system(f'cp {args.inFileName} {outFileName}')

outFile = TFile.Open(outFileName, 'update')
for var in FONLLVar:
    hPromptLcpK0s[var].Write(hPromptLcpK0s[var].GetName(), TObject.kOverwrite)
    hPromptLcpKpi[var].Write(hPromptLcpKpi[var].GetName(), TObject.kOverwrite)
outFile.Close()

# test that Lc/D0 is the expected one
cLcOverD0 = TCanvas('cLcOverD0', '', 800, 800)
hRatioLcpK0s = hPromptD0['central'].Clone('hRatioLcpK0s')
hRatioLcpKpi = hPromptD0['central'].Clone('hRatioLcpKpi')
hRatioLcpK0s.Divide(hPromptLcpK0s['central'], hPromptD0['central'], 1./args.BRpK0s, 1., '')
hRatioLcpKpi.Divide(hPromptLcpKpi['central'], hPromptD0['central'], 1./args.BRpKpi, 1., '')
hRatioLcpK0s.SetStats(0)
hRatioLcpKpi.SetStats(0)
SetObjectStyle(hRatioLcpKpi, color=kBlack, markerstyle=kFullSquare)
SetObjectStyle(hRatioLcpK0s, color=kRed, markerstyle=kOpenCircle)
cLcOverD0.DrawFrame(0., 0., 50., 1., ';#it{p}_{T} (GeV/#it{c});#Lambda_{c}^{+} / D^{0}')
hRatioLcpKpi.Draw('same')
hRatioLcpK0s.Draw('same')
cLcOverD0.Modified()
cLcOverD0.Update()

input('Press enter to exit')
