'''
python script for the computation of the fractions of prompt and feed-down D for a given cut set
run: python ComputeDataDrivenFraction.py effAccFile.root fracFile.root outFile.root
'''

import argparse
import numpy as np
from ROOT import TFile, TCanvas, TLegend, TGraphErrors  # pylint: disable=import-error,no-name-in-module
from utils.StyleFormatter import SetGlobalStyle

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('effAccFileName', metavar='text', default='effAccFile.root',
                    help='root file with efficiency and acceptance')
parser.add_argument('fracFileName', metavar='text', default='fracFile.root',
                    help='root file with prompt (FD) fraction')
parser.add_argument('outFileName', metavar='text', default='outFile.root', help='root output file name')
args = parser.parse_args()

# load input file
effAccFile = TFile.Open(args.effAccFileName)
hEffAccPrompt = effAccFile.Get('hAccEffPrompt')
hEffAccFD = effAccFile.Get('hAccEffFD')

fracFile = TFile.Open(args.fracFileName)
hCorrYieldPrompt = fracFile.Get('hCorrYieldPrompt')
hCorrYieldFD = fracFile.Get('hCorrYieldFD')
hCovPromptPrompt = fracFile.Get('hCovPromptPrompt')
hCovPromptFD = fracFile.Get('hCovPromptFD')
hCovFDFD = fracFile.Get('hCovFDFD')

# systematic uncertainties --> total taken as sum in quadrature (uncorrelated)
hPromptFrac = hEffAccPrompt.Clone('hPromptFrac')
hFDFrac = hEffAccFD.Clone('hFDFrac')
hPromptFrac.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{f}_{prompt}')
hFDFrac.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{f}_{FD}')

for iPt in range(hEffAccPrompt.GetNbinsX()):
    ptMin = hEffAccPrompt.GetBinLowEdge(iPt+1)
    ptMax = ptMin+hEffAccPrompt.GetBinWidth(iPt+1)
    ptCent = hEffAccPrompt.GetBinCenter(iPt+1)
    effAccPrompt = hEffAccPrompt.GetBinContent(iPt+1)
    effAccFD = hEffAccFD.GetBinContent(iPt+1)
    effAccPromptUnc = hEffAccPrompt.GetBinError(iPt+1)
    effAccFDUnc = hEffAccFD.GetBinError(iPt+1)

    # ingredients for (prompt or FD) fraction computation
    corrYieldPrompt = hCorrYieldPrompt.GetBinContent(iPt+1)
    corrYieldFD = hCorrYieldFD.GetBinContent(iPt+1)
    covPromptPrompt = hCovPromptPrompt.GetBinContent(iPt+1)
    covPromptFD = hCovPromptFD.GetBinContent(iPt+1)
    covFDFD = hCovFDFD.GetBinContent(iPt+1)

    # prompt fraction
    fPrompt = effAccPrompt * corrYieldPrompt / (effAccPrompt * corrYieldPrompt + effAccFD * corrYieldFD)
    defPdeNP = (effAccPrompt * (effAccPrompt * corrYieldPrompt + effAccFD * corrYieldFD) - effAccPrompt**2
                * corrYieldPrompt) / (effAccPrompt * corrYieldPrompt + effAccFD * corrYieldFD)**2
    defPdeNF = - effAccFD * effAccPrompt * corrYieldPrompt / \
        (effAccPrompt * corrYieldPrompt + effAccFD * corrYieldFD)**2
    fPromptUnc = np.sqrt(defPdeNP**2 * covPromptPrompt + defPdeNF**2 * covFDFD + 2 * defPdeNP * defPdeNF * covPromptFD)

    # feed-down fraction
    fFD = effAccFD * corrYieldFD / (effAccPrompt * corrYieldPrompt + effAccFD * corrYieldFD)
    defFdeNF = (effAccFD * (effAccFD * corrYieldFD + effAccPrompt * corrYieldPrompt) - effAccFD**2
                * corrYieldFD) / (effAccPrompt * corrYieldPrompt + effAccFD * corrYieldFD)**2
    defFdeNP = - effAccFD * effAccPrompt * corrYieldFD / \
        (effAccPrompt * corrYieldPrompt + effAccFD * corrYieldFD)**2
    fFDUnc = np.sqrt(defFdeNF**2 * covFDFD + defFdeNP**2 * covPromptPrompt + 2 * defFdeNF * defFdeNP * covPromptFD)

    hPromptFrac.SetBinContent(iPt+1, fPrompt)
    hPromptFrac.SetBinError(iPt+1, fPromptUnc)
    hFDFrac.SetBinContent(iPt+1, fFD)
    hFDFrac.SetBinError(iPt+1, fFDUnc)

SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14)

legFrac = TLegend(0.6, 0.7, 0.9, 0.9)
legFrac.SetBorderSize(0)
legFrac.SetFillStyle(0)
legFrac.SetTextSize(0.045)
legFrac.AddEntry(hPromptFrac, 'Prompt', 'p')
legFrac.AddEntry(hFDFrac, 'Non-prompt', 'p')

legEff = legFrac.Clone('legEff')
legEff.SetY1(0.2)
legEff.SetY2(0.4)

cFrac = TCanvas('cFrac', '', 800, 800)
cFrac.DrawFrame(hPromptFrac.GetBinLowEdge(1), 0., ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); fraction')
hPromptFrac.Draw('same')
hFDFrac.Draw('same')
legFrac.Draw()
cFrac.Update()

cEff = TCanvas('cEff', '', 800, 800)
cEff.DrawFrame(hPromptFrac.GetBinLowEdge(1), 1.e-4, ptMax, 1., ';#it{p}_{T} (GeV/#it{c}); (Acc#times#font[152]{e})')
cEff.SetLogy()
hEffAccPrompt.Draw('same')
hEffAccFD.Draw('same')
legEff.Draw()
cEff.Update()

outFile = TFile(args.outFileName, 'recreate')
hEffAccPrompt.Write()
hEffAccFD.Write()
hPromptFrac.Write()
hFDFrac.Write()
cFrac.Write()
cEff.Write()
outFile.Close()

input('Press enter to exit')
