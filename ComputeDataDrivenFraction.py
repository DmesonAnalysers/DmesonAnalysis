'''
python script for the computation of the fractions of prompt and feed-down D for a given cut set
run: python ComputeDataDrivenFraction.py effAccFile.root fracFile.root outFile.root
'''

import argparse
from ROOT import TFile, TCanvas, TLegend  # pylint: disable=import-error,no-name-in-module
from utils.AnalysisUtils import GetPromptFDFractionCutSet
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
hEffAccPrompt = effAccFile.Get('hEffPrompt')
hEffAccFD = effAccFile.Get('hEffFD')

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
hPromptFracCorr = hEffAccPrompt.Clone('hPromptFracCorr')
hFDFracCorr = hEffAccFD.Clone('hFDFracCorr')
hPromptFracCorr.SetTitle(';#it{p}_{T} (GeV/#it{c}); corrected #it{f}_{prompt}')
hFDFracCorr.SetTitle(';#it{p}_{T} (GeV/#it{c}); corrected #it{f}_{FD}')

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

    # prompt and FD and fractions
    fracPromptFD, uncFracPromptFD = GetPromptFDFractionCutSet(effAccPrompt, effAccFD, corrYieldPrompt, corrYieldFD,
                                                              covPromptPrompt, covFDFD, covPromptFD)

    fracPromptFDcorr, uncFracPromptFDcorr = GetPromptFDFractionCutSet(1., 1., corrYieldPrompt, corrYieldFD,
                                                                      covPromptPrompt, covFDFD, covPromptFD)

    hPromptFrac.SetBinContent(iPt+1, fracPromptFD[0])
    hPromptFrac.SetBinError(iPt+1, uncFracPromptFD[0])
    hFDFrac.SetBinContent(iPt+1, fracPromptFD[1])
    hFDFrac.SetBinError(iPt+1, uncFracPromptFD[1])
    hPromptFracCorr.SetBinContent(iPt+1, fracPromptFDcorr[0])
    hPromptFracCorr.SetBinError(iPt+1, uncFracPromptFDcorr[0])
    hFDFracCorr.SetBinContent(iPt+1, fracPromptFDcorr[1])
    hFDFracCorr.SetBinError(iPt+1, uncFracPromptFDcorr[1])

SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14)

legFrac = TLegend(0.2, 0.84, 0.4, 0.94)
legFrac.SetBorderSize(0)
legFrac.SetFillStyle(0)
legFrac.SetTextSize(0.045)
legFrac.AddEntry(hPromptFrac, 'Prompt', 'p')
legFrac.AddEntry(hFDFrac, 'Non-prompt', 'p')

legEff = legFrac.Clone('legEff')
legEff.SetY1(0.2)
legEff.SetY2(0.4)

ptMin = hPromptFrac.GetBinLowEdge(1)
cFrac = TCanvas('cFrac', '', 800, 800)
cFrac.DrawFrame(ptMin, 0., ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); fraction')
hPromptFrac.Draw('same')
hFDFrac.Draw('same')
legFrac.Draw()
cFrac.Update()

cFracCorrFrac = TCanvas('cFracCorrFrac', '', 800, 800)
cFracCorrFrac.DrawFrame(ptMin, 0., ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); corrected fraction')
hPromptFracCorr.Draw('same')
hFDFracCorr.Draw('same')
legFrac.Draw()
cFracCorrFrac.Update()

cEff = TCanvas('cEff', '', 800, 800)
cEff.DrawFrame(ptMin, 1.e-4, ptMax, 1., ';#it{p}_{T} (GeV/#it{c}); (Acc#times#font[152]{e})')
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
hPromptFracCorr.Write()
hFDFracCorr.Write()
cFrac.Write()
cFracCorrFrac.Write()
cEff.Write()
outFile.Close()

input('Press enter to exit')
