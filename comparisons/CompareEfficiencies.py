import sys
from os.path import join
from ROOT import TCanvas, TFile, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kAzure, kFullCircle, kFullSquare, kFullDiamond # pylint: disable=import-error,no-name-in-module,unused-import
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error,no-name-in-module

inDir = '~/Desktop/Analyses/pp5TeV/Ds_wML_mult/outputs/100320'
inFileNames = ['eff/Efficiency_Ds_FDen_pt2_12_conservative.root',
               'systematics/eff_noImprover/Efficiency_Ds_FDen_conservative_pt2_12_noImprover.root']
histoNames = ['hEff', 'hEff']
colors = [kRed, kBlue]
markers = [kFullSquare, kFullCircle]
legNames = ['Improver', 'No Improver (same ML models)']
outDir = '~/Desktop/Analyses/pp5TeV/Ds_wML_mult/outputs/100320/systematics/eff_noImprover'
outSuffix = 'Improver_sameModels'
showUnc = False

SetGlobalStyle(padleftmargin=0.16, padtopmargin=0.05, padbottommargin=0.14,
               titleoffsety=1.5, titlesize=0.05, labelsize=0.045)

hEffPrompt, hEffFD, hEffPromptRatio, hEffFDRatio = ([] for _ in range(4))

leg = TLegend(0.2, 0.18, 0.9, 0.33)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile, inFileName in enumerate(inFileNames):
    inFileName = join(inDir, inFileName)
    inFile = TFile.Open(inFileName)
    hEffPrompt.append(inFile.Get(f'{histoNames[iFile]}Prompt'))
    hEffFD.append(inFile.Get(f'{histoNames[iFile]}FD'))
    hEffPrompt[iFile].SetDirectory(0)
    hEffFD[iFile].SetDirectory(0)
    SetObjectStyle(hEffPrompt[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hEffFD[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile], linestyle=1)
    hEffPrompt[iFile].GetYaxis().SetRangeUser(1e-3, 1.)
    hEffFD[iFile].GetYaxis().SetRangeUser(1e-3, 1.)
    hEffPromptRatio.append(hEffPrompt[iFile].Clone(f'hEffPromptRatio{iFile}'))
    hEffPromptRatio[iFile].SetDirectory(0)
    hEffPromptRatio[iFile].Divide(hEffPrompt[iFile], hEffPrompt[0])
    hEffPromptRatio[iFile].GetYaxis().SetRangeUser(0.8, 1.2)
    hEffFDRatio.append(hEffFD[iFile].Clone(f'hEffFDRatio{iFile}'))
    hEffFDRatio[iFile].SetDirectory(0)
    hEffFDRatio[iFile].Divide(hEffFD[iFile], hEffFD[0])
    hEffFDRatio[iFile].GetYaxis().SetRangeUser(0.8, 1.2)
    leg.AddEntry(hEffFD[iFile], legNames[iFile], 'p')
    if not showUnc:
        for iBin in range(hEffPromptRatio[iFile].GetNbinsX()):
            hEffPromptRatio[iFile].SetBinError(iBin+1, 1.e-20)
            hEffFDRatio[iFile].SetBinError(iBin+1, 1.e-20)

ptmin = hEffPrompt[0].GetBinLowEdge(1)
ptmax = hEffPrompt[0].GetBinLowEdge(hEffPrompt[0].GetNbinsX())+hEffPrompt[0].GetBinWidth(hEffPrompt[0].GetNbinsX())

cPrompt = TCanvas('cPrompt', '', 1000, 500)
cPrompt.Divide(2, 1)
cPrompt.cd(1).DrawFrame(ptmin, 1.e-4, ptmax, 1., ';#it{p}_{T} (GeV/#it{c}); Prompt efficiency')
cPrompt.cd(1).SetLogy()
for iFile in range(len(inFileNames)):
    hEffPrompt[iFile].Draw('same')
leg.Draw()
cPrompt.cd(2).DrawFrame(ptmin, 0.5, ptmax, 1.5, ';#it{p}_{T} (GeV/#it{c}); Prompt efficiency ratio')
for iFile in range(len(inFileNames)):
    if iFile == 0:
        continue
    hEffPromptRatio[iFile].Draw('same')

cFD = TCanvas('cFD', '', 1000, 500)
cFD.Divide(2, 1)
cFD.cd(1).DrawFrame(ptmin, 1.e-4, ptmax, 1., ';#it{p}_{T} (GeV/#it{c}); Feed-down efficiency')
cFD.cd(1).SetLogy()
for iFile in range(len(inFileNames)):
    hEffFD[iFile].Draw('same')
leg.Draw()
cFD.cd(2).DrawFrame(ptmin, 0.5, ptmax, 1.5, ';#it{p}_{T} (GeV/#it{c}); Feed-down efficiency ratio')
for iFile in range(len(inFileNames)):
    if iFile == 0:
        continue
    hEffFDRatio[iFile].Draw('same')

cPrompt.SaveAs(f'{outDir}/PromptEfficiencyComparison_{outSuffix}.pdf')
cFD.SaveAs(f'{outDir}/FDEfficiencyComparison_{outSuffix}.pdf')

input('Press enter to exit')
