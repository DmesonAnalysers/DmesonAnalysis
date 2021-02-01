import sys
from os.path import join
from ROOT import TCanvas, TFile, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kAzure, kOrange, kFullCircle, kFullSquare, kFullDiamond # pylint: disable=import-error,no-name-in-module,unused-import
sys.path.append('..')
from utils.AnalysisUtils import ComputeRatioDiffBins #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

inDir = 'inputdir'
inFileNames = ['file1.root', 'file2.root']
histoNames = ['hAccEff', 'hAccEff']
colors = [kRed+1, kAzure+4]
markers = [kFullSquare, kFullCircle]
legNames = ['title1', 'title2']
outDir = 'outputdir'
outSuffix = 'suffix'
showUnc = True

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
    hEffPrompt[iFile].GetYaxis().SetRangeUser(1e-4, 1.)
    hEffFD[iFile].GetYaxis().SetRangeUser(1e-4, 1.)
    hEffPromptRatio.append(ComputeRatioDiffBins(hEffPrompt[iFile], hEffPrompt[0], 'B'))
    hEffPromptRatio[iFile].SetDirectory(0)
    hEffPromptRatio[iFile].SetName(f'hEffPromptRatio{iFile}')
    hEffFDRatio.append(ComputeRatioDiffBins(hEffFD[iFile], hEffFD[0], 'B'))
    hEffFDRatio[iFile].SetDirectory(0)
    hEffFDRatio[iFile].SetName(f'hEffFDRatio{iFile}')
    SetObjectStyle(hEffPrompt[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hEffFD[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile], linestyle=1)
    SetObjectStyle(hEffPromptRatio[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile], linestyle=1)
    SetObjectStyle(hEffFDRatio[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    leg.AddEntry(hEffFD[iFile], legNames[iFile], 'p')
    if not showUnc:
        for iBin in range(hEffPromptRatio[iFile].GetNbinsX()):
            hEffPromptRatio[iFile].SetBinError(iBin+1, 1.e-20)
            hEffFDRatio[iFile].SetBinError(iBin+1, 1.e-20)

PtMin = hEffPrompt[0].GetBinLowEdge(1)
PtMax = hEffPrompt[0].GetBinLowEdge(hEffPrompt[0].GetNbinsX())+hEffPrompt[0].GetBinWidth(hEffPrompt[0].GetNbinsX())
for histo in hEffPrompt:
    if histo.GetBinLowEdge(1) < PtMin:
        PtMin = histo.GetBinLowEdge(1)
    if histo.GetBinLowEdge(histo.GetNbinsX())+histo.GetBinWidth(histo.GetNbinsX()) > PtMax:
        PtMax = histo.GetBinLowEdge(histo.GetNbinsX())+histo.GetBinWidth(histo.GetNbinsX())

cPrompt = TCanvas('cPrompt', '', 1000, 500)
cPrompt.Divide(2, 1)
cPrompt.cd(1).DrawFrame(PtMin, hEffPrompt[0].GetMinimum()/5, PtMax, 1.,
                        ';#it{p}_{T} (GeV/#it{c}); Prompt (Acc #times #epsilon)')
cPrompt.cd(1).SetLogy()
for iFile in range(len(inFileNames)):
    hEffPrompt[iFile].Draw('same')
leg.Draw()
cPrompt.cd(2).DrawFrame(PtMin, hEffPromptRatio[1].GetMinimum()/2, PtMax, hEffPromptRatio[1].GetMaximum()*1.5,
                        ';#it{p}_{T} (GeV/#it{c}); Prompt (Acc #times #epsilon) ratio')
for iFile in range(len(inFileNames)):
    if iFile == 0:
        continue
    hEffPromptRatio[iFile].Draw('same')

cFD = TCanvas('cFD', '', 1000, 500)
cFD.Divide(2, 1)
cFD.cd(1).DrawFrame(PtMin, hEffFD[0].GetMinimum()/5, PtMax, 1.,
                    ';#it{p}_{T} (GeV/#it{c}); Feed-down (Acc #times #epsilon)')
cFD.cd(1).SetLogy()
for iFile in range(len(inFileNames)):
    hEffFD[iFile].Draw('same')
leg.Draw()
cFD.cd(2).DrawFrame(PtMin, hEffFDRatio[1].GetMinimum()/2, PtMax, hEffFDRatio[1].GetMaximum()*2,
                    ';#it{p}_{T} (GeV/#it{c}); Feed-down (Acc #times #epsilon) ratio')
for iFile in range(len(inFileNames)):
    if iFile == 0:
        continue
    hEffFDRatio[iFile].Draw('same')

cPrompt.SaveAs(f'{outDir}/PromptEfficiencyComparison_{outSuffix}.pdf')
cFD.SaveAs(f'{outDir}/FDEfficiencyComparison_{outSuffix}.pdf')

input('Press enter to exit')
