import sys
from ROOT import TCanvas, TFile, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kGreen, kFullCircle, kFullSquare, kFullDiamond # pylint: disable=import-error,no-name-in-module,unused-import
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error,no-name-in-module

inputdir = 'outputs/efficiency'
inputfilenames = ['Efficiency_Ds_010_central_18q.root', 'Efficiency_Ds_010_central_18r.root']
histonames = ['hEff', 'hEff']
colors = [kRed, kBlue]
markers = [kFullSquare, kFullCircle]
legendnames = ['LHC18q', 'LHC18r']
outputsuffix = 'LHC18q_LHC18r'

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, titlesize=0.045, labelsize=0.04)

hEffPrompt, hEffFD, hEffPromptRatio, hEffFDRatio = ([] for _ in range(4))

leg = TLegend(0.4, 0.38, 0.9, 0.53)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile, _ in enumerate(inputfilenames):
    inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
    hEffPrompt.append(inputfile.Get('%sPrompt' % histonames[iFile]))
    hEffFD.append(inputfile.Get('%sFD' % histonames[iFile]))
    hEffPrompt[iFile].SetDirectory(0)
    hEffFD[iFile].SetDirectory(0)
    SetObjectStyle(hEffPrompt[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hEffFD[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile], linestyle=1)
    hEffPrompt[iFile].GetYaxis().SetRangeUser(1e-3, 1.)
    hEffFD[iFile].GetYaxis().SetRangeUser(1e-3, 1.)
    hEffPromptRatio.append(hEffPrompt[iFile].Clone('hEffPromptRatio%d' % iFile))
    hEffPromptRatio[iFile].SetDirectory(0)
    hEffPromptRatio[iFile].Divide(hEffPrompt[iFile], hEffPrompt[0])
    hEffPromptRatio[iFile].GetYaxis().SetRangeUser(0.8, 1.2)
    hEffFDRatio.append(hEffFD[iFile].Clone('hEffFDRatio%d' % iFile))
    hEffFDRatio[iFile].SetDirectory(0)
    hEffFDRatio[iFile].Divide(hEffFD[iFile], hEffFD[0])
    hEffFDRatio[iFile].GetYaxis().SetRangeUser(0.8, 1.2)
    leg.AddEntry(hEffFD[iFile], legendnames[iFile], 'p')
    #for iBin in range(hEffPromptRatio[iFile].GetNbinsX()):
    #  hEffPromptRatio[iFile].SetBinError(iBin+1, 1.e-20)
    #  hEffFDRatio[iFile].SetBinError(iBin+1, 1.e-20)

ptmin = hEffPrompt[0].GetBinLowEdge(1)
ptmax = hEffPrompt[0].GetBinLowEdge(hEffPrompt[0].GetNbinsX())+hEffPrompt[0].GetBinWidth(hEffPrompt[0].GetNbinsX())

cPrompt = TCanvas('cPrompt', '', 1000, 500)
cPrompt.Divide(2, 1)
cPrompt.cd(1).DrawFrame(ptmin, 1.e-4, ptmax, 1., ';#it{p}_{T} (GeV/#it{c}); Prompt efficiency')
cPrompt.cd(1).SetLogy()
for iFile in range(len(inputfilenames)):
    hEffPrompt[iFile].Draw('same')
leg.Draw()
cPrompt.cd(2).DrawFrame(ptmin, 0.5, ptmax, 1.5, ';#it{p}_{T} (GeV/#it{c}); Prompt efficiency ratio')
for iFile in range(len(inputfilenames)):
    if iFile == 0:
        continue
    hEffPromptRatio[iFile].Draw('same')

cFD = TCanvas('cFD', '', 1000, 500)
cFD.Divide(2, 1)
cFD.cd(1).DrawFrame(ptmin, 1.e-4, ptmax, 1., ';#it{p}_{T} (GeV/#it{c}); Feed-down efficiency')
cFD.cd(1).SetLogy()
for iFile in range(len(inputfilenames)):
    hEffFD[iFile].Draw('same')
leg.Draw()
cFD.cd(2).DrawFrame(ptmin, 0.5, ptmax, 1.5, ';#it{p}_{T} (GeV/#it{c}); Feed-down efficiency ratio')
for iFile in range(len(inputfilenames)):
    if iFile == 0:
        continue
    hEffFDRatio[iFile].Draw('same')

cPrompt.SaveAs('%s/PromptEfficiencyComparison_%s.pdf' % (inputdir, outputsuffix))
cFD.SaveAs('%s/FDEfficiencyComparison_%s.pdf' % (inputdir, outputsuffix))

input('Press enter to exit')
