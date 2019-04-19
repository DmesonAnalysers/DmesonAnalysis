import sys
from ROOT import TCanvas, TFile, TH1F, TList, TLine, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import gStyle, kRed, kBlue, kBlack, kGreen, kAzure, kOrange, kFullCircle, kFullSquare, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCross # pylint: disable=import-error,no-name-in-module

#TO BE SET
############################################################
inFilePtWeightsName = 'ptweights/PtWeigths_LHC19c3a.root'
inFileEffNames = ['outputs/genptshape/Eff_Ds_010_Pythia.root', 'outputs/genptshape/Eff_Ds_010_FONLL.root', 'outputs/genptshape/Eff_Ds_010_FONLL_times_TAMU.root', 'outputs/genptshape/Eff_Ds_010_FONLL_times_PHSD.root', 'outputs/genptshape/Eff_Ds_010_FONLL_times_Catania.root', 'outputs/genptshape/Eff_Ds_010_FONLL_times_Gossiaux.root']

outFilePtWeightsName = 'outputs/genptshape/PtWeights_010_LHC19c3a.pdf'
outFileEffName = 'outputs/genptshape/SystPtWeights_010_LHC19c3a.pdf'

PtShapes = ['Pythia', 'FONLLcent', 'FONLLtimesTAMUcent', 'FONLLtimesPHSDcent', 'FONLLtimesCataniacent', 'FONLLtimesGossiauxcent'] 
colors = [kBlack, kBlue+1, kOrange+7, kGreen+2, kRed+2, kAzure+4]
markers = [kFullCircle, kFullSquare, kFullSquare, kFullTriangleUp, kFullTriangleDown, kFullCross]
############################################################

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadTopMargin(0.035)
gStyle.SetTitleSize(0.045,'xy')
gStyle.SetLabelSize(0.040,'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)

hPtShapes, hPtWeights, hEffPrompt, hEffPromptRatio = ([] for iList in range(4))

inFilePtWeights = TFile.Open(inFilePtWeightsName)
for iShape in PtShapes:
  hPtShapes.append(inFilePtWeights.Get('hPt%s' % iShape))
  if iShape is not 'Pythia':
    hPtWeights.append(inFilePtWeights.Get('hPtWeights%s' % iShape))

if len(inFileEffNames)!=0 and len(inFileEffNames)==len(PtShapes):
  for iFile, filename in enumerate(inFileEffNames):
    print(filename)
    inFileEff = TFile.Open(filename)
    hEffPrompt.append(inFileEff.Get('hEffPrompt'))
    hEffPrompt[iFile].SetName('hEffPrompt%s' % PtShapes[iFile])
    hEffPrompt[iFile].SetDirectory(0)

  for iFile, filename in enumerate(inFileEffNames):
    hEffPromptRatio.append(hEffPrompt[iFile].Clone('hEffPromptRatio%s' % PtShapes[iFile]))
    hEffPromptRatio[iFile].SetDirectory(0)
    hEffPromptRatio[iFile].Divide(hEffPrompt[iFile],hEffPrompt[0])
    for iBin in range(hEffPromptRatio[iFile].GetNbinsX()):
      hEffPromptRatio[iFile].SetBinError(iBin+1,1.e-20)

leg = TLegend(0.3,0.7,0.7,0.9)
leg.SetFillStyle(0)
leg.SetTextSize(0.04)
leg.AddEntry(hPtShapes[0],'Pythia','l')
leg.AddEntry(hPtShapes[1],'FONLL','l')
leg.AddEntry(hPtShapes[2],'FONLL #times TAMU','l')
leg.AddEntry(hPtShapes[3],'FONLL #times PHSD','l')
leg.AddEntry(hPtShapes[4],'FONLL #times Catania','l')
leg.AddEntry(hPtShapes[5],'FONLL #times MC@sHQ-EPOS2','l')

if len(inFileEffNames)!=0 and len(inFileEffNames)==len(PtShapes):
  legEff = TLegend(0.3,0.2,0.7,0.4)
  legEff.SetFillStyle(0)
  legEff.SetTextSize(0.04)
  legEff.AddEntry(hEffPrompt[0],'Pythia','p')
  legEff.AddEntry(hEffPrompt[1],'FONLL','p')
  legEff.AddEntry(hEffPrompt[2],'FONLL #times TAMU','p')
  legEff.AddEntry(hEffPrompt[3],'FONLL #times PHSD','p')
  legEff.AddEntry(hEffPrompt[4],'FONLL #times Catania','p')
  legEff.AddEntry(hEffPrompt[5],'FONLL #times MC@sHQ-EPOS2','p')

lineatone = TLine(0.,1.,36.,1.)
lineatone.SetLineWidth(2)
lineatone.SetLineStyle(9)
lineatone.SetLineColor(kBlack)

cShapes = TCanvas('cShapes','',1000,500)
cShapes.Divide(2,1)
cShapes.cd(1).DrawFrame(0,1.e-8,36,1.,';#it{p}_{T} GeV/#it{c}; d#it{N}/d#it{p}_{T} (a. u.)')
cShapes.cd(1).SetLogy()
for iShape in range(len(hPtShapes)):
  hPtShapes[iShape].SetLineWidth(2)
  hPtShapes[iShape].SetLineColor(colors[iShape])
  hPtShapes[iShape].Draw('chistsame')
leg.Draw()
cShapes.cd(2).DrawFrame(0,1.e-3,36,1.e+1,';#it{p}_{T} GeV/#it{c}; #it{p}_{T} weights')
cShapes.cd(2).SetLogy()
lineatone.Draw('same')
for iShape in range(len(hPtWeights)):
  hPtWeights[iShape].Draw('chistsame')
  hPtWeights[iShape].SetLineWidth(2)
  hPtWeights[iShape].SetLineColor(colors[iShape+1])

cShapes.SaveAs(outFilePtWeightsName)

if len(inFileEffNames)!=0 and len(inFileEffNames)==len(PtShapes):
  ptmin = hEffPrompt[0].GetBinLowEdge(1)
  ptmax = hEffPrompt[0].GetBinLowEdge(hEffPrompt[0].GetNbinsX())+hEffPrompt[0].GetBinWidth(hEffPrompt[0].GetNbinsX())
  cEff = TCanvas('cEff','',1000,500)
  cEff.Divide(2,1)
  cEff.cd(1).DrawFrame(ptmin,1.e-4,ptmax,1.,';#it{p}_{T} GeV/#it{c}; prompt efficiency')
  cEff.cd(1).SetLogy()
  for iEff in range(len(hEffPrompt)):
    hEffPrompt[iEff].SetLineWidth(2)
    hEffPrompt[iEff].SetLineColor(colors[iEff])
    hEffPrompt[iEff].SetMarkerColor(colors[iEff])
    hEffPrompt[iEff].SetMarkerSize(1)
    hEffPrompt[iEff].Draw('same')
  legEff.Draw()
  cEff.cd(2).DrawFrame(ptmin,0.8,ptmax,1.2,';#it{p}_{T} GeV/#it{c}; prompt efficiency ratio')
  for iEff in range(len(hEffPromptRatio)):
    hEffPromptRatio[iEff].SetLineWidth(2)
    hEffPromptRatio[iEff].SetLineColor(colors[iEff])
    hEffPromptRatio[iEff].SetMarkerColor(colors[iEff])
    hEffPromptRatio[iEff].SetMarkerSize(1)
    hEffPromptRatio[iEff].Draw('same')
  cEff.SaveAs(outFileEffName)

raw_input()