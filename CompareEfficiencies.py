from ROOT import TCanvas, TFile, TLegend, TLine
from ROOT import gStyle, kRed, kBlack, kBlue, kFullCircle, kFullSquare, kFullDiamond 
import math

inputdir = 'outputs'
inputfilenames = [ '2015results/AccEffDs_3050_AOD198.root', '3_24bin_merge/eff/EffAcc_Ds_3050_2015cuts_Improver.root']
histonames = [ 'hEffAcc', 'hAccEff' ]
colors = [ kRed, kBlue ]
markers = [ kFullSquare, kFullCircle]
legendnames = [ 'LHC16i2b (2015 cuts)', 'LHC19c3b(2) (2015 cuts)']
outputsuffix = 'LHC15o_LHC18qr'

hEffPrompt, hEffFD, hEffPromptRatio, hEffFDRatio = ([] for iList in range(4))

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadTopMargin(0.05)
gStyle.SetTitleSize(0.045,'xy')
gStyle.SetLabelSize(0.040,'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)
gStyle.SetOptStat(0)

leg = TLegend(0.2,0.78,0.8,0.93)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile in range(len(inputfilenames)):
  inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
  hEffPrompt.append(inputfile.Get('%sPrompt' % histonames[iFile]))
  hEffFD.append(inputfile.Get('%sFD' % histonames[iFile]))
  hEffPrompt[iFile].SetDirectory(0)
  hEffPrompt[iFile].SetLineColor(colors[iFile])
  hEffPrompt[iFile].SetLineWidth(2)
  hEffPrompt[iFile].SetLineStyle(1)
  hEffPrompt[iFile].SetMarkerColor(colors[iFile])
  hEffPrompt[iFile].SetMarkerStyle(markers[iFile])
  hEffPrompt[iFile].SetMarkerSize(1)
  hEffPrompt[iFile].GetYaxis().SetRangeUser(1e-3,1.)
  hEffFD[iFile].SetDirectory(0)
  hEffFD[iFile].SetLineColor(colors[iFile])
  hEffFD[iFile].SetLineWidth(2)
  hEffFD[iFile].SetLineStyle(1)
  hEffFD[iFile].SetMarkerColor(colors[iFile])
  hEffFD[iFile].SetMarkerStyle(markers[iFile])
  hEffFD[iFile].SetMarkerSize(1)
  hEffFD[iFile].GetYaxis().SetRangeUser(1e-3,1.)
  hEffPromptRatio.append(hEffPrompt[iFile].Clone('hEffPromptRatio%d' % iFile))
  hEffPromptRatio[iFile].SetDirectory(0)
  hEffPromptRatio[iFile].Divide(hEffPrompt[iFile],hEffPrompt[0])
  hEffPromptRatio[iFile].GetYaxis().SetRangeUser(0.8, 1.2)
  hEffFDRatio.append(hEffFD[iFile].Clone('hEffFDRatio%d' % iFile))
  hEffFDRatio[iFile].SetDirectory(0)
  hEffFDRatio[iFile].Divide(hEffFD[iFile],hEffFD[0])
  hEffFDRatio[iFile].GetYaxis().SetRangeUser(0.8, 1.2)
  leg.AddEntry(hEffFD[iFile],legendnames[iFile],'p')
  #for iBin in range(hEffPromptRatio[iFile].GetNbinsX()):
  #  hEffPromptRatio[iFile].SetBinError(iBin+1,1.e-20)
  #  hEffFDRatio[iFile].SetBinError(iBin+1,1.e-20)

cPrompt = TCanvas('cPrompt','',1000,500)
cPrompt.Divide(2,1)
cPrompt.cd(1).DrawFrame(4,1.e-4,16,1.,';#it{p}_{T} (GeV/#it{c}); Prompt efficiency')
cPrompt.cd(1).SetLogy()
for iFile in range(len(inputfilenames)):
  hEffPrompt[iFile].Draw('same')
leg.Draw()
cPrompt.cd(2).DrawFrame(4,0.5,16,1.5,';#it{p}_{T} (GeV/#it{c}); Prompt efficiency ratio')
for iFile in range(len(inputfilenames)):
  if iFile==0: 
    continue
  hEffPromptRatio[iFile].Draw('same')

cFD = TCanvas('cFD','',1000,500)
cFD.Divide(2,1)
cFD.cd(1).DrawFrame(4,1.e-4,16,1.,';#it{p}_{T} (GeV/#it{c}); Feed-down efficiency')
cFD.cd(1).SetLogy()
for iFile in range(len(inputfilenames)):
  hEffFD[iFile].Draw('same')
leg.Draw()
cFD.cd(2).DrawFrame(4,0.5,16,1.5,';#it{p}_{T} (GeV/#it{c}); Feed-down efficiency ratio')
for iFile in range(len(inputfilenames)):
  if iFile==0: 
    continue
  hEffFDRatio[iFile].Draw('same')

cPrompt.SaveAs('%s/PromptEfficiencyComparison_%s.pdf' % (inputdir,outputsuffix))
cFD.SaveAs('%s/FDEfficiencyComparison_%s.pdf' % (inputdir,outputsuffix))

input()