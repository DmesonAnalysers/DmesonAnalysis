from ROOT import TCanvas, TFile, TLegend
from ROOT import gStyle, kRed, kBlack, kBlue, kGreen, kFullCircle, kFullSquare, kFullDiamond 
import math, six

inputdir = 'outputs/'
inputfilenames = [ '2015results/RawYieldDs_3050.root', '3_24bin_merge/raw_yields/RawYieldsDs_3050_2015cuts.root']
outputdir = 'outputs/3_24bin_merge/raw_yields/2015_comp/'
outputsuffix = 'LHC15o_LHC18qr'
signalhistonames = [ 'hSignal', 'hRawYields']
bkghistonames = [ 'hBackground', 'hRawYieldsBkg' ]
SoverBhistonames = [ 'hSignal', 'hRawYieldsSoverB' ]
signifhistonames = [ 'hSignificance', 'hRawYieldsSignificance' ]
evhistonames = [ 'hNEvents', 'hEvForNorm' ]
colors = [ kRed, kBlue ]
markers = [ kFullSquare, kFullCircle]
legendnames = [ 'LHC15o (2015 cuts)', 'LHC18qr (2015 cuts)']

# inputfilenames = [ 'RawYieldsDs_010_centralcuts_LHC18qr.root', 'RawYieldsDs_010_centralcuts_LHC18qr_old.root' ]
# signalhistonames = [ 'hRawYields', 'hRawYields' ]
# bkghistonames = [ 'hRawYieldsBkg', 'hRawYieldsBkg' ]
# SoverBhistonames = [ 'hRawYieldsSoverB', 'hRawYieldsSoverB' ]
# signifhistonames = [ 'hRawYieldsSignificance', 'hRawYieldsSignificance' ]
# evhistonames = [ 'hEvForNorm', 'hEvForNorm' ]
# colors = [ kRed, kBlack, kBlue ]
# markers = [ kFullCircle, kFullSquare, kFullDiamond ]
# legendnames = [ 'LHC18qr (PbPb splines, 88M events)', 'LHC18qr (pp splines, 40M events)' ]
# outputsuffix = '2018_diffsplines_2018cuts'

hSignal, hRatioSignal, hBackground, hRatioBkg, hSoverB, hSignif, hRatioSignif, hEv = ([] for iList in range(8))

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadTopMargin(0.05)
gStyle.SetTitleSize(0.045,'xy')
gStyle.SetLabelSize(0.040,'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)
gStyle.SetOptStat(0)

leg = TLegend(0.5,0.73,0.8,0.93)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile in range(len(inputfilenames)):
  inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
  hSignal.append(inputfile.Get(signalhistonames[iFile]))
  hBackground.append(inputfile.Get(bkghistonames[iFile]))
  hSignif.append(inputfile.Get(signifhistonames[iFile]))
  hSoverB.append(inputfile.Get(SoverBhistonames[iFile]))
  hEv.append(inputfile.Get(evhistonames[iFile]))
  hSignal[iFile].SetDirectory(0)
  hBackground[iFile].SetDirectory(0)
  hSoverB[iFile].SetDirectory(0)
  hSignif[iFile].SetDirectory(0)
  hEv[iFile].SetDirectory(0)
  hSignal[iFile].SetLineColor(colors[iFile])
  hSignal[iFile].SetLineWidth(2)
  hSignal[iFile].SetMarkerColor(colors[iFile])
  hSignal[iFile].SetMarkerStyle(markers[iFile])
  hSignal[iFile].Scale(1./hEv[iFile].GetBinContent(1))
  hBackground[iFile].SetLineColor(colors[iFile])
  hBackground[iFile].SetLineWidth(2)
  hBackground[iFile].SetMarkerColor(colors[iFile])
  hBackground[iFile].SetMarkerStyle(markers[iFile])
  hBackground[iFile].Scale(1./hEv[iFile].GetBinContent(1))
  hSoverB[iFile].SetLineColor(colors[iFile])
  hSoverB[iFile].SetLineWidth(2)
  hSoverB[iFile].SetMarkerColor(colors[iFile])
  hSoverB[iFile].SetMarkerStyle(markers[iFile])
  hSignif[iFile].SetLineColor(colors[iFile])
  hSignif[iFile].SetLineWidth(2)
  hSignif[iFile].SetMarkerColor(colors[iFile])
  hSignif[iFile].SetMarkerStyle(markers[iFile])
  hSignif[iFile].Scale(1./math.sqrt(hEv[iFile].GetBinContent(1)))
  leg.AddEntry(hSignal[iFile],legendnames[iFile],'lp')
  hRatioSignal.append(hSignal[iFile].Clone("hRatioSignal%d" % iFile))
  hRatioSignal[iFile].SetDirectory(0)
  hRatioSignal[iFile].Divide(hSignal[iFile],hSignal[0],1.,1.,"")
  #for iBin in range(hRatioSignal[iFile].GetNbinsX()):
  #  hRatioSignal[iFile].SetBinError(iBin+1,1.e-20)
  hRatioBkg.append(hBackground[iFile].Clone("hRatioBkg%d" % iFile))
  hRatioBkg[iFile].SetDirectory(0)
  hRatioBkg[iFile].Divide(hBackground[iFile],hBackground[0],1.,1.,"")
  #for iBin in range(hRatioBkg[iFile].GetNbinsX()):
  #  hRatioBkg[iFile].SetBinError(iBin+1,1.e-20)
  hRatioSignif.append(hSignif[iFile].Clone("hRatioSignif%d" % iFile))
  hRatioSignif[iFile].SetDirectory(0)
  hRatioSignif[iFile].Divide(hSignif[iFile],hSignif[0],1.,1.,"")
  #for iBin in range(hRatioSignif[iFile].GetNbinsX()):
  #  hRatioSignif[iFile].SetBinError(iBin+1,1.e-20)

PtMin = hSignal[0].GetBinLowEdge(1)
PtMax = hSignal[0].GetBinLowEdge(hSignal[0].GetNbinsX())+hSignal[0].GetBinWidth(hSignal[0].GetNbinsX())

cSignal = TCanvas('cSignal','',1000,500)
cSignal.Divide(2,1)
cSignal.cd(1).DrawFrame(PtMin,0.,PtMax,3.e-5,';#it{p}_{T} (GeV/#it{c}); raw yields / #it{N}_{events}')
for iFile in range(len(inputfilenames)):
  hSignal[iFile].Draw('same')
leg.Draw()
cSignal.cd(2).DrawFrame(PtMin,0.,PtMax,2.,';#it{p}_{T} (GeV/#it{c}); ratio of raw yields / #it{N}_{events}')
for iFile in range(len(inputfilenames)):
  if iFile==0:
    continue
  hRatioSignal[iFile].Draw('same')

cRatioSignal = TCanvas('cRatioSignal','',800,800)
cRatioSignal.DrawFrame(PtMin,0.5,PtMax,2.,';#it{p}_{T} (GeV/#it{c}); ratio of raw yields / #it{N}_{events}')
for iFile in range(len(inputfilenames)):
  hRatioSignal[iFile].Draw('same')
leg.Draw()

cBkg = TCanvas('cBkg','',1000,500)
cBkg.Divide(2,1)
cBkg.cd(1).DrawFrame(PtMin,1.e-7,PtMax,1.e-2,';#it{p}_{T} (GeV/#it{c}); background(3#sigma) / #it{N}_{events}')
cBkg.cd(1).SetLogy()
for iFile in range(len(inputfilenames)):
  hBackground[iFile].Draw('same')
leg.Draw()
cBkg.cd(2).DrawFrame(PtMin,0.,PtMax,2.,';#it{p}_{T} (GeV/#it{c}); ratio of background(3#sigma) / #it{N}_{events}')
for iFile in range(len(inputfilenames)):
  if iFile==0:
    continue
  hRatioBkg[iFile].Draw('same')

cSoverB = TCanvas('cSoverB','',800,800)
cSoverB.DrawFrame(PtMin,0.01,PtMax,10.,';#it{p}_{T} (GeV/#it{c}); S/B(3#sigma)')
cSoverB.SetLogy()
for iFile in range(len(inputfilenames)):
  hSoverB[iFile].Draw('same')
leg.Draw()

cSignificance = TCanvas('cSignificance','',1000,500)
cSignificance.Divide(2,1)
cSignificance.cd(1).DrawFrame(PtMin,0.,PtMax,3.e-3,';#it{p}_{T} (GeV/#it{c}); significance / #sqrt{#it{N}_{events}}')
for iFile in range(len(inputfilenames)):
  hSignif[iFile].Draw('same')
leg.Draw()
cSignificance.cd(2).DrawFrame(PtMin,0.,PtMax,2.,';#it{p}_{T} (GeV/#it{c}); ratio of significance / #sqrt{#it{N}_{events}}')
for iFile in range(len(inputfilenames)):
  if iFile==0:
    continue
  hRatioSignif[iFile].Draw('same')

cSignal.SaveAs('%s/SignalPerEventComparison_%s.eps' % (outputdir,outputsuffix))
cRatioSignal.SaveAs('%s/SignalPerEventRatio_%s.eps' % (outputdir,outputsuffix))
cBkg.SaveAs('%s/BkgPerEventComparison_%s.eps' % (outputdir,outputsuffix))
cSoverB.SaveAs('%s/SoverB_%s.eps' % (outputdir,outputsuffix))
cSignificance.SaveAs('%s/SignificancePerEventComparison_%s.eps' % (outputdir,outputsuffix))

if six.PY2:
    raw_input('Press enter to exit')
elif six.PY3:
    input('Press enter to exit')