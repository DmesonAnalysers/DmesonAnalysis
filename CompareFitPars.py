from ROOT import TCanvas, TFile, TLegend, TLine # pylint: disable=import-error,no-name-in-module
from ROOT import gStyle, kRed, kBlack, kBlue, kOrange, kGreen, kFullCircle, kOpenCircle, kFullSquare, kFullDiamond, kFullCross, kOpenCross # pylint: disable=import-error,no-name-in-module
import math
import six

inputdir = 'outputs/improver'
inputfilenames = [ 'RawYieldsDs_MC_010_improver.root', 'RawYieldsDs_MC_010_noimprover.root' ]
colors = [ kRed, kBlue, kBlack ]
markers = [ kFullCircle, kFullSquare, kFullDiamond ]
legendnames = [ 'w/ improver', 'w/o improver' ]
outputsuffix = 'Improver'

hMean, hSigma = ([] for iList in range(2))

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.18)
gStyle.SetPadTopMargin(0.05)
gStyle.SetTitleSize(0.045,'xy')
gStyle.SetLabelSize(0.040,'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)
gStyle.SetOptStat(0)

legSigma = TLegend(0.2,0.78,0.8,0.93)
legSigma.SetFillStyle(0)
legSigma.SetBorderSize(0)
legSigma.SetTextSize(0.04)

legMean = TLegend(0.5,0.73,0.8,0.93)
legMean.SetFillStyle(0)
legMean.SetBorderSize(0)
legMean.SetTextSize(0.04)

lineMass = TLine(3.,1.96850,36,1.96850)
lineMass.SetLineWidth(2)
lineMass.SetLineColor(kBlack)
lineMass.SetLineStyle(9)

for iFile in range(len(inputfilenames)):
  inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
  hMean.append(inputfile.Get('hRawYieldsMean'))
  hSigma.append(inputfile.Get('hRawYieldsSigma'))
  hMean[iFile].SetDirectory(0)
  hMean[iFile].SetLineColor(colors[iFile])
  hMean[iFile].SetLineWidth(2)
  hMean[iFile].SetMarkerColor(colors[iFile])
  hMean[iFile].SetMarkerStyle(markers[iFile])
  hSigma[iFile].SetDirectory(0)
  hSigma[iFile].SetLineColor(colors[iFile])
  hSigma[iFile].SetLineWidth(2)
  hSigma[iFile].SetMarkerColor(colors[iFile])
  hSigma[iFile].SetMarkerStyle(markers[iFile])
  legMean.AddEntry(hSigma[iFile],legendnames[iFile],'p')
  legSigma.AddEntry(hSigma[iFile],legendnames[iFile],'p')

legMean.AddEntry(lineMass,"PDG",'l')

cMean = TCanvas('cMean','',800,800)
cMean.DrawFrame(3,1.9641,36,1.9759,';#it{p}_{T} (GeV/#it{c}); peak mean (GeV/#it{c}^{2})')
lineMass.Draw("same")
for iFile in range(len(inputfilenames)):
  hMean[iFile].Draw('same')
legMean.Draw()

cSigma = TCanvas('cSigma','',800,800)
cSigma.DrawFrame(3,0.,36,0.025,';#it{p}_{T} (GeV/#it{c}); peak width (GeV/#it{c}^{2})')
for iFile in range(len(inputfilenames)):
  hSigma[iFile].Draw('same')
legSigma.Draw()

cMean.SaveAs('%s/Mean_%s.pdf' % (inputdir,outputsuffix))
cSigma.SaveAs('%s/Sigma_%s.pdf' % (inputdir,outputsuffix))

if six.PY2:
    raw_input('Press enter to exit')
elif six.PY3:
    input('Press enter to exit')