from ROOT import TCanvas, TFile, TLegend, TLine  # pylint: disable=import-error,no-name-in-module
from ROOT import gStyle, kRed, kBlack, kBlue, kOrange, kGreen, kFullCircle, kOpenCircle, kFullSquare, kFullDiamond, kFullCross, kOpenCross  # pylint: disable=import-error,no-name-in-module
import math
import six

inputdir = 'DplusNclsTPC/'

inputfilenames = ['HFPtSpectrumRaaDplus_010_Ncls50_pt2_50.root',
                  'HFPtSpectrumRaaDplus_010_Ncls70_pt2_50.root']
histonames = ['hRABvsPt', 'hRABvsPt']
graphnames = ['gRAB_GlobalSystematics', 'gRAB_GlobalSystematics']
colors = [kRed, kBlue]
linecolors = [kRed, kBlue]
markers = [kFullSquare, kFullCircle]
legendnames = [
    'N_{cls}^{TPC} > 50 (filterbit 4)', 'N_{cls}^{TPC} > 70']
outputsuffix = 'DplusNTPCcls'

hRaa, gRaa, hRaaRatio = ([] for i in range(3))

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.12)
gStyle.SetPadTopMargin(0.05)
gStyle.SetTitleSize(0.045, 'xy')
gStyle.SetLabelSize(0.040, 'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)
gStyle.SetOptStat(0)

leg = TLegend(0.3, 0.68, 0.75, 0.83)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile in range(len(inputfilenames)):
    inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
    hRaa.append(inputfile.Get(histonames[iFile]))
    gRaa.append(inputfile.Get(graphnames[iFile]))
    hRaa[iFile].SetDirectory(0)
    hRaa[iFile].SetLineColor(linecolors[iFile])
    hRaa[iFile].SetLineWidth(2)
    hRaa[iFile].SetLineStyle(1)
    hRaa[iFile].SetMarkerSize(1.5)
    hRaa[iFile].SetMarkerColor(colors[iFile])
    hRaa[iFile].SetMarkerStyle(markers[iFile])
    gRaa[iFile].SetLineColor(colors[iFile])
    gRaa[iFile].SetLineWidth(2)
    gRaa[iFile].SetFillStyle(0)
    leg.AddEntry(hRaa[iFile], legendnames[iFile], 'p')
    hRaaRatio.append(hRaa[iFile].Clone("hRaa%d" % iFile))
    hRaaRatio[iFile].SetDirectory(0)
    hRaaRatio[iFile].Divide(hRaa[iFile], hRaa[0], 1, 1, "B")

ptmin = hRaa[0].GetBinLowEdge(1)
ptmax = hRaa[0].GetBinLowEdge(hRaa[0].GetNbinsX()) + \
    hRaa[0].GetBinWidth(hRaa[0].GetNbinsX())

lineatone = TLine(ptmin, 1., ptmax, 1.)
lineatone.SetLineWidth(1)
lineatone.SetLineColor(kBlack)
lineatone.SetLineStyle(9)

cRaa = TCanvas('cRaa', '', 1000, 500)
cRaa.Divide(2, 1)
cRaa.cd(1).DrawFrame(ptmin, 0., ptmax, 2.,
                     ';#it{p}_{T} (GeV/#it{c}); #it{R}_{AA}')
# cRaa.SetLogx()
lineatone.Draw('same')
for iFile in range(len(inputfilenames)):
    # gRaa[iFile].Draw('2')
    hRaa[iFile].Draw('same')
leg.Draw()
cRaa.cd(2).DrawFrame(ptmin, 0.9, ptmax, 1.1,
                     ';#it{p}_{T} (GeV/#it{c}); ratio #it{R}_{AA}')
for iFile in range(len(inputfilenames)):
    if iFile == 0:
        continue
    # gRaa[iFile].Draw('2')
    hRaaRatio[iFile].Draw('same')

cRaa.SaveAs('%s/RaaComparison_%s.pdf' % (inputdir, outputsuffix))

if six.PY2:
    raw_input('Press enter to exit')
elif six.PY3:
    input('Press enter to exit')
