import sys
from ROOT import TCanvas, TFile, TLegend, TLine  # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kFullCircle, kOpenCircle, kFullSquare  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

inputdir = 'outputs/raa'

inputfilenames = ['HFPtSpectrumRaaDs_010_RightPredFile_2015.root',
                  'HFPtSpectrumRaaDs_010_central_2015.root', 'HFPtSpectrumRaaDs_010_central_2015_oldppref.root']
histonames = ['hRABvsPt', 'hRABvsPt', 'hRABvsPt']
graphnames = ['gRAB_GlobalSystematics',
              'gRAB_GlobalSystematics', 'gRAB_GlobalSystematics']
colors = [kRed, kBlue, kBlack, kBlue, kRed]
linecolors = [kRed, kBlue, kBlack, kBlue, kRed]
markers = [kFullSquare, kFullCircle, kOpenCircle]
legendnames = ['D_{s}^{+} 2015 (2015 cuts)', 'D_{s}^{+} 2018 (2015 cuts)',
               'D_{s}^{+} 2018 (2015 cuts, 2015 pp ref)']
outputsuffix = 'Ds2015_Ds2018_samecuts'

SetGlobalStyle(padleftmargin=0.15, padtopmargin=0.05, titlesize=0.045, labelsize=0.04)

hRaa, gRaa, hRaaRatio = ([] for i in range(3))

leg = TLegend(0.3, 0.68, 0.75, 0.83)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile, _ in enumerate(inputfilenames):
    inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
    hRaa.append(inputfile.Get(histonames[iFile]))
    gRaa.append(inputfile.Get(graphnames[iFile]))
    hRaa[iFile].SetDirectory(0)
    SetObjectStyle(hRaa[iFile], linecolor=linecolors[iFile], markercolor=colors[iFile],
                   markersize=1.5, markerstyle=markers[iFile])
    SetObjectStyle(gRaa[iFile], linecolor=colors[iFile], fillstyle=0)
    leg.AddEntry(hRaa[iFile], legendnames[iFile], 'p')
    hRaaRatio.append(hRaa[iFile].Clone("hRaa%d" % iFile))
    hRaaRatio[iFile].SetDirectory(0)
    hRaaRatio[iFile].Divide(hRaa[iFile], hRaa[0])

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
cRaa.cd(2).DrawFrame(ptmin, 0., ptmax, 2.,
                     ';#it{p}_{T} (GeV/#it{c}); ratio #it{R}_{AA}')
for iFile in range(len(inputfilenames)):
    if iFile == 0:
        continue
    # gRaa[iFile].Draw('2')
    hRaaRatio[iFile].Draw('same')

cRaa.SaveAs('%s/RaaComparison_%s.pdf' % (inputdir, outputsuffix))

input('Press enter to exit')
