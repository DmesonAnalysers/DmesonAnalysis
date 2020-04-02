import sys
from ROOT import TCanvas, TFile, TLegend, TLine  # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kFullCircle, kOpenCircle, kFullSquare  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error,no-name-in-module

inputdir = 'outputs/crosssec'

inputfilenames = ['DsYield_method2_fd1_br1_010_2015_paper.root',
                  'DsYield_method2_fd1_br1_010_central_2015.root']
histonames = ['hAAC', 'hAAC']
graphnames = ['gaaCsystTot', 'gaaCsystTot']
colors = [kRed, kBlue, kBlack, kBlue, kRed]
linecolors = [kRed, kBlue, kBlack, kBlue, kRed]
markers = [kFullSquare, kFullCircle, kOpenCircle]
legendnames = ['D_{s}^{+} 2015 (2015 cuts)', 'D_{s}^{+} 2018 (2015 cuts)']
outputsuffix = 'Ds2015_Ds2018_samecuts'

hCorrYield, gCorrYield, hCorrYieldRatio = ([] for _ in range(3))

SetGlobalStyle(padleftmargin=0.15, padtopmargin=0.05, titlesize=0.045, labelsize=0.04)

leg = TLegend(0.3, 0.68, 0.75, 0.83)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile, _ in enumerate(inputfilenames):
    inputfile = TFile('%s/%s' % (inputdir, inputfilenames[iFile]))
    hCorrYield.append(inputfile.Get(histonames[iFile]))
    gCorrYield.append(inputfile.Get(graphnames[iFile]))
    hCorrYield[iFile].SetDirectory(0)
    SetObjectStyle(hCorrYield[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile], markersize=1.5)
    SetObjectStyle(gCorrYield[iFile], linecolor=colors[iFile], fillstyle=0)
    leg.AddEntry(hCorrYield[iFile], legendnames[iFile], 'p')
    hCorrYieldRatio.append(hCorrYield[iFile].Clone("hCorrYield%d" % iFile))
    hCorrYieldRatio[iFile].SetDirectory(0)
    hCorrYieldRatio[iFile].Divide(hCorrYield[iFile], hCorrYield[0])

ptmin = hCorrYield[0].GetBinLowEdge(1)+2
ptmax = hCorrYield[0].GetBinLowEdge(hCorrYield[0].GetNbinsX())+hCorrYield[0].GetBinWidth(hCorrYield[0].GetNbinsX())

lineatone = TLine(ptmin, 1., ptmax, 1.)
lineatone.SetLineWidth(1)
lineatone.SetLineColor(kBlack)
lineatone.SetLineStyle(9)

cCorrYield = TCanvas('cCorrYield', '', 1000, 500)
cCorrYield.Divide(2, 1)
cCorrYield.cd(1).DrawFrame(ptmin, 1.e-4, ptmax, 1.,
                           ';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T} (#it{c} GeV^{-1})')
cCorrYield.cd(1).SetLogy()
lineatone.Draw('same')
for iFile in range(len(inputfilenames)):
    # gCorrYield[iFile].Draw('2')
    hCorrYield[iFile].Draw('same')
leg.Draw()
cCorrYield.cd(2).DrawFrame(ptmin, 0., ptmax, 2.,
                           ';#it{p}_{T} (GeV/#it{c}); d#it{N}/d#it{p}_{T} ratio')
for iFile in range(len(inputfilenames)):
    if iFile == 0:
        continue
    # gCorrYield[iFile].Draw('2')
    hCorrYieldRatio[iFile].Draw('same')

cCorrYield.SaveAs('%s/CorrYieldComparison_%s.pdf' % (inputdir, outputsuffix))

input('Press enter to exit')
