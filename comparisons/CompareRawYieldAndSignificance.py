import sys
import math
from ROOT import TCanvas, TFile, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlue # pylint: disable=import-error,no-name-in-module
from ROOT import kFullCircle, kFullSquare # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error,no-name-in-module

inputdir = 'outputs/'
inputfilenames = ['2015results/RawYieldDs_3050.root', '3_24bin_merge/raw_yields/RawYieldsDs_3050_2015cuts.root']
outputdir = 'outputs/3_24bin_merge/raw_yields/2015_comp/'
outputsuffix = 'LHC15o_LHC18qr'
signalhistonames = ['hSignal', 'hRawYields']
bkghistonames = ['hBackground', 'hRawYieldsBkg']
SoverBhistonames = ['hSignal', 'hRawYieldsSoverB']
signifhistonames = ['hSignificance', 'hRawYieldsSignificance']
evhistonames = ['hNEvents', 'hEvForNorm']
colors = [kRed, kBlue]
markers = [kFullSquare, kFullCircle]
legendnames = ['LHC15o (2015 cuts)', 'LHC18qr (2015 cuts)']

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, titlesize=0.045, labelsize=0.04)

hSignal, hRatioSignal, hBackground, hRatioBkg, hSoverB, hSignif, hRatioSignif, hEv = [], [], [], [], [], [], [], []

leg = TLegend(0.5, 0.73, 0.8, 0.93)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile, filename in enumerate(inputfilenames):
    inputfile = TFile(f'{inputdir}/{filename}')
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
    hSignal[iFile].Scale(1./hEv[iFile].GetBinContent(1))
    hBackground[iFile].Scale(1./hEv[iFile].GetBinContent(1))
    hSignif[iFile].Scale(1./math.sqrt(hEv[iFile].GetBinContent(1)))
    SetObjectStyle(hSignal[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hBackground[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hSoverB[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hSignif[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    hRatioSignal.append(hSignal[iFile].Clone(f'hRatioSignal{iFile}'))
    hRatioSignal[iFile].SetDirectory(0)
    hRatioSignal[iFile].Divide(hSignal[iFile], hSignal[0], 1., 1., "")
    hRatioBkg.append(hBackground[iFile].Clone(f'hRatioBkg{iFile}'))
    hRatioBkg[iFile].SetDirectory(0)
    hRatioBkg[iFile].Divide(hBackground[iFile], hBackground[0], 1., 1., "")
    hRatioSignif.append(hSignif[iFile].Clone(f'hRatioSignif{iFile}'))
    hRatioSignif[iFile].SetDirectory(0)
    hRatioSignif[iFile].Divide(hSignif[iFile], hSignif[0], 1., 1., "")
    leg.AddEntry(hSignal[iFile], legendnames[iFile], 'lp')

PtMin = hSignal[0].GetBinLowEdge(1)
PtMax = hSignal[0].GetBinLowEdge(hSignal[0].GetNbinsX())+hSignal[0].GetBinWidth(hSignal[0].GetNbinsX())

cSignal = TCanvas('cSignal', '', 1000, 500)
cSignal.Divide(2, 1)
cSignal.cd(1).DrawFrame(PtMin, 0., PtMax, 1.e-6, ';#it{p}_{T} (GeV/#it{c}); raw yields / #it{N}_{events}')
for iFile, _ in enumerate(inputfilenames):
    hSignal[iFile].Draw('same')
leg.Draw()
cSignal.cd(2).DrawFrame(PtMin, 0., PtMax, 3., ';#it{p}_{T} (GeV/#it{c}); ratio of raw yields / #it{N}_{events}')
for iFile, _ in enumerate(inputfilenames):
    if iFile == 0:
        continue
    hRatioSignal[iFile].Draw('same')

cRatioSignal = TCanvas('cRatioSignal', '', 800, 800)
cRatioSignal.DrawFrame(PtMin, 0.5, PtMax, 3., ';#it{p}_{T} (GeV/#it{c}); ratio of raw yields / #it{N}_{events}')
for iFile, _ in enumerate(inputfilenames):
    hRatioSignal[iFile].Draw('same')
leg.Draw()

cBkg = TCanvas('cBkg', '', 1000, 500)
cBkg.Divide(2, 1)
cBkg.cd(1).DrawFrame(PtMin, 1.e-9, PtMax, 5.e-6, ';#it{p}_{T} (GeV/#it{c}); background(3#sigma) / #it{N}_{events}')
cBkg.cd(1).SetLogy()
for iFile, _ in enumerate(inputfilenames):
    hBackground[iFile].Draw('same')
leg.Draw()
cBkg.cd(2).DrawFrame(PtMin, 0., PtMax, 3., ';#it{p}_{T} (GeV/#it{c}); ratio of background(3#sigma) / #it{N}_{events}')
for iFile, _ in enumerate(inputfilenames):
    if iFile == 0:
        continue
    hRatioBkg[iFile].Draw('same')

cSoverB = TCanvas('cSoverB', '', 800, 800)
cSoverB.DrawFrame(PtMin, 0.1, PtMax, 15., ';#it{p}_{T} (GeV/#it{c}); S/B(3#sigma)')
cSoverB.SetLogy()
for iFile, _ in enumerate(inputfilenames):
    hSoverB[iFile].Draw('same')
leg.Draw()

cSignificance = TCanvas('cSignificance', '', 1000, 500)
cSignificance.Divide(2, 1)
cSignificance.cd(1).DrawFrame(PtMin, 0., PtMax, 1.e-3,
                              ';#it{p}_{T} (GeV/#it{c}); significance / #sqrt{#it{N}_{events}}')
for iFile, _ in enumerate(inputfilenames):
    hSignif[iFile].Draw('same')
leg.Draw()
cSignificance.cd(2).DrawFrame(PtMin, 0., PtMax, 2.,
                              ';#it{p}_{T} (GeV/#it{c}); ratio of significance / #sqrt{#it{N}_{events}}')
for iFile, _ in enumerate(inputfilenames):
    if iFile == 0:
        continue
    hRatioSignif[iFile].Draw('same')

cSignal.SaveAs(f'{outputdir}/SignalPerEventComparison_{outputsuffix}.pdf')
cRatioSignal.SaveAs(f'{outputdir}/SignalPerEventRatio_{outputsuffix}.pdf')
cBkg.SaveAs(f'{outputdir}/BkgPerEventComparison_{outputsuffix}.pdf')
cSoverB.SaveAs(f'{outputdir}/SoverB_{outputsuffix}.pdf')
cSignificance.SaveAs(f'{outputdir}/SignificancePerEventComparison_{outputsuffix}.pdf')

input('Press enter to exit')
