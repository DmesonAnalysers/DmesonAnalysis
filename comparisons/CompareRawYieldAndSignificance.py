import sys
import math
from ROOT import TCanvas, TFile, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kAzure # pylint: disable=import-error,no-name-in-module
from ROOT import kFullCircle, kFullSquare # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.AnalysisUtils import ComputeRatioDiffBins #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

inputdir = 'inputdir'
inputfilenames = ['file1.root', 'file2.root']
outputdir = 'outputdir'
outputsuffix = 'suffix'
signalhistonames = ['hRawYields', 'hRawYields']
bkghistonames = ['hRawYieldsBkg', 'hRawYieldsBkg']
SoverBhistonames = ['hRawYieldsSoverB', 'hRawYieldsSoverB']
signifhistonames = ['hRawYieldsSignificance', 'hRawYieldsSignificance']
evhistonames = ['hEvForNorm', 'hEvForNorm']
colors = [kRed+1, kAzure+4]
markers = [kFullSquare, kFullCircle]
legendnames = ['title1', 'title2']

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14,
               titleoffsety=1.8, titlesize=0.045, labelsize=0.04, maxdigits=2)

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
    hRatioSignal.append(ComputeRatioDiffBins(hSignal[iFile], hSignal[0]))
    hRatioSignal[iFile].SetDirectory(0)
    hRatioSignal[iFile].SetName(f'hRatioSignal{iFile}')
    hRatioBkg.append(ComputeRatioDiffBins(hBackground[iFile], hBackground[0]))
    hRatioBkg[iFile].SetDirectory(0)
    hRatioBkg[iFile].SetName(f'hRatioBkg{iFile}')
    hRatioSignif.append(ComputeRatioDiffBins(hSignif[iFile], hSignif[0]))
    hRatioSignif[iFile].SetDirectory(0)
    hRatioSignif[iFile].SetName(f'hRatioSignif{iFile}')
    SetObjectStyle(hSignal[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hBackground[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hSoverB[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hSignif[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hRatioSignal[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hRatioBkg[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hRatioSignif[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    leg.AddEntry(hSignal[iFile], legendnames[iFile], 'lp')

PtMin = hSignal[0].GetBinLowEdge(1)
PtMax = hSignal[0].GetBinLowEdge(hSignal[0].GetNbinsX())+hSignal[0].GetBinWidth(hSignal[0].GetNbinsX())
for histo in hSignal:
    if histo.GetBinLowEdge(1) < PtMin:
        PtMin = histo.GetBinLowEdge(1)
    if histo.GetBinLowEdge(histo.GetNbinsX())+histo.GetBinWidth(histo.GetNbinsX()) > PtMax:
        PtMax = histo.GetBinLowEdge(histo.GetNbinsX())+histo.GetBinWidth(histo.GetNbinsX())

cSignal = TCanvas('cSignal', '', 1000, 500)
cSignal.Divide(2, 1)
cSignal.cd(1).DrawFrame(PtMin, 0., PtMax, hSignal[1].GetMaximum()*2,
                        ';#it{p}_{T} (GeV/#it{c}); raw yields / #it{N}_{events}')
for iFile, _ in enumerate(inputfilenames):
    hSignal[iFile].Draw('same')
leg.Draw()
cSignal.cd(2).DrawFrame(PtMin, 0., PtMax, hRatioSignal[1].GetMaximum()*2,
                        ';#it{p}_{T} (GeV/#it{c}); ratio of raw yields / #it{N}_{events}')
for iFile, _ in enumerate(inputfilenames):
    if iFile == 0:
        continue
    hRatioSignal[iFile].Draw('same')

cRatioSignal = TCanvas('cRatioSignal', '', 800, 800)
cRatioSignal.DrawFrame(PtMin, 0.5, PtMax, hRatioSignal[1].GetMaximum()*2,
                       ';#it{p}_{T} (GeV/#it{c}); ratio of raw yields / #it{N}_{events}')
for iFile, _ in enumerate(inputfilenames):
    hRatioSignal[iFile].Draw('same')
leg.Draw()

cBkg = TCanvas('cBkg', '', 1000, 500)
cBkg.Divide(2, 1)
cBkg.cd(1).DrawFrame(PtMin, hBackground[0].GetMinimum()*0.5, PtMax, hBackground[0].GetMaximum()*2,
                     ';#it{p}_{T} (GeV/#it{c}); background(3#sigma) / #it{N}_{events}')
cBkg.cd(1).SetLogy()
for iFile, _ in enumerate(inputfilenames):
    hBackground[iFile].Draw('same')
leg.Draw()
cBkg.cd(2).DrawFrame(PtMin, hRatioBkg[1].GetMinimum()*0.5, PtMax, hRatioBkg[1].GetMaximum()*2,
                     ';#it{p}_{T} (GeV/#it{c}); ratio of background(3#sigma) / #it{N}_{events}')
for iFile, _ in enumerate(inputfilenames):
    if iFile == 0:
        continue
    hRatioBkg[iFile].Draw('same')

cSoverB = TCanvas('cSoverB', '', 800, 800)
cSoverB.DrawFrame(PtMin, hSoverB[0].GetMinimum()*0.5, PtMax,
                  hSoverB[0].GetMaximum()*2, ';#it{p}_{T} (GeV/#it{c}); S/B(3#sigma)')
cSoverB.SetLogy()
for iFile, _ in enumerate(inputfilenames):
    hSoverB[iFile].Draw('same')
leg.Draw()

cSignificance = TCanvas('cSignificance', '', 1000, 500)
cSignificance.Divide(2, 1)
cSignificance.cd(1).DrawFrame(PtMin, 0., PtMax, hSignif[0].GetMaximum()*2,
                              ';#it{p}_{T} (GeV/#it{c}); significance / #sqrt{#it{N}_{events}}')
for iFile, _ in enumerate(inputfilenames):
    hSignif[iFile].Draw('same')
leg.Draw()
cSignificance.cd(2).DrawFrame(PtMin, 0., PtMax, hRatioSignif[1].GetMaximum()*2,
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
