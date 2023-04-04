import sys
from os.path import join
import math
from ROOT import TCanvas, TFile, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue,kGreen, kAzure, kOrange,kFullCircle, kOpenCircle, kFullSquare, kFullDiamond,kOpenDiamond, kFullCross, kOpenCross, kOpenSquare # pylint: disable=import-error,no-name-in-module,unused-import
sys.path.append('..')
from utils.AnalysisUtils import ComputeRatioDiffBins #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

inputdir = '/home/fchinu/Ds_pp_13TeV/output_analysis'
inputfilenames = ['/final/RawYieldDs_data_pp13TeV.root',
                '/stefano/RawYieldsDs_Ds_pp13TeV_PromptEn_22112021.root',
                '/stefano/RawYieldsDs_CentralConsPID.root',
                '/stefano/RawYieldsDs_CentralStrongPID.root']
mergeFileNames=['/stefano/RawYieldsDs_CentralConsPID.root',
                '/stefano/RawYieldsDs_CentralStrongPID.root']
outputdir = '/home/fchinu/Ds_pp_13TeV/output_analysis/final/thesis'
outputsuffix = '_final_binary_multiclass_standard'
signalhistonames = ['hRawYields', 'hRawYields', 'hRawYields', 'hRawYields']
bkghistonames = ['hRawYieldsBkg', 'hRawYieldsBkg', 'hRawYieldsBkg', 'hRawYieldsBkg']
SoverBhistonames = ['hRawYieldsSoverB', 'hRawYieldsSoverB', 'hRawYieldsSoverB', 'hRawYieldsSoverB']
signifhistonames = ['hRawYieldsSignificance', 'hRawYieldsSignificance', 'hRawYieldsSignificance', 'hRawYieldsSignificance']
evhistonames = ['hEvForNorm', 'hEvForNorm', 'hEvForNorm', 'hEvForNorm']
colors = [kAzure+3, kGreen-2, kOrange-3]
markers = [kFullCircle, kFullDiamond, kFullSquare,  kFullCross]
legendnames = ['Binary', 'Multiclass', 'Standard','Std (Strong PID)']

SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14,
               titleoffsety=1.8, titlesize=0.045, labelsize=0.04, maxdigits=2)

def StdMerge(hCons, hStrong, ptSwitch = 6):
    '''
    Helper function to merge the standard histograms
    '''
    hMerged = hCons.Clone(f'hMergedStd')
    for iBin in range(hMerged.GetNbinsX()):
        if (hStrong.GetBinLowEdge(iBin+1)) >= ptSwitch:        
            hMerged.SetBinContent(iBin+1, hCons.GetBinContent(iBin+1))
            hMerged.SetBinError(iBin+1, hCons.GetBinError(iBin+1))
        else:
            hMerged.SetBinContent(iBin+1, hStrong.GetBinContent(iBin+1))
            hMerged.SetBinError(iBin+1, hStrong.GetBinError(iBin+1))
        
    return  hMerged

temphSignal, temphRatioSignal, temphBackground, temphRatioBkg, temphSoverB, temphSignif, temphRatioSignif, temphEv = [], [], [], [], [], [], [], []

for iFile, inFileName in enumerate(mergeFileNames):
    inFileName = inputdir+ inFileName
    print(inFileName)
    inFile = TFile.Open(inFileName)
    temphSignal.append(inFile.Get('hRawYields'))
    temphBackground.append(inFile.Get('hRawYieldsBkg'))
    temphSoverB.append(inFile.Get('hRawYieldsSoverB'))
    temphSignif.append(inFile.Get('hRawYieldsSignificance'))
    temphEv.append(inFile.Get('hEvForNorm'))
    temphSignal[iFile].SetDirectory(0)
    temphBackground[iFile].SetDirectory(0)
    temphSoverB[iFile].SetDirectory(0)
    temphSignif[iFile].SetDirectory(0)
    temphEv[iFile].SetDirectory(0)
    inFile.Close()

file = TFile.Open('/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/mergedrawyields.root','recreate')
StdMerge(temphSignal[0],temphSignal[1]).Write('hRawYields')
StdMerge(temphBackground[0],temphBackground[1]).Write('hRawYieldsBkg')
StdMerge(temphSoverB[0],temphSoverB[1]).Write('hRawYieldsSoverB')
StdMerge(temphSignif[0],temphSignif[1]).Write('hRawYieldsSignificance')
StdMerge(temphEv[0],temphEv[1]).Write('hEvForNorm')
file.Close()

inputfilenames.pop()
inputfilenames.pop()
inputfilenames.append('/stefano/mergedrawyields.root')


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

PtMax=36

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

cSignifOnly = TCanvas('cSignifOnly', '', 800, 800)
cSignifOnly.DrawFrame(PtMin, 0., PtMax, hSignif[0].GetMaximum()*1.5,
                      ';#it{p}_{T} (GeV/#it{c}); significance / #sqrt{#it{N}_{events}}')
for iFile, _ in enumerate(inputfilenames):
    hSignif[iFile].Draw('same')
leg.Draw()

cSignal.SaveAs(f'{outputdir}/SignalPerEventComparison_{outputsuffix}.png')
cRatioSignal.SaveAs(f'{outputdir}/SignalPerEventRatio_{outputsuffix}.png')
cBkg.SaveAs(f'{outputdir}/BkgPerEventComparison_{outputsuffix}.png')
cSoverB.SaveAs(f'{outputdir}/SoverB_{outputsuffix}.png')
cSignificance.SaveAs(f'{outputdir}/SignificancePerEventComparison_{outputsuffix}.png')
cSignifOnly.SaveAs(f'{outputdir}/SignifPerEventComp_NoRatio_{outputsuffix}.png')

input('Press enter to exit')
