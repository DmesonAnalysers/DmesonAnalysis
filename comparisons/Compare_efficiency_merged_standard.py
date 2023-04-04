import sys
from os.path import join
from ROOT import TCanvas, TFile, TLegend, TH1F, TH1 # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue,kGreen, kAzure, kOrange,kFullCircle, kOpenCircle, kFullSquare, kFullDiamond,kOpenDiamond, kFullCross, kOpenCross, kOpenSquare # pylint: disable=import-error,no-name-in-module,unused-import
sys.path.append('..')
from utils.AnalysisUtils import ComputeRatioDiffBins #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

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

inDir = ''
inFileNames = ['/home/fchinu/Ds_pp_13TeV/output_analysis/final/AccTimesEfficiencyDs_pp13TeV.root',
    '/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/Eff_times_Acc_Ds_Ds_pp13TeV_PromptEn_22112021.root',
    '/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/AccxEff_Ds_CentralConsPID_LHC20f4_ptwhisto.root',
    '/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/AccxEff_Ds_CentralStrongPID_LHC20f4_ptwhisto.root']
mergeFileNames=['/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/AccxEff_Ds_CentralConsPID_LHC20f4_ptwhisto.root',
    '/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/AccxEff_Ds_CentralStrongPID_LHC20f4_ptwhisto.root']
histoNames = ['hAccEff', 'hAccEff', 'hAxe', 'hAxe']
colors= [kAzure+3, kGreen-2, kOrange-3]
markers= [kFullCircle, kFullDiamond, kFullSquare,  kFullCross]
legNames = ['Binary', 'Multiclass', 'Standard','Std (Strong PID)']
outDir = '/home/fchinu/Ds_pp_13TeV/output_analysis/final/thesis'
outSuffix = 'Comp_binary_multiclass_standard'
showUnc = True

SetGlobalStyle(padleftmargin=0.16, padtopmargin=0.05, padbottommargin=0.14,
               titleoffsety=1.5, titlesize=0.05, labelsize=0.045)

tempEffPrompt = [] 
tempEffFD= []

for iFile, inFileName in enumerate(mergeFileNames):
    inFileName = join(inDir, inFileName)
    print(inFileName)
    inFile = TFile.Open(inFileName)
    tempEffPrompt.append(inFile.Get('hAxePrompt'))
    tempEffFD.append(inFile.Get('hAxeFeeddw'))
    tempEffPrompt[iFile].SetDirectory(0)
    tempEffFD[iFile].SetDirectory(0)
    inFile.Close()

file = TFile.Open('/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/mergedacc.root','recreate')
StdMerge(tempEffPrompt[0],tempEffPrompt[1]).Write('hAxePrompt')
StdMerge(tempEffFD[0],tempEffFD[1]).Write('hAxeFD')
file.Close()

inFileNames.pop()
inFileNames.pop()
inFileNames.append('/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/mergedacc.root')


hEffPrompt, hEffFD, hEffPromptRatio, hEffFDRatio = ([] for _ in range(4))

leg = TLegend(0.2, 0.18, 0.9, 0.33)
leg.SetFillStyle(0)
leg.SetBorderSize(0)
leg.SetTextSize(0.04)

for iFile, inFileName in enumerate(inFileNames):
    inFileName = join(inDir, inFileName)
    inFile = TFile.Open(inFileName)
    if not (inFileName=='/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/AccxEff_Ds_CentralConsPID_LHC20f4_ptwhisto.root' or inFileName=='/home/fchinu/Ds_pp_13TeV/output_analysis/stefano/AccxEff_Ds_CentralStrongPID_LHC20f4_ptwhisto.root'):
        hEffPrompt.append(inFile.Get(f'{histoNames[iFile]}Prompt'))
        hEffFD.append(inFile.Get(f'{histoNames[iFile]}FD'))
    else:
        hEffPrompt.append(inFile.Get(f'{histoNames[iFile]}Prompt'))
        hEffFD.append(inFile.Get(f'{histoNames[iFile]}Feeddw'))
    hEffPrompt[iFile].SetDirectory(0)
    hEffFD[iFile].SetDirectory(0)
    hEffPrompt[iFile].GetYaxis().SetRangeUser(1e-4, 1.)
    hEffFD[iFile].GetYaxis().SetRangeUser(1e-4, 1.)
    hEffPrompt[iFile].GetXaxis().SetRangeUser(1, 36.)
    hEffFD[iFile].GetXaxis().SetRangeUser(1, 36)
    hEffPromptRatio.append(ComputeRatioDiffBins(hEffPrompt[iFile], hEffPrompt[0], 'B'))
    hEffPromptRatio[iFile].SetDirectory(0)
    hEffPromptRatio[iFile].SetName(f'hEffPromptRatio{iFile}')
    hEffFDRatio.append(ComputeRatioDiffBins(hEffFD[iFile], hEffFD[0], 'B'))
    hEffFDRatio[iFile].SetDirectory(0)
    hEffFDRatio[iFile].SetName(f'hEffFDRatio{iFile}')
    SetObjectStyle(hEffPrompt[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    SetObjectStyle(hEffFD[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile], linestyle=1)
    SetObjectStyle(hEffPromptRatio[iFile], linecolor=colors[iFile], markercolor=colors[iFile],
                   markerstyle=markers[iFile], linestyle=1)
    SetObjectStyle(hEffFDRatio[iFile], linecolor=colors[iFile], markercolor=colors[iFile], markerstyle=markers[iFile])
    leg.AddEntry(hEffFD[iFile], legNames[iFile], 'p')
    if not showUnc:
        for iBin in range(hEffPromptRatio[iFile].GetNbinsX()):
            hEffPromptRatio[iFile].SetBinError(iBin+1, 1.e-20)
            hEffFDRatio[iFile].SetBinError(iBin+1, 1.e-20)

PtMin = hEffPrompt[0].GetBinLowEdge(1)
PtMax = hEffPrompt[0].GetBinLowEdge(hEffPrompt[0].GetNbinsX())+hEffPrompt[0].GetBinWidth(hEffPrompt[0].GetNbinsX())
for histo in hEffPrompt:
    if histo.GetBinLowEdge(1) < PtMin:
        PtMin = histo.GetBinLowEdge(1)
    if histo.GetBinLowEdge(histo.GetNbinsX())+histo.GetBinWidth(histo.GetNbinsX()) > PtMax:
        PtMax = histo.GetBinLowEdge(histo.GetNbinsX())+histo.GetBinWidth(histo.GetNbinsX())
PtMax =36.
cPrompt = TCanvas('cPrompt', '', 1000, 500)
cPrompt.Divide(2, 1)
cPrompt.cd(1).DrawFrame(PtMin, hEffPrompt[0].GetMinimum()/5, PtMax, 1.,
                        ';#it{p}_{T} (GeV/#it{c}); Prompt (Acc #times #epsilon)')
cPrompt.cd(1).SetLogy()
for iFile in range(len(inFileNames)):
    hEffPrompt[iFile].Draw('same')
leg.Draw()
cPrompt.cd(2).DrawFrame(PtMin, hEffPromptRatio[1].GetMinimum()/2, PtMax, hEffPromptRatio[1].GetMaximum()*1.5,
                        ';#it{p}_{T} (GeV/#it{c}); Prompt (Acc #times #epsilon) ratio')
for iFile in range(len(inFileNames)):
    if iFile == 0:
        continue
    hEffPromptRatio[iFile].Draw('same')

cFD = TCanvas('cFD', '', 1000, 500)
cFD.Divide(2, 1)
cFD.cd(1).DrawFrame(PtMin, hEffFD[0].GetMinimum()/5, PtMax, 1.,
                    ';#it{p}_{T} (GeV/#it{c}); Feed-down (Acc #times #epsilon)')
cFD.cd(1).SetLogy()
for iFile in range(len(inFileNames)):
    hEffFD[iFile].Draw('same')
leg.Draw()
cFD.cd(2).DrawFrame(PtMin, hEffFDRatio[1].GetMinimum()/2, PtMax, hEffFDRatio[1].GetMaximum()*2,
                    ';#it{p}_{T} (GeV/#it{c}); Feed-down (Acc #times #epsilon) ratio')
for iFile in range(len(inFileNames)):
    if iFile == 0:
        continue
    hEffFDRatio[iFile].Draw('same')

cPrompt.SaveAs(f'{outDir}/PromptEfficiencyComparison_{outSuffix}.png')
cFD.SaveAs(f'{outDir}/FDEfficiencyComparison_{outSuffix}.png')

input('Press enter to exit')


