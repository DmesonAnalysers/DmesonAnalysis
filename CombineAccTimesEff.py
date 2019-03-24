#******************************************************************************************#
# python script for the combination of efficiency and acceptance factors                   #
# run: python CombineAccTimesEff.py effFileName.root accFileName.root outFileName.root     #
# author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino                           #
#******************************************************************************************#

from ROOT import TFile, TCanvas, TH1F, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import gStyle, kBlack, kFullCircle, kOpenSquare, kFullDiamond # pylint: disable=import-error,no-name-in-module
import yaml, sys, array, math, string

effFileName = sys.argv[1]
accFileName = sys.argv[2]
outFileName = sys.argv[3]

effFile = TFile.Open(effFileName)
hEffPrompt = effFile.Get('hEffPrompt')
hEffFD = effFile.Get('hEffFD')

accFile = TFile.Open(accFileName)
hAccNum = accFile.Get('hPtGenAcc')
hAccDen = accFile.Get('hPtGenLimAcc')
hAcc = hEffPrompt.Clone('hAcc')
hAcc.SetMarkerColor(kBlack)
hAcc.SetLineColor(kBlack)
hAcc.SetMarkerStyle(kFullDiamond)

hTmpNum = TH1F('hTmpNum','',1,0,1)
hTmpDen = TH1F('hTmpDen','',1,0,1)
nPtBins = hEffPrompt.GetNbinsX()
for iPt in range(nPtBins):
  PtMin = hEffPrompt.GetBinLowEdge(iPt+1)
  PtMax = PtMin+hEffPrompt.GetBinWidth(iPt+1)
  PtBinMin = hAccNum.GetXaxis().FindBin(PtMin*1.0001)
  PtBinMax = hAccNum.GetXaxis().FindBin(PtMax*0.9999)
  
  hTmpNum.SetBinContent(1,hAccNum.Integral(PtBinMin,PtBinMax))  
  hTmpNum.SetBinError(1,math.sqrt(hAccNum.Integral(PtBinMin,PtBinMax)))  
  hTmpDen.SetBinContent(1,hAccDen.Integral(PtBinMin,PtBinMax))
  hTmpDen.SetBinError(1,math.sqrt(hAccDen.Integral(PtBinMin,PtBinMax)))  

  hTmpNum.Divide(hTmpNum,hTmpDen,1.,1.,'B')
  hAcc.SetBinContent(iPt+1,hTmpNum.GetBinContent(1))
  hAcc.SetBinError(iPt+1,hTmpNum.GetBinError(1))

hAccEffPrompt = hEffPrompt.Clone('hAccEffPrompt')
hAccEffPrompt.Multiply(hAcc)
hAccEffFD = hEffFD.Clone('hAccEffFD')
hAccEffFD.Multiply(hAcc)

hAccEffPrompt.SetTitle(';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon);')
hAccEffFD.SetTitle(';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon);')

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadTopMargin(0.035)
gStyle.SetTitleSize(0.045,'xy')
gStyle.SetLabelSize(0.040,'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)
gStyle.SetOptStat(0)

cAccEff = TCanvas('cAccEff','',800,800)
cAccEff.DrawFrame(hEffPrompt.GetBinLowEdge(1),1.e-5,hEffPrompt.GetBinLowEdge(nPtBins)+hEffPrompt.GetBinWidth(nPtBins),1.,';#it{p}_{T} (GeV/#it{c});(Acc #times #epsilon);')
cAccEff.SetLogy()
hAccEffPrompt.Draw('same')
hAccEffFD.Draw('same')

outFile = TFile(outFileName,'recreate')
cAccEff.Write()
hEffPrompt.Write()
hEffFD.Write()
hAccEffPrompt.Write()
hAccEffFD.Write()
outFile.Close()

outFileNamePDF = string.replace(outFileName,'.root','.pdf')
cAccEff.SaveAs(outFileNamePDF)

raw_input('Press enter to exit')