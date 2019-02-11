from ROOT import TFile, TCanvas, TH1F, TF1, TNtuple, TGraph, TSpline3
from ROOT import kRed, kBlack, kBlue, kGreen, kOrange
import sys, yaml

cfgFileName = sys.argv[1]
inputFileName = sys.argv[2]
PtMin = sys.argv[3]
PtMax = sys.argv[4]
minSignif = sys.argv[5]
minEffPrompt = sys.argv[6]

with open(cfgFileName, 'r') as ymlCfgFile:
  inputCfg = yaml.load(ymlCfgFile)

cutVars = inputCfg['cutvars']

infile = TFile(inputFileName)
ntuple = infile.Get('tSignif')

cDist = TCanvas('cDist','',1920,1080)
cDist.Divide((len(cutVars)+2)/2,2)
counter = 0
for iVar in cutVars :
  counter += 1
  cDist.cd(counter)
  ntuple.Draw(iVar,'PtMin>=%f && PtMax<=%f && Signif>%f && EffPrompt>%f' % (float(PtMin), float(PtMax), float(minSignif), float(minEffPrompt)))
cDist.cd(counter+1)
ntuple.Draw('EffPrompt','PtMin>=%f && PtMax<=%f && Signif>%f && EffPrompt>%f' % (float(PtMin), float(PtMax), float(minSignif), float(minEffPrompt)))
cDist.cd(counter+2)
ntuple.Draw('Signif','PtMin>=%f && PtMax<=%f && Signif>%f && EffPrompt>%f' % (float(PtMin), float(PtMax), float(minSignif), float(minEffPrompt)))

cDist.SaveAs('DsNtupleProj_pt_%f_%f.pdf')
