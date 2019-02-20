from ROOT import TFile, TCanvas, TH1F, TF1, TNtuple, TGraph, TSpline3 # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kGreen, kOrange # pylint: disable=import-error,no-name-in-module
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
cDist.Divide((len(cutVars)+3)/2,2)
counter = 0
sel_string = 'PtMin>=%f && PtMax<=%f && Signif>%f && EffPrompt>%f' % (float(PtMin), float(PtMax), float(minSignif), float(minEffPrompt))
for iVar in cutVars :
  counter += 1
  cDist.cd(counter)
  ntuple.Draw(iVar, sel_string)
cDist.cd(counter+1)
ntuple.Draw('EffPrompt', sel_string)
cDist.cd(counter+2)
ntuple.Draw('Signif', sel_string)
cDist.cd(counter+3)
ntuple.Draw('EffPrompt:Signif', sel_string, 'colz')

cDist.SaveAs('DsNtupleProj_pt_%d_%d.pdf' % (int(PtMin), int(PtMax)))
raw_input("Press enter to exit")