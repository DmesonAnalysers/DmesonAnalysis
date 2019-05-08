from ROOT import TFile, TCanvas, TH1F, TF1, TNtuple, TGraph, TSpline3 # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kBlack, kBlue, kGreen, kOrange # pylint: disable=import-error,no-name-in-module
import sys, yaml, six

cfgFileName = sys.argv[1]
inputFileName = sys.argv[2]
PtMin = sys.argv[3]
PtMax = sys.argv[4]
SignifMin = sys.argv[5]
SignifMax = sys.argv[6]
EffPromptMin = sys.argv[7]
EffPromptMax = sys.argv[8]

with open(cfgFileName, 'r') as ymlCfgFile:
  inputCfg = yaml.load(ymlCfgFile)

cutVars = inputCfg['cutvars']

infile = TFile(inputFileName)
ntuple = infile.Get('tSignif')

cDist = TCanvas('cDist','',1920,1080)
numCol = int(round((len(cutVars) + 4)/3.))
cDist.Divide(numCol,3)
counter = 0
sel_string = 'PtMin>=%f && PtMax<=%f && Signif>%f && Signif<%f && EffPrompt>%f && EffPrompt<%f' % (float(PtMin),
             float(PtMax), float(SignifMin), float(SignifMax), float(EffPromptMin), float(EffPromptMax))
for iVar in cutVars :
  counter += 1
  cDist.cd(counter)
  ntuple.Draw(iVar, sel_string)
cDist.cd(counter+1)
ntuple.Draw('EffPrompt', sel_string)
cDist.cd(counter+2)
ntuple.Draw('Signif', sel_string)
cDist.cd(counter+3)
ntuple.Draw('SoverB', sel_string)
cDist.cd(counter+4)
ntuple.Draw('Signif:EffPrompt', sel_string, 'colz')
cDist.cd(counter+5)
ntuple.Draw('Signif:SoverB', sel_string, 'colz')

cDist.SaveAs('DsNtupleProj_pt_%d_%d.pdf' % (int(PtMin), int(PtMax)))
if six.PY2:
    raw_input('Press enter to exit')
elif six.PY3:
    input('Press enter to exit')