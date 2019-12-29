'''
python script for the number of events from Ds or Dplus task
run: python GetNumberOfEvents.py cfgFileName.yml outFileName.root
'''

import argparse
import yaml
from ROOT import TH1F, TFile  # pylint: disable=import-error,no-name-in-module
from utils.TaskFileLoader import LoadNormObjFromTask

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('outFileName', metavar='text', default='outFileName.root',
                    help='output root file name')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

infilenames = inputCfg['filename']
if not isinstance(infilenames, list):
    infilenames = [infilenames]

for iFile, infilename in enumerate(infilenames):
    if iFile == 0:
        hEv, normCounter = LoadNormObjFromTask(infilename, inputCfg)
    else:
        hEvPart, normCounterPart = LoadNormObjFromTask(infilename, inputCfg)
        hEv.Add(hEvPart)
        normCounterPart.Add(normCounterPart)

hEvForNorm = TH1F("hEvForNorm", ";;Number of events", 2, 0., 2.)
hEvForNorm.GetXaxis().SetBinLabel(1, "norm counter")
hEvForNorm.GetXaxis().SetBinLabel(2, "accepted events")
hEvForNorm.SetBinContent(1, normCounter.GetNEventsForNorm())
for iBin in range(1, hEv.GetNbinsX()+1):
    binLabel = hEv.GetXaxis().GetBinLabel(iBin)
    if binLabel.Contains('isEvSelected') or binLabel.Contains('accepted'):
        hEvForNorm.SetBinContent(2, hEv.GetBinContent(iBin))
        break

outFile = TFile(args.outFileName, 'recreate')
hEvForNorm.Write()
outFile.Close()
