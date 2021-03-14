'''
Script for the computation of multiplicity weights
run: python3 ComputeMultWeights.py inConfigMC.yml inConfigData.yml outFile.root
'''

import sys
import argparse
import yaml
from ROOT import TFile, TCanvas, TLegend # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, GetROOTColor, GetROOTMarker

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('inConfigMC', metavar='text', default='config_MC.yml',
                    help='config file name with root input files for MC')
parser.add_argument('inConfigData', metavar='text', default='config_data.yml',
                    help='config file name with root input files for data')
parser.add_argument('outFileName', metavar='text', default='outFile.root',
                    help='output root file name')
args = parser.parse_args()

with open(args.inConfigMC, 'r') as ymlCfgFile:
    inputCfgMC = yaml.load(ymlCfgFile, yaml.FullLoader)

with open(args.inConfigData, 'r') as ymlCfgFile:
    inputCfgData = yaml.load(ymlCfgFile, yaml.FullLoader)

inFileNames = {'MC': inputCfgMC['filename'], 'data': inputCfgData['filename']}
inDirNames = {'MC': inputCfgMC['dirname'], 'data': inputCfgData['dirname']}
inListNames = {'MC': inputCfgMC['listname'], 'data': inputCfgData['listname']}
colors = {'MC': GetROOTColor('kAzure+4'), 'data': GetROOTColor('kRed+1')}
markers = {'MC': GetROOTMarker('kFullCircle'), 'data': GetROOTMarker('kFullSquare')}
hNtrkl, hNtrklCand, hNtrklCandInMass = ({} for _ in range(3))
for inpt in inFileNames:
    if not isinstance(inFileNames[inpt], list):
        inFileNames[inpt] = [inFileNames[inpt]]

    for iFile, fileName in enumerate(inFileNames[inpt]):
        inFile = TFile.Open(fileName)
        inDir = inFile.Get(inDirNames[inpt])
        inList = inDir.Get(inListNames[inpt])
        if iFile == 0:
            hNtrkl[inpt] = inList.FindObject('hSPDMult')
            hNtrklCand[inpt] = inList.FindObject('hSPDMultCand')
            hNtrklCandInMass[inpt] = inList.FindObject('hSPDMultCandInMass')
            if not hNtrkl[inpt] or not hNtrklCand[inpt] or not hNtrklCandInMass[inpt]:
                print('ERROR: sparse without multiplicity histograms! Exit')
                sys.exit()
        else:
            hNtrkl[inpt].Add(inList.FindObject('hSPDMult'))
            hNtrklCand[inpt].Add(inList.FindObject('hSPDMultCand'))
            hNtrklCandInMass[inpt].Add(inList.FindObject('hSPDMultCandInMass'))

    hNtrkl[inpt].SetName(f'hNtrkl{inpt}')
    hNtrklCand[inpt].SetName(f'hNtrklCand{inpt}')
    hNtrklCandInMass[inpt].SetName(f'hNtrklCandInMass{inpt}')
    hNtrkl[inpt].Sumw2()
    hNtrklCand[inpt].Sumw2()
    hNtrklCandInMass[inpt].Sumw2()
    hNtrkl[inpt].Scale(1./hNtrkl[inpt].Integral())
    hNtrklCand[inpt].Scale(1./hNtrklCand[inpt].Integral())
    hNtrklCandInMass[inpt].Scale(1./hNtrklCandInMass[inpt].Integral())
    SetObjectStyle(hNtrkl[inpt], fillstyle=0, color=colors[inpt], markerstyle=markers[inpt])
    SetObjectStyle(hNtrklCand[inpt], fillstyle=0, color=colors[inpt], markerstyle=markers[inpt])
    SetObjectStyle(hNtrklCandInMass[inpt], fillstyle=0, color=colors[inpt], markerstyle=markers[inpt])

hNtrklWeights = hNtrkl['data'].Clone('hNtrklWeights')
hNtrklWeights.Divide(hNtrkl['MC'])
hNtrklWeightsCand = hNtrklCand['data'].Clone('hNtrklWeightsCand')
hNtrklWeightsCand.Divide(hNtrklCand['MC'])
hNtrklWeightsCandInMass = hNtrklCandInMass['data'].Clone('hNtrklWeightsCandInMass')
hNtrklWeightsCandInMass.Divide(hNtrklCandInMass['MC'])

SetGlobalStyle(padbottommargin=0.12, padleftmargin=0.14, padtopmargin=0.075, titleoffsety=1.4, opttitle=1)

leg = TLegend(0.7, 0.8, 0.9, 0.9)
leg.SetTextSize(0.045)
leg.SetBorderSize(0)
leg.SetFillStyle(0)
leg.AddEntry(hNtrkl['data'], 'data', 'pl')
leg.AddEntry(hNtrkl['MC'], 'MC', 'pl')

cNtrkl = TCanvas('cNtrkl', '', 1500, 500)
cNtrkl.Divide(3, 1)
cNtrkl.cd(1).DrawFrame(0., 1.e-10, hNtrkl[inpt].GetXaxis().GetXmax(), 1.,
                       'all events;#it{N}_{tracklets} (|#it{#eta}|<1);Normalised counts')
cNtrkl.cd(1).SetLogy()
for inpt in hNtrkl:
    hNtrkl[inpt].Draw('same')
leg.Draw()
cNtrkl.cd(2).DrawFrame(0., 1.e-10, hNtrkl[inpt].GetXaxis().GetXmax(), 1.,
                       'events with cand;#it{N}_{tracklets} (|#it{#eta}|<1);Normalised counts')
cNtrkl.cd(2).SetLogy()
for inpt in hNtrklCand:
    hNtrklCand[inpt].Draw('same')
leg.Draw()
cNtrkl.cd(3).DrawFrame(0., 1.e-10, hNtrkl[inpt].GetXaxis().GetXmax(), 1.,
                       'events with cand in mass;#it{N}_{tracklets} (|#it{#eta}|<1);Normalised counts')
cNtrkl.cd(3).SetLogy()
for inpt in hNtrklCandInMass:
    hNtrklCandInMass[inpt].Draw('same')
leg.Draw()
cNtrkl.Modified()
cNtrkl.Update()

cWeights = TCanvas('cWeights', '', 1500, 500)
cWeights.Divide(3, 1)
minW = hNtrklWeights.GetMinimum() if hNtrklWeights.GetMinimum() > 0 else 1.e-6
cWeights.cd(1).DrawFrame(0., minW, hNtrklWeights.GetXaxis().GetXmax(),
                         hNtrklWeights.GetMaximum()*2,
                         'all events;#it{N}_{tracklets} (|#it{#eta}|<1);weights')
cWeights.cd(1).SetLogy()
hNtrklWeights.Draw('same')
minW = hNtrklWeightsCand.GetMinimum() if hNtrklWeightsCand.GetMinimum() > 0 else 1.e-6
cWeights.cd(2).DrawFrame(0., minW, hNtrklWeightsCand.GetXaxis().GetXmax(),
                         hNtrklWeightsCand.GetMaximum()*2,
                         'events with cand;#it{N}_{tracklets} (|#it{#eta}|<1);weights')
cWeights.cd(2).SetLogy()
hNtrklWeightsCand.Draw('same')
minW = hNtrklWeightsCandInMass.GetMinimum() if hNtrklWeightsCandInMass.GetMinimum() > 0 else 1.e-6
cWeights.cd(3).DrawFrame(0., minW, hNtrklWeightsCandInMass.GetXaxis().GetXmax(),
                         hNtrklWeightsCandInMass.GetMaximum()*2,
                         'events with cand in mass;#it{N}_{tracklets} (|#it{#eta}|<1);weights')
cWeights.cd(3).SetLogy()
hNtrklWeightsCandInMass.Draw('same')
cWeights.Modified()
cWeights.Update()

outFile = TFile.Open(args.outFileName, 'recreate')
for inpt in hNtrklCand:
    hNtrkl[inpt].Write()
    hNtrklCand[inpt].Write()
    hNtrklCandInMass[inpt].Write()
hNtrklWeights.Write()
hNtrklWeightsCand.Write()
hNtrklWeightsCandInMass.Write()
outFile.Close()

input('Press enter to exit')
