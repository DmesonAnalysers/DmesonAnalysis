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
parser.add_argument('--centMin', type=float, default=-1,
                    help='minimum centrality (percentile)')
parser.add_argument('--centMax', type=float, default=101,
                    help='maximum centrality (percentile)')
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
        if inList.FindObject('hSPDMult'):
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
        else:
            if iFile == 0:
                hTmp = inList.FindObject('hSPDMultVsCent')
                if 'MC' not in inpt:
                    binCentMin = hTmp.GetXaxis().FindBin(args.centMin)
                    binCentMax = hTmp.GetXaxis().FindBin(args.centMax)
                else:
                    binCentMin = -1
                    binCentMax = -1
                hNtrkl[inpt] = hTmp.ProjectionY('hSPDMult', binCentMin, binCentMax).Clone()
                hTmp = inList.FindObject('hSPDMultVsCentCand')
                hNtrklCand[inpt] = hTmp.ProjectionY('hSPDMultCand', binCentMin, binCentMax).Clone()
                hTmp = inList.FindObject('hSPDMultVsCentCandInMass')
                hNtrklCandInMass[inpt] = hTmp.ProjectionY('hSPDMultCandInMass', binCentMin, binCentMax).Clone()
                if not hNtrkl[inpt] or not hNtrklCand[inpt] or not hNtrklCandInMass[inpt]:
                    print('ERROR: sparse without multiplicity histograms! Exit')
                    sys.exit()
                hNtrkl[inpt].SetDirectory(0)
                hNtrklCand[inpt].SetDirectory(0)
                hNtrklCandInMass[inpt].SetDirectory(0)
            else:
                hTmp = inList.FindObject('hSPDMultVsCent')
                hNtrkl[inpt].Add(hTmp.ProjectionY('hSPDMult', binCentMin, binCentMax).Clone())
                hTmp = inList.FindObject('hSPDMultVsCentCand')
                hNtrklCand[inpt].Add(hTmp.ProjectionY('hSPDMultCand', binCentMin, binCentMax).Clone())
                hTmp = inList.FindObject('hSPDMultVsCentCandInMass')
                hNtrklCandInMass[inpt].Add(hTmp.ProjectionY('hSPDMultCandInMass', binCentMin, binCentMax).Clone())

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

lNtrkl = TLegend(0.2, 0.8, 0.7, 0.9)
lNtrkl.SetTextSize(0.045)
lNtrkl.SetBorderSize(0)
lNtrkl.SetFillStyle(0)
lNtrkl.AddEntry(hNtrkl['data'],
                f'data, mean = {hNtrkl["data"].GetMean():0.1f}, RMS = {hNtrkl["data"].GetRMS():0.1f}', 'pl')
lNtrkl.AddEntry(hNtrkl['MC'],
                f'MC, mean = {hNtrkl["MC"].GetMean():0.1f}, RMS = {hNtrkl["MC"].GetRMS():0.1f}', 'pl')

lNtrklCand = TLegend(0.2, 0.8, 0.7, 0.9)
lNtrklCand.SetTextSize(0.045)
lNtrklCand.SetBorderSize(0)
lNtrklCand.SetFillStyle(0)
lNtrklCand.AddEntry(hNtrklCand['data'],
                    f'data, mean = {hNtrklCand["data"].GetMean():0.1f},'
                    f' RMS = {hNtrklCand["data"].GetRMS():0.1f}', 'pl')
lNtrklCand.AddEntry(hNtrklCand['MC'],
                    f'MC, mean = {hNtrklCand["MC"].GetMean():0.1f},'
                    f' RMS = {hNtrklCand["MC"].GetRMS():0.1f}', 'pl')

lNtrklCandInMass = TLegend(0.2, 0.8, 0.7, 0.9)
lNtrklCandInMass.SetTextSize(0.045)
lNtrklCandInMass.SetBorderSize(0)
lNtrklCandInMass.SetFillStyle(0)
lNtrklCandInMass.AddEntry(hNtrklCandInMass['data'],
                          f'data, mean = {hNtrklCandInMass["data"].GetMean():0.1f},'
                          f' RMS = {hNtrklCandInMass["data"].GetRMS():0.1f}', 'pl')
lNtrklCandInMass.AddEntry(hNtrklCandInMass['MC'],
                          f'MC, mean = {hNtrklCandInMass["MC"].GetMean():0.1f},'
                          f' RMS = {hNtrklCandInMass["MC"].GetRMS():0.1f}', 'pl')

cNtrkl = TCanvas('cNtrkl', '', 1500, 500)
cNtrkl.Divide(3, 1)
cNtrkl.cd(1).DrawFrame(0., 1.e-10, hNtrkl[inpt].GetXaxis().GetXmax(), 10.,
                       'all events;#it{N}_{tracklets} (|#it{#eta}|<1);Normalised counts')
cNtrkl.cd(1).SetLogy()
for inpt in hNtrkl:
    hNtrkl[inpt].Draw('same')
lNtrkl.Draw()
cNtrkl.cd(2).DrawFrame(0., 1.e-10, hNtrkl[inpt].GetXaxis().GetXmax(), 10.,
                       'events with cand;#it{N}_{tracklets} (|#it{#eta}|<1);Normalised counts')
cNtrkl.cd(2).SetLogy()
for inpt in hNtrklCand:
    hNtrklCand[inpt].Draw('same')
lNtrklCand.Draw()
cNtrkl.cd(3).DrawFrame(0., 1.e-10, hNtrkl[inpt].GetXaxis().GetXmax(), 10.,
                       'events with cand in mass;#it{N}_{tracklets} (|#it{#eta}|<1);Normalised counts')
cNtrkl.cd(3).SetLogy()
for inpt in hNtrklCandInMass:
    hNtrklCandInMass[inpt].Draw('same')
lNtrklCandInMass.Draw()
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
