'''
Script to get the non-prompt D cross-section from Pythia8 (to run after SimulateDecays.cc)
'''

import sys
import argparse
import yaml
from ROOT import TFile, TCanvas, TLegend, TH1F #pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure, kOrange, kGreen #pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle  #pylint: disable=wrong-import-position,import-error
from utils.DfUtils import LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
inKineFileName = inputCfg['inKineFileName']
Bhadrons = inputCfg['pdgCodeB']
FFbtoB = inputCfg['FFbtoB']
Dhadrons = inputCfg['pdgCodeD']
BRDhadrons = inputCfg['BRD']
partPlusAntiPart = inputCfg['partPlusAntiPart']
outFileName = inputCfg['outFileName']

kineDf = LoadDfFromRootOrParquet(inKineFileName, None, 'fTreeDecays')

SetGlobalStyle(padleftmargin=0.14, padbottommargin=0.12, titleoffsety=1.3, optstat=0)
Bcolors = [kRed+1, kAzure+4, kOrange+7, kGreen+2]
Bnames = {511:'B^{+}', 521:'B^{0}', 531:'B_{s}^{0}', 5122:'#Lambda_{b}^{0}'}
Dnames = {411:'D^{+}', 421:'D^{0}', 431:'D_{s}^{+}', 4122:'#Lambda_{c}^{+}'}

leg = TLegend(0.6, 0.65, 0.95, 0.95)
leg.SetTextSize(0.045)
leg.SetBorderSize(0)
leg.SetFillStyle(0)

# get normalisation
norm = (kineDf['norm'].values)[0]

# get BRs of B -> D from sim
BRBhadronsToD = {}
if inputCfg['BRB']['fromPythia']:
    for B in Bhadrons:
        BRBhadronsToD[B] = {}
        kineDfSelB = kineDf.query(f'pdgB == {B}')
        # remove double-counted B mesons (i.e. those decaying to Ds+ and Ds- simultaneously)
        kineDfSelBForDen = kineDfSelB.drop_duplicates(subset=['ptB', 'pB', 'yB'], keep='first')
        for D in Dhadrons:
            if partPlusAntiPart:
                kineDfSelBSelD = kineDfSelB.query(f'abs(pdgD) == {D}')
            else:
                kineDfSelBSelD = kineDfSelB.query(f'pdgD == {D}')
            BRBhadronsToD[B][D] = float(len(kineDfSelBSelD)/len(kineDfSelBForDen))
            if partPlusAntiPart:
                print(f'BR of {B} -> +/-{abs(D)} + X: {BRBhadronsToD[B][D]:.3f}')
            else:
                print(f'BR of {B} -> {D} + X: {BRBhadronsToD[B][D]:.3f}')
else:
    BRcustom = inputCfg['BRB']['customValues']
    for iB, B in enumerate(Bhadrons):
        BRBhadronsToD[B] = {}
        for iD, D in enumerate(Dhadrons):
            BRBhadronsToD[B][D] = BRcustom[iB][iD]
            if partPlusAntiPart:
                print(f'BR of {B} -> +/-{abs(D)} + X: {BRBhadronsToD[B][D]:.3f}')
            else:
                print(f'BR of {B} -> {D} + X: {BRBhadronsToD[B][D]:.3f}')

hPtDFromB = {}
cPtDFromB = {}
for D, BRD in zip(Dhadrons, BRDhadrons):
    if partPlusAntiPart:
        kineDfSelD = kineDf.query(f'abs(pdgD) == {D}')
    else:
        kineDfSelD = kineDf.query(f'pdgD == {D}')

    if inputCfg['accFactor']['apply']: # acceptance cut
        kineDfSelDAcc = kineDfSelD.query(inputCfg['accFactor']['selection'])
        acc = float(len(kineDfSelDAcc)/len(kineDfSelD))
    else:
        kineDfSelDAcc = kineDfSelD
        acc = 1.

    hPtDFromB[D] = {}
    cPtDFromB[D] = TCanvas('cPtDFromB', '', 800, 800)
    cPtDFromB[D].DrawFrame(0., 1.e-7, 50.05, 2.,
                           ';d#it{p}_{T} (GeV/#it{c});d#sigma/d#it{p}_{T} #times BR (#mub GeV^{-1} #it{c});')
    cPtDFromB[D].SetLogy()
    for iB, (B, FFb) in enumerate(zip(Bhadrons, FFbtoB)):
        kineDfSelDSelB = kineDfSelDAcc.query(f'pdgB == {B}')
        hPtDFromB[D][B] = TH1F(f'hPt{D}From{B}', '', 1001, 0., 50.05)
        for pt in kineDfSelDSelB['ptD'].to_numpy():
            hPtDFromB[D][B].Fill(pt)
        hPtDFromB[D][B].Sumw2()
        hPtDFromB[D][B].Scale(1.e-6 * BRBhadronsToD[B][D] * FFb * norm * acc * BRD / hPtDFromB[D][B].Integral())
        SetObjectStyle(hPtDFromB[D][B], linecolor=Bcolors[iB], markercolor=Bcolors[iB], linewidth=1, markerstyle=0)
        hPtDFromB[D][B].Draw('same')
        leg.AddEntry(hPtDFromB[D][B], f'{Bnames[B]} #rightarrow {Dnames[D]}', 'le')

    hPtDFromB[D]['sum'] = hPtDFromB[D][B].Clone(f'hPt{D}FromAllB')
    SetObjectStyle(hPtDFromB[D]['sum'], linecolor=kBlack, markercolor=kBlack, linewidth=1, markerstyle=0)
    hPtDFromB[D]['sum'].Reset()
    for B in Bhadrons:
        hPtDFromB[D]['sum'].Add(hPtDFromB[D][B])
    hPtDFromB[D]['sum'].Draw('same')
    leg.AddEntry(hPtDFromB[D]['sum'], f'H_{{b}} #rightarrow {Dnames[D]}', 'le')
    leg.Draw()
    outFileNamePDF = outFileName.replace('.root', f'_{D}.pdf')
    cPtDFromB[D].SaveAs(outFileNamePDF)

outFile = TFile(outFileName, 'recreate')
for D in Dhadrons:
    cPtDFromB[D].Write()
    for B in Bhadrons:
        hPtDFromB[D][B].Write()
    hPtDFromB[D]['sum'].Write()
outFile.Close()

input('Press enter to exit')
