'''
python script for the check of the mass shaping effect
run: python CheckMassShaping.py cfgFileName.yml PtMin PtMax MLoutMin MLoutMax MLoutStep outFileName
'''

import argparse
import six
import yaml
from ROOT import TCanvas, TLegend, kFullCircle, kRed  # pylint: disable=import-error,no-name-in-module
from utils.TaskFileLoader import LoadSingleSparseFromTask
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle

SetGlobalStyle()

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('PtMin', type=float, default=2., help='minimum pT')
parser.add_argument('PtMax', type=float, default=3., help='maximum pT')
parser.add_argument('MLoutMin', type=float, default=0.99, help='minimum value of ML output to be checked')
parser.add_argument('MLoutMax', type=float, default=0.999, help='maximum value of ML output to be checked')
parser.add_argument('MLoutStep', type=float, default=0.001, help='step of ML output to be checked')
parser.add_argument('outFileName', metavar='text', default='outFileName', help='output file name w/o extension')
parser.add_argument('--rebin', type=int, required=False, default=1, help='mass rebin (optional)')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

infilenames = inputCfg['filename']
if not isinstance(infilenames, list):
    infilenames = [infilenames]

for iFile, infilename in enumerate(infilenames):
    if iFile == 0:
        sparseBkg = LoadSingleSparseFromTask(infilename, inputCfg, 'sparsenameBkg')
    else:
        sparseBkgPart = LoadSingleSparseFromTask(infilename, inputCfg, 'sparsenameBkg')

#draw mass vs. ML output
PtBinMin = sparseBkg.GetAxis(1).FindBin(args.PtMin*1.0001)
PtBinMax = sparseBkg.GetAxis(1).FindBin(args.PtMax*0.9999)
sparseBkg.GetAxis(1).SetRange(PtBinMin, PtBinMax)

hMassVsML = sparseBkg.Projection(0, 13)
nMassbins = hMassVsML.GetYaxis().GetNbins()

cMassVsML = TCanvas('cMassVsML', '', 800, 800)
hMassVsML.GetXaxis().SetNdivisions(505)
hMassVsML.GetYaxis().SetTitle('#it{M}(KK#pi) (GeV/#it{c}^{2})')
hMassVsML.GetYaxis().SetDecimals()
hMassVsML.Draw('colz')

#draw mass for several ML outputs
MLoutBinMin = sparseBkg.GetAxis(13).FindBin(args.MLoutMin*1.0001)
MLoutBinMax = sparseBkg.GetAxis(13).FindBin(args.MLoutMax*0.9999)
MLoutStepBin = int(round(args.MLoutStep / sparseBkg.GetAxis(13).GetBinWidth(1), 0))
nMLbins = sparseBkg.GetAxis(13).GetNbins()

if MLoutStepBin == 0:
    print('ERROR: ML step passed is smaller than THnSparse binning! Exit')
    exit()

leg = TLegend(0.45, 0.6, 0.85, 0.9)
leg.SetTextSize(0.04)

hMass = []
cMassMLSel = TCanvas('cMassMLSel', '', 800, 800)
cMassMLSel.cd().DrawFrame(hMassVsML.GetYaxis().GetBinLowEdge(1), 0., \
    hMassVsML.GetYaxis().GetBinUpEdge(nMassbins), hMassVsML.ProjectionY().GetMaximum()*1.5*args.rebin, \
        ';#it{M}(KK#pi) (GeV/#it{c}^{2}); Counts')

for iBin, MLbin in enumerate(range(MLoutBinMin, MLoutBinMax+MLoutStepBin, MLoutStepBin)):
    MLoutBinMax = sparseBkg.GetAxis(13).SetRange(MLbin, nMLbins+1)
    hMass.append(sparseBkg.Projection(0))
    SetObjectStyle(hMass[iBin], color=kRed-5+iBin, markerstyle=kFullCircle, linewidth=2)
    if args.rebin:
        hMass[iBin].Rebin(args.rebin)
    hMass[iBin].Draw('Esame')
    leg.AddEntry(hMass[iBin], 'ML output > {:0.4f}'.format(sparseBkg.GetAxis(13).GetBinLowEdge(MLbin)), 'p')

leg.Draw()

cMassVsML.SaveAs('{0}_MassVsML.pdf'.format(args.outFileName))
cMassMLSel.SaveAs('{0}_MassDiffSel.pdf'.format(args.outFileName))

if six.PY2:
    raw_input('Press Enter to exit')
elif six.PY3:
    input('Press Enter to exit')
