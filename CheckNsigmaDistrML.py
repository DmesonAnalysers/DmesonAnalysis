'''
python script for the projection of Nsigma distributions from TH3 with PID var vs. ML output vs. pT
run: python CheckNsigmaDistrML.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import argparse
import copy
import string
import six
import yaml
from ROOT import TFile, TCanvas, TH1F, TDirectoryFile, TLegend, gStyle  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlue, kRed, kFullCircle, kOpenCircle  # pylint: disable=import-error,no-name-in-module
from TaskFileLoader import LoadPIDTH3


def SetHistoStyle(histo, color, alpha, marker, markersize=1.5, linewidth=2):
    '''
    method to set histo style
    '''
    histo.SetMarkerColor(color)
    histo.SetFillColorAlpha(color, alpha)
    histo.SetLineColorAlpha(color, alpha)
    histo.SetLineWidth(linewidth)
    histo.SetMarkerStyle(marker)
    histo.SetMarkerSize(markersize)


# main function
gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadBottomMargin(0.14)
gStyle.SetPadTopMargin(0.1)
gStyle.SetTitleSize(0.055, 'xyz')
gStyle.SetLabelSize(0.050, 'xyz')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)
gStyle.SetOptStat(0)

PARSER = argparse.ArgumentParser(description='Arguments to pass')
PARSER.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
PARSER.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                    help='input file with cut set')
PARSER.add_argument('outFileName', metavar='text', default='outFileName.root',
                    help='output root file name')
ARGS = PARSER.parse_args()

if six.PY2:
    outFileNamePDF = string.replace(ARGS.outFileName, '.root', '.pdf')
elif six.PY3:
    outFileNamePDF = ARGS.outFileName.replace('.root', '.pdf')

with open(ARGS.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

infilenames = inputCfg['filename']
if not isinstance(infilenames, list):
    infilenames = [infilenames]

hNsigmaVsPtVsML = {}
for iFile, infilename in enumerate(infilenames):
    hNsigmaTmp = LoadPIDTH3(infilename, inputCfg)
    if iFile == 0:
        hNsigmaVsPtVsML = copy.deepcopy(hNsigmaTmp)

    for det in hNsigmaTmp:
        for spe in hNsigmaTmp[det]:
            for prong in hNsigmaTmp[det][spe]:
                hNsigmaVsPtVsML[det][spe][prong].Add(
                    hNsigmaTmp[det][spe][prong])

with open(ARGS.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']
MLmin, MLmax = [], []
for var in cutVars:
    if 'BDT' in var or 'ML' in var:
        MLmin = cutVars[var]['min']
        MLmax = cutVars[var]['max']
        break

if MLmin == [] or MLmax == []:
    print('ERROR: You did not set any ML selection! Exit')
    exit()

hNsigma, hNsigmaSel = {}, {}
for det in hNsigmaVsPtVsML:
    hNsigma[det], hNsigmaSel[det] = {}, {}
    for spe in hNsigmaVsPtVsML[det]:
        hNsigma[det][spe], hNsigmaSel[det][spe] = {}, {}
        for prong in hNsigmaVsPtVsML[det][spe]:
            hNsigma[det][spe][prong], hNsigmaSel[det][spe][prong] = {}, {}
            for iPt, (ptmin, ptmax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):

                ptbinmin = hNsigmaVsPtVsML[det][spe][prong].GetXaxis().FindBin(ptmin*1.0001)
                ptbinmax = hNsigmaVsPtVsML[det][spe][prong].GetXaxis().FindBin(ptmax*0.9999)

                hNsigma[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                    hNsigmaVsPtVsML[det][spe][prong].ProjectionZ(
                        'hNsigma{:s}_{:s}_Prong{:s}_Pt{:.0f}_{:.0f}'.format(
                            det, spe, prong, ptmin, ptmax), ptbinmin, ptbinmax, 0, 100000)
                SetHistoStyle(hNsigma[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)], \
                    kBlue, 0.25, kFullCircle, 0.5)

                MLbinmin = hNsigmaVsPtVsML[det][spe][prong].GetYaxis().FindBin(MLmin[iPt]*1.0001)
                MLbimax = hNsigmaVsPtVsML[det][spe][prong].GetYaxis().FindBin(MLmax[iPt]*1.0001)

                hNsigmaSel[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                    hNsigmaVsPtVsML[det][spe][prong].ProjectionZ(
                        'hNsigmaSel{:s}_{:s}_Prong{:s}_Pt{:.0f}_{:.0f}'.format(
                            det, spe, prong, ptmin, ptmax), ptbinmin, ptbinmax, MLbinmin, MLbimax)
                SetHistoStyle(hNsigmaSel[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)], \
                    kRed, 0.25, kOpenCircle, 0.5)

leg = TLegend(0.55, 0.76, 0.89, 0.89)
leg.SetTextSize(0.05)
leg.SetFillStyle(0)
leg.AddEntry(list(hNsigma[det][spe][prong].values())[0], 'w/o ML selection', 'fp')
leg.AddEntry(list(hNsigmaSel[det][spe][prong].values())[0], 'w/ ML selection', 'fp')

#plot on canvases and save on root file
outfile = TFile(ARGS.outFileName, 'recreate')
cNsigma, dirNsigma = {}, {}
for iDet, det in enumerate(hNsigmaVsPtVsML):
    cNsigma[det] = {}
    outfile.cd()
    dirNsigma[det] = TDirectoryFile(det, det)
    dirNsigma[det].Write()
    dirNsigma[det].cd()
    for iPt, (ptmin, ptmax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
            TCanvas('cNsigma{:s}_Pt{:.0f}_{:.0f}'.format(
                det, ptmin, ptmax), '', 1920, 1080)
        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Divide(3, 2)
        for iSpe, spe in enumerate(hNsigmaVsPtVsML[det]):
            for iProng, prong in enumerate(hNsigmaVsPtVsML[det][spe]):
                cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).SetLogy()
                his = hNsigma[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)]
                hissel = hNsigmaSel[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)]
                if his.Integral() > 0:
                    his.Scale(1./his.Integral())
                if hissel.Integral() > 0:
                    hissel.Scale(1./hissel.Integral())
                his.SetTitle(' %0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{N}_{#sigma}^{%s} (%s);Normalised entries' % (
                    ptmin, ptmax, det, spe))
                hissel.SetTitle(' %0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{N}_{#sigma}^{%s} (%s);Normalised entries' % (
                    ptmin, ptmax, det, spe))
                his.GetYaxis().SetRangeUser(1.e-5, 1.)
                his.GetXaxis().SetTitleSize(0.055)
                hissel.GetXaxis().SetTitleSize(0.055)
                his.GetXaxis().SetLabelSize(0.050)
                hissel.GetXaxis().SetLabelSize(0.050)
                his.DrawCopy('hist')
                his.DrawCopy('Esame')
                hissel.DrawCopy('histsame')
                hissel.DrawCopy('Esame')
                his.Write()
                hissel.Write()
                leg.Draw('same')

                cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).Modified()
                cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).Update()

        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Write()

        #save also pdf
        if iDet == 0 and iPt == 0:
            cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Print('{0}['.format(outFileNamePDF))
        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Print(outFileNamePDF)
        if iDet == len(hNsigmaVsPtVsML)-1 and iPt == len(cutVars['Pt']['min'])-1:
            cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Print('{0}]'.format(outFileNamePDF))

outfile.Close()

if six.PY2:
    raw_input('Press enter to exit')
elif six.PY3:
    input('Press enter to exit')
