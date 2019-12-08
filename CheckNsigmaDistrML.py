'''
python script for the projection of Nsigma distributions from THnSparses with PID var vs. ML output vs. pT
run: python CheckNsigmaDistrML.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import argparse
import copy
import string
import os
import six
import numpy as np
import yaml
from PIL import Image
from ROOT import TFile, TCanvas, TDirectoryFile, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlue, kRed, kFullCircle, kOpenCircle  # pylint: disable=import-error,no-name-in-module
from utils.TaskFileLoader import LoadPIDSparses
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle


# main function
SetGlobalStyle()
PARSER = argparse.ArgumentParser(description='Arguments to pass')
PARSER.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
PARSER.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml',
                    help='input file with cut set')
PARSER.add_argument('outFileName', metavar='text', default='outFileName.root',
                    help='output root file name')
ARGS = PARSER.parse_args()

if six.PY2:
    outFileNameWoExt = string.replace(ARGS.outFileName, '.root', '')
elif six.PY3:
    outFileNameWoExt = ARGS.outFileName.replace('.root', '')

with open(ARGS.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

infilenames = inputCfg['filename']
if not isinstance(infilenames, list):
    infilenames = [infilenames]

for iFile, infilename in enumerate(infilenames):
    sPIDNsigmaTmp, sPIDNsigmaCombTmp, axes = LoadPIDSparses(infilename, inputCfg)
    if iFile == 0:
        sPIDNsigma = copy.deepcopy(sPIDNsigmaTmp)
        sPIDNsigmaComb = copy.deepcopy(sPIDNsigmaCombTmp)
    else:
        sPIDNsigma.Add(sPIDNsigmaTmp)
        sPIDNsigmaComb.Add(sPIDNsigmaCombTmp)

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
#1D distributions
for det in axes:
    hNsigma[det], hNsigmaSel[det] = {}, {}
    for spe in axes[det]:
        hNsigma[det][spe], hNsigmaSel[det][spe] = {}, {}
        for prong in axes[det][spe]:
            hNsigma[det][spe][prong], hNsigmaSel[det][spe][prong] = {}, {}
            for iPt, (ptmin, ptmax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):

                ptbinmin = sPIDNsigma.GetAxis(0).FindBin(ptmin*1.0001)
                ptbinmax = sPIDNsigma.GetAxis(0).FindBin(ptmax*0.9999)

                MLbinmin = sPIDNsigma.GetAxis(1).FindBin(MLmin[iPt]*1.0001)
                MLbinmax = sPIDNsigma.GetAxis(1).FindBin(MLmax[iPt]*0.9999)


                if det != 'Comb':
                    sparse = sPIDNsigma
                else:
                    sparse = sPIDNsigmaComb

                sparse.GetAxis(0).SetRange(ptbinmin, ptbinmax)

                hNsigma[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                    sparse.Projection(axes[det][spe][prong])
                hNsigma[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SetName(\
                    'hNsigma{0}_{1}_Prong{2}_Pt{3:.0f}_{4:.0f}'.format(det, spe, prong, ptmin, ptmax))
                SetObjectStyle(hNsigma[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(\
                    ptmin, ptmax)], color=kBlue, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.5, linewidth=1)

                sparse.GetAxis(1).SetRange(MLbinmin, MLbinmax)

                hNsigmaSel[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                    sparse.Projection(axes[det][spe][prong])
                hNsigmaSel[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SetName(\
                    'hNsigmaSel{0}_{1}_Prong{2}_Pt{3:.0f}_{4:.0f}'.format(det, spe, prong, ptmin, ptmax))
                SetObjectStyle(hNsigmaSel[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(\
                    ptmin, ptmax)], color=kRed, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.5, linewidth=1)

                sparse.GetAxis(0).SetRange(-1, -1)
                sparse.GetAxis(1).SetRange(-1, -1)

#2D distributions (TPC vs TOF)
hNsigma['TPCTOF'], hNsigmaSel['TPCTOF'] = {}, {}
for spe in axes[det]:
    hNsigma['TPCTOF'][spe], hNsigmaSel['TPCTOF'][spe] = {}, {}
    for prong in axes['TPC'][spe]:
        hNsigma['TPCTOF'][spe][prong], hNsigmaSel['TPCTOF'][spe][prong] = {}, {}
        for iPt, (ptmin, ptmax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):

            ptbinmin = sPIDNsigma.GetAxis(0).FindBin(ptmin*1.0001)
            ptbinmax = sPIDNsigma.GetAxis(0).FindBin(ptmax*0.9999)

            MLbinmin = sPIDNsigma.GetAxis(1).FindBin(MLmin[iPt]*1.0001)
            MLbinmax = sPIDNsigma.GetAxis(1).FindBin(MLmax[iPt]*0.9999)

            sPIDNsigma.GetAxis(0).SetRange(ptbinmin, ptbinmax)

            hNsigma['TPCTOF'][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                sPIDNsigma.Projection(axes['TPC'][spe][prong], axes['TOF'][spe][prong])
            hNsigma['TPCTOF'][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SetName(\
                'hNsigmaTPCTOF_{0}_Prong{1}_Pt{2:.0f}_{3:.0f}'.format(spe, prong, ptmin, ptmax))
            SetObjectStyle(hNsigma['TPCTOF'][spe][prong]['Pt{:.0f}_{:.0f}'.format(\
                ptmin, ptmax)], color=kBlue, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.3, linewidth=1)

            sPIDNsigma.GetAxis(1).SetRange(MLbinmin, MLbinmax)

            hNsigmaSel['TPCTOF'][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                sPIDNsigma.Projection(axes['TPC'][spe][prong], axes['TOF'][spe][prong])
            hNsigmaSel['TPCTOF'][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SetName(\
                'hNsigmaSelTPCTOF_{0}_Prong{1}_Pt{2:.0f}_{3:.0f}'.format(spe, prong, ptmin, ptmax))
            SetObjectStyle(hNsigmaSel['TPCTOF'][spe][prong]['Pt{:.0f}_{:.0f}'.format(\
                ptmin, ptmax)], color=kRed, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.3, linewidth=1)

            sPIDNsigma.GetAxis(0).SetRange(-1, -1)
            sPIDNsigma.GetAxis(1).SetRange(-1, -1)

#2D distributions (prong0 vs prong2)
for det in axes:
    for spe in axes[det]:
        hNsigma[det][spe]['0-2'], hNsigmaSel[det][spe]['0-2'] = {}, {}
        for iPt, (ptmin, ptmax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):

            ptbinmin = sPIDNsigma.GetAxis(0).FindBin(ptmin*1.0001)
            ptbinmax = sPIDNsigma.GetAxis(0).FindBin(ptmax*0.9999)

            MLbinmin = sPIDNsigma.GetAxis(1).FindBin(MLmin[iPt]*1.0001)
            MLbinmax = sPIDNsigma.GetAxis(1).FindBin(MLmax[iPt]*0.9999)


            if det != 'Comb':
                sparse = sPIDNsigma
            else:
                sparse = sPIDNsigmaComb

            sparse.GetAxis(0).SetRange(ptbinmin, ptbinmax)

            hNsigma[det][spe]['0-2']['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                sparse.Projection(axes[det][spe]['0'], axes[det][spe]['2'])
            hNsigma[det][spe]['0-2']['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SetName(\
                'hNsigma{0}_{1}_Pront02_Pt{2:.0f}_{3:.0f}'.format(det, spe, ptmin, ptmax))
            SetObjectStyle(hNsigma[det][spe]['0-2']['Pt{:.0f}_{:.0f}'.format(\
                ptmin, ptmax)], color=kBlue, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.3, linewidth=1)

            sparse.GetAxis(1).SetRange(MLbinmin, MLbinmax)

            hNsigmaSel[det][spe]['0-2']['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                sparse.Projection(axes[det][spe]['0'], axes[det][spe]['2'])
            hNsigmaSel[det][spe]['0-2']['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SetName(\
                'hNsigmaSel{:s}_{:s}_Prong02_Pt{:.0f}_{:.0f}'.format(det, spe, ptmin, ptmax))
            SetObjectStyle(hNsigmaSel[det][spe]['0-2']['Pt{:.0f}_{:.0f}'.format(\
                ptmin, ptmax)], color=kRed, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.3, linewidth=1)

            sparse.GetAxis(0).SetRange(-1, -1)
            sparse.GetAxis(1).SetRange(-1, -1)


leg = TLegend(0.65, 0.76, 0.99, 0.89)
leg.SetTextSize(0.05)
leg.SetFillStyle(0)
leg.AddEntry(list(hNsigma['TPC']['Pi']['0'].values())[0], 'w/o ML sel', 'fp')
leg.AddEntry(list(hNsigmaSel['TPC']['Pi']['0'].values())[0], 'w/ ML sel', 'fp')

leg2D = TLegend(0.18, 0.76, 0.52, 0.89)
leg2D.SetTextSize(0.045)
leg2D.SetFillStyle(0)
leg2D.AddEntry(list(hNsigma['TPC']['Pi']['0'].values())[0], 'w/o ML sel', 'fp')
leg2D.AddEntry(list(hNsigmaSel['TPC']['Pi']['0'].values())[0], 'w/ ML sel', 'fp')

#plot on canvases and save on root file
listoffigures = []
outfile = TFile(ARGS.outFileName, 'recreate')
cNsigma, cNsigma02, cNsigmaNorm, dirNsigma = {}, {}, {}, {}
for iDet, det in enumerate(hNsigma):
    cNsigma[det] = {}
    cNsigmaNorm[det] = {}
    cNsigma02[det] = {}
    outfile.cd()
    dirNsigma[det] = TDirectoryFile(det, det)
    dirNsigma[det].Write()
    dirNsigma[det].cd()
    for iPt, (ptmin, ptmax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
            TCanvas('cNsigma{:s}_Pt{:.0f}_{:.0f}'.format(
                det, ptmin, ptmax), '', 1920, 1080)
        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Divide(3, 2)
        cNsigmaNorm[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
            TCanvas('cNsigmaNorm{:s}_Pt{:.0f}_{:.0f}'.format(
                det, ptmin, ptmax), '', 1920, 1080)
        cNsigmaNorm[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Divide(3, 2)
        if det != 'TPCTOF':
            cNsigma02[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
                TCanvas('cNsigma02{:s}_Pt{:.0f}_{:.0f}'.format(
                    det, ptmin, ptmax), '', 1920, 1080)
            cNsigma02[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Divide(2, 1)

        for iSpe, spe in enumerate(hNsigma[det]):
            for iProng, prong in enumerate(hNsigma[det][spe]):
                if prong != '0-2':
                    his = hNsigma[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)]
                    hissel = hNsigmaSel[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)]
                    his.GetXaxis().SetTitleSize(0.055)
                    hissel.GetXaxis().SetTitleSize(0.055)
                    his.GetXaxis().SetLabelSize(0.050)
                    hissel.GetXaxis().SetLabelSize(0.050)

                    if det != 'TPCTOF':
                        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).SetLogy()
                        his.GetYaxis().SetRangeUser(0.8, his.GetMaximum()*10)
                        his.SetTitle(\
                            ' %0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{N}_{#sigma}^{%s} (%s) prong%d;Entries'\
                             % (ptmin, ptmax, det, spe, iProng))
                        hissel.SetTitle(\
                            ' %0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{N}_{#sigma}^{%s} (%s) prong%d;Entries'\
                             % (ptmin, ptmax, det, spe, iProng))
                        his.DrawCopy('hist')
                        his.DrawCopy('Esame')
                        hissel.DrawCopy('histsame')
                        hissel.DrawCopy('Esame')
                        leg.Draw('same')
                        his.Write()
                        hissel.Write()
                    else:
                        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1)
                        his.SetTitle(\
                            ' %0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{N}_{#sigma}^{TOF} (%s) prong%d;#it{N}_{#sigma}^{TPC} (%s) prong%d'\
                             % (ptmin, ptmax, spe, iProng, spe, iProng))
                        hissel.SetTitle(\
                            ' %0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{N}_{#sigma}^{TOF} (%s) prong%d;#it{N}_{#sigma}^{TPC} (%s) prong%d'\
                             % (ptmin, ptmax, spe, iProng, spe, iProng))
                        his.DrawCopy()
                        hissel.DrawCopy('same')
                        leg2D.Draw('same')
                        his.Write()
                        hissel.Write()
                        continue

                    hisnorm = his.Clone(his.GetName() + '_Norm')
                    hisselnorm = hissel.Clone(hissel.GetName() + '_Norm')

                    if hisnorm.Integral() > 0:
                        hisnorm.Scale(1./hisnorm.Integral())
                        hisnorm.GetYaxis().SetTitle('Normalised entries')
                        hisnorm.GetYaxis().SetRangeUser(1.e-5, 1.)

                    if hisselnorm.Integral() > 0:
                        hisselnorm.Scale(1./hisselnorm.Integral())
                        hisselnorm.GetYaxis().SetTitle('Normalised entries')
                        hisselnorm.GetYaxis().SetRangeUser(1.e-5, 1.)

                    cNsigmaNorm[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).SetLogy()
                    hisnorm.DrawCopy('hist')
                    hisnorm.DrawCopy('Esame')
                    hisselnorm.DrawCopy('histsame')
                    hisselnorm.DrawCopy('Esame')
                    leg.Draw('same')

                    cNsigmaNorm[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).Modified()
                    cNsigmaNorm[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).Update()
                    hisnorm.Write()
                    hisselnorm.Write()

                cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).Modified()
                cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(3*iSpe+iProng+1).Update()

            if det != 'TPCTOF':
                his = hNsigma[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)]
                hissel = hNsigmaSel[det][spe][prong]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)]
                his.GetXaxis().SetTitleSize(0.055)
                hissel.GetXaxis().SetTitleSize(0.055)
                his.GetXaxis().SetLabelSize(0.050)
                hissel.GetXaxis().SetLabelSize(0.050)

                cNsigma02[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iSpe+1)
                his.SetTitle(\
                    ' %0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{N}_{#sigma}^{%s} (%s) prong0;#it{N}_{#sigma}^{%s} (%s) prong2'\
                    % (ptmin, ptmax, det, spe, det, spe))
                hissel.SetTitle(\
                    ' %0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{N}_{#sigma}^{%s} (%s) prong0;#it{N}_{#sigma}^{%s} (%s) prong2'\
                    % (ptmin, ptmax, det, spe, det, spe))
                his.DrawCopy()
                hissel.DrawCopy('same')
                leg2D.Draw('same')

        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Write()
        if det != 'TPCTOF':
            cNsigmaNorm[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Write()
            cNsigma02[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Write()

        #save also pdf

        cNsigma[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SaveAs('{:s}{:s}{:d}.png'.format(outFileNameWoExt, det, iPt))
        listoffigures.append('{:s}{:s}{:d}.png'.format(outFileNameWoExt, det, iPt))

        if det != 'TPCTOF':
            cNsigmaNorm[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SaveAs('{:s}{:s}{:d}Norm.png'.format(outFileNameWoExt, det, iPt))
            cNsigma02[det]['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SaveAs('{:s}{:s}{:d}02.png'.format(outFileNameWoExt, det, iPt))
            listoffigures.append('{:s}{:s}{:d}Norm.png'.format(outFileNameWoExt, det, iPt))
            listoffigures.append('{:s}{:s}{:d}02.png'.format(outFileNameWoExt, det, iPt))

outfile.Close()

# combine png files
figures = [Image.open(figures) for figures in listoffigures]
min_shape = sorted([(np.sum(fig.size), fig.size) for fig in figures])[0][1]
imgs_comb = np.vstack(figures)
imgs_comb = Image.fromarray(imgs_comb)
imgs_comb.save('{:s}.png'.format(outFileNameWoExt))

for fig in listoffigures:
    os.remove(fig)

if six.PY2:
    raw_input('Press enter to exit')
elif six.PY3:
    input('Press enter to exit')
