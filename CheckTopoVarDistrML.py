'''
python script for the projection of distributions of topological variables from THnSparses vs. ML output vs. pT
run: python CheckTopoVarDistrML.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import argparse
import string
import os
import six
import numpy as np
import yaml
from PIL import Image
from ROOT import TFile, TCanvas, TDirectoryFile, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlue, kRed, kFullCircle, kOpenCircle  # pylint: disable=import-error,no-name-in-module
from TaskFileLoader import LoadSparseFromTask
from StyleFormatter import SetObjectStyle, SetGlobalStyle

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
    if iFile == 0:
        sparseReco, sparseGen = LoadSparseFromTask(infilename, inputCfg)
    else:
        sparseRecoPart, sparseGenPart = LoadSparseFromTask(
            infilename, inputCfg)
        for sparsetype in sparseRecoPart:
            sparseReco[sparsetype].Add(sparseRecoPart[sparsetype])

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

topovars = {1: 'pt', 2: 'deltamKK', 3: 'decL', 5: 'normdecLxy',
            7: 'cospxy', 8: 'sigvtx', 10: 'abscospiKphi3', 11: 'd0d0exp'}

vartitles = {1: '#it{p}_{T} (GeV/#it{c})', 2: '#Delta#it{M}_{KK} (MeV/#it{c})',
             3: 'decay length (#mum #times 0.1)', 4: 'decay length xy (#mum #times 0.1)',
             5: 'normalised decay length xy', 6: 'cos(#theta_{p}) #times 100',
             7: 'cos(#theta_{p}^{xy}) #times 100', 8: 'sigma vertex (#mum #times 0.1)',
             10: '|cos^{3}(#piK#phi)|', 11: 'max norm #it{d}_{0}-#it{d}_{0}^{exp}', 12: 'impact parameter xy'}

outfile = TFile(ARGS.outFileName, 'recreate')
dir1D = TDirectoryFile('Distr1D', 'Distr1D')
dir1D.Write()
outfile.cd()
dirVsPt = TDirectoryFile('DistrVsPt', 'DistrVsPt')
dirVsPt.Write()
outfile.cd()
dirVsML = TDirectoryFile('DistrVsML', 'DistrVsML')
dirVsML.Write()

cVars, cVarsNorm, cVarsVsPt, cVarsVsML = {}, {}, {}, {}
hVars, hVarsSel, hVarsNorm, hVarsSelNorm, hVarsVsPt, hVarsSelVsPt, hVarsVsML, hVarsSelVsML = (
    {} for iDic in range(8))

listoffigures = []
for iPt, (ptmin, ptmax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
    cVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
        TCanvas('cVarsPt{:.0f}_{:.0f}'.format(ptmin, ptmax), '', 1920, 1080)
    cVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Divide(4, 2)
    cVarsNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
        TCanvas('cVarsNormPt{:.0f}_{:.0f}'.format(
            ptmin, ptmax), '', 1920, 1080)
    cVarsNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Divide(4, 2)
    cVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
        TCanvas('cVarsVsPtPt{:.0f}_{:.0f}'.format(
            ptmin, ptmax), '', 1920, 1080)
    cVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Divide(4, 2)
    cVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = \
        TCanvas('cVarsVsMLPt{:.0f}_{:.0f}'.format(
            ptmin, ptmax), '', 1920, 1080)
    cVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Divide(4, 2)

    hVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = {}
    hVarsNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = {}
    hVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = {}
    hVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = {}
    hVarsSel['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = {}
    hVarsSelNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = {}
    hVarsSelVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = {}
    hVarsSelVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)] = {}

    ptbinmin = sparseReco['RecoPrompt'].GetAxis(1).FindBin(ptmin*1.0001)
    ptbinmax = sparseReco['RecoPrompt'].GetAxis(1).FindBin(ptmax*0.9999)

    MLbinmin = sparseReco['RecoPrompt'].GetAxis(13).FindBin(MLmin[iPt]*1.0001)
    MLbinmax = sparseReco['RecoPrompt'].GetAxis(13).FindBin(MLmax[iPt]*0.9999)

    sparseReco['RecoPrompt'].GetAxis(1).SetRange(ptbinmin, ptbinmax)

    leg = TLegend(0.18, 0.76, 0.52, 0.89)
    leg.SetTextSize(0.05)
    leg.SetFillStyle(0)

    leg2D = TLegend(0.18, 0.17, 0.52, 0.4)
    leg2D.SetTextSize(0.05)
    leg2D.SetFillStyle(0)

    for iVar, axis in enumerate(topovars):

        hVars['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis)
        hVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetName(
            'hVarsSel_{:s}_Pt{:.0f}_{:.0f}'.format(topovars[axis], ptmin, ptmax))
        SetObjectStyle(hVars['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]], color=kBlue, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.5, linewidth=1)

        hVarsVsPt['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis, 1)
        hVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetName(
            'hVarsSel_{:s}_Pt{:.0f}_{:.0f}'.format(topovars[axis], ptmin, ptmax))
        SetObjectStyle(hVarsVsPt['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]], color=kBlue, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.3, linewidth=1)

        hVarsVsML['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis, 13)
        hVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetName(
            'hVarsSel_{:s}_Pt{:.0f}_{:.0f}'.format(topovars[axis], ptmin, ptmax))
        SetObjectStyle(hVarsVsML['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]], color=kBlue, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.3, linewidth=1)

        sparseReco['RecoPrompt'].GetAxis(13).SetRange(MLbinmin, MLbinmax)

        hVarsSel['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis)
        hVarsSel['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetName(
            'hVarsSel_{:s}_Pt{:.0f}_{:.0f}'.format(topovars[axis], ptmin, ptmax))
        SetObjectStyle(hVarsSel['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]], color=kRed, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.5, linewidth=1)

        hVarsSelVsPt['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis, 1)
        hVarsSelVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetName(
            'hVarsSelVsPt_{:s}_Pt{:.0f}_{:.0f}'.format(topovars[axis], ptmin, ptmax))
        SetObjectStyle(hVarsSelVsPt['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]], color=kRed, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.3, linewidth=1)

        hVarsSelVsML['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis, 13)
        hVarsSelVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetName(
            'hVarsSelVsML_{:s}_Pt{:.0f}_{:.0f}'.format(topovars[axis], ptmin, ptmax))
        SetObjectStyle(hVarsSelVsML['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]], color=kRed, linealpha=0.25, fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.3, linewidth=1)

        sparseReco['RecoPrompt'].GetAxis(13).SetRange(-1, -1)

        his = hVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]]
        hissel = hVarsSel['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]]
        his.SetTitle('%0.f < #it{p}_{T} < %0.f GeV/#it{c};%s;Entries' %
                     (ptmin, ptmax, vartitles[axis]))
        hissel.SetTitle(
            '%0.f < #it{p}_{T} < %0.f GeV/#it{c};%s;Entries' % (ptmin, ptmax, vartitles[axis]))
        his.GetXaxis().SetTitleSize(0.055)
        hissel.GetXaxis().SetTitleSize(0.055)
        his.GetXaxis().SetLabelSize(0.050)
        hissel.GetXaxis().SetLabelSize(0.050)
        his.GetYaxis().SetRangeUser(1., his.GetMaximum()*10)

        if iVar == 0:
            leg.AddEntry(hVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]], 'w/o ML sel', 'fp')
            leg.AddEntry(hVarsSel['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]], 'w/ ML sel', 'fp')
            leg2D.AddEntry(hVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]], 'w/o ML sel', 'fp')
            leg2D.AddEntry(hVarsSel['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]], 'w/ ML sel', 'fp')

        cVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).SetLogy()
        his.DrawCopy('hist')
        his.DrawCopy('Esame')
        hissel.DrawCopy('histsame')
        hissel.DrawCopy('Esame')
        leg.Draw('same')

        hisnorm = his.Clone(his.GetName() + '_Norm')
        if hisnorm.Integral() > 0:
            hisnorm.Scale(1./hisnorm.Integral())
            hisnorm.GetYaxis().SetRangeUser(1.e-5, 5.)
            hisnorm.GetYaxis().SetTitle('Normalised entries')
        hisselnorm = hissel.Clone(hissel.GetName() + '_Norm')
        if hisselnorm.Integral() > 0:
            hisselnorm.Scale(1./hisselnorm.Integral())
            hisselnorm.GetYaxis().SetTitle('Normalised entries')

        cVarsNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).SetLogy()
        hisnorm.DrawCopy('hist')
        hisnorm.DrawCopy('Esame')
        hisselnorm.DrawCopy('histsame')
        hisselnorm.DrawCopy('Esame')
        leg.Draw('same')

        dir1D.cd()
        cVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Write()
        cVarsNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Write()
        hVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].Write()
        hVarsSel['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].Write()
        hisnorm.Write()
        hisselnorm.Write()

        hVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetTitle(
            '%0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{p}_{T} (GeV/#it{c});%s' % (ptmin, ptmax, vartitles[axis]))
        hVarsSelVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetTitle(
            '%0.f < #it{p}_{T} < %0.f GeV/#it{c};#it{p}_{T} (GeV/#it{c});%s' % (ptmin, ptmax, vartitles[axis]))

        cVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1)
        hVarsVsPt['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].DrawCopy()
        hVarsSelVsPt['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].DrawCopy('same')
        if topovars[axis] == 'cospxy':
            leg2D.Draw('same')

        dirVsPt.cd()
        cVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Write()
        hVarsVsPt['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].Write()
        hVarsSelVsPt['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].Write()

        hVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetTitle(
            '%0.f < #it{p}_{T} < %0.f GeV/#it{c};ML output;%s' % (ptmin, ptmax, vartitles[axis]))
        hVarsSelVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].SetTitle(
            '%0.f < #it{p}_{T} < %0.f GeV/#it{c};ML output;%s' % (ptmin, ptmax, vartitles[axis]))
        hVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)][topovars[axis]].GetXaxis().SetNdivisions(505)

        cVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1)
        hVarsVsML['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].DrawCopy()
        hVarsSelVsML['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].DrawCopy('same')
        if topovars[axis] == 'cospxy':
            leg2D.Draw('same')

        dirVsML.cd()
        cVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].Write()
        hVarsVsML['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].Write()
        hVarsSelVsML['Pt{:.0f}_{:.0f}'.format(
            ptmin, ptmax)][topovars[axis]].Write()

        cVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).Modified()
        cVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).Update()
        cVarsNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).Modified()
        cVarsNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).Update()
        cVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).Modified()
        cVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).Update()
        cVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).Modified()
        cVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].cd(iVar+1).Update()

    #save PNG also
    cVars['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SaveAs('{:s}{:d}.png'.format(outFileNameWoExt, iPt))
    cVarsNorm['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SaveAs('{:s}Norm{:d}.png'.format(outFileNameWoExt, iPt))
    cVarsVsPt['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SaveAs('{:s}VsPt{:d}.png'.format(outFileNameWoExt, iPt))
    cVarsVsML['Pt{:.0f}_{:.0f}'.format(ptmin, ptmax)].SaveAs('{:s}VsML{:d}.png'.format(outFileNameWoExt, iPt))
    listoffigures.append('{:s}{:d}.png'.format(outFileNameWoExt, iPt))
    listoffigures.append('{:s}Norm{:d}.png'.format(outFileNameWoExt, iPt))
    listoffigures.append('{:s}VsPt{:d}.png'.format(outFileNameWoExt, iPt))
    listoffigures.append('{:s}VsML{:d}.png'.format(outFileNameWoExt, iPt))

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
