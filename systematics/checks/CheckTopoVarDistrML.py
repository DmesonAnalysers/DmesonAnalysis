'''
python script for the projection of distributions of topological variables from THnSparses vs. ML output vs. pT
run: python CheckTopoVarDistrML.py cfgFileName.yml cutSetFileName.yml outFileName.root
'''

import sys
import argparse
import os
import numpy as np
import yaml
from PIL import Image
from ROOT import TFile, TCanvas, TDirectoryFile, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlue, kRed, kFullCircle, kOpenCircle  # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.TaskFileLoader import LoadSparseFromTask  #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle  #pylint: disable=wrong-import-position,import-error

# main function
SetGlobalStyle()

PARSER = argparse.ArgumentParser(description='Arguments to pass')
PARSER.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
PARSER.add_argument('cutSetFileName', metavar='text', default='cutSetFileName.yml', help='input file with cut set')
PARSER.add_argument('outFileName', metavar='text', default='outFileName.root', help='output root file name')
ARGS = PARSER.parse_args()
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
    sys.exit()

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
hVars, hVarsSel, hVarsNorm, hVarsSelNorm, hVarsVsPt, hVarsSelVsPt, hVarsVsML, hVarsSelVsML = ({} for iDic in range(8))

listoffigures = []
for iPt, (ptmin, ptmax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
    cVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = TCanvas(f'cVarsPt{ptmin:.0f}_{ptmax:.0f}', '', 1920, 1080)
    cVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'].Divide(4, 2)
    cVarsNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = TCanvas(f'cVarsNormPt{ptmin:.0f}_{ptmax:.0f}', '', 1920, 1080)
    cVarsNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'].Divide(4, 2)
    cVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = TCanvas(f'cVarsVsPtPt{ptmin:.0f}_{ptmax:.0f}', '', 1920, 1080)
    cVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'].Divide(4, 2)
    cVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = TCanvas(f'cVarsVsMLPt{ptmin:.0f}_{ptmax:.0f}', '', 1920, 1080)
    cVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'].Divide(4, 2)

    hVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = {}
    hVarsNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = {}
    hVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = {}
    hVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = {}
    hVarsSel[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = {}
    hVarsSelNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = {}
    hVarsSelVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = {}
    hVarsSelVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'] = {}

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

        hVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis)
        hVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetName(
            f'hVarsSel_{topovars[axis]}_Pt{ptmin:.0f}_{ptmax:.0f}')
        SetObjectStyle(hVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], color=kBlue, linealpha=0.25,
                       fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.5, linewidth=1)

        hVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis, 1)
        hVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetName(
            f'hVarsSel_{topovars[axis]}_Pt{ptmin:.0f}_{ptmax:.0f}')
        SetObjectStyle(hVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], color=kBlue, linealpha=0.25,
                       fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.3, linewidth=1)

        hVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis, 13)
        hVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetName(
            f'hVarsSel_{topovars[axis]}_Pt{ptmin:.0f}_{ptmax:.0f}')
        SetObjectStyle(hVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], color=kBlue, linealpha=0.25,
                       fillalpha=0.25, markeralpha=1, markerstyle=kFullCircle, markersize=0.3, linewidth=1)

        sparseReco['RecoPrompt'].GetAxis(13).SetRange(MLbinmin, MLbinmax)

        hVarsSel[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis)
        hVarsSel[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetName(
            f'hVarsSel_{topovars[axis]}_Pt{ptmin:.0f}_{ptmax:.0f}')
        SetObjectStyle(hVarsSel[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], color=kRed, linealpha=0.25,
                       fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.5, linewidth=1)

        hVarsSelVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis, 1)
        hVarsSelVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetName(
            f'hVarsSelVsPt_{topovars[axis]}_Pt{ptmin:.0f}_{ptmax:.0f}')
        SetObjectStyle(hVarsSelVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], color=kRed, linealpha=0.25,
                       fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.3, linewidth=1)

        hVarsSelVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]] = sparseReco['RecoPrompt'].Projection(axis, 13)
        hVarsSelVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetName(
            f'hVarsSelVsML_{topovars[axis]}_Pt{ptmin:.0f}_{ptmax:.0f}')
        SetObjectStyle(hVarsSelVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], color=kRed, linealpha=0.25,
                       fillalpha=0.25, markeralpha=1, markerstyle=kOpenCircle, markersize=0.3, linewidth=1)

        sparseReco['RecoPrompt'].GetAxis(13).SetRange(-1, -1)

        his = hVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]]
        hissel = hVarsSel[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]]
        his.SetTitle('%0.f < #it{p}_{T} < %0.f GeV/#it{c};%s;Entries' % (ptmin, ptmax, vartitles[axis]))
        hissel.SetTitle('%0.f < #it{p}_{T} < %0.f GeV/#it{c};%s;Entries' % (ptmin, ptmax, vartitles[axis]))
        his.GetXaxis().SetTitleSize(0.055)
        hissel.GetXaxis().SetTitleSize(0.055)
        his.GetXaxis().SetLabelSize(0.050)
        hissel.GetXaxis().SetLabelSize(0.050)
        his.GetYaxis().SetRangeUser(1., his.GetMaximum()*10)

        if iVar == 0:
            leg.AddEntry(hVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], 'w/o ML sel', 'fp')
            leg.AddEntry(hVarsSel[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], 'w/ ML sel', 'fp')
            leg2D.AddEntry(hVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], 'w/o ML sel', 'fp')
            leg2D.AddEntry(hVarsSel[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]], 'w/ ML sel', 'fp')

        cVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).SetLogy()
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

        cVarsNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).SetLogy()
        hisnorm.DrawCopy('hist')
        hisnorm.DrawCopy('Esame')
        hisselnorm.DrawCopy('histsame')
        hisselnorm.DrawCopy('Esame')
        leg.Draw('same')

        dir1D.cd()
        cVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'].Write()
        cVarsNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'].Write()
        hVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].Write()
        hVarsSel[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].Write()
        hisnorm.Write()
        hisselnorm.Write()

        hVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetTitle(
            f'{ptmin} < #it{{p}}_{{T}} < {ptmax} GeV/#it{{c}};#it{{p}}_{{T}} (GeV/#it{{c}});{vartitles[axis]}')
        hVarsSelVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetTitle(
            f'{ptmin} < #it{{p}}_{{T}} < {ptmax} GeV/#it{{c}};#it{{p}}_{{T}} (GeV/#it{{c}});{vartitles[axis]}')

        cVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1)
        hVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].DrawCopy()
        hVarsSelVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].DrawCopy('same')
        if topovars[axis] == 'cospxy':
            leg2D.Draw('same')

        dirVsPt.cd()
        cVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'].Write()
        hVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].Write()
        hVarsSelVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].Write()

        hVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetTitle(
            f'{ptmin} < #it{{p}}_{{T}} < {ptmax} GeV/#it{{c}};ML output;{vartitles[axis]}')
        hVarsSelVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].SetTitle(
            f'{ptmin} < #it{{p}}_{{T}} < {ptmax} GeV/#it{{c}};ML output;{vartitles[axis]}')
        hVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].GetXaxis().SetNdivisions(505)

        cVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1)
        hVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].DrawCopy()
        hVarsSelVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].DrawCopy('same')
        if topovars[axis] == 'cospxy':
            leg2D.Draw('same')

        dirVsML.cd()
        cVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'].Write()
        hVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].Write()
        hVarsSelVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'][topovars[axis]].Write()

        cVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).Modified()
        cVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).Update()
        cVarsNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).Modified()
        cVarsNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).Update()
        cVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).Modified()
        cVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).Update()
        cVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).Modified()
        cVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'].cd(iVar+1).Update()

    #save PNG also
    cVars[f'Pt{ptmin:.0f}_{ptmax:.0f}'].SaveAs(f'{outFileNameWoExt}{iPt}.png')
    cVarsNorm[f'Pt{ptmin:.0f}_{ptmax:.0f}'].SaveAs(f'{outFileNameWoExt}Norm{iPt}.png')
    cVarsVsPt[f'Pt{ptmin:.0f}_{ptmax:.0f}'].SaveAs(f'{outFileNameWoExt}VsPt{iPt}.png')
    cVarsVsML[f'Pt{ptmin:.0f}_{ptmax:.0f}'].SaveAs(f'{outFileNameWoExt}VsML{iPt}.png')
    listoffigures.append(f'{outFileNameWoExt}{iPt}.png')
    listoffigures.append(f'{outFileNameWoExt}Norm{iPt}.png')
    listoffigures.append(f'{outFileNameWoExt}VsPt{iPt}.png')
    listoffigures.append(f'{outFileNameWoExt}VsML{iPt}.png')

outfile.Close()

# combine png files
figures = [Image.open(figures) for figures in listoffigures]
min_shape = sorted([(np.sum(fig.size), fig.size) for fig in figures])[0][1]
imgs_comb = np.vstack(figures)
imgs_comb = Image.fromarray(imgs_comb)
imgs_comb.save(f'{outFileNameWoExt}.png')

for fig in listoffigures:
    os.remove(fig)

input('Press enter to exit')
