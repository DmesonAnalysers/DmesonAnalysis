'''
Script for the computation of the sysematic uncertainty
on the nonprompt fraction and nonprompt fraction ratio vs mult
run: python ComputeDataDrivenFracVsMultSyst.py config.yml
'''

import sys
import os
import argparse
import numpy as np
import yaml
from ROOT import TCanvas, TFile, TH1F, TGraphErrors, TLine, TBox, TLatex, TLegend, kRainBow, kRed, kBlack # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle, DivideCanvas #pylint: disable=wrong-import-position,import-error,no-name-in-module

SetGlobalStyle(padbottommargin=0.2, padleftmargin=0.2,
               padtopmargin=0.1, opttitle=1, titleoffset=1.6, maxdigits=3)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('config', metavar='text', default='config.yml',
                    help='yaml config file')
args = parser.parse_args()

with open(args.config, 'r') as ymlCfg: # pylint: disable=unspecified-encoding
    cfg = yaml.load(ymlCfg, yaml.FullLoader)

multClasses = cfg['inputs'].keys()
if 'MB' not in multClasses:
    print('ERROR: MB files must be included! Exit.')
    sys.exit()

hFracPrompt, hFracNonPrompt, hRatioFracPrompt, hRatioFracNonPrompt, corr, label = ({} for _ in range(6))
hDistrFracPrompt, hDistrFracNonPrompt, hDistrRatioFracPrompt, hDistrRatioFracNonPrompt = ({} for _ in range(4))
gFracPromptVsTrial, gFracNonPromptVsTrial, gRatioFracPromptVsTrial, gRatioFracNonPromptVsTrial = ({} for _ in range(4))
rmsPrompt, rmsNonPrompt, rmsRatioPrompt, rmsRatioNonPrompt = ({} for _ in range(4))
meanPrompt, meanNonPrompt, meanRatioPrompt, meanRatioNonPrompt = ({} for _ in range(4))
sysPrompt, sysNonPrompt, sysRatioPrompt, sysRatioNonPrompt = ({} for _ in range(4))
hSystPrompt, hSystNonPrompt, hSystRatioPrompt, hSystRatioNonPrompt = ({} for _ in range(4))
for iMult, mult in enumerate(multClasses):
    hFracPrompt[mult] = []
    hFracNonPrompt[mult] = []
    hDistrFracPrompt[mult] = []
    hDistrFracNonPrompt[mult] = []
    gFracPromptVsTrial[mult] = []
    gFracNonPromptVsTrial[mult] = []
    rmsPrompt[mult] = []
    rmsNonPrompt[mult] = []
    meanPrompt[mult] = []
    meanNonPrompt[mult] = []
    sysPrompt[mult] = []
    sysNonPrompt[mult] = []
    label[mult] = cfg['inputs'][mult]['label']
    if mult != 'MB':
        hRatioFracPrompt[mult] = []
        hRatioFracNonPrompt[mult] = []
        corr[mult] = cfg['inputs'][mult]['corr']
        hDistrRatioFracPrompt[mult] = []
        hDistrRatioFracNonPrompt[mult] = []
        gRatioFracPromptVsTrial[mult] = []
        gRatioFracNonPromptVsTrial[mult] = []
        rmsRatioPrompt[mult] = []
        rmsRatioNonPrompt[mult] = []
        meanRatioPrompt[mult] = []
        meanRatioNonPrompt[mult] = []
        sysRatioPrompt[mult] = []
        sysRatioNonPrompt[mult] = []
    inFileNames = cfg['inputs'][mult]['files']
    defFile = cfg['inputs'][mult]['default'].rsplit('/', 1)[1]
    defDir = cfg['inputs'][mult]['default'].rsplit('/', 1)[0]
    inFileNames.insert(0, defFile)
    for iVar, file in enumerate(inFileNames):
        if iVar == 0:
            inDir = defDir
        else:
            inDir = cfg['inputs'][mult]['directory']
        inFile = TFile.Open(os.path.join(inDir, file))
        hFracNonPrompt[mult].append(inFile.Get('hFDFracCorr'))
        hFracNonPrompt[mult][iVar].SetDirectory(0)
        hFracNonPrompt[mult][iVar].SetName(f'hFracNonPrompt_{mult}_{iVar}')
        hFracPrompt[mult].append(inFile.Get('hPromptFracCorr'))
        hFracPrompt[mult][iVar].SetDirectory(0)
        hFracPrompt[mult][iVar].SetName(f'hFracPrompt_{mult}_{iVar}')
        stepColor = int(np.round(40/len(inFileNames)))
        SetObjectStyle(hFracNonPrompt[mult][-1], color=kRainBow+stepColor*iVar)
        SetObjectStyle(hFracPrompt[mult][-1], color=kRainBow+stepColor*iVar)
        inFile.Close()
        if mult != 'MB':
            hRatioFracPrompt[mult].append(hFracPrompt[mult][iVar].Clone(f'hRatioFracPrompt_{mult}_{iVar}'))
            hRatioFracNonPrompt[mult].append(hFracNonPrompt[mult][iVar].Clone(f'hRatioFracNonPrompt_{mult}_{iVar}'))
            hRatioFracPrompt[mult][iVar].Divide(hFracPrompt['MB'][iVar])
            hRatioFracNonPrompt[mult][iVar].Divide(hFracNonPrompt['MB'][iVar])
            if corr[mult]:
                for iPt in range(1, hRatioFracPrompt[mult][iVar].GetNbinsX()+1):
                    uncMult = hFracPrompt[mult][iVar].GetBinError(iPt)
                    uncMB = hFracPrompt['MB'][iVar].GetBinError(iPt)
                    valmult = hFracPrompt[mult][iVar].GetBinContent(iPt)
                    valMB = hFracPrompt['MB'][iVar].GetBinContent(iPt)
                    rho = uncMult / uncMB
                    if rho > 1.:
                        rho = 1. / rho
                    uncRatio = np.sqrt(
                        uncMult**2 / valMB**2 + valmult**2 / valMB**4 * uncMB**2 -
                        2 * rho * valmult / valMB**3 * uncMult * uncMB
                    )
                    hRatioFracPrompt[mult][iVar].SetBinError(iPt, uncRatio)

                    uncMult = hFracNonPrompt[mult][iVar].GetBinError(iPt)
                    uncMB = hFracNonPrompt['MB'][iVar].GetBinError(iPt)
                    valmult = hFracNonPrompt[mult][iVar].GetBinContent(iPt)
                    valMB = hFracNonPrompt['MB'][iVar].GetBinContent(iPt)
                    rho = uncMult / uncMB
                    if rho > 1.:
                        rho = 1. / rho
                    uncRatio = np.sqrt(
                        uncMult**2 / valMB**2 + valmult**2 / valMB**4 * uncMB**2 -
                        2 * rho * valmult / valMB**3 * uncMult * uncMB
                    )
                    hRatioFracNonPrompt[mult][iVar].SetBinError(iPt, uncRatio)

    hSystPrompt[mult] = hFracPrompt[mult][0].Clone(f'hSystPrompt_mult{mult}')
    hSystNonPrompt[mult] = hFracNonPrompt[mult][0].Clone(f'hSystPrompt_mult{mult}')
    stepColor = int(np.round(50/len(multClasses)))
    SetObjectStyle(hSystPrompt[mult], color=kRainBow+stepColor*iMult, fillstyle=0)
    SetObjectStyle(hSystNonPrompt[mult], color=kRainBow+stepColor*iMult, fillstyle=0)
    if mult != 'MB':
        hSystRatioPrompt[mult] = hRatioFracPrompt[mult][0].Clone(f'hSystRatioPrompt_mult{mult}')
        hSystRatioNonPrompt[mult] = hRatioFracNonPrompt[mult][0].Clone(f'hSystRatioNonPrompt_mult{mult}')
        SetObjectStyle(hSystRatioPrompt[mult], color=kRainBow+stepColor*iMult, fillstyle=0)
        SetObjectStyle(hSystRatioNonPrompt[mult], color=kRainBow+stepColor*iMult, fillstyle=0)

    for iPt in range(hFracPrompt[mult][0].GetNbinsX()):
        ptMin = hFracPrompt[mult][0].GetBinLowEdge(iPt+1)
        ptMax = hFracPrompt[mult][0].GetXaxis().GetBinUpEdge(iPt+1)
        hDistrFracPrompt[mult].append(
            TH1F(
                f'hDistrFracPrompt_mult{mult}_ptbin{iPt}',
                f'{label[mult]}, {ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};'
                '#it{f}_{prompt} / #it{f}_{prompt}^{ default};',
                100, 0.95, 1.05
                )
            )
        hDistrFracNonPrompt[mult].append(
            TH1F(
                f'hDistrFracNonPrompt_mult{mult}_ptbin{iPt}',
                f'{label[mult]}, {ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};'
                '#it{f}_{non-prompt} / #it{f}_{non-prompt}^{ default};',
                100, 0.5, 1.5
                )
            )
        gFracPromptVsTrial[mult].append(TGraphErrors(0))
        gFracNonPromptVsTrial[mult].append(TGraphErrors(0))
        gFracPromptVsTrial[mult][-1].SetNameTitle(
            f'gFracPromptVsTrial_mult{mult}_ptbin{iPt}',
            f'{label[mult]}, {ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};trial;'
            '#it{f}_{prompt} / #it{f}_{prompt}^{ default}'
        )
        gFracNonPromptVsTrial[mult][-1].SetNameTitle(
            f'gFracNonPromptVsTrial_mult{mult}_ptbin{iPt}',
            f'{label[mult]}, {ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};trial;'
            '#it{f}_{non-prompt} / #it{f}_{non-prompt}^{ default}'
        )
        SetObjectStyle(hDistrFracPrompt[mult][-1])
        SetObjectStyle(hDistrFracNonPrompt[mult][-1])
        SetObjectStyle(gFracPromptVsTrial[mult][-1])
        SetObjectStyle(gFracNonPromptVsTrial[mult][-1])
        hDistrFracPrompt[mult][-1].GetXaxis().SetNdivisions(505)
        hDistrFracNonPrompt[mult][-1].GetXaxis().SetNdivisions(505)
        gFracPromptVsTrial[mult][-1].GetYaxis().SetNdivisions(505)
        gFracNonPromptVsTrial[mult][-1].GetYaxis().SetNdivisions(505)
        if mult != 'MB':
            hDistrRatioFracPrompt[mult].append(
                TH1F(
                    f'hDistrRatioFracPrompt_mult{mult}_ptbin{iPt}',
                    f'{label[mult]}, {ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};'
                    '#frac{#it{f}_{prompt}}{#it{f}_{prompt}^{MB}} / #frac{#it{f}_{prompt}}{#it{f}_{prompt}^{MB}}|_{default};',
                    100, 0.95, 1.05
                    )
                )
            hDistrRatioFracNonPrompt[mult].append(
                TH1F(
                    f'hDistrRatioFracNonPrompt_mult{mult}_ptbin{iPt}',
                    f'{label[mult]}, {ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};'
                    '#frac{#it{f}_{non-prompt}}{#it{f}_{non-prompt}^{MB}} / #frac{#it{f}_{non-prompt}}{#it{f}_{non-prompt}^{MB}}|_{default};',
                    100, 0.5, 1.5
                    )
                )
            gRatioFracPromptVsTrial[mult].append(TGraphErrors(0))
            gRatioFracNonPromptVsTrial[mult].append(TGraphErrors(0))
            gRatioFracPromptVsTrial[mult][-1].SetNameTitle(
                f'gRatioFracPromptVsTrial_mult{mult}_ptbin{iPt}',
                f'{label[mult]}, {ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};trial;'
                '#frac{#it{f}_{prompt}}{#it{f}_{prompt}^{MB}} / #frac{#it{f}_{prompt}}{#it{f}_{prompt}^{MB}}|_{default}'
            )
            gRatioFracNonPromptVsTrial[mult][-1].SetNameTitle(
                f'gRatioFracNonPromptVsTrial_mult{mult}_ptbin{iPt}',
                f'{label[mult]}, {ptMin} < #it{{p}}_{{T}} < {ptMax} GeV/#it{{c}};trial;'
                '#frac{#it{f}_{non-prompt}}{#it{f}_{non-prompt}^{MB}} / #frac{#it{f}_{non-prompt}}{#it{f}_{non-prompt}^{MB}}|_{default}'
            )
            SetObjectStyle(hDistrRatioFracPrompt[mult][-1])
            SetObjectStyle(hDistrRatioFracNonPrompt[mult][-1])
            SetObjectStyle(gRatioFracPromptVsTrial[mult][-1])
            SetObjectStyle(gRatioFracNonPromptVsTrial[mult][-1])
            hDistrRatioFracPrompt[mult][-1].GetXaxis().SetNdivisions(505)
            hDistrRatioFracNonPrompt[mult][-1].GetXaxis().SetNdivisions(505)
            gRatioFracPromptVsTrial[mult][-1].GetYaxis().SetNdivisions(505)
            gRatioFracNonPromptVsTrial[mult][-1].GetYaxis().SetNdivisions(505)

        for iVar, _ in enumerate(hFracPrompt[mult]):
            ratio = hFracPrompt[mult][iVar].GetBinContent(iPt+1) / hFracPrompt[mult][0].GetBinContent(iPt+1)
            hDistrFracPrompt[mult][iPt].Fill(ratio)
            gFracPromptVsTrial[mult][-1].SetPoint(iVar, iVar, ratio)
            ratio = hFracNonPrompt[mult][iVar].GetBinContent(iPt+1) / hFracNonPrompt[mult][0].GetBinContent(iPt+1)
            hDistrFracNonPrompt[mult][iPt].Fill(ratio)
            gFracNonPromptVsTrial[mult][-1].SetPoint(iVar, iVar, ratio)

            if mult != 'MB':
                ratio = hRatioFracPrompt[mult][iVar].GetBinContent(iPt+1) / hRatioFracPrompt[mult][0].GetBinContent(iPt+1)
                hDistrRatioFracPrompt[mult][iPt].Fill(ratio)
                gRatioFracPromptVsTrial[mult][iPt].SetPoint(iVar, iVar, ratio)
                ratio = hRatioFracNonPrompt[mult][iVar].GetBinContent(iPt+1) / hRatioFracNonPrompt[mult][0].GetBinContent(iPt+1)
                hDistrRatioFracNonPrompt[mult][iPt].Fill(ratio)
                gRatioFracNonPromptVsTrial[mult][iPt].SetPoint(iVar, iVar, ratio)

        rmsPrompt[mult].append(hDistrFracPrompt[mult][iPt].GetRMS())
        rmsNonPrompt[mult].append(hDistrFracNonPrompt[mult][iPt].GetRMS())
        meanPrompt[mult].append(abs(1-hDistrFracPrompt[mult][iPt].GetMean()))
        meanNonPrompt[mult].append(abs(1-hDistrFracNonPrompt[mult][iPt].GetMean()))
        sysPrompt[mult].append(np.sqrt(rmsPrompt[mult][iPt]**2 + meanPrompt[mult][iPt]**2))
        sysNonPrompt[mult].append(np.sqrt(rmsNonPrompt[mult][iPt]**2 + meanNonPrompt[mult][iPt]**2))
        hSystPrompt[mult].SetBinContent(iPt+1, sysPrompt[mult][iPt])
        hSystPrompt[mult].SetBinError(iPt+1, 0)
        hSystNonPrompt[mult].SetBinContent(iPt+1, sysNonPrompt[mult][iPt])
        hSystNonPrompt[mult].SetBinError(iPt+1, 0)
        if mult != 'MB':
            rmsRatioPrompt[mult].append(hDistrRatioFracPrompt[mult][iPt].GetRMS())
            rmsRatioNonPrompt[mult].append(hDistrRatioFracNonPrompt[mult][iPt].GetRMS())
            meanRatioPrompt[mult].append(abs(1-hDistrRatioFracPrompt[mult][iPt].GetMean()))
            meanRatioNonPrompt[mult].append(abs(1-hDistrRatioFracNonPrompt[mult][iPt].GetMean()))
            sysRatioPrompt[mult].append(np.sqrt(rmsRatioPrompt[mult][iPt]**2 + meanRatioPrompt[mult][iPt]**2))
            sysRatioNonPrompt[mult].append(np.sqrt(rmsRatioNonPrompt[mult][iPt]**2 + meanRatioNonPrompt[mult][iPt]**2))
            hSystRatioPrompt[mult].SetBinContent(iPt+1, sysRatioPrompt[mult][iPt])
            hSystRatioPrompt[mult].SetBinError(iPt+1, 0)
            hSystRatioNonPrompt[mult].SetBinContent(iPt+1, sysRatioNonPrompt[mult][iPt])
            hSystRatioNonPrompt[mult].SetBinError(iPt+1, 0)

# plots
lat = TLatex()
lat.SetNDC()
lat.SetTextSize(0.045)
lat.SetTextFont(42)
lat.SetTextColor(kBlack)

# plot fractions
cNonPromptFracs = TCanvas('cNonPromptFracs', '', 800, 800)
DivideCanvas(cNonPromptFracs, len(multClasses))
for iMult, mult in enumerate(multClasses):
    cNonPromptFracs.cd(iMult+1).DrawFrame(
        0.,
        0.,
        hFracNonPrompt[mult][-1].GetXaxis().GetBinUpEdge(hFracNonPrompt[mult][-1].GetNbinsX())+1,
        hFracNonPrompt[mult][-1].GetMaximum()*2,
        ';#it{p}_{T} (GeV/#it{c}); #it{f}_{non-prompt}'
    )
    for iVar, _ in enumerate(hFracNonPrompt[mult]):
        hFracNonPrompt[mult][iVar].Draw('same')
    lat.DrawLatex(0.25, 0.8, label[mult])
cNonPromptFracs.Modified()
cNonPromptFracs.Update()

cPromptFracs = TCanvas('cPromptFracs', '', 800, 800)
DivideCanvas(cPromptFracs, len(multClasses))
for iMult, mult in enumerate(multClasses):
    cPromptFracs.cd(iMult+1).DrawFrame(
        0.,
        hFracPrompt[mult][-1].GetMinimum()*0.75,
        hFracPrompt[mult][-1].GetXaxis().GetBinUpEdge(hFracPrompt[mult][-1].GetNbinsX())+1,
        1,
        ';#it{p}_{T} (GeV/#it{c}); #it{f}_{prompt}'
    )
    for iVar, _ in enumerate(hFracPrompt[mult]):
        hFracPrompt[mult][iVar].DrawCopy('same')
    lat.DrawLatex(0.25, 0.25, label[mult])
cPromptFracs.Modified()
cPromptFracs.Update()

cRatioNonPromptFracs = TCanvas('cRatioNonPromptFracs', '', 1920, 1080)
DivideCanvas(cRatioNonPromptFracs, len(multClasses)-1)
for iMult, mult in enumerate(list(multClasses)[1:]):
    cRatioNonPromptFracs.cd(iMult+1).DrawFrame(
        0.,
        0.,
        hRatioFracNonPrompt[mult][-1].GetXaxis().GetBinUpEdge(hRatioFracNonPrompt[mult][-1].GetNbinsX()) + 1,
        hRatioFracNonPrompt[mult][-1].GetMaximum()*2,
        ';#it{p}_{T} (GeV/#it{c}); #it{f}_{non-prompt} / #it{f}_{non-prompt}^{MB}'
    )
    for iVar, _ in enumerate(hRatioFracNonPrompt[mult]):
        hRatioFracNonPrompt[mult][iVar].DrawCopy('same')
    lat.DrawLatex(0.25, 0.8, label[mult])
cRatioNonPromptFracs.Modified()
cRatioNonPromptFracs.Update()

cRatioPromptFracs = TCanvas('cRatioPromptFracs', '', 1920, 1080)
DivideCanvas(cRatioPromptFracs, len(multClasses)-1)
for iMult, mult in enumerate(list(multClasses)[1:]):
    cRatioPromptFracs.cd(iMult+1).DrawFrame(
        0.,
        0.,
        hRatioFracPrompt[mult][-1].GetXaxis().GetBinUpEdge(hRatioFracPrompt[mult][-1].GetNbinsX()) + 1,
        hRatioFracPrompt[mult][-1].GetMaximum()*2,
        ';#it{p}_{T} (GeV/#it{c}); #it{f}_{non-prompt} / #it{f}_{non-prompt}^{MB}'
    )
    for iVar, _ in enumerate(hRatioFracPrompt[mult]):
        hRatioFracPrompt[mult][iVar].DrawCopy('same')
    lat.DrawLatex(0.25, 0.8, label[mult])
cRatioPromptFracs.Modified()
cRatioPromptFracs.Update()

# plot fraction distributions
cDistrFracPrompt, cDistrFracNonPrompt, cDistrRatioFracPrompt, cDistrRatioFracNonPrompt = ({} for _ in range(4))
lineVertPrompt, lineVertNonPrompt, lineVertRatioPrompt, lineVertRatioNonPrompt = ({} for _ in range(4))
boxVertPrompt, boxVertNonPrompt, boxVertRatioPrompt, boxVertRatioNonPrompt = ({} for _ in range(4))
for mult in hDistrFracPrompt:
    cDistrFracPrompt[mult] = TCanvas(f'cDistrFracPrompt_mult{mult}', '', 1920, 1080)
    cDistrFracNonPrompt[mult] = TCanvas(f'cDistrFracNonPrompt_mult{mult}', '', 1920, 1080)
    nPtBins = hFracPrompt[mult][0].GetNbinsX()
    DivideCanvas(cDistrFracPrompt[mult], nPtBins)
    DivideCanvas(cDistrFracNonPrompt[mult], nPtBins)
    lineVertPrompt[mult] = []
    lineVertNonPrompt[mult] = []
    boxVertPrompt[mult] = []
    boxVertNonPrompt[mult] = []
    if mult != 'MB':
        cDistrRatioFracPrompt[mult] = TCanvas(f'cDistrRatioFracPrompt_mult{mult}', '', 1920, 1080)
        cDistrRatioFracNonPrompt[mult] = TCanvas(f'cDistrRatioFracNonPrompt_mult{mult}', '', 1920, 1080)
        DivideCanvas(cDistrRatioFracPrompt[mult], nPtBins)
        DivideCanvas(cDistrRatioFracNonPrompt[mult], nPtBins)
        lineVertRatioPrompt[mult] = []
        lineVertRatioNonPrompt[mult] = []
        boxVertRatioPrompt[mult] = []
        boxVertRatioNonPrompt[mult] = []
    for iPt in range(nPtBins):
        cDistrFracPrompt[mult].cd(iPt+1)
        hDistrFracPrompt[mult][iPt].DrawCopy()
        maxEntries = hDistrFracPrompt[mult][iPt].GetMaximum()
        lineVertPrompt[mult].append(TLine(1., 0., 1., maxEntries))
        lineVertPrompt[mult][iPt].SetLineColor(kRed+1)
        lineVertPrompt[mult][iPt].SetLineWidth(2)
        lineVertPrompt[mult][iPt].Draw()
        boxVertPrompt[mult].append(
            TBox(1.-sysPrompt[mult][iPt], 0., 1.+sysPrompt[mult][iPt], maxEntries)
        )
        boxVertPrompt[mult][iPt].SetFillColorAlpha(kRed+1, 0.3)
        boxVertPrompt[mult][iPt].SetLineWidth(2)
        boxVertPrompt[mult][iPt].Draw()
        cDistrFracNonPrompt[mult].cd(iPt+1)
        hDistrFracNonPrompt[mult][iPt].DrawCopy()
        maxEntries = hDistrFracNonPrompt[mult][iPt].GetMaximum()
        lineVertNonPrompt[mult].append(TLine(1., 0., 1., maxEntries))
        lineVertNonPrompt[mult][iPt].SetLineColor(kRed+1)
        lineVertNonPrompt[mult][iPt].SetLineWidth(2)
        lineVertNonPrompt[mult][iPt].Draw()
        boxVertNonPrompt[mult].append(
            TBox(1.-sysNonPrompt[mult][iPt], 0., 1.+sysNonPrompt[mult][iPt], maxEntries)
        )
        boxVertNonPrompt[mult][iPt].SetFillColorAlpha(kRed+1, 0.3)
        boxVertNonPrompt[mult][iPt].SetLineWidth(2)
        boxVertNonPrompt[mult][iPt].Draw()
        if mult != 'MB':
            cDistrRatioFracPrompt[mult].cd(iPt+1)
            hDistrRatioFracPrompt[mult][iPt].DrawCopy()
            maxEntries = hDistrRatioFracPrompt[mult][iPt].GetMaximum()
            lineVertRatioPrompt[mult].append(TLine(1., 0., 1., maxEntries))
            lineVertRatioPrompt[mult][iPt].SetLineColor(kRed+1)
            lineVertRatioPrompt[mult][iPt].SetLineWidth(2)
            lineVertRatioPrompt[mult][iPt].Draw()
            boxVertRatioPrompt[mult].append(
                TBox(1.-sysRatioPrompt[mult][iPt], 0., 1.+sysRatioPrompt[mult][iPt], maxEntries)
            )
            boxVertRatioPrompt[mult][iPt].SetFillColorAlpha(kRed+1, 0.3)
            boxVertRatioPrompt[mult][iPt].SetLineWidth(2)
            boxVertRatioPrompt[mult][iPt].Draw()
            cDistrRatioFracNonPrompt[mult].cd(iPt+1)
            hDistrRatioFracNonPrompt[mult][iPt].DrawCopy()
            maxEntries = hDistrRatioFracNonPrompt[mult][iPt].GetMaximum()
            lineVertRatioNonPrompt[mult].append(TLine(1., 0., 1., maxEntries))
            lineVertRatioNonPrompt[mult][iPt].SetLineColor(kRed+1)
            lineVertRatioNonPrompt[mult][iPt].SetLineWidth(2)
            lineVertRatioNonPrompt[mult][iPt].Draw()
            boxVertRatioNonPrompt[mult].append(
                TBox(1.-sysRatioNonPrompt[mult][iPt], 0., 1.+sysRatioNonPrompt[mult][iPt], maxEntries)
            )
            boxVertRatioNonPrompt[mult][iPt].SetFillColorAlpha(kRed+1, 0.3)
            boxVertRatioNonPrompt[mult][iPt].SetLineWidth(2)
            boxVertRatioNonPrompt[mult][iPt].Draw()
    cDistrFracPrompt[mult].Modified()
    cDistrFracPrompt[mult].Update()
    cDistrFracNonPrompt[mult].Modified()
    cDistrFracNonPrompt[mult].Update()
    if mult != 'MB':
        cDistrRatioFracPrompt[mult].Modified()
        cDistrRatioFracPrompt[mult].Update()
        cDistrRatioFracNonPrompt[mult].Modified()
        cDistrRatioFracNonPrompt[mult].Update()

# plot fractions vs trial
cFracPromptVsTrial, cFracNonPromptVsTrial, cRatioFracPromptVsTrial, cRatioFracNonPromptVsTrial = ({} for _ in range(4))
lineHorPrompt, lineHorNonPrompt, lineHorRatioPrompt, lineHorRatioNonPrompt = ({} for _ in range(4))
boxHorPrompt, boxHorNonPrompt, boxHorRatioPrompt, boxHorRatioNonPrompt = ({} for _ in range(4))
for mult in gFracPromptVsTrial:
    cFracPromptVsTrial[mult] = TCanvas(f'cFracPromptVsTrial_mult{mult}', '', 1920, 1080)
    cFracNonPromptVsTrial[mult] = TCanvas(f'cFracNonPromptVsTrial_mult{mult}', '', 1920, 1080)
    nPtBins = hFracPrompt[mult][0].GetNbinsX()
    DivideCanvas(cFracPromptVsTrial[mult], nPtBins)
    DivideCanvas(cFracNonPromptVsTrial[mult], nPtBins)
    lineHorPrompt[mult] = []
    lineHorNonPrompt[mult] = []
    boxHorPrompt[mult] = []
    boxHorNonPrompt[mult] = []
    if mult != 'MB':
        cRatioFracPromptVsTrial[mult] = TCanvas(f'cRatioFracPromptVsTrial_mult{mult}', '', 1920, 1080)
        cRatioFracNonPromptVsTrial[mult] = TCanvas(f'cRatioFracNonPromptVsTrial_mult{mult}', '', 1920, 1080)
        DivideCanvas(cRatioFracPromptVsTrial[mult], nPtBins)
        DivideCanvas(cRatioFracNonPromptVsTrial[mult], nPtBins)
        lineHorRatioPrompt[mult] = []
        lineHorRatioNonPrompt[mult] = []
        boxHorRatioPrompt[mult] = []
        boxHorRatioNonPrompt[mult] = []
    for iPt in range(nPtBins):
        cFracPromptVsTrial[mult].cd(iPt+1)
        gFracPromptVsTrial[mult][iPt].GetYaxis().SetRangeUser(0.95, 1.05)
        gFracPromptVsTrial[mult][iPt].Draw('AP')
        nTrials = gFracPromptVsTrial[mult][iPt].GetN()+1
        lineHorPrompt[mult].append(TLine(0., 1., nTrials, 1.))
        lineHorPrompt[mult][iPt].SetLineColor(kRed+1)
        lineHorPrompt[mult][iPt].SetLineWidth(2)
        lineHorPrompt[mult][iPt].Draw()
        boxHorPrompt[mult].append(
            TBox(0., 1.-sysPrompt[mult][iPt], nTrials, 1.+sysPrompt[mult][iPt])
        )
        boxHorPrompt[mult][iPt].SetFillColorAlpha(kRed+1, 0.3)
        boxHorPrompt[mult][iPt].SetLineWidth(2)
        boxHorPrompt[mult][iPt].Draw()
        cFracNonPromptVsTrial[mult].cd(iPt+1)
        gFracNonPromptVsTrial[mult][iPt].GetYaxis().SetRangeUser(0.5, 1.5)
        gFracNonPromptVsTrial[mult][iPt].Draw('AP')
        lineHorNonPrompt[mult].append(TLine(0., 1., nTrials, 1.))
        lineHorNonPrompt[mult][iPt].SetLineColor(kRed+1)
        lineHorNonPrompt[mult][iPt].SetLineWidth(2)
        lineHorNonPrompt[mult][iPt].Draw()
        boxHorNonPrompt[mult].append(
            TBox(0., 1.-sysNonPrompt[mult][iPt], nTrials, 1.+sysNonPrompt[mult][iPt])
        )
        boxHorNonPrompt[mult][iPt].SetFillColorAlpha(kRed+1, 0.3)
        boxHorNonPrompt[mult][iPt].SetLineWidth(2)
        boxHorNonPrompt[mult][iPt].Draw()
        if mult != 'MB':
            cRatioFracPromptVsTrial[mult].cd(iPt+1)
            gRatioFracPromptVsTrial[mult][iPt].Draw('AP')
            gRatioFracPromptVsTrial[mult][iPt].GetYaxis().SetRangeUser(0.95, 1.05)
            lineHorRatioPrompt[mult].append(TLine(0., 1., nTrials, 1.))
            lineHorRatioPrompt[mult][iPt].SetLineColor(kRed+1)
            lineHorRatioPrompt[mult][iPt].SetLineWidth(2)
            lineHorRatioPrompt[mult][iPt].Draw()
            boxHorRatioPrompt[mult].append(
                TBox(0., 1.-sysRatioPrompt[mult][iPt], nTrials, 1.+sysRatioPrompt[mult][iPt])
            )
            boxHorRatioPrompt[mult][iPt].SetFillColorAlpha(kRed+1, 0.3)
            boxHorRatioPrompt[mult][iPt].SetLineWidth(2)
            boxHorRatioPrompt[mult][iPt].Draw()
            cRatioFracNonPromptVsTrial[mult].cd(iPt+1)
            gRatioFracNonPromptVsTrial[mult][iPt].GetYaxis().SetRangeUser(0.5, 1.5)
            gRatioFracNonPromptVsTrial[mult][iPt].Draw('AP')
            lineHorRatioNonPrompt[mult].append(TLine(0., 1., nTrials, 1.))
            lineHorRatioNonPrompt[mult][iPt].SetLineColor(kRed+1)
            lineHorRatioNonPrompt[mult][iPt].SetLineWidth(2)
            lineHorRatioNonPrompt[mult][iPt].Draw()
            boxHorRatioNonPrompt[mult].append(
                TBox(0., 1.-sysRatioNonPrompt[mult][iPt], nTrials, 1.+sysRatioNonPrompt[mult][iPt])
            )
            boxHorRatioNonPrompt[mult][iPt].SetFillColorAlpha(kRed+1, 0.3)
            boxHorRatioNonPrompt[mult][iPt].SetLineWidth(2)
            boxHorRatioNonPrompt[mult][iPt].Draw()
    cFracPromptVsTrial[mult].Modified()
    cFracPromptVsTrial[mult].Update()
    cFracNonPromptVsTrial[mult].Modified()
    cFracNonPromptVsTrial[mult].Update()
    if mult != 'MB':
        cRatioFracPromptVsTrial[mult].Modified()
        cRatioFracPromptVsTrial[mult].Update()
        cRatioFracNonPromptVsTrial[mult].Modified()
        cRatioFracNonPromptVsTrial[mult].Update()

# systematics
leg = TLegend(0.5, 0.85-0.04*len(multClasses), 0.85, 0.85)
leg.SetTextSize(0.04)
leg.SetFillStyle(0)
leg.SetFillStyle(0)

cSyst = TCanvas('cSyst', '', 800, 800)
cSyst.Divide(2, 2)
for iMult, mult in enumerate(multClasses):
    leg.AddEntry(hSystPrompt[mult], label[mult], 'l')
    if iMult == 0:
        cSyst.cd(1).DrawFrame(
            0.,
            0.,
            hSystPrompt[mult].GetXaxis().GetBinUpEdge(hSystPrompt[mult].GetNbinsX())+1,
            0.02,
            ';#it{p}_{T} (GeV/#it{c}); #sqrt{RMS^{2} + shift^{2}} #it{f}_{prompt}'
        )
    else:
        cSyst.cd(1)
    hSystPrompt[mult].Draw('hist same')

    if iMult == 0:
        cSyst.cd(2).DrawFrame(
            0.,
            0.,
            hSystNonPrompt[mult].GetXaxis().GetBinUpEdge(hSystNonPrompt[mult].GetNbinsX())+1,
            0.2,
            ';#it{p}_{T} (GeV/#it{c}); #sqrt{RMS^{2} + shift^{2}} #it{f}_{non-prompt}'
        )
    else:
        cSyst.cd(2)
    hSystNonPrompt[mult].Draw('hist same')
    if mult != 'MB':
        if iMult == 1:
            cSyst.cd(3).DrawFrame(
                0.,
                0.,
                hSystRatioPrompt[mult].GetXaxis().GetBinUpEdge(hSystRatioPrompt[mult].GetNbinsX())+1,
                0.02,
                ';#it{p}_{T} (GeV/#it{c}); #sqrt{RMS^{2} + shift^{2}} #it{f}_{prompt} / #it{f}_{prompt}^{ MB}'
            )
        else:
            cSyst.cd(3)
        hSystRatioPrompt[mult].Draw('hist same')

        if iMult == 1:
            cSyst.cd(4).DrawFrame(
                0.,
                0.,
                hSystRatioNonPrompt[mult].GetXaxis().GetBinUpEdge(hSystRatioNonPrompt[mult].GetNbinsX())+1,
                0.2,
                ';#it{p}_{T} (GeV/#it{c}); #sqrt{RMS^{2} + shift^{2}} #it{f}_{non-prompt} / #it{f}_{non-prompt}^{ MB}'
            )
        else:
            cSyst.cd(4)
        hSystRatioNonPrompt[mult].Draw('hist same')
cSyst.cd(1)
leg.Draw()

cSyst.Modified()
cSyst.Update()

outFileName = os.path.join(cfg['output']['directory'], cfg['output']['filewoext'])
cSyst.SaveAs(outFileName + '.pdf[')
cSyst.SaveAs(outFileName + '.pdf')
cNonPromptFracs.SaveAs(outFileName + '.pdf')
cPromptFracs.SaveAs(outFileName + '.pdf')
cRatioNonPromptFracs.SaveAs(outFileName + '.pdf')
cRatioPromptFracs.SaveAs(outFileName + '.pdf')
for iMult, mult in enumerate(multClasses):
    cFracPromptVsTrial[mult].SaveAs(outFileName + '.pdf')
    cFracNonPromptVsTrial[mult].SaveAs(outFileName + '.pdf')
    cDistrFracPrompt[mult].SaveAs(outFileName + '.pdf')
    cDistrFracNonPrompt[mult].SaveAs(outFileName + '.pdf')
    if mult != 'MB':
        cRatioFracPromptVsTrial[mult].SaveAs(outFileName + '.pdf')
        cRatioFracNonPromptVsTrial[mult].SaveAs(outFileName + '.pdf')
        cDistrRatioFracPrompt[mult].SaveAs(outFileName + '.pdf')
        cDistrRatioFracNonPrompt[mult].SaveAs(outFileName + '.pdf')
        if iMult == len(multClasses)-1:
            cDistrRatioFracNonPrompt[mult].SaveAs(outFileName + '.pdf]')

outFile = TFile.Open(outFileName + '.root', 'recreate')
cSyst.Write()
cNonPromptFracs.Write()
cPromptFracs.Write()
cRatioNonPromptFracs.Write()
cRatioPromptFracs.Write()
for mult in multClasses:
    cFracPromptVsTrial[mult].Write()
    cFracNonPromptVsTrial[mult].Write()
    cDistrFracPrompt[mult].Write()
    cDistrFracNonPrompt[mult].Write()
    if mult != 'MB':
        cRatioFracPromptVsTrial[mult].Write()
        cRatioFracNonPromptVsTrial[mult].Write()
        cDistrRatioFracPrompt[mult].Write()
        cDistrRatioFracNonPrompt[mult].Write()
    for iVar, _ in enumerate(hFracPrompt[mult]):
        hFracPrompt[mult][iVar].Write()
        hFracNonPrompt[mult][iVar].Write()
        if mult != 'MB':
            hRatioFracPrompt[mult][iVar].Write()
            hRatioFracNonPrompt[mult][iVar].Write()
    for iPt in range(hFracPrompt[mult][0].GetNbinsX()):
        hDistrFracPrompt[mult][iPt].Write()
        hDistrFracNonPrompt[mult][iPt].Write()
        gFracPromptVsTrial[mult][iPt].Write()
        gFracNonPromptVsTrial[mult][iPt].Write()
        if mult != 'MB':
            hDistrRatioFracPrompt[mult][iPt].Write()
            hDistrRatioFracNonPrompt[mult][iPt].Write()
            gRatioFracPromptVsTrial[mult][iPt].Write()
            gRatioFracNonPromptVsTrial[mult][iPt].Write()
outFile.Close()

input('Press enter to exit')
