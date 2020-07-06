'''
python script for the check of the mass shaping effect
run: python CheckMassShaping.py cfgFileName.yml cutSetFileName.yml outputPath
                                [--tree] [--Dspecie specie] [--rebin rebin] [--fitfunc func]
'''

import os
import sys
import argparse
import yaml
from root_numpy import fill_hist
from ROOT import TCanvas, TLegend, TDatabasePDG, TH1F, TH2F, TF1, TLine, TGraph # pylint: disable=import-error,no-name-in-module
from ROOT import kFullCircle, kFullSquare, kRainBow, kRed, kAzure, kGray, kBlack # pylint: disable=import-error,no-name-in-module
from ROOT import RooStats # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
#pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.TaskFileLoader import LoadSingleSparseFromTask
from utils.DfUtils import LoadDfFromRootOrParquet
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle, DivideCanvas

SetGlobalStyle(palette=kRainBow, padbottommargin=0.14, padrightmargin=0.14,
               padleftmargin=0.14, padtopmargin=0.075, titleoffsety=1.4, opttitle=1)

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
parser.add_argument('cutSetFileName', metavar='text', default='cutset.yml',
                    help='config file name with root input files')
parser.add_argument('outputPath', metavar='text', default='outputPath', help='output path')
parser.add_argument('--tree', action='store_true', default=False,
                    help='flag for imput tree/dataframe instead of sparse')
parser.add_argument('--Dspecie', metavar='text', default='Ds',
                    help='D-meson species, options: Ds, Dplus')
parser.add_argument('--rebin', type=int, required=False, default=1, help='mass rebin (optional)')
parser.add_argument('--fitfunc', metavar='text', default='expo',
                    help='fit function for bkg distribution, options: expo and polN, with N integer number')
args = parser.parse_args()

if args.Dspecie == 'Ds':
    mD = TDatabasePDG.Instance().GetParticle(431).Mass()
    massTitle = '#it{M}(KK#pi) (GeV/#it{c}^{2})'
    mesonName = 'Ds'
    mesonTitle = 'D_{s}^{+}'
elif args.Dspecie == 'Dplus':
    mD = TDatabasePDG.Instance().GetParticle(411).Mass()
    massTitle = '#it{M}(K#pi#pi) (GeV/#it{c}^{2})'
    mesonName = 'Dplus'
    mesonTitle = 'D^{+}'
else:
    print('ERROR: the D-meson specie must be Ds or Dplus! Exit')
    sys.exit()

if args.fitfunc != 'expo' and 'pol' not in args.fitfunc:
    print('ERROR: the fit function must be expo or polN, with N integer number! Exit')
    sys.exit()

with open(args.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
inFileNames = inputCfg['filename']
if not isinstance(inFileNames, list):
    inFileNames = [inFileNames]

with open(args.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']

hMassVsML, hMassNoSel, hMassSel = ([] for _ in range(3))
if not args.tree: # data from sparse
    for iFile, infilename in enumerate(inFileNames):
        if iFile == 0:
            sparseBkg = LoadSingleSparseFromTask(infilename, inputCfg, 'sparsenameBkg')
        else:
            sparseBkg.Add(LoadSingleSparseFromTask(infilename, inputCfg, 'sparsenameBkg'))

    for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
        print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')
        ptBinMin = sparseBkg.GetAxis(cutVars['Pt']['axisnum']).FindBin(ptMin * 1.0001)
        ptBinMax = sparseBkg.GetAxis(cutVars['Pt']['axisnum']).FindBin(ptMax * 0.9999)
        sparseBkg.GetAxis(cutVars['Pt']['axisnum']).SetRange(ptBinMin, ptBinMax)
        hMassNoSel.append(sparseBkg.Projection(0))
        hMassNoSel[iPt].SetNameTitle(f'hMassNoSelPt{ptMin:.0f}_{ptMax:.0f}',
                                     f'{ptMin} < #it{{p}}_{{T}} < {ptMax} (GeV/#it{{c}});{massTitle};Counts')

        hMassVsML.append({})
        hMassVsML[iPt]['ML_output'] = sparseBkg.Projection(13, 0)
        hMassVsML[iPt].SetNameTitle(f'hMassVsMLPt{ptMin:.0f}_{ptMax:.0f}',
                                    f'{ptMin} < #it{{p}}_{{T}} < {ptMax} (GeV/#it{{c}});{massTitle};ML_output')

        for var in cutVars:
            if 'InvMass' in var:
                continue
            binMin = sparseBkg.GetAxis(cutVars[var]['axisnum']).FindBin(cutVars[var]['min'][iPt] * 1.0001)
            binMax = sparseBkg.GetAxis(cutVars[var]['axisnum']).FindBin(cutVars[var]['max'][iPt] * 0.9999)
            sparseBkg.GetAxis(cutVars[var]['axisnum']).SetRange(binMin, binMax)
        hMassSel.append(sparseBkg.Projection(0))
        hMassSel[iPt].SetNameTitle(f'hMassSelPt{ptMin:.0f}_{ptMax:.0f}',
                                   f'{ptMin} < #it{{p}}_{{T}} < {ptMax} (GeV/#it{{c}});{massTitle};Counts')

        for iAxis in range(sparseBkg.GetNdimensions()):
            sparseBkg.GetAxis(iAxis).SetRange(-1, -1)

else: # data from tree/dataframe
    dataFrameBkg = LoadDfFromRootOrParquet(inputCfg['tree']['filenameBkg'], inputCfg['tree']['dirname'],
                                           inputCfg['tree']['treename'])

    massBins = 500
    massLimLow = min(dataFrameBkg['inv_mass'])
    massLimHigh = max(dataFrameBkg['inv_mass'])

    # selections to be applied
    selToApply = []
    for iPt, _ in enumerate(cutVars['Pt']['min']):
        selToApply.append('')
        for varName in cutVars:
            if varName == 'InvMass':
                continue
            if selToApply[iPt] != '':
                selToApply[iPt] += ' & '
            selToApply[iPt] += \
                f"{cutVars[varName]['min'][iPt]}<{cutVars[varName]['name']}<{cutVars[varName]['max'][iPt]}"

    for iPt, (cuts, ptMin, ptMax) in enumerate(zip(selToApply, cutVars['Pt']['min'], cutVars['Pt']['max'])):
        print(f'Projecting distributions for {ptMin:.1f} < pT < {ptMax:.1f} GeV/c')

        hMassNoSel.append(TH1F(f'hMassNoSelPt{ptMin:.0f}_{ptMax:.0f}',
                               f'{ptMin} < #it{{p}}_{{T}} < {ptMax} (GeV/#it{{c}});{massTitle};Counts',
                               massBins, massLimLow, massLimHigh))
        hMassSel.append(TH1F(f'hMassSelPt{ptMin:.0f}_{ptMax:.0f}',
                             f'{ptMin} < #it{{p}}_{{T}} < {ptMax} (GeV/#it{{c}});{massTitle};Counts',
                             massBins, massLimLow, massLimHigh))

        hMassVsML.append({})
        for var in dataFrameBkg.columns:
            if 'ML' in var:
                hMassVsML[iPt][var] = TH2F(f'hMassVs{var}Pt{ptMin:.0f}_{ptMax:.0f}', f';{massTitle};{var}',
                                           massBins, massLimLow, massLimHigh,
                                           100, min(dataFrameBkg[var]), max(dataFrameBkg[var]))

        dataFrameBkgPtSel = dataFrameBkg.astype(float).query(f'{ptMin} < pt_cand < {ptMax}')
        fill_hist(hMassNoSel[iPt], dataFrameBkgPtSel['inv_mass'].values)
        for var in hMassVsML[iPt]:
            fill_hist(hMassVsML[iPt][var],
                      list(zip(dataFrameBkgPtSel['inv_mass'].values, dataFrameBkgPtSel[var].values)))
        dataFrameBkgSel = dataFrameBkg.astype(float).query(cuts)
        fill_hist(hMassSel[iPt], dataFrameBkgSel['inv_mass'].values)

hMassSel[iPt].Rebin(args.rebin)
hMassNoSel[iPt].Rebin(args.rebin)
SetObjectStyle(hMassNoSel[iPt], color=kAzure+4, markers=kFullSquare)
SetObjectStyle(hMassSel[iPt], color=kRed+1, markers=kFullCircle)

line = TLine(mD, 0, mD, hMassSel[0].GetMaximum())
line.SetLineColor(kGray+1)
line.SetLineWidth(1)
line.SetLineStyle(9)

lineVsML = {}
for var in hMassVsML[0]:
    lineVsML[var] = TLine(mD, hMassVsML[0][var].GetYaxis().GetBinLowEdge(1),
                          mD, hMassVsML[0][var].GetYaxis().GetBinUpEdge(hMassVsML[0][var].GetYaxis().GetNbins()))
    lineVsML[var].SetLineColor(kBlack)
    lineVsML[var].SetLineWidth(2)
    lineVsML[var].SetLineStyle(9)

leg = TLegend(0.18, 0.73, 0.4, 0.88)
leg.SetTextSize(0.05)
leg.SetFillStyle(0)

legFit = TLegend(0.18, 0.73, 0.78, 0.88)
legFit.SetTextSize(0.05)
legFit.SetFillStyle(0)

cMass, cMassFit, cMassVsML, fMass, gMassPValue, gMassNsigma = ([] for _ in range(6))
for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):

    minMass = hMassSel[iPt].GetBinLowEdge(1)
    maxMass = hMassSel[iPt].GetXaxis().GetBinUpEdge(hMassSel[iPt].GetNbinsX())
    hMassNoSelNorm = hMassNoSel[iPt].Clone()
    hMassSelNorm = hMassSel[iPt].Clone()
    hMassNoSelNorm.Scale(1./hMassNoSelNorm.Integral())
    hMassSelNorm.Scale(1./hMassSelNorm.Integral())
    lineNorm = line.Clone()
    lineNorm.SetY2(hMassSelNorm.GetBinContent(int(hMassSelNorm.GetNbinsX()/2))*1.5)

    if iPt == 0:
        leg.AddEntry(hMassNoSel[iPt], 'no selection', 'lp')
        leg.AddEntry(hMassSel[iPt], 'ML selection', 'lp')
        leg.AddEntry(lineNorm, f'{mesonTitle} mass', 'l')

    cMass.append(TCanvas(f'cMassPt{ptMin*10:.0f}_{ptMax*10:.0f}', '', 500, 500))
    cMass[iPt].cd().DrawFrame(minMass, 0., maxMass, hMassSelNorm.GetMaximum()*1.5, f';{massTitle}; Normalised counts')
    hMassNoSelNorm.DrawCopy('esame')
    hMassSelNorm.DrawCopy('esame')
    leg.Draw()
    lineNorm.Draw()
    cMass[iPt].Update()
    cMass[iPt].Modified()

    if len(hMassVsML[iPt]) > 1:
        cMassVsML.append(TCanvas(f'cMassVsML{ptMin*10:.0f}_{ptMax*10:.0f}', '', 1200, 400))
    else:
        cMassVsML.append(TCanvas(f'cMassVsML{ptMin*10:.0f}_{ptMax*10:.0f}', '', 500, 500))

    if len(hMassVsML[iPt]) > 1:
        DivideCanvas(cMassVsML[iPt], len(hMassVsML[iPt]))
        for iVar, var in enumerate(hMassVsML[iPt]):
            cMassVsML[iPt].cd(iVar+1).SetLogz()
            hMassVsML[iPt][var].DrawCopy('colz')
            lineVsML[var].Draw('same')
    else:
        cMassVsML[iPt].cd().SetLogz()
        hMassVsML[iPt][list(hMassVsML[iPt].keys())[0]].DrawCopy('colz')
    cMassVsML[iPt].Update()
    cMassVsML[iPt].Modified()

    fMass.append(TF1(f'fMassPt{ptMin*10:.0f}_{ptMax*10:.0f}', args.fitfunc, minMass, maxMass))
    fMass[iPt].SetLineColor(kBlack)

    if iPt == 0:
        legFit.AddEntry(hMassSel[iPt], 'MC bkg ML selected', 'lp')
        legFit.AddEntry(fMass[iPt], f'{args.fitfunc} fit', 'lp')
        legFit.AddEntry(line, f'{mesonTitle} mass', 'l')

    gMassPValue.append(TGraph(0))
    gMassNsigma.append(TGraph(0))
    SetObjectStyle(gMassPValue[iPt], color=kRed+1)
    SetObjectStyle(gMassNsigma[iPt], color=kRed+1)

    cMassFit.append(TCanvas(f'cMassFitPt{ptMin*10:.0f}_{ptMax*10:.0f}', '', 1200, 400))
    cMassFit[iPt].Divide(3, 1)
    hMassToFit = hMassSel[iPt].Clone()
    for iBin in range(hMassToFit.GetNbinsX()): #remove signal region
        if mD-0.02 < hMassToFit.GetBinCenter(iBin+1) < mD + 0.02:
            hMassToFit.SetBinContent(iBin+1, 0)
    cMassFit[iPt].cd(1).DrawFrame(minMass, 0., maxMass, hMassSel[iPt].GetMaximum()*1.5,
                                  f';{massTitle};Counts')
    hMassToFit.Fit(fMass[iPt], 'RE0')
    hMassSel[iPt].DrawCopy('esame')
    fMass[iPt].Draw('same')
    line.Draw()
    legFit.Draw()

    linePval = line.Clone()
    linePval.SetY2(1.)
    lineNsigma = line.Clone()
    lineNsigma.SetY1(-5.)
    lineNsigma.SetY2(5.)

    lineOneSigma = TLine(minMass, 1, maxMass, 1)
    lineOneSigma.SetLineColor(kBlack)
    lineOneSigma.SetLineWidth(1)
    lineOneSigma.SetLineStyle(2)

    lineMinusOneSigma = TLine(minMass, -1, maxMass, -1)
    lineMinusOneSigma.SetLineColor(kBlack)
    lineMinusOneSigma.SetLineWidth(1)
    lineMinusOneSigma.SetLineStyle(2)

    lineProbOneSigma = TLine(minMass, 0.1585, maxMass, 0.1585)
    lineProbOneSigma.SetLineColor(kBlack)
    lineProbOneSigma.SetLineWidth(1)
    lineProbOneSigma.SetLineStyle(2)

    lineProbTwoSigma = TLine(minMass, 0.0225, maxMass, 0.0225)
    lineProbTwoSigma.SetLineColor(kBlack)
    lineProbTwoSigma.SetLineWidth(1)
    lineProbTwoSigma.SetLineStyle(2)

    lineProbThreeSigma = TLine(minMass, 0.0015, maxMass, 0.0015)
    lineProbThreeSigma.SetLineColor(kBlack)
    lineProbThreeSigma.SetLineWidth(1)
    lineProbThreeSigma.SetLineStyle(2)

    lineAtZero = TLine(minMass, 0, maxMass, 0)
    lineAtZero.SetLineColor(kBlack)
    lineAtZero.SetLineWidth(1)
    lineAtZero.SetLineStyle(2)

    for iBin in range(hMassSel[iPt].GetNbinsX()):
        massCentre = hMassSel[iPt].GetBinCenter(iBin+1)
        obs = hMassSel[iPt].GetBinContent(iBin+1)
        exp = fMass[iPt].Eval(massCentre)
        unc = hMassSel[iPt].GetBinError(iBin+1)
        pval = RooStats.NumberCountingUtils.BinomialObsP(obs, exp, unc/obs)
        nsigma = (obs-exp)/unc
        gMassPValue[iPt].SetPoint(iBin, massCentre, pval)
        gMassNsigma[iPt].SetPoint(iBin, massCentre, nsigma)

    cMassFit[iPt].cd(2).DrawFrame(minMass, -5., maxMass, 5., f';{massTitle}; n#sigma')
    lineAtZero.Draw()
    lineOneSigma.Draw()
    lineMinusOneSigma.Draw()
    gMassNsigma[iPt].Draw('c')
    lineNsigma.Draw()

    cMassFit[iPt].cd(3).DrawFrame(minMass, gMassPValue[iPt].GetMinimum(), maxMass, 1., f';{massTitle}; p-value')
    cMassFit[iPt].cd(3).SetLogy()
    gMassPValue[iPt].Draw('c')
    linePval.Draw()
    lineProbOneSigma.Draw()
    lineProbTwoSigma.Draw()
    lineProbThreeSigma.Draw()
    cMassFit[iPt].Update()
    cMassFit[iPt].Modified()

# save outputs
if not os.path.isdir(args.outputPath):
    os.mkdir(args.outputPath)
for iPt, (ptMin, ptMax) in enumerate(zip(cutVars['Pt']['min'], cutVars['Pt']['max'])):
    cMass[iPt].SaveAs(os.path.join(args.outputPath, f'{mesonName}MassDistr_w_wo_selection_pT{ptMin}-{ptMax}.pdf'))
    cMassVsML[iPt].SaveAs(os.path.join(args.outputPath, f'{mesonName}Mass_vs_MLoutput_pT{ptMin}-{ptMax}.pdf'))
    cMassFit[iPt].SaveAs(os.path.join(args.outputPath,
                                      f'{mesonName}MassDistr_sel_fit_{args.fitfunc}_pT{ptMin}-{ptMax}.pdf'))

input('Press Enter to exit')
