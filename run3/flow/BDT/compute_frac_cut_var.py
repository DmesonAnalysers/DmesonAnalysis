import argparse
import yaml
import sys
import os
import numpy as np
from itertools import product
import ROOT
from ROOT import TFile, TCanvas, TLegend, TLatex, gROOT
from ROOT import TFile, TH1F, TH2F, TCanvas, TLegend, TGraphAsymmErrors, TLatex, gRandom, TF1  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure, kGreen, kRainBow # pylint: disable=import-error,no-name-in-module
from ROOT import kFullCircle, kFullSquare, kOpenSquare, kOpenCircle, kOpenCross, kOpenDiamond # pylint: disable=import-error,no-name-in-module
from os.path import exists
sys.path.append('../../../')
sys.path.append('..')
from flow_analysis_utils import get_cut_sets_config
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.AnalysisUtils import GetPromptFDYieldsAnalyticMinimisation, ApplyVariationToList

def compute_frac_cut_var(config, inputdir, outputdir, suffix, batch=False):

    gROOT.SetBatch(batch)

    CutSets, _, _, _, _ = get_cut_sets_config(config)
    nCutSets = max(CutSets)
    with open(config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)

    if os.path.exists(f'{inputdir}/eff'):
        effFiles = [f'{inputdir}/eff/{file}'
                    for file in os.listdir(f'{inputdir}/eff') if file.endswith('.root') and suffix in file]
    else:
        raise ValueError(f'No eff fodel found in {inputdir}')
    
    if os.path.exists(f'{inputdir}/ry'):
        rawYieldFiles = [f'{inputdir}/ry/{file}' 
                        for file in os.listdir(f'{inputdir}/ry') if file.endswith('.root') and suffix in file]
    else:
        raise ValueError(f'No ry folder found in {inputdir}')
    
    effFiles.sort()
    rawYieldFiles.sort()

    # load configuration
    ptmins = config['ptmins']
    ptmaxs = config['ptmaxs']

    #TODO: apply smearing

    hRawYields, hEffPrompt, hEffFD = [], [], []

    # load inputs raw yields and efficiencies
    hCrossSecPrompt, hCrossSecFD = [], []

    for inFileNameRawYield, inFileNameEff in zip(rawYieldFiles, effFiles):

        inFileRawYield = ROOT.TFile.Open(inFileNameRawYield)
        hRawYields.append(inFileRawYield.Get('hRawYieldsSimFit'))
        hRawYields[-1].SetDirectory(0)
        
        inFileEff = TFile.Open(inFileNameEff)
        hEffPrompt.append(inFileEff.Get('hEffPrompt'))
        hEffFD.append(inFileEff.Get('hEffFD'))
        hEffPrompt[-1].SetDirectory(0)
        hEffFD[-1].SetDirectory(0)

    # create output histograms
    hCorrYieldPrompt = hRawYields[0].Clone('hCorrYieldPrompt')
    hCorrYieldPrompt.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{N}_{prompt}')
    SetObjectStyle(hCorrYieldPrompt, color=kRed+1, fillcolor=kRed+1, markerstyle=kFullCircle)

    hCorrYieldFD = hRawYields[0].Clone('hCorrYieldFD')
    hCorrYieldFD.SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{N}_{non-prompt}')
    SetObjectStyle(hCorrYieldFD, color=kAzure+4, fillcolor=kAzure+4, markerstyle=kFullSquare)

    hCovCorrYields = [[hRawYields[0].Clone('hCovPromptPrompt'), hRawYields[0].Clone('hCovPromptFD')],
                    [hRawYields[0].Clone('hCovFDPrompt'), hRawYields[0].Clone('hCovFDFD')]]
    
    for iRow, row in enumerate(hCovCorrYields):
        for iCol, hCov in enumerate(row):
            SetObjectStyle(hCov, linecolor=kBlack)
            if iRow == 0:
                rowName = '#it{N}_{prompt}'
            else:
                rowName = '#it{N}_{non-prompt}'
            if iCol == 0:
                colName = '#it{N}_{prompt}'
            else:
                colName = '#it{N}_{non-prompt}'
            hCov.SetTitle(f';#it{{p}}_{{T}} (GeV/#it{{c}}); #sigma({rowName}, {colName})')

    SetGlobalStyle(padleftmargin=0.15, padtopmargin=0.08, titleoffsetx=1.,
                titleoffsety=1.4, opttitle=1, palette=kRainBow, maxdigits=2)

    legDistr = TLegend(0.45, 0.69, 0.75, 0.89)
    legDistr.SetFillStyle(0)
    legDistr.SetBorderSize(0)
    legDistr.SetTextSize(0.045)

    legEff = TLegend(0.2, 0.2, 0.4, 0.4)
    legEff.SetFillStyle(0)
    legEff.SetBorderSize(0)
    legEff.SetTextSize(0.045)

    legFrac = TLegend(0.2, 0.79, 0.4, 0.89)
    legFrac.SetFillStyle(0)
    legFrac.SetBorderSize(0)
    legFrac.SetTextSize(0.045)

    latInfo = TLatex()
    latInfo.SetNDC()
    latInfo.SetTextSize(0.045)
    latInfo.SetTextFont(42)
    latInfo.SetTextColor(1)

    hRawYieldsVsCut, hRawYieldsVsCutReSum, hRawYieldPromptVsCut, hRawYieldFDVsCut, cDistr = [], [], [], [], []
    hEffPromptVsCut, hEffFDVsCut, cEff = [], [], []
    hPromptFracVsCut, hFDFracVsCut, gPromptFracFcVsCut, gFDFracFcVsCut, gPromptFracNbVsCut, \
    gFDFracNbVsCut, cFrac = [], [], [], [], [], [], []
    hCorrMatrixCutSets, cCorrMatrix = [], []

    for iPt, (ptMin, ptMax) in enumerate(zip(ptmins, ptmaxs)):
        nSets = CutSets[iPt]

        listRawYield = [hRaw.GetBinContent(iPt+1) for hRaw in hRawYields[:nSets]]
        listRawYieldUnc = [hRaw.GetBinError(iPt+1) for hRaw in hRawYields[:nSets]]
        listEffPrompt = [hEffP.GetBinContent(iPt+1) for hEffP in hEffPrompt[:nSets]]
        listEffPromptUnc = [hEffP.GetBinError(iPt+1) for hEffP in hEffPrompt[:nSets]]
        listEffFD = [hEffF.GetBinContent(iPt+1) for hEffF in hEffFD[:nSets]]
        listEffFDUnc = [hEffF.GetBinError(iPt+1) for hEffF in hEffFD[:nSets]]

        print(f'Pt: {ptMin:.1f}-{ptMax:.1f}')
        for i in range(len(listEffPrompt)):
            print(f'Eff Prompt: {listEffPrompt[i]:.3f}    Eff FD: {listEffFD[i]:.3f}    Raw Yield: {listRawYield[i]:.2f}')

        corrYields, covMatrixCorrYields, chiSquare, matrices = \
            GetPromptFDYieldsAnalyticMinimisation(listEffPrompt, listEffFD, listRawYield, listEffPromptUnc, listEffFDUnc,
                                                listRawYieldUnc, config['minimisation']['correlated'])

        hCorrYieldPrompt.SetBinContent(iPt+1, corrYields.item(0))
        hCorrYieldPrompt.SetBinError(iPt+1, np.sqrt(covMatrixCorrYields.item(0, 0)))
        hCorrYieldFD.SetBinContent(iPt+1, corrYields.item(1))
        hCorrYieldFD.SetBinError(iPt+1, np.sqrt(covMatrixCorrYields.item(1, 1)))
        for covElem in product(range(2), range(2)):
            hCovCorrYields[covElem[0]][covElem[1]].SetBinContent(iPt+1, covMatrixCorrYields.item(covElem))
            hCovCorrYields[covElem[0]][covElem[1]].SetBinError(iPt+1, 0.)

        ptString = f'pT{ptMin:.1f}_{ptMax:.1f}'
        commonString = f'{ptMin:.1f} < #it{{p}}_{{T}} < {ptMax:.1f}  GeV/#it{{c}};cut set'
        hRawYieldsVsCut.append(TH1F(f'hRawYieldsVsCutPt_{ptString}', f'{commonString};raw yield', nSets, 0.5, nSets + 0.5))
        hRawYieldsVsCutReSum.append(TH1F(f'hRawYieldsVsCutReSum_{ptString}', f'{commonString};raw yield',
                                        nSets, 0.5, nSets + 0.5))
        hRawYieldPromptVsCut.append(TH1F(f'hRawYieldPromptVsCut_{ptString}', f'{commonString};raw yield',
                                        nSets, 0.5, nSets + 0.5))
        hRawYieldFDVsCut.append(TH1F(f'hRawYieldFDVsCut_{ptString}', f'{commonString};raw yield', nSets, 0.5, nSets + 0.5))
        hEffPromptVsCut.append(TH1F(f'hEffPromptVsCut_{ptString}', f'{commonString};efficiency', nSets, 0.5, nSets + 0.5))
        hEffFDVsCut.append(TH1F(f'hEffFDVsCut_{ptString}', f'{commonString};efficiency', nSets, 0.5, nSets + 0.5))
        hPromptFracVsCut.append(TH1F(f'hPromptFracVsCut_{ptString}', f'{commonString};#it{{f}}_{{prompt}}',
                                    nSets, 0.5, nSets + 0.5))
        hFDFracVsCut.append(TH1F(f'hFDFracVsCut_{ptString}', f'{commonString};#it{{f}}_{{FD}}', nSets, 0.5, nSets + 0.5))

        SetObjectStyle(hRawYieldsVsCut[iPt], linecolor=kBlack, markercolor=kBlack, markerstyle=kFullCircle)
        SetObjectStyle(hRawYieldPromptVsCut[iPt], color=kRed+1, fillcolor=kRed+1, markerstyle=kOpenCircle, fillalpha=0.3)
        SetObjectStyle(hRawYieldFDVsCut[iPt], color=kAzure+4, fillcolor=kAzure+4, markerstyle=kOpenSquare, fillalpha=0.3)
        SetObjectStyle(hRawYieldsVsCutReSum[iPt], linecolor=kGreen+2)
        SetObjectStyle(hEffPromptVsCut[iPt], color=kRed+1, markerstyle=kFullCircle)
        SetObjectStyle(hEffFDVsCut[iPt], color=kAzure+4, markerstyle=kFullSquare)
        SetObjectStyle(hPromptFracVsCut[iPt], color=kRed+1, markerstyle=kFullCircle)
        SetObjectStyle(hFDFracVsCut[iPt], color=kAzure+4, markerstyle=kFullSquare)

        hCorrMatrixCutSets.append(TH2F(f'hCorrMatrixCutSets_{ptString}', f'{commonString};cut set',
                                    nSets, 0.5, nSets + 0.5, nSets, 0.5, nSets + 0.5))
        for mEl in product(range(nSets), range(nSets)):
            hCorrMatrixCutSets[iPt].SetBinContent(mEl[0]+1, mEl[1]+1, matrices['corrMatrix'].item(mEl[0], mEl[1]))

        for iCutSet, (rawY, effP, effF, rawYunc, effPunc, effFunc) in enumerate(zip(listRawYield, listEffPrompt, listEffFD,
                                                                                listRawYieldUnc, listEffPromptUnc,
                                                                                listEffFDUnc)):
            # efficiency
            hEffPromptVsCut[iPt].SetBinContent(iCutSet+1, effP)
            hEffPromptVsCut[iPt].SetBinError(iCutSet+1, effPunc)
            hEffFDVsCut[iPt].SetBinContent(iCutSet+1, effF)
            hEffFDVsCut[iPt].SetBinError(iCutSet+1, effFunc)

            # raw yields (including prompt and non-prompt raw yields)
            hRawYieldsVsCut[iPt].SetBinContent(iCutSet+1, rawY)
            hRawYieldsVsCut[iPt].SetBinError(iCutSet+1, rawYunc)
            hRawYieldPromptVsCut[iPt].SetBinContent(iCutSet+1, corrYields.item(0) * effP)
            hRawYieldPromptVsCut[iPt].SetBinError(iCutSet+1, np.sqrt(covMatrixCorrYields.item(0, 0)) * effP)
            hRawYieldFDVsCut[iPt].SetBinContent(iCutSet+1, corrYields.item(1) * effF)
            hRawYieldFDVsCut[iPt].SetBinError(iCutSet+1, np.sqrt(covMatrixCorrYields.item(1, 1)) * effF)
            hRawYieldsVsCutReSum[iPt].SetBinContent(iCutSet+1, hRawYieldPromptVsCut[iPt].GetBinContent(iCutSet+1) +
                                                    hRawYieldFDVsCut[iPt].GetBinContent(iCutSet+1))

            # prompt fraction
            fPrompt = effP * corrYields.item(0) / (effP * corrYields.item(0) + effF * corrYields.item(1))
            defPdeNP = (effP * (effP * corrYields.item(0) + effF * corrYields.item(1)) - effP**2
                        * corrYields.item(0)) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
            defPdeNF = - effF * effP * corrYields.item(0) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
            fPromptUnc = np.sqrt(defPdeNP**2 * covMatrixCorrYields.item(0, 0) +
                                defPdeNF**2 * covMatrixCorrYields.item(1, 1) +
                                2 * defPdeNP * defPdeNF * covMatrixCorrYields.item(1, 0))

            # feed-down fraction
            fFD = effF * corrYields.item(1) / (effP * corrYields.item(0) + effF * corrYields.item(1))
            defFdeNF = (effF * (effF * corrYields.item(1) + effP * corrYields.item(0)) - effF**2
                        * corrYields.item(1)) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
            defFdeNP = - effF * effP * corrYields.item(1) / (effP * corrYields.item(0) + effF * corrYields.item(1))**2
            fFDUnc = np.sqrt(defFdeNF**2 * covMatrixCorrYields.item(1, 1) +
                            defFdeNP**2 * covMatrixCorrYields.item(0, 0) +
                            2 * defFdeNF * defFdeNP * covMatrixCorrYields.item(1, 0))

            hPromptFracVsCut[iPt].SetBinContent(iCutSet+1, fPrompt)
            hPromptFracVsCut[iPt].SetBinError(iCutSet+1, fPromptUnc)
            hFDFracVsCut[iPt].SetBinContent(iCutSet+1, fFD)
            hFDFracVsCut[iPt].SetBinError(iCutSet+1, fFDUnc)

        if iPt == 0:
            legDistr.AddEntry(hRawYieldsVsCut[iPt], 'Measured raw yield', 'lpe')
            legDistr.AddEntry(hRawYieldPromptVsCut[iPt], 'Prompt', 'f')
            legDistr.AddEntry(hRawYieldFDVsCut[iPt], 'Non-prompt', 'f')
            legDistr.AddEntry(hRawYieldsVsCutReSum[iPt], 'Prompt + non-prompt', 'l')
            legEff.AddEntry(hEffPromptVsCut[iPt], 'Prompt', 'lpe')
            legEff.AddEntry(hEffFDVsCut[iPt], 'Non-prompt', 'lpe')
            legFrac.AddEntry(hPromptFracVsCut[iPt], 'Prompt', 'lpe')
            legFrac.AddEntry(hFDFracVsCut[iPt], 'Non-prompt', 'lpe')

        cEff.append(TCanvas(f'cEff_{ptString}', '', 800, 800))
        cEff[iPt].DrawFrame(0.5, hEffPromptVsCut[iPt].GetMinimum()/5, nSets + 0.5, 1., f'{commonString};efficiency')
        cEff[iPt].SetLogy()
        hEffPromptVsCut[iPt].DrawCopy('same')
        hEffFDVsCut[iPt].DrawCopy('same')
        legEff.Draw()

        cDistr.append(TCanvas(f'cDistr_{ptString}', '', 800, 800))
        hFrameDistr = cDistr[iPt].DrawFrame(0.5, 0., nSets + 0.5, hRawYieldsVsCut[iPt].GetMaximum() * 1.2,
                                            f'{commonString};raw yield')
        hFrameDistr.GetYaxis().SetDecimals()
        hRawYieldsVsCut[iPt].Draw('same')
        hRawYieldPromptVsCut[iPt].DrawCopy('histsame')
        hRawYieldFDVsCut[iPt].DrawCopy('histsame')
        hRawYieldsVsCutReSum[iPt].Draw('same')
        legDistr.Draw()
        latInfo.DrawLatex(0.47, 0.65, f'#chi^{{2}} / ndf = {chiSquare:.3f}')

        cFrac.append(TCanvas(f'cFrac_{ptString}', '', 800, 800))
        cFrac[iPt].DrawFrame(0.5, 0., nSets + 0.5, 1.8, f'{commonString};fraction')
        hPromptFracVsCut[iPt].DrawCopy('Esame')
        hFDFracVsCut[iPt].DrawCopy('Esame')
        legFrac.Draw()

        cCorrMatrix.append(TCanvas(f'cCorrMatrix_{ptString}', '', 800, 800))
        cCorrMatrix[-1].cd().SetRightMargin(0.14)
        hCorrMatrixCutSets[iPt].Draw('colz')

    nPtBins = hCorrYieldPrompt.GetNbinsX()
    cCorrYield = TCanvas('cCorrYield', '', 800, 800)
    cCorrYield.DrawFrame(hCorrYieldPrompt.GetBinLowEdge(1), 1.,
                        hCorrYieldPrompt.GetBinLowEdge(nPtBins) + hCorrYieldPrompt.GetBinWidth(nPtBins),
                        hCorrYieldPrompt.GetMaximum() * 1.2, ';#it{p}_{T} (GeV/#it{c});corrected yield')
    cCorrYield.SetLogy()
    hCorrYieldPrompt.Draw('same')
    hCorrYieldFD.Draw('same')
    legEff.Draw()

    os.makedirs(f'{outputdir}/CutVarFrac', exist_ok=True)
    outFileName = f'{outputdir}/CutVarFrac/CutVarFrac_{suffix}.root'
    outFile = TFile(outFileName, 'recreate')
    cCorrYield.Write()
    hCorrYieldPrompt.Write()
    hCorrYieldFD.Write()
    for covElem in product(range(2), range(2)):
        hCovCorrYields[covElem[0]][covElem[1]].Write()
    for iPt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs)):
        outFile.mkdir(f"pt{ptmin:.1f}_{ptmax:.1f}")
        outFile.cd(f"pt{ptmin:.1f}_{ptmax:.1f}")
        print(f"Writing to dir pt{ptmin:.1f}_{ptmax:.1f} of file {outFile}")
        cDistr[iPt].Write()
        cEff[iPt].Write()
        cFrac[iPt].Write()
        hRawYieldsVsCut[iPt].Write()
        hRawYieldPromptVsCut[iPt].Write()
        hRawYieldFDVsCut[iPt].Write()
        hRawYieldsVsCutReSum[iPt].Write()
        hEffPromptVsCut[iPt].Write()
        hEffFDVsCut[iPt].Write()
        hPromptFracVsCut[iPt].Write()
        hFDFracVsCut[iPt].Write()
        hCorrMatrixCutSets[iPt].Write()
    outFile.Close()

    outFileNameEffPDF = outFileName.replace('.root', '_Eff.pdf')
    outFileNameDistrPDF = outFileName.replace('.root', '_Distr.pdf')
    outFileNameFracPDF = outFileName.replace('.root', '_Frac.pdf')
    outFileNameCorrMatrixPDF = outFileName.replace('.root', '_CorrMatrix.pdf')
    for iPt in range(len(ptmins)):
        if iPt == 0:
            cEff[iPt].SaveAs(f'{outFileNameEffPDF}[')
            cDistr[iPt].SaveAs(f'{outFileNameDistrPDF}[')
            cFrac[iPt].SaveAs(f'{outFileNameFracPDF}[')
            cCorrMatrix[iPt].SaveAs(f'{outFileNameCorrMatrixPDF}[')
        cEff[iPt].SaveAs(outFileNameEffPDF)
        cDistr[iPt].SaveAs(outFileNameDistrPDF)
        cFrac[iPt].SaveAs(outFileNameFracPDF)
        cCorrMatrix[iPt].SaveAs(outFileNameCorrMatrixPDF)
        if iPt == hRawYields[0].GetNbinsX() - 1:
            cEff[iPt].SaveAs(f'{outFileNameEffPDF}]')
            cDistr[iPt].SaveAs(f'{outFileNameDistrPDF}]')
            cFrac[iPt].SaveAs(f'{outFileNameFracPDF}]')
            cCorrMatrix[iPt].SaveAs(f'{outFileNameCorrMatrixPDF}]')
    
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="flow configuration file")
    parser.add_argument('inputdir', metavar='text',
                        default='path/to/eff/proj_mc', help='input path')
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--batch", "-b", action="store_true",
                        help="run in batch mode")
    args = parser.parse_args()

    compute_frac_cut_var(args.config, args.inputdir, args.outputdir, args.suffix, args.batch)