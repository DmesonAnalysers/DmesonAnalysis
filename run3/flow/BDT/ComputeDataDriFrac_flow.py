'''
python script for the computation of the fractions of prompt and feed-down D for all cutset
run: python ComputeDataDrivenFraction.py --inputdir path/to/input --outputdir path/to/output --suffix text
'''

import argparse
import os
import sys
sys.path.append('../../../')
from ROOT import TFile, TCanvas, TLegend, gROOT, kRed, kBlue  # pylint: disable=import-error,no-name-in-module
from utils.AnalysisUtils import GetPromptFDFractionCutSet
from utils.StyleFormatter import SetGlobalStyle

def data_driven_frac(outputdir, suffix, iFile, hEffPrompt, hEffFD, \
                        hPromptFrac, hFDFrac, hPromptFracCorr, hFDFracCorr, \
                        hCorrYieldPrompt, hCorrYieldFD, hCovPromptPrompt, hCovPromptFD, hCovFDFD):
    for iPt in range(hEffPrompt.GetNbinsX()):
        ptMin = hEffPrompt.GetBinLowEdge(iPt+1)
        ptMax = ptMin+hEffPrompt.GetBinWidth(iPt+1)
        ptCent = hEffPrompt.GetBinCenter(iPt+1)
        effAccPrompt = hEffPrompt.GetBinContent(iPt+1)
        effAccFD = hEffFD.GetBinContent(iPt+1)
        effAccPromptUnc = hEffPrompt.GetBinError(iPt+1)
        effAccFDUnc = hEffFD.GetBinError(iPt+1)

        corrYieldPrompt = hCorrYieldPrompt.GetBinContent(iPt+1)
        corrYieldFD = hCorrYieldFD.GetBinContent(iPt+1)
        covPromptPrompt = hCovPromptPrompt.GetBinContent(iPt+1)
        covPromptFD = hCovPromptFD.GetBinContent(iPt+1)
        covFDFD = hCovFDFD.GetBinContent(iPt+1)

        fracPromptFD, uncFracPromptFD = GetPromptFDFractionCutSet(effAccPrompt, effAccFD, corrYieldPrompt, corrYieldFD,
                                                              covPromptPrompt, covFDFD, covPromptFD)

        fracPromptFDcorr, uncFracPromptFDcorr = GetPromptFDFractionCutSet(1., 1., corrYieldPrompt, corrYieldFD,
                                                                      covPromptPrompt, covFDFD, covPromptFD)

        hPromptFrac.SetBinContent(iPt+1, fracPromptFD[0])
        hPromptFrac.SetBinError(iPt+1, uncFracPromptFD[0])
        hFDFrac.SetBinContent(iPt+1, fracPromptFD[1])
        hFDFrac.SetBinError(iPt+1, uncFracPromptFD[1])
        hPromptFracCorr.SetBinContent(iPt+1, fracPromptFDcorr[0])
        hPromptFracCorr.SetBinError(iPt+1, uncFracPromptFDcorr[0])
        hFDFracCorr.SetBinContent(iPt+1, fracPromptFDcorr[1])
        hFDFracCorr.SetBinError(iPt+1, uncFracPromptFDcorr[1])

    SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14)

    legFrac = TLegend(0.2, 0.84, 0.4, 0.94)
    legFrac.SetBorderSize(0)
    legFrac.SetFillStyle(0)
    legFrac.SetTextSize(0.045)
    legFrac.AddEntry(hPromptFrac, 'Prompt', 'p')
    legFrac.AddEntry(hFDFrac, 'Non-prompt', 'p')

    legEff = legFrac.Clone('legEff')
    legEff.SetY1(0.2)
    legEff.SetY2(0.4)

    ptMin = hPromptFrac.GetBinLowEdge(1)
    cFrac = TCanvas('cFrac', '', 800, 800)
    cFrac.DrawFrame(ptMin, 0., ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); fraction')
    hPromptFrac.SetLineColor(kRed)
    hPromptFrac.Draw('same')
    hFDFrac.SetLineColor(kBlue)
    hFDFrac.Draw('same')
    legFrac.Draw()
    cFrac.Update()

    cFracCorrFrac = TCanvas('cFracCorrFrac', '', 800, 800)
    cFracCorrFrac.DrawFrame(ptMin, 0., ptMax, 1.2, ';#it{p}_{T} (GeV/#it{c}); corrected fraction')
    hPromptFracCorr.SetLineColor(kRed)
    hPromptFracCorr.Draw('same')
    hFDFracCorr.SetLineColor(kBlue)
    hFDFracCorr.Draw('same')
    legFrac.Draw()
    cFracCorrFrac.Update()

    cEff = TCanvas('cEff', '', 800, 800)
    cEff.DrawFrame(ptMin, 1.e-4, ptMax, 1., ';#it{p}_{T} (GeV/#it{c}); (Acc#times#font[152]{e})')
    cEff.SetLogy()
    hEffPrompt.SetLineColor(kRed)
    hEffPrompt.Draw('same')
    hEffFD.SetLineColor(kBlue)
    hEffFD.Draw('same')
    legEff.Draw()
    cEff.Update()
        
    os.makedirs(outputdir + '/DataDrivenFrac', exist_ok=True)
    outFile = TFile(outputdir + '/DataDrivenFrac/' + f'DataDrivenFrac_{suffix}_{iFile:02}.root', 'recreate')
    hEffPrompt.Write()
    hEffFD.Write()
    hPromptFrac.Write()
    hFDFrac.Write()
    hPromptFracCorr.Write()
    hFDFracCorr.Write()
    cFrac.Write()
    cFracCorrFrac.Write()
    cEff.Write()
    outFile.Close()

def load_cutVarFrac_histos(cutVarFracFiles):
    cutVarFracFile = TFile.Open(cutVarFracFiles[0])
    hCorrYieldPrompt = cutVarFracFile.Get('hCorrYieldPrompt')
    hCorrYieldPrompt.SetDirectory(0)
    hCorrYieldFD = cutVarFracFile.Get('hCorrYieldFD')
    hCorrYieldFD.SetDirectory(0)
    hCovPromptPrompt = cutVarFracFile.Get('hCovPromptPrompt')
    hCovPromptPrompt.SetDirectory(0)
    hCovPromptFD = cutVarFracFile.Get('hCovPromptFD')
    hCovPromptFD.SetDirectory(0)
    hCovFDFD = cutVarFracFile.Get('hCovFDFD')
    hCovFDFD.SetDirectory(0)
    return hCorrYieldPrompt, hCorrYieldFD, hCovPromptPrompt, hCovPromptFD, hCovFDFD

def load_eff_histos(effFiles):
    
    hEffPrompts, hEffFDs, hPromptFracs, hFDFracs, hPromptFracCorrs, hFDFracCorrs = [], [], [], [], [], []
    for effFile in effFiles:
        effFile = TFile.Open(effFile)
        hEffPrompt = effFile.Get('hEffPrompt')
        hEffPrompt.SetDirectory(0)
        hEffFD = effFile.Get('hEffFD')
        hEffFD.SetDirectory(0)
                             
        hEffPrompts.append(hEffPrompt)
        hEffFDs.append(hEffFD)
        
        hPromptFracs.append(hEffPrompt.Clone('hPromptFrac'))
        hFDFracs.append(hEffFD.Clone('hFDFrac'))
        hPromptFracs[-1].SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{f}_{prompt}')
        hFDFracs[-1].SetTitle(';#it{p}_{T} (GeV/#it{c}); #it{f}_{FD}')
        hPromptFracs[-1].SetDirectory(0)
        hFDFracs[-1].SetDirectory(0)
        
        hPromptFracCorrs.append(hEffPrompt.Clone('hPromptFracCorr'))
        hFDFracCorrs.append(hEffFD.Clone('hFDFracCorr'))
        hPromptFracCorrs[-1].SetTitle(';#it{p}_{T} (GeV/#it{c}); corrected #it{f}_{prompt}')
        hFDFracCorrs[-1].SetTitle(';#it{p}_{T} (GeV/#it{c}); corrected #it{f}_{FD}')
        hPromptFracCorrs[-1].SetDirectory(0)
        hFDFracCorrs[-1].SetDirectory(0)
    return hEffPrompts, hEffFDs, hPromptFracs, hFDFracs, hPromptFracCorrs, hFDFracCorrs

def load_eff_files(inputdir):
    if os.path.exists(f'{inputdir}/eff'):
        print(f'Loading {inputdir}/eff')
        effFiles = [f'{inputdir}/eff/{file}'
                        for file in os.listdir(f'{inputdir}/eff') if file.endswith('.root')]
        effFiles.sort()
    else:
        raise ValueError(f'No eff folder found in {inputdir}')
    return effFiles

def load_cutVarFrac_files(inputdir):
    if os.path.exists(f'{inputdir}/CutVarFrac'):
        print(f'Loading {inputdir}/CutVarFrac')
        cutVarFracFiles = [f'{inputdir}/CutVarFrac/{file}'
                        for file in os.listdir(f'{inputdir}/CutVarFrac') if file.endswith('.root')]
        cutVarFracFiles.sort()
    else:
        raise ValueError(f'No CutVarFrac folder found in {inputdir}')
    return cutVarFracFiles

def main_data_driven_frac(inputdir, outputdir, suffix, batch, combined=False, correlatedCutVarPath="", outputdir_combined=""):

    if batch:
        gROOT.SetBatch()

    effFiles = load_eff_files(inputdir)
    hEffPrompts, hEffFDs, hPromptFracs, hFDFracs, hPromptFracCorrs, hFDFracCorrs = load_eff_histos(effFiles)

    cutVarFracFiles = load_cutVarFrac_files(inputdir)
    hCorrYieldPrompt, hCorrYieldFD, hCovPromptPrompt, hCovPromptFD, hCovFDFD = load_cutVarFrac_histos(cutVarFracFiles)

    for iFile in range(len(effFiles)):
        data_driven_frac(
            outputdir,
            suffix,
            iFile,
            hEffPrompts[iFile],
            hEffFDs[iFile],
            hPromptFracs[iFile],
            hFDFracs[iFile],
            hPromptFracCorrs[iFile],
            hFDFracCorrs[iFile],
            hCorrYieldPrompt,
            hCorrYieldFD,
            hCovPromptPrompt,
            hCovPromptFD,
            hCovFDFD
        )
    if combined:
        print('Computing combined DataDriFrac')
        cutVarFracFiles = load_cutVarFrac_files(correlatedCutVarPath)
        hCorrYieldPrompt, hCorrYieldFD, hCovPromptPrompt, hCovPromptFD, hCovFDFD = load_cutVarFrac_histos(cutVarFracFiles)
        
        for iFile in range(len(effFiles)):
            data_driven_frac(
                outputdir_combined,
                suffix,
                iFile,
                hEffPrompts[iFile],
                hEffFDs[iFile],
                hPromptFracs[iFile],
                hFDFracs[iFile],
                hPromptFracCorrs[iFile],
                hFDFracCorrs[iFile],
                hCorrYieldPrompt,
                hCorrYieldFD,
                hCovPromptPrompt,
                hCovPromptFD,
                hCovFDFD
            )

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('--inputdir', '-i', metavar='text',
                        default='.', help='input directory containing eff and frac files')
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--batch", '-b',action='store_true', help="run in batch mode")
    parser.add_argument("--combined", '-comb', default=False, required=False,
                        action='store_true', help="combined method")
    parser.add_argument("--correlatedCutVarPath", "-cc", metavar="text", required=False,
                        default="", help="path to the correlated cut method")
    parser.add_argument("--outputdir_combined", "-oc", metavar="text", required=False,
                        default="", help="output directory for the combined method")
    args = parser.parse_args()
    
    main_data_driven_frac(
        args.inputdir,
        args.outputdir,
        args.suffix,
        args.batch,
        combined=args.combined if args.combined else False,
        correlatedCutVarPath=args.correlatedCutVarPath if args.correlatedCutVarPath else "",
        outputdir_combined=args.outputdir_combined if args.outputdir_combined else ""
        )
