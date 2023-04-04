"""
Utility script to "cook" the FONLL prediction for the non-prompt Ds
"""
import sys
import pandas as pd
import numpy as np
from ROOT import TFile, TCanvas, TLegend, TH1D, gPad #pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kAzure, kFullCircle, kFullSquare #pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

def main(): #pylint: disable=too-many-statements
    """
    main function of script
    """
    SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14, titleoffsety=1.6, maxdigits=2, optstat=0)

    BR_dsPhitoKKpi = 0.0224 # (PDG 2020)
    BR_BtoDs = np.array([0.09, 0.117, 0.93, 0.011]) # B+, B0, Bs, Lb (PDG 2020)
    FF_btoB = np.array([0.344, 0.344, 0.115, 0.198]) # B+, B0, Bs, Lb (PDG 2020 from ppbar collisions)
    Ds_from_B = np.inner(BR_BtoDs, FF_btoB)
    cooked_factor = Ds_from_B * BR_dsPhitoKKpi

    n_bins = 1001
    bin_width = 0.05 # GeV
    bin_lims = [i * bin_width for i in range(0, n_bins + 1)]
    hFONLLCentral = TH1D('hDsPhipitoKkpifromBpred_central_corr',
                         ';#it{p}_{T} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (pb GeV^{-1} #it{c})',
                         n_bins, np.asarray(bin_lims, 'd'))
    hFONLLMin = TH1D('hDsPhipitoKkpifromBpred_min_corr',
                     ';#it{p}_{T} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (pb GeV^{-1} #it{c})',
                     n_bins, np.asarray(bin_lims, 'd'))
    hFONLLMax = TH1D('hDsPhipitoKkpifromBpred_max_corr',
                     ';#it{p}_{T} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (pb GeV^{-1} #it{c})',
                     n_bins, np.asarray(bin_lims, 'd'))

    fonll_df = pd.read_csv("fonll/FONLL_DfromB_pp5_y05_toCook.txt", sep=" ", header=14).astype('float64')

    for pt, cen, low, up in zip(fonll_df['pt'].to_numpy(), fonll_df['central'].to_numpy(),
                                fonll_df['min'].to_numpy(), fonll_df['max'].to_numpy()):
        hFONLLCentral.Fill(pt, cen * cooked_factor)
        hFONLLMin.Fill(pt, low * cooked_factor)
        hFONLLMax.Fill(pt, up * cooked_factor)

    for iBin in range(hFONLLCentral.GetNbinsX()):
        hFONLLCentral.SetBinError(iBin + 1, 1.e-3)
        hFONLLMin.SetBinError(iBin + 1, 1.e-3)
        hFONLLMax.SetBinError(iBin + 1, 1.e-3)

    stdFONLLfile = TFile.Open("fonll/feeddown/D0DplusDstarPredictions_502TeV_y05_noYShift_all_191017_BDShapeCorrected.root")
    hStdFONLLCentral = stdFONLLfile.Get('hDsPhipitoKkpifromBpred_central_corr')
    hStdFONLLMin = stdFONLLfile.Get('hDsPhipitoKkpifromBpred_min_corr')
    hStdFONLLMax = stdFONLLfile.Get('hDsPhipitoKkpifromBpred_max_corr')
    hFONLLPromptCentral = stdFONLLfile.Get('hDsPhipitoKkpipred_central')
    hFONLLPromptMin = stdFONLLfile.Get('hDsPhipitoKkpipred_min')
    hFONLLPromptMax = stdFONLLfile.Get('hDsPhipitoKkpipred_max')
    hStdFONLLCentral.SetDirectory(0)
    hStdFONLLMin.SetDirectory(0)
    hStdFONLLMax.SetDirectory(0)
    hFONLLPromptCentral.SetDirectory(0)
    hFONLLPromptMin.SetDirectory(0)
    hFONLLPromptMax.SetDirectory(0)
    hStdFONLLCentral.SetStats(0)
    hStdFONLLMin.SetStats(0)
    hStdFONLLMax.SetStats(0)
    hFONLLPromptCentral.SetStats(0)
    hFONLLPromptMin.SetStats(0)
    hFONLLPromptMax.SetStats(0)
    stdFONLLfile.Close()

    outFile = TFile('fonll/feeddown/NonPromptDsPredictions_502TeV_y05_cooked.root', 'recreate')
    hFONLLCentral.Write()
    hFONLLMin.Write()
    hFONLLMax.Write()
    hFONLLPromptCentral.Write()
    hFONLLPromptMin.Write()
    hFONLLPromptMax.Write()
    outFile.Close()

    hFONLL = [hFONLLCentral, hFONLLMin, hFONLLMax]
    hStdFONLL = [hStdFONLLCentral, hStdFONLLMin, hStdFONLLMax]
    labels = ['Central', 'Min', 'Max']
    canvas = []

    for label in labels:
        canvas.append(TCanvas(f'cComp{label}', '', 2000, 1000))

    for (histo, histo_std, label, canv) in zip(hFONLL, hStdFONLL, labels, canvas):
        pt_lims = np.array(histo.GetXaxis().GetXbins(), 'd')

        pt_min = list(pt_lims)[0]
        pt_max = list(pt_lims)[-1]
        histo_std = histo_std.Rebin(histo.GetNbinsX(), f'hStdFONLL{label}', pt_lims)
        histo.Rebin(5)
        histo_std.Rebin(5)
        histo.Scale(1.e-6)
        histo_std.Scale(1.e-6)
        histo_ratio = histo.Clone(f'hRatioCookedOverStd{label}')
        histo_ratio.Divide(histo, histo_std)
        histo_ratio.GetYaxis().SetTitle('Cooked FONLL / "Std" FONLL')

        print(f'\n"STD" FONLL {label} integral (ub): {histo_std.Integral("width"):.4e}')
        print(f'Cooked FONLL {label} integral (ub): {histo.Integral("width"):.4e}\n')

        canv.Divide(2, 1)
        canv.cd(1).DrawFrame(pt_min, 1e-7, pt_max, 2.,
                             ';#it{p}_{T} (GeV/#it{c});#frac{d#sigma}{d#it{p}_{T}} (#mub GeV^{-1} #it{c})')
        leg = TLegend(0.4, 0.7, 0.8, 0.9)
        leg.SetTextSize(0.045)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        SetObjectStyle(histo, color=kAzure+4, markerstyle=kFullSquare)
        SetObjectStyle(histo_std, color=kRed, markerstyle=kFullCircle)
        histo.Draw('same')
        histo_std.Draw('same')
        leg.AddEntry(histo, f'Cooked FONLL {label}', 'p')
        leg.AddEntry(histo_std, f'"Std" FONLL {label}', 'p')
        gPad.SetLogy()
        leg.Draw()

        hframe_ratio = canv.cd(2).DrawFrame(pt_min, 0., pt_max, 3.,
                                            ';#it{p}_{T} (GeV/#it{c}); Cooked FONLL / "Std" FONLL')
        hframe_ratio.GetYaxis().SetDecimals()
        SetObjectStyle(histo_ratio, color=kAzure+4, markerstyle=kFullSquare)
        histo_ratio.Draw('same')
        canv.Update()
        canv.SaveAs(f'CompFONLL_CookedVsStd_{label}.pdf')

    input('Press enter to exit')

main()
