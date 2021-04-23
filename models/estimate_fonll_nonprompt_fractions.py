"""
Utility script to estimate the "natural" non-prompt fraction from FONLL + PYTHIA8
"""

import sys
import numpy as np
from ROOT import TFile, TCanvas, TLegend #pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kOrange, kGreen #pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

def main(): #pylint: disable=too-many-statements
    SetGlobalStyle(padleftmargin=0.18, padbottommargin=0.14, titleoffsety=1.6, maxdigits=2, optstat=0)

    inFileFONLL = TFile.Open('fonll/feeddown/DmesonLcPredictions_502TeV_y05_FFee_BRpythia8_SepContr_PDG2020.root')
    hFONLLPromptD0 = inFileFONLL.Get('hD0Kpipred_max')
    hFONLLPromptDplus = inFileFONLL.Get('hDpluskpipipred_max')
    hFONLLPromptDs = inFileFONLL.Get('hDsPhipitoKkpipred_max')
    hFONLLNonPromptD0 = inFileFONLL.Get('hD0KpifromBpred_central_corr')
    hFONLLNonPromptDplus = inFileFONLL.Get('hDpluskpipifromBpred_central_corr')
    hFONLLNonPromptDs = inFileFONLL.Get('hDsPhipitoKkpifromBpred_central_corr')

    binLims = np.array([0., 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 16., 24.], 'd')
    nBins = 12

    hFONLLPromptD0 = hFONLLPromptD0.Rebin(nBins, 'hFONLLPromptD0', binLims)
    hFONLLPromptDplus = hFONLLPromptDplus.Rebin(nBins, 'hFONLLPromptDplus', binLims)
    hFONLLPromptDs = hFONLLPromptDs.Rebin(nBins, 'hFONLLPromptDs', binLims)
    hFONLLNonPromptD0 = hFONLLNonPromptD0.Rebin(nBins, 'hFONLLNonPromptD0', binLims)
    hFONLLNonPromptDplus = hFONLLNonPromptDplus.Rebin(nBins, 'hFONLLNonPromptDplus', binLims)
    hFONLLNonPromptDs = hFONLLNonPromptDs.Rebin(nBins, 'hFONLLNonPromptDs', binLims)

    hFONLLNPfracD0 = hFONLLPromptD0.Clone('hFONLLNPfracD0')
    hFONLLNPfracD0.Divide(hFONLLNonPromptD0, hFONLLNPfracD0)
    hFONLLNPfracDplus = hFONLLPromptDplus.Clone('hFONLLNPfracDplus')
    hFONLLNPfracDplus.Divide(hFONLLNonPromptDplus, hFONLLNPfracDplus)
    hFONLLNPfracDs = hFONLLPromptDs.Clone('hFONLLNPfracDs')
    hFONLLNPfracDs.Divide(hFONLLNonPromptDs, hFONLLNPfracDs)

    SetObjectStyle(hFONLLNPfracD0, color=kRed, linewitdh=2, markerstyle=0, fillstyle=0)
    SetObjectStyle(hFONLLNPfracDplus, color=kGreen, linewitdh=2, markerstyle=0, fillstyle=0)
    SetObjectStyle(hFONLLNPfracDs, color=kOrange, linewitdh=2, markerstyle=0, fillstyle=0)

    mesonList = ['D0', 'D+', 'Ds']
    histoList = [hFONLLNPfracD0, hFONLLNPfracDplus, hFONLLNPfracDs]

    print(*binLims, sep=' ')
    for meson, histo in zip(mesonList, histoList):
        valueList = []
        for iBin in range(1, nBins + 1):
            valueList.append(f'{histo.GetBinContent(iBin):.3f}')
        print(meson)
        print(*valueList, sep=' ')

    legMeson = TLegend(0.65, 0.6, 0.85, 0.8)
    legMeson.SetTextSize(0.04)
    legMeson.SetFillStyle(0)
    legMeson.AddEntry(hFONLLNPfracD0, 'D^{0}', 'l')
    legMeson.AddEntry(hFONLLNPfracDplus, 'D^{+}', 'l')
    legMeson.AddEntry(hFONLLNPfracDs, 'D_{s}^{+}', 'l')

    cFrac = TCanvas('cFrac', '', 500, 500)
    hFrame = cFrac.DrawFrame(0., 0.01, 24., 1.05, ';#it{p}_{T} (GeV/#it{c}); #it{f}_{non-prompt}')
    hFrame.GetYaxis().SetDecimals()
    hFrame.Draw('axis same')
    hFONLLNPfracD0.Draw('same')
    hFONLLNPfracDplus.Draw('same')
    hFONLLNPfracDs.Draw('same')
    legMeson.Draw('same')
    cFrac.Update()

    TFile('./NaturalNonPromptFrac_FONLL_FFee_BRpythia8.root', 'recreate')
    hFONLLNPfracD0.Write()
    hFONLLNPfracDplus.Write()
    hFONLLNPfracDs.Write()
    cFrac.Write()

    input('Press enter to exit')

main()
