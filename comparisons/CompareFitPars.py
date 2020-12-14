import sys
import argparse
from ROOT import TCanvas, TFile, TLegend, TLine, TDatabasePDG # pylint: disable=import-error,no-name-in-module
from ROOT import gStyle, kRed, kAzure, kBlack, kBlue, kOrange, kGreen # pylint: disable=import-error,no-name-in-module,unused-import
from ROOT import kFullCircle, kOpenCircle, kFullSquare, kFullDiamond, kFullCross, kOpenCross, kOpenSquare # pylint: disable=import-error,no-name-in-module,unused-import
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error
from utils.AnalysisUtils import ComputeRatioDiffBins #pylint: disable=wrong-import-position,import-error

def comp_fit_pars(do_ratio=False, meson='Ds'): #pylint: disable-msg=too-many-statements,too-many-locals
    inputdir = '../../AnalysisNonPromptDpp2017/Dplus/outputs/rawyields'
    input_files = ['RawYieldsDplus_pp5TeV_prompt_central.root', 'RawYieldsDplus_pp5TeV_FD_central_freesigma.root']
    input_files_MC = ['RawYieldsDplusMC_pp5TeV_prompt_central.root', 'RawYieldsDplusMC_pp5TeV_FD_central.root']
    colors = [kOrange+7, kAzure+2, kRed+1, kAzure+4]
    markers = [kOpenCircle, kOpenSquare, kFullCircle, kFullSquare]
    legendnames = ['MC - prompt enhanced', 'MC - FD enhanced', 'data - prompt enhanced', 'data - FD enhanced']
    suffix = 'CompMCData'
    min_pt = 2.
    max_pt = 16.

    SetGlobalStyle(padleftmargin=0.18, padtopmargin=0.05, padbottommargin=0.14,
                   titleoffsety=1.6, titlesize=0.045, labelsize=0.04)
    hMean, hSigma = [], []
    input_files = input_files_MC + input_files

    if meson == 'Ds':
        massD = TDatabasePDG.Instance().GetParticle(431).Mass()
    elif meson == 'Dplus':
        massD = TDatabasePDG.Instance().GetParticle(411).Mass()

    lineMass = TLine(min_pt, massD, max_pt, massD)
    lineMass.SetLineWidth(2)
    lineMass.SetLineColor(kBlack)
    lineMass.SetLineStyle(9)

    legSigma = TLegend(0.2, 0.78, 0.8, 0.93)
    legSigma.SetFillStyle(0)
    legSigma.SetBorderSize(0)
    legSigma.SetTextSize(0.04)

    legMean = TLegend(0.4, 0.73, 0.7, 0.93)
    legMean.SetFillStyle(0)
    legMean.SetBorderSize(0)
    legMean.SetTextSize(0.04)
    legMean.AddEntry(lineMass, "PDG", 'l')

    for file_path, color, marker, legend_name in zip(input_files, colors, markers, legendnames):
        input_file = TFile(f'{inputdir}/{file_path}')
        histo_mean = input_file.Get('hRawYieldsMean')
        histo_sigma = input_file.Get('hRawYieldsSigma')
        histo_mean.SetDirectory(0)
        histo_sigma.SetDirectory(0)
        SetObjectStyle(histo_mean, linecolor=color, markercolor=color, markerstyle=marker)
        SetObjectStyle(histo_sigma, linecolor=color, markercolor=color, markerstyle=marker)
        legMean.AddEntry(histo_mean, legend_name, 'p')
        legSigma.AddEntry(histo_sigma, legend_name, 'p')
        hMean.append(histo_mean)
        hSigma.append(histo_sigma)

    if do_ratio:
        mean_num_list = list(hMean)
        mean_num_list.pop(0)
        mean_den_hist = hMean[0]
        mean_ratio_list = []
        for histo in mean_num_list:
            mean_ratio_list.append(ComputeRatioDiffBins(histo, mean_den_hist))
        cMeanRatio = TCanvas('cMeanRatio', '', 800, 800)
        cMeanRatio.DrawFrame(min_pt, 0.99, max_pt, 1.01, ';#it{p}_{T} (GeV/#it{c}); peak mean / peak mean MC')
        lineRatio = TLine(min_pt, 1., max_pt, 1.)
        lineRatio.SetLineWidth(2)
        lineRatio.SetLineColor(kBlack)
        lineRatio.SetLineStyle(9)
        legMeanRatio = TLegend(0.4, 0.73, 0.7, 0.93)
        legMeanRatio.SetFillStyle(0)
        legMeanRatio.SetBorderSize(0)
        legMeanRatio.SetTextSize(0.04)
        for hist, color, marker, legend_name in zip(mean_ratio_list, colors[1:],
                                                    markers[1:], legendnames[1:]):
            SetObjectStyle(hist, linecolor=color, markercolor=color, markerstyle=marker)
            hist.Draw('same')
            legMeanRatio.AddEntry(hist, legend_name, 'p')
        lineRatio.Draw()
        legMeanRatio.Draw()

        sigma_num_list = list(hSigma)
        sigma_num_list.pop(0)
        sigma_den_hist = hSigma[0]
        sigma_ratio_list = []
        for histo in sigma_num_list:
            sigma_ratio_list.append(ComputeRatioDiffBins(histo, sigma_den_hist))
        cSigmaRatio = TCanvas('cSigmaRatio', '', 800, 800)
        hFrameSigma = cSigmaRatio.DrawFrame(min_pt, 0.6, max_pt, 1.8,
                                            ';#it{p}_{T} (GeV/#it{c}); peak width / peak width MC')
        hFrameSigma.GetYaxis().SetDecimals()
        legSigmaRatio = TLegend(0.4, 0.73, 0.7, 0.93)
        legSigmaRatio.SetFillStyle(0)
        legSigmaRatio.SetBorderSize(0)
        legSigmaRatio.SetTextSize(0.04)
        for hist, color, marker, legend_name in zip(sigma_ratio_list, colors[1:],
                                                    markers[1:], legendnames[1:]):
            SetObjectStyle(hist, linecolor=color, markercolor=color, markerstyle=marker)
            hist.Draw('same')
            legSigmaRatio.AddEntry(hist, legend_name, 'p')
        lineRatio.Draw()
        legSigmaRatio.Draw()
        cMeanRatio.SaveAs(f'{inputdir}/MeanRatio_{suffix}.pdf')
        cSigmaRatio.SaveAs(f'{inputdir}/SigmaRatio_{suffix}.pdf')

    cMean = TCanvas('cMean', '', 800, 800)
    hFrameMean = cMean.DrawFrame(min_pt, massD*0.995, max_pt, massD*1.01,
                                 ';#it{p}_{T} (GeV/#it{c}); peak mean (GeV/#it{c}^{2})')
    hFrameMean.GetYaxis().SetDecimals()
    lineMass.Draw("same")
    for histo_mean in hMean:
        histo_mean.Draw('same')
    legMean.Draw()

    cSigma = TCanvas('cSigma', '', 800, 800)
    cSigma.DrawFrame(min_pt, 0., max_pt, 0.025, ';#it{p}_{T} (GeV/#it{c}); peak width (GeV/#it{c}^{2})')
    for histo_sigma in hSigma:
        histo_sigma.Draw('same')
    legSigma.Draw()

    cMean.SaveAs(f'{inputdir}/Mean_{suffix}.pdf')
    cSigma.SaveAs(f'{inputdir}/Sigma_{suffix}.pdf')

def main():
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument("--ratio", help="make ratio to MC", action="store_true")
    parser.add_argument("--Dspecie", help="Dplus or Ds mesons", metavar="text", default="Ds")
    args = parser.parse_args()
    comp_fit_pars(args.ratio, args.Dspecie)
    input('Press enter to exit')

main()
