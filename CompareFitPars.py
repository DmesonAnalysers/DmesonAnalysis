import argparse
from ROOT import TCanvas, TFile, TLegend, TLine # pylint: disable=import-error,no-name-in-module
from ROOT import gStyle, kRed, kBlack, kBlue, kOrange, kGreen # pylint: disable=import-error,no-name-in-module
from ROOT import kFullCircle, kOpenCircle, kFullSquare, kFullDiamond, kFullCross, kOpenCross # pylint: disable=import-error,no-name-in-module

def set_style():
    gStyle.SetPadRightMargin(0.035)
    gStyle.SetPadLeftMargin(0.18)
    gStyle.SetPadTopMargin(0.05)
    gStyle.SetTitleSize(0.045, 'xy')
    gStyle.SetLabelSize(0.040, 'xy')
    gStyle.SetPadTickX(1)
    gStyle.SetPadTickY(1)
    gStyle.SetLegendBorderSize(0)
    gStyle.SetOptStat(0)

def create_ratio_hist(num_hist_list, den_hist):
    ratio_list = []
    for num_hist in num_hist_list:
        ratio_hist = den_hist.Clone()
        ratio_hist.Reset()
        n_bin_orig = num_hist.GetNbinsX()
        for ipt in range(ratio_hist.GetNbinsX()+1):
            pt_cent = ratio_hist.GetBinCenter(ipt)
            pt_bin = num_hist.FindBin(pt_cent)
            if 0 < pt_bin <= n_bin_orig:
                ratio_hist.SetBinContent(ipt, num_hist.GetBinContent(pt_bin))
                ratio_hist.SetBinError(ipt, num_hist.GetBinError(pt_bin))
            else:
                ratio_hist.SetBinContent(ipt, 0.)
                ratio_hist.SetBinError(ipt, 0.)
        ratio_hist.Divide(den_hist)
        ratio_list.append(ratio_hist)
    return ratio_list

def comp_fit_pars(do_ratio=False): #pylint: disable-msg=too-many-statements
    inputdir = './'
    input_files = ['../ML_DsAnalysis/QM_prel/outputs/3050/raw_yields/RawYieldsDs_3050_pt2_36.root',
                   'outputs/3_24bin_merge/raw_yields/RawYieldsDs_3050_100419cuts.root']
    input_files_MC = ['../ML_DsAnalysis/QM_prel/outputs/3050/raw_yields/RawYieldsDs_3050_MC_pt2_36.root']
    colors = [kBlack, kGreen+3, kRed]
    markers = [kFullCircle, kFullDiamond, kFullSquare]
    legendnames = ['ML', 'Std', 'ML - MC']
    suffix = 'compMCStd'
    min_pt = 2.
    max_pt = 36.

    set_style()
    hMean, hSigma = [], []
    input_files = input_files + input_files_MC

    lineMass = TLine(min_pt, 1.96850, max_pt, 1.96850)
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
        histo_mean.SetLineColor(color)
        histo_mean.SetLineWidth(2)
        histo_mean.SetMarkerColor(color)
        histo_mean.SetMarkerStyle(marker)
        histo_sigma.SetDirectory(0)
        histo_sigma.SetLineColor(color)
        histo_sigma.SetLineWidth(2)
        histo_sigma.SetMarkerColor(color)
        histo_sigma.SetMarkerStyle(marker)
        legMean.AddEntry(histo_mean, legend_name, 'p')
        legSigma.AddEntry(histo_sigma, legend_name, 'p')
        hMean.append(histo_mean)
        hSigma.append(histo_sigma)

    if do_ratio:
        mean_num_list = list(hMean)
        mean_num_list.pop()
        mean_den_hist = hMean[-1]
        mean_ratio_list = create_ratio_hist(mean_num_list, mean_den_hist)
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
        for hist, color, marker, legend_name in zip(mean_ratio_list, colors[:-1], 
                                                    markers[:-1], legendnames[:-1]):
            hist.SetLineColor(color)
            hist.SetLineWidth(2)
            hist.SetMarkerColor(color)
            hist.SetMarkerStyle(marker)
            hist.Draw('same')
            legMeanRatio.AddEntry(hist, legend_name, 'p')
        lineRatio.Draw()
        legMeanRatio.Draw()

        sigma_num_list = list(hSigma)
        sigma_num_list.pop()
        sigma_den_hist = hSigma[-1]
        sigma_ratio_list = create_ratio_hist(sigma_num_list, sigma_den_hist)
        cSigmaRatio = TCanvas('cSigmaRatio', '', 800, 800)
        cSigmaRatio.DrawFrame(min_pt, 0.6, max_pt, 1.8, ';#it{p}_{T} (GeV/#it{c}); peak width / peak width MC')
        legSigmaRatio = TLegend(0.4, 0.73, 0.7, 0.93)
        legSigmaRatio.SetFillStyle(0)
        legSigmaRatio.SetBorderSize(0)
        legSigmaRatio.SetTextSize(0.04)
        for hist, color, marker, legend_name in zip(sigma_ratio_list, colors[:-1], 
                                                    markers[:-1], legendnames[:-1]):
            hist.SetLineColor(color)
            hist.SetLineWidth(2)
            hist.SetMarkerColor(color)
            hist.SetMarkerStyle(marker)
            hist.Draw('same')
            legSigmaRatio.AddEntry(hist, legend_name, 'p')
        lineRatio.Draw()
        legSigmaRatio.Draw()
        cMeanRatio.SaveAs(f'{inputdir}/MeanRatio_{suffix}.pdf')
        cSigmaRatio.SaveAs(f'{inputdir}/SigmaRatio_{suffix}.pdf')

    cMean = TCanvas('cMean', '', 800, 800)
    cMean.DrawFrame(min_pt, 1.9641, max_pt, 1.9759, ';#it{p}_{T} (GeV/#it{c}); peak mean (GeV/#it{c}^{2})')
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
    args = parser.parse_args()
    comp_fit_pars(args.ratio)
    input('Press enter to exit')

main()
