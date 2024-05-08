import sys
import argparse
import ROOT
from flow_analysis_utils import get_resolution, get_centrality_bins, getListOfHisots
sys.path.append('../../')
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle
SetGlobalStyle(padleftmargin=0.15, padbottommargin=0.15,
               padrightmargin=0.15, titleoffsety=1.1, maxdigits=3, titlesizex=0.03,
               labelsizey=0.04, setoptstat=0, setopttitle=0, palette=ROOT.kGreyScale)

ROOT.gROOT.SetBatch(False)

# TODO: move this to the StyleFormatter
def SetFrameStyle(hFrame, xtitle, ytitle, ytitleoffset, ytitlesize, ylabelsize,
                  ylabeloffset, xticklength, yticklength, xtitlesize, xlabelsize,
                  xtitleoffset, xlabeloffset, ydivisions, xmoreloglabels, ycentertitle, ymaxdigits):
    hFrame.GetXaxis().SetTitle(xtitle)
    hFrame.GetYaxis().SetTitle(ytitle)
    hFrame.GetYaxis().SetTitleOffset(ytitleoffset)
    hFrame.GetYaxis().SetTitleSize(ytitlesize)
    hFrame.GetYaxis().SetLabelSize(ylabelsize)
    hFrame.GetYaxis().SetLabelOffset(ylabeloffset)
    hFrame.GetXaxis().SetTickLength(xticklength)
    hFrame.GetYaxis().SetTickLength(yticklength)
    hFrame.GetXaxis().SetTitleSize(xtitlesize)
    hFrame.GetXaxis().SetLabelSize(xlabelsize)
    hFrame.GetXaxis().SetTitleOffset(xtitleoffset)
    hFrame.GetXaxis().SetLabelOffset(xlabeloffset)
    hFrame.GetYaxis().SetNdivisions(ydivisions)
    hFrame.GetXaxis().SetMoreLogLabels(xmoreloglabels)
    hFrame.GetYaxis().CenterTitle(ycentertitle)
    hFrame.GetYaxis().SetMaxDigits(ymaxdigits)

def compute_reso(an_res_file, vn_method,
                 centClass, wagon_id, outputdir, suffix):

    _, cent_min_max = get_centrality_bins(centClass)
    histos_triplets, histos_triplets_lables = getListOfHisots(an_res_file, wagon_id, vn_method)

    # prepare output file
    if vn_method == 'sp':
        ytitle = 'Q^{A} Q^{B}'
        outfile_name = f'{outputdir}resoSP{suffix}.root'
    elif vn_method == 'ep' or vn_method == 'deltaphi':
        ytitle = 'cos(2(#Psi^{A}-#Psi^{B}))'
        outfile_name = f'{outputdir}resoEP{suffix}.root'
    else:
        sys.exit('\033[91mFATAL: Invalid vn_method. Only sp, ep, deltaphi implemented. Exit!\033[0m')
    outfile = ROOT.TFile(outfile_name, 'RECREATE')

    # loop over all possible combinations of detectors
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextSize(0.05)
    for i, (histo_triplet, histo_triplet_label) in enumerate(zip(histos_triplets, histos_triplets_lables)):
        histos_mean, histos_mean_deltacent, histo_reso, histo_reso_deltacent = get_resolution(histo_triplet,
                                                                                              histo_triplet_label,
                                                                                              cent_min_max)
        detA_label = histo_triplet_label[0]
        detB_label = histo_triplet_label[1]
        detC_label = histo_triplet_label[2]
        outfile.cd()
        outfile.mkdir(f'{detA_label}_{detB_label}_{detC_label}')
        outfile.cd(f'{detA_label}_{detB_label}_{detC_label}')
        canvas = ROOT.TCanvas(f'canvas_{detA_label}_{detB_label}_{detC_label}',
                              f'canvas_{detA_label}_{detB_label}_{detC_label}',
                              2400, 800)
        canvas.Divide(3, 1)
        leg = ROOT.TLegend(0.2, 0.2, 0.5, 0.3)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.03)
        for i, (hist_det, hist_mean, histo_mean_deltacent) in enumerate(zip(histo_triplet,
                                                                            histos_mean,
                                                                            histos_mean_deltacent)):
            SetObjectStyle(hist_mean, color=ROOT.kRed, markerstyle=ROOT.kFullCircle,
                           markersize=1, fillstyle=0, linewidth=2)
            SetObjectStyle(histo_mean_deltacent, color=ROOT.kBlue, markerstyle=ROOT.kOpenCircle,
                           markersize=1, fillstyle=0, linestyle=2, linewidth=3)
            canvas.cd(i+1)
            canvas.cd(i+1).SetLogz()
            hFrame = canvas.cd(i+1).DrawFrame(0, -2, 100, 2)
            SetFrameStyle(hFrame,
                          xtitle='Cent. FT0c (%)',
                          ytitle=ytitle,
                          ytitleoffset=1.15,
                          ytitlesize=0.05,
                          ylabelsize=0.04,
                          ylabeloffset=0.01,
                          xticklength=0.04,
                          yticklength=0.03,
                          xtitlesize=0.05,
                          xlabelsize=0.04,
                          xtitleoffset=1.1,
                          xlabeloffset=0.020,
                          ydivisions=406,
                          xmoreloglabels=True,
                          ycentertitle=True,
                          ymaxdigits=5)
            hist_det.Draw('same colz')
            histo_mean_deltacent.Draw('same pl')
            hist_mean.Draw('same pl')
            if i == 0:
                leg.AddEntry(hist_mean, 'Average 1% centrality', 'lp')
                leg.AddEntry(histo_mean_deltacent,
                             f'Average {cent_min_max[1]-cent_min_max[0]}% centrality', 'lp')
                leg.Draw()
                latex.DrawLatex(0.2, 0.85, f'A: {detA_label}, B: {detB_label}')
            elif i == 1:
                latex.DrawLatex(0.2, 0.85, f'A: {detA_label}, B: {detC_label}')
            else:
                latex.DrawLatex(0.2, 0.85, f'A: {detB_label}, B: {detC_label}')
            histo_mean_deltacent.Write()
            hist_mean.Write()
            hist_det.Write()
        canvas.Update()
        canvas.Write()
        histo_reso.SetDirectory(outfile)
        histo_reso_deltacent.SetDirectory(outfile)
        histo_reso.Write()
        histo_reso_deltacent.Write()
        outfile.cd('..')

    input('Resolutions computed. Press any key to continue')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("an_res_file", metavar="text",
                        default="an_res.root", help="input ROOT file with anres")
    parser.add_argument('--centClass', '-c', metavar='text', default='k0100')
    parser.add_argument('--vn_method', '-vn', metavar='text', default='sp')
    parser.add_argument("--wagon_id", "-w", metavar="text",
                        default="", help="wagon ID", required=False)
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    compute_reso(
        an_res_file=args.an_res_file,
        vn_method=args.vn_method,
        centClass=args.centClass,
        wagon_id=args.wagon_id,
        outputdir=args.outputdir,
        suffix=args.suffix
        )
