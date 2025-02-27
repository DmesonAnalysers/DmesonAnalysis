import sys
import os
import numpy as np
import argparse
import yaml
import array
from ROOT import TFile, TCanvas, kFullSquare, kBlack, kOrange, kAzure, kGray, kRed, TLegend, kCyan, kSpring, kGreen, kBlue, kMagenta, kFullCircle, TGraphErrors, TH1F, TMultiGraph
sys.path.append('../../../')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

SetGlobalStyle(titleoffsety=1.4, maxdigits=3, topmargin=0.1, bottommargin=0.4, leftmargin=1, rightmargin=0.15,
               labelsizey=0.04, setoptstat=0, setopttitle=0, setdecimals=True,titleoffsetx=0.94)
cols = [kBlack, kOrange+1, kAzure+4, kRed+1, kCyan+2, kGreen+2, kBlue-4, kMagenta+2, kSpring+2]
labels = ['default', 'random', 'even', 'odd', 'step1', 'step2', 'minus3low', 'minus3high', 'minus3']

def get_hdeltav(href, hsyst):
    hdelta = hsyst.Clone()
    hdelta.Reset()
    hdelta.SetDirectory(0)
    for ibin in range(1, href.GetNbinsX()+1):
        hdelta.SetBinContent(ibin, hsyst.GetBinContent(ibin)-href.GetBinContent(ibin))
        hdelta.SetBinError(ibin, 1.e-9)

    return hdelta

def get_rms_shift(histsdeltav2):
    gsyst = TGraphErrors(-1)
    SetObjectStyle(gsyst, color=kAzure+2, fillcolor=kAzure+2, linewidth=0, fillalpha=0.2, linestyle=8, fillstyle=3154)
    for ibin in range(1, histsdeltav2[0].GetNbinsX()+1):
        hsyst = TH1F('', '', 20000, -0.1, 0.1)
        for hist in histsdeltav2[1:]: # skipping default
            hsyst.Fill(hist.GetBinContent(ibin))
        #hsyst.Draw()
        gsyst.SetPoint(ibin, histsdeltav2[0].GetBinCenter(ibin), 0)
        gsyst.SetPointError(ibin, histsdeltav2[0].GetBinWidth(ibin)/2, np.sqrt(hsyst.GetRMS()**2 + hsyst.GetMean()**2))
        #print(f'(pT cent: {histsdeltav2[0].GetBinCenter(ibin)}) Syst = {np.sqrt(hsyst.GetRMS()**2 + hsyst.GetMean()**2)}, rms = {hsyst.GetRMS()}, mu = {hsyst.GetMean()}')

    return gsyst


def compute_syst(infiles, outputdir, suffix):
    #______________________________________________________________________________________
    # Collect all .root files
    hv2_prompt = []
    hv2_fd = []
    hdeltav2_prompt = []
    hdeltav2_fd = []
    gv2vsfrac = {}
    for ifile, file in enumerate(infiles):
        print(f'INFO: Found {file}')
        tfile = TFile().Open(file)
        # check if the file has the right keys
        if not tfile.GetListOfKeys().Contains('hV2VsPtFD') or not tfile.GetListOfKeys().Contains('hV2VsPtPrompt'):
            print(f'ERROR: {file} does not contain hV2VsPtFD or hV2VsPtPrompt')
            continue
        hv2_prompt.append(tfile.Get('hV2VsPtPrompt'))
        hv2_fd.append(tfile.Get('hV2VsPtFD'))

        gv2vsfrac[labels[ifile]] = []
        for ibin in range(1, hv2_prompt[0].GetNbinsX()+1):
            ptmin = hv2_prompt[0].GetXaxis().GetBinLowEdge(ibin)
            ptmax = hv2_prompt[0].GetXaxis().GetBinLowEdge(ibin) + hv2_prompt[0].GetXaxis().GetBinWidth(ibin)
            print(f'pt_{ptmin*10:.0f}_{ptmax*10:.0f}/cFrac_{ptmin:.0f}_{ptmax:.0f}')
            gv2vsfrac[labels[ifile]].append(tfile.Get(f'pt_{ptmin*10:.0f}_{ptmax*10:.0f}/gV2VsFrac'))
            SetObjectStyle(gv2vsfrac[labels[ifile]][-1], color=cols[ifile], fillcolor=cols[ifile], markerstyle=kFullCircle, markersize=2, linewidth=2, fillalpha=0.2)

        hv2_prompt[-1].SetDirectory(0)
        hv2_fd[-1].SetDirectory(0)
        SetObjectStyle(hv2_prompt[-1], color=cols[ifile], fillcolor=cols[ifile], markerstyle=kFullCircle, markersize=2, linewidth=2, fillalpha=0.2)
        SetObjectStyle(hv2_fd[-1], color=cols[ifile], fillcolor=cols[ifile], markerstyle=kFullSquare, markersize=2, linewidth=2, fillalpha=0.2)

        hdeltav2_prompt.append(get_hdeltav(hv2_prompt[0], hv2_prompt[ifile]))
        hdeltav2_fd.append(get_hdeltav(hv2_fd[0], hv2_fd[ifile]))
    gsyst_prompt = get_rms_shift(hdeltav2_prompt)
    gsyst_fd = get_rms_shift(hdeltav2_fd)

    canv_fFD = TCanvas('canv', 'canv', 1600, 1600)
    canv_fFD.Divide(4, 4)
    for igs, gfracs in enumerate(gv2vsfrac):
        print(gfracs)
        for ig, gfrac in enumerate(gv2vsfrac[gfracs]):
            canv_fFD.cd(ig+1)
            if igs == 0:
                gfrac.Draw('pea')
            else:
                gfrac.Draw('same pe')
    canv_fFD.Update()

    #______________________________________________________________________________________
    # Plot the systematics
    canv = TCanvas('canv', 'canv', 1600, 1600)
    canv.Divide(2, 2)
    
    legv2 = TLegend(0.6, 0.6, 0.9, 0.9)
    legv2.SetBorderSize(0)
    legv2.SetFillStyle(0)
    legv2.SetTextSize(0.03)
    legv2.SetTextFont(42)
    
    legsyst = TLegend(0.6, 0.79, 0.9, 0.89)
    legsyst.SetBorderSize(0)
    legsyst.SetFillStyle(0)
    legsyst.SetTextSize(0.03)
    legsyst.SetTextFont(42)
    
    canv.cd(1).SetLeftMargin(0.16)
    canv.cd(1)
    for ihist, hist in enumerate(hv2_prompt):
        hist.SetStats(0)
        hist.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
        hist.GetYaxis().SetTitle('prompt #it{v}_{2}')
        hist.GetXaxis().SetRangeUser(0, 25)
        hist.GetYaxis().SetRangeUser(-0.1, 0.35)
        legv2.AddEntry(hist, labels[ihist], 'lp')
        hist.Draw('same')
    legv2.Draw()
    canv.cd(2).SetLeftMargin(0.16)
    for _, hist in enumerate(hv2_fd):
        hist.SetStats(0)
        hist.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
        hist.GetYaxis().SetTitle('non-prompt #it{v}_{2}')
        hist.GetXaxis().SetRangeUser(0, 25)
        hist.GetYaxis().SetRangeUser(-0.2, 0.35)
        hist.Draw('same')
    canv.cd(3).SetGridy()
    canv.cd(3).SetLeftMargin(0.16)
    for _, hist in enumerate(hdeltav2_prompt):
        hist.SetStats(0)
        hist.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
        hist.GetYaxis().SetTitle('prompt #it{v}_{2}^{syst.} #minus #it{v}_{2}^{ref.}')
        hist.GetXaxis().SetRangeUser(0, 25)
        hist.GetYaxis().SetRangeUser(-0.05, 0.05)
        hist.Draw('samepe')
    gsyst_prompt.Draw('same5')
    legsyst.AddEntry(gsyst_prompt, f'#sqrt{{shift^{{2}} + rms^{{2}}}}', 'f')
    legsyst.Draw()
    canv.cd(4).SetGridy()
    canv.cd(4).SetLeftMargin(0.16)
    for _, hist in enumerate(hdeltav2_fd):
        hist.SetStats(0)
        hist.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
        hist.GetYaxis().SetTitle('non-prompt #it{v}_{2}^{syst.} #minus #it{v}_{2}^{ref.}')
        hist.GetXaxis().SetRangeUser(0, 25)
        hist.GetYaxis().SetRangeUser(-0.05, 0.05)
        hist.Draw('samepe')
    gsyst_fd.Draw('same5')

    canv.Update()
    input('Press enter to exit.')

    #______________________________________________________________________________________
    # Save output
    outputfile = os.path.join(outputdir, f'syst_fFD_{suffix}.root')
    outfile = TFile(outputfile, 'RECREATE')
    canv.Write()
    outfile.Close()
    print(f'INFO: Output saved in {outputfile}')

    canv.SaveAs(f'{outputdir}/SystfFD_{suffix}.pdf')
    canv.SaveAs(f'{outputdir}/SystfFD_{suffix}.png')

def check_frac_variation_range(infiles, outputdir):
    
    # get pt bins from first file
    tfile_forptbins = TFile.Open(infiles[0])
    ptmins = []
    ptmaxs = []
    hRef = tfile_forptbins.Get('hV2VsPtPrompt')
    for ibin in range(1, hRef.GetNbinsX()+1):
        ptmins.append(int(hRef.GetXaxis().GetBinLowEdge(ibin) * 10))  # Convert to integer
        ptmaxs.append(int((hRef.GetXaxis().GetBinLowEdge(ibin) + hRef.GetXaxis().GetBinWidth(ibin)) * 10))

    cAllPtFracs = []
    gSinglePtAllFracs = [] 
    legends = []
    
    for iPt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs)):
        cAllPtFracs.append(TCanvas(f"cFrac_pt_{ptmin}_{ptmax};cutset;fFD", "", 800, 600))
        gSinglePtAllFracs.append(TMultiGraph(f"gFrac_pt_{ptmin}_{ptmax};cutset;fFD", ""))
        legends.append(TLegend(0.7, 0.7, 0.9, 0.9))  # Legend position

    for ifile, file in enumerate(infiles):
        tfile = TFile.Open(file)
        if not tfile or tfile.IsZombie():
            print(f"Error: Cannot open file {file}")
            continue
        
        for iPt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs)):
            gFracPt = tfile.Get(f'pt_{ptmin}_{ptmax}/gV2VsFrac')
            x_values = gFracPt.GetX()
            frac_values = [x_values[i] for i in range(gFracPt.GetN())]

            graph = TGraphErrors(gFracPt.GetN(), array.array('d', range(len(frac_values))), array.array('d', frac_values))
            graph.SetMarkerColor(ifile + 1)  # Different color for each file
            graph.SetLineColor(ifile + 1)
            graph.SetLineWidth(2)
            graph.SetMarkerStyle(20)

            gSinglePtAllFracs[iPt].Add(graph, "LP")
            legends[iPt].AddEntry(graph, os.path.basename(file).replace("V2VsFrac_", "").replace(".root", ""), "LP")

    # Save to ROOT file
    FracScanFile = TFile(f'{outputdir}/FracVariations.root', 'recreate')
    for iPt, (ptmin, ptmax) in enumerate(zip(ptmins, ptmaxs)):
        FracScanFile.mkdir(f"pt_{ptmin}_{ptmax}")
        FracScanFile.cd(f"pt_{ptmin}_{ptmax}")
        cAllPtFracs[iPt].cd()
        gSinglePtAllFracs[iPt].Draw("A")
        legends[iPt].Draw()
        gSinglePtAllFracs[iPt].Write()
        cAllPtFracs[iPt].SetLogy()
        cAllPtFracs[iPt].Write()

    FracScanFile.Close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('infiles', metavar='text', nargs='*', default='path/to/syst/files')
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    args = parser.parse_args()

    compute_syst(args.infiles,
                 args.outputdir,
                 args.suffix)

    check_frac_variation_range(args.infiles, args.outputdir)
