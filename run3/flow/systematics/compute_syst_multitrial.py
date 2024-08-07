import sys
import os
import numpy as np
import argparse
from ROOT import TFile, TCanvas, TH1F, TGraphAsymmErrors, TLegend, kOrange, kAzure, kGray
sys.path.append('../../../')
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle

SetGlobalStyle(titleoffsety=1.1, maxdigits=3, topmargin=0.1, bottommargin=0.4, leftmargin=0.3, rightmargin=0.15,
               labelsizey=0.04, setoptstat=0, setopttitle=0, setdecimals=True,titleoffsetx=0.71)

def compute_syst_multitrial(rypathsyst, ry_default, outputdir):
    #______________________________________________________________________________________
    # Collect all .root files from rypathsyst
    root_files = []
    for file in os.listdir(rypathsyst):
        if file.endswith('.root'):
            root_files.append(file)
    if not root_files:
        sys.exit(f'ERROR: No .root files found in {rypathsyst}')

    #______________________________________________________________________________________
    # Get vn vs mass from default .root file
    tfile = TFile().Open(ry_default)
    gvn_vs_mass_default = tfile.Get('gvnSimFit')

    #______________________________________________________________________________________
    # Get vn vs mass from all .root files
    gvn_vs_mass = []
    hchi2 = []
    hsignificance = []
    for file in root_files:
        tfile = TFile().Open(os.path.join(rypathsyst, file))
        # check if the file has the right keys
        if not tfile.GetListOfKeys().Contains('gvnSimFit'):
            print(f'ERROR: {file} does not contain gvnSimFit')
            continue
        gvn_vs_mass.append(tfile.Get('gvnSimFit'))
        hchi2.append(tfile.Get('hRedChi2SimFit'))
        hchi2[-1].SetDirectory(0)
        hsignificance.append(tfile.Get('hRawYieldsSignificanceSimFit'))
        hsignificance[-1].SetDirectory(0)
        tfile.Close()

    #______________________________________________________________________________________
    # Compute systematics for each pt bin
    hvn, hsyst, hchi2_vs_trial, hsignificance_vs_trial = [], [], [], []
    gref, gref_two, gsyst = [], [], []
    hsyst_final = hchi2[0].Clone('hsyst_final')
    hsyst_final.SetTitle('hsyst_final;#it{p}_{T} (GeV/#it{c});syst (%)')
    hsyst_final.Reset()
    canvas, leg = [], []
    for i in range(0, gvn_vs_mass[0].GetN()): # loop over pt bins
        canvas.append(TCanvas(f'syst_multitrail_pt{gvn_vs_mass[0].GetX()[i]}', f'syst_multitrail_pt{gvn_vs_mass[0].GetX()[i] - gvn_vs_mass[0].GetEXlow()[i]}_{gvn_vs_mass[0].GetX()[i] + gvn_vs_mass[0].GetEXhigh()[i]}', 800, 800))
        canvas[-1].cd().Divide(2, 2)
        hvn.append(TH1F('hvn', 'hvn;trial;#it{v}_{n}', len(gvn_vs_mass), 0, len(gvn_vs_mass)+1))
        hsignificance_vs_trial.append(TH1F('hsignificance_vs_trial', 'hsignificance_vs_trial;trial;significance', len(gvn_vs_mass), 0, len(gvn_vs_mass)+1))
        hchi2_vs_trial.append(TH1F('hchi2_vs_trial', 'hchi2_vs_trial;trial;#chi^{2}', len(gvn_vs_mass), 0, len(gvn_vs_mass)+1))
        hsyst.append(TH1F('hsyst', 'hsyst;#it{v}_{n};', 40, 0.8, 1.2))
        
        SetObjectStyle(hvn[-1], markerstyle=20, markercolor=kGray+1, markersize=1., linecolor=kGray+1)
        SetObjectStyle(hsignificance_vs_trial[-1], markerstyle=20, markercolor=kGray+1, markersize=1., linecolor=kGray+1)
        SetObjectStyle(hchi2_vs_trial[-1], markerstyle=20, markercolor=kGray+1, markersize=1., linecolor=kGray+1)
        SetObjectStyle(hsyst[-1], markerstyle=20, markercolor=kGray+1, markersize=1., linecolor=kGray+1)
        ipt = i+1
        for j, _ in enumerate(gvn_vs_mass):    # loop over trials
            chi2 = hchi2[j].GetBinContent(ipt)
            hchi2_vs_trial[-1].SetBinContent(j, chi2)
            hchi2_vs_trial[-1].SetBinError(j, hchi2[j].GetBinError(ipt))
            significance = hsignificance[j].GetBinContent(ipt)
            hsignificance_vs_trial[-1].SetBinContent(j, significance)
            hsignificance_vs_trial[-1].SetBinError(j, hsignificance[j].GetBinError(ipt))
            hvn[-1].SetBinContent(j, gvn_vs_mass[j].GetY()[i])
            hvn[-1].SetBinError(j, gvn_vs_mass[j].GetEYlow()[i])

            # Skip chi2 higher than 5 and significance lower than 3
            if (chi2 > 3 and chi2 != 0) or (significance < 3 and significance != 0):
                print(f'Skipping trial {j} for pt bin {ipt} with chi2 = {chi2} and significance = {significance}')
                continue
            hsyst[-1].Fill(gvn_vs_mass[j].GetY()[i] / gvn_vs_mass_default.GetY()[i])

        # Compute systematic uncertainty
        rms = hsyst[-1].GetRMS()
        mean = hsyst[-1].GetMean() - 1
        syst = np.sqrt(rms**2 + mean**2)
        print(f'pt bin {ipt}, mean = {mean}, rms = {rms}, syst = {syst}')
        gsyst.append(TGraphAsymmErrors())
        gsyst[-1].SetPoint(0, 1, hsyst[-1].GetMaximum()*0.5)
        gsyst[-1].SetPointError(0, syst, syst,
                                hsyst[-1].GetMaximum()*0.5, hsyst[-1].GetMaximum()*0.5)
        SetObjectStyle(gsyst[-1], markerstyle=20, markercolor=kOrange+2,
                       markersize=1, linecolor=kOrange+2,
                          linewidth=2, fillcolor=kOrange+2, fillstyle=3153,
                          fillalpha=0.5, linestyle=9)

        # Pad 1: vn vs trial
        canvas[-1].cd(1).SetGrid()
        # Define reference line
        gref.append(TGraphAsymmErrors())
        gref[-1].SetPoint(0, 0, gvn_vs_mass_default.GetY()[i])
        gref[-1].SetPointError(0, 0, 0, gvn_vs_mass_default.GetEYlow()[i],
                           gvn_vs_mass_default.GetEYhigh()[i])
        gref[-1].SetPoint(1, len(gvn_vs_mass), gvn_vs_mass_default.GetY()[i])
        gref[-1].SetPointError(1, 0, 0, gvn_vs_mass_default.GetEYlow()[i],
                           gvn_vs_mass_default.GetEYhigh()[i])
        SetObjectStyle(gref[-1], markerstyle=20, markercolor=kAzure+2,
                       markersize=0, linecolor=kAzure+2,
                       linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        hvn[-1].Draw('same')
        gref[-1].Draw('c3 same')
        canvas[-1].cd(2)

        # Legend
        leg.append(TLegend(0.6, 0.6, 0.9, 0.9))
        leg[-1].SetBorderSize(0)
        leg[-1].SetFillStyle(0)
        leg[-1].SetTextSize(0.03)
        leg[-1].SetHeader(f'{gvn_vs_mass[0].GetX()[i] - gvn_vs_mass[0].GetEXlow()[i]} < #it{{p}}_{{T}} < {gvn_vs_mass[0].GetX()[i] + gvn_vs_mass[0].GetEXhigh()[i]} GeV/#it{{c}}')
        leg[-1].AddEntry(hvn[-1], 'vn vs trial', 'l')
        leg[-1].AddEntry(gsyst[-1], f'#sqrt{{shift^{{2}} + rms^{{2}}}}) = {syst:.2f}', 'f')
        leg[-1].AddEntry(gref[-1], 'vn stat. unc.', 'f')

        # Define reference vertical line at 1
        gref_two.append(gref[-1].Clone())
        gref_two[-1].SetPoint(0, 1, hsyst[-1].GetMaximum()*0.5)
        gref_two[-1].SetPointError(0, gvn_vs_mass_default.GetEYlow()[i],
                                   gvn_vs_mass_default.GetEYhigh()[i],
                                   hsyst[-1].GetMaximum()*0.5, hsyst[-1].GetMaximum()*0.5)
        SetObjectStyle(gref_two[-1], markerstyle=20, markercolor=kAzure+2,
                       markersize=1, linecolor=kAzure+2,
                       linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        hsyst[-1].GetYaxis().SetRangeUser(0, hsyst[-1].GetMaximum()*1.8)
        hsyst[-1].Draw('same')
        gsyst[-1].Draw('2')
        gref_two[-1].Draw('2')
        leg[-1].Draw()
        # Pad 3: chi2 vs trial
        canvas[-1].cd(3)
        hchi2_vs_trial[-1].Draw('same')
        # Pad 4: significance vs trial
        canvas[-1].cd(4)
        hsignificance_vs_trial[-1].Draw('same')
        hsyst_final.SetBinContent(ipt, syst)
        canvas[-1].Update()

    canvsyst = TCanvas('canvsyst', 'canvsyst', 800, 800)
    canvsyst.cd()
    canvsyst.SetGrid()
    SetObjectStyle(hsyst_final, markerstyle=20, markercolor=kAzure+2,
                   markersize=1.,linecolor=kOrange+2,
                   linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hsyst_final.Draw('same e')

    for icanv, canv in enumerate(canvas):
        if icanv == 0:
            suffix_pdf = '('
        elif icanv == len(canvas) - 1:
            suffix_pdf = ')'
        else:
            suffix_pdf = ''
        canv.SaveAs(f'{outputdir}/SystMultitrial.pdf{suffix_pdf}')
    canvsyst.SaveAs(f'{outputdir}/SystMultitrial_vs_pt.pdf)')
    input()

    outdir = os.path.join(outputdir, 'syst_multitrial')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfilename = os.path.join(outdir, 'syst_multitrial.root')
    outfile = TFile(outfilename, 'recreate')
    for h in hvn:
        h.Write()
    for h in hsyst:
        h.Write()
    hsyst_final.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('rypathsyst', metavar='text', default='path to the directory containing the .root files from multitrial')
    parser.add_argument('ry_default', metavar='text', default='default .root file')
    parser.add_argument("--outputdir", "-o", metavar="text", default=".", help="output directory")
    args = parser.parse_args()

    compute_syst_multitrial(args.rypathsyst,
                            args.ry_default,
                            args.outputdir)
