import sys
import os
import numpy as np
import argparse
from ROOT import TFile, TCanvas, TH1F, TGraphAsymmErrors, TLegend, kOrange, kAzure, kBlack
sys.path.append('../../../')
from utils.StyleFormatter import SetObjectStyle, SetGlobalStyle

class SystMultitrial:
    prompt = False
    ptmin = 0
    ptmax = 0

class Trails:
    def __init__(self, trail_path, default=False):
        self.trail_path = trail_path
        self.name = trail_path.split('_')[-1]
        self.v2vsfracfile = self._get_v2_vs_frac_file(trail_path)
        self.ryfiles = self._get_ry_files(trail_path)
        self.hv2vsptprompt = self._get_hv2vsptprompt(self.v2vsfracfile)
        self.gvnsimfits, self.hchi2s, self.hsignificances = self._get_fit_results(self.ryfiles, default)
        # self.gvnsimfits = self._get_gvnSimFit(self.ryfiles)
        if default:
            self.ncut = len(self.ryfiles)
    
    def _get_v2_vs_frac_file(self, trail_path):
        try:
            for file in os.listdir(path=f'{trail_path}/V2VsFrac'):
                if file.endswith('.root'):
                    return os.path.join(trail_path, 'V2VsFrac', file)
        except FileNotFoundError:
            print(f'No V2VsFrac directory found in {trail_path}')
            self.v2vsfracfile = None
            return None
        
    def _get_ry_files(self, trail_path):
        files = []
        try:
            for raw_file in os.listdir(path=f'{trail_path}/ry'):
                if raw_file.endswith('.root'):
                    files.append(os.path.join(trail_path, 'ry', raw_file))
                # make sure the files are sorted
                files.sort()
        except FileNotFoundError:
            print(f'No ry directory found in {trail_path}')
            return []
        return files
    
    def _get_hv2vsptprompt(self, v2vsfracfile):# -> None | Any:
        if v2vsfracfile is None:
            return None
        tfile = TFile().Open(v2vsfracfile)
        if not tfile.GetListOfKeys().Contains('hV2VsPtPrompt'):
            print(f'Warning: {v2vsfracfile} does not contain hV2VsPtPrompt')
            return None
        hv2vsptprompt = tfile.Get('hV2VsPtPrompt')
        hv2vsptprompt.SetDirectory(0)
        tfile.Close()
        return hv2vsptprompt
    
    def _get_fit_results(self, ryfiles, default):
        gvnsimfits = []
        hchi2s = []
        hsignificances = []
        for file in ryfiles:
            tfile = TFile().Open(file)
            if not tfile.GetListOfKeys().Contains('gvnSimFit'):
                print(f'Warning: {file} does not contain gvnSimFit')
                gvnsimfits.append(None)
                hchi2s.append(None)
                hsignificances.append(None)
                tfile.Close()
                continue
            gvnsimfits.append(tfile.Get('gvnSimFit'))
            hchi2s.append(tfile.Get('hRedChi2SimFit'))
            hchi2s[-1].SetDirectory(0)
            hsignificances.append(tfile.Get('hRawYieldsSignificanceSimFit'))
            hsignificances[-1].SetDirectory(0)
            tfile.Close()
        return gvnsimfits, hchi2s, hsignificances

SetGlobalStyle(titleoffsety=1.1, maxdigits=3, topmargin=0.1, bottommargin=0.4, leftmargin=0.3, rightmargin=0.15,
               labelsizey=0.04, setoptstat=0, setopttitle=0, setdecimals=True,titleoffsetx=0.74)

def compute_syst_multitrial(rypathsyst, ry_default, outputdir):
    trails_path = os.listdir(rypathsyst)
    trails = []
    # collect all trails
    for trail_path in trails_path:
        trails.append(Trails(os.path.join(rypathsyst, trail_path)))
        
    # default trail i.e. central value
    default_trail = Trails(ry_default, default=True)
    
    outputdir = os.path.join(outputdir, 'syst_multitrial')
    os.makedirs(outputdir, exist_ok=True)

    if SystMultitrial.prompt:
        # for iFile in range(len(default_trail.ryfiles)):
        #     compute_systematics(outputdir, default_trail.gvnsimfits[iFile], \
        #         [trail.gvnsimfits[iFile] for trail in trails if len(trail.gvnsimfits) > iFile], \
        #                 [trail.hchi2s[iFile] for trail in trails if len(trail.hchi2s) > iFile], \
        #         [trail.hsignificances[iFile] for trail in trails if len(trail.hsignificances) > iFile], \
        #             trails, iFile)
        compute_systematics_prompt(outputdir, default_trail, trails)
    else:
        compute_systematics(outputdir, default_trail.gvnsimfits, \
                [trail.gvnsimfits for trail in trails], \
                        [trail.hchi2s for trail in trails], \
                [trail.hsignificances for trail in trails])

def compute_systematics(outputdir, gvn_vs_mass_default, gvn_vs_mass, hchi2, hsignificance, trails, iFile=0):
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
        hsyst.append(TH1F('hsyst', 'hsyst;#it{v}_{n}^{trial} - #it{v}_{n}^{default};"', 200, -0.2, 0.2))
        
        SetObjectStyle(hvn[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
        SetObjectStyle(hsignificance_vs_trial[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
        SetObjectStyle(hchi2_vs_trial[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
        SetObjectStyle(hsyst[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
        ipt = i+1
        for j, _ in enumerate(gvn_vs_mass):    # loop over trials
            if trails[j].gvnsimfits == None:
                print(f'No prompt v2 vs pt found for trial {j}: {trails[j].trail_path}')
                continue
            significance = hsignificance[j].GetBinContent(ipt)
            chi2 = hchi2[j].GetBinContent(ipt)

            if len(trails[j].ryfiles) <= iFile:
                print(f'No ry file found for trial {j}: {trails[j].trail_path}')
                continue
            # Skip chi2 higher than 5 and significance lower than 3
            if (chi2 > 3 and chi2 != 0) or (significance < 6 or significance > 600):
                print(f'Skipping trial {j}: {trails[j].ryfiles[iFile]} for pt bin {ipt} with chi2 = {chi2} and significance = {significance}')
                continue
            hchi2_vs_trial[-1].SetBinContent(j, chi2)
            hchi2_vs_trial[-1].SetBinError(j, hchi2[j].GetBinError(ipt))
            hsignificance_vs_trial[-1].SetBinContent(j, significance)
            hsignificance_vs_trial[-1].SetBinError(j, hsignificance[j].GetBinError(ipt))
            hvn[-1].SetBinContent(j, gvn_vs_mass[j].GetY()[i])
            hvn[-1].SetBinError(j, gvn_vs_mass[j].GetEYlow()[i])
            hsyst[-1].Fill(gvn_vs_mass[j].GetY()[i] - gvn_vs_mass_default.GetY()[i])
        #input(f'pt bin {ipt} done. Number of trials skipped = {counter}')

        # Compute systematic uncertainty
        rms = hsyst[-1].GetRMS()
        mean = hsyst[-1].GetMean()
        syst = np.sqrt(rms**2 + mean**2)
        print(f'pt bin {ipt}, mean = {mean}, rms = {rms}, syst = {syst}')
        gsyst.append(TGraphAsymmErrors())
        gsyst[-1].SetPoint(0, 0, hsyst[-1].GetMaximum()*0.5)
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
        leg[-1].AddEntry(hvn[-1], f'vn vs trial ({hvn[-1].GetEntries()} trials)', 'p')
        leg[-1].AddEntry(gsyst[-1], f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst:.3f}', 'f')
        leg[-1].AddEntry(gref[-1], 'vn stat. unc.', 'f')

        # Define reference vertical line at 1
        gref_two.append(gref[-1].Clone())
        gref_two[-1].SetPoint(0, 0, hsyst[-1].GetMaximum()*0.5)
        gref_two[-1].SetPointError(0, gvn_vs_mass_default.GetEYlow()[i],
                                   gvn_vs_mass_default.GetEYhigh()[i],
                                   hsyst[-1].GetMaximum()*0.5, hsyst[-1].GetMaximum()*0.5)
        SetObjectStyle(gref_two[-1], markerstyle=20, markercolor=kAzure+2,
                       markersize=1, linecolor=kAzure+2,
                       linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        hsyst[-1].GetYaxis().SetRangeUser(0, hsyst[-1].GetMaximum()*1.8)
        hsyst[-1].GetXaxis().SetRangeUser(hsyst[-1].GetXaxis().GetXmin(), hsyst[-1].GetXaxis().GetXmax())
        hsyst[-1].Draw('same')
        gsyst[-1].Draw('2')
        gref_two[-1].Draw('2')
        gsyst[-1].Draw('2')
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
        if len(canvas) == 1:
            suffix_pdf = ''
        canv.SaveAs(f'{outputdir}/SystRy_{iFile}.pdf{suffix_pdf}')
    canvsyst.SaveAs(f'{outputdir}/SystRy_vs_pt_{iFile}.pdf)')

    outdir = os.path.join(outputdir, 'syst_multitrial')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfilename = os.path.join(outdir, f'syst_ry_{iFile}.root')
    outfile = TFile(outfilename, 'recreate')
    for h in hvn:
        h.Write()
    for h in hsyst:
        h.Write()
    hsyst_final.Write()
    outfile.Close()

def compute_systematics_prompt(outputdir, default_trail, trails):

    hvn, hsyst, hchi2_vs_trial, hsignificance_vs_trial = [], [], [], []
    gref, gref_two, gsyst = [], [], []
    hsyst_final = trails[0].hchi2s[0].Clone('hsyst_final')
    hsyst_final.SetTitle('hsyst_final;#it{p}_{T} (GeV/#it{c});syst (%)')
    hsyst_final.Reset()
    canvas, leg = [], []
    for i in range(0, default_trail.gvnsimfits[0].GetN()): # loop over pt bins
        canvas.append(TCanvas(f'syst_multitrail_pt{default_trail.gvnsimfits[0].GetX()[i]}', f'syst_multitrail_pt{default_trail.gvnsimfits[0].GetX()[i] - default_trail.gvnsimfits[0].GetEXlow()[i]}_{default_trail.gvnsimfits[0].GetX()[i] + default_trail.gvnsimfits[0].GetEXhigh()[i]}', 800, 800))
        canvas[-1].cd().Divide(1, 2)
        hvn.append(TH1F('hvn', 'hv2;trial;prompt #it{v}_{2}', len(trails), 0, len(trails)+1))
        hsyst.append(TH1F('hsyst', 'hsyst;prompt #it{v}_{2}^{trial} - prompt #it{v}_{2}^{default};"', 200, -0.2, 0.2))
        
        SetObjectStyle(hvn[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
        SetObjectStyle(hsyst[-1], markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)

        ipt = i+1
        for j, _ in enumerate(trails):
            if trails[j].hv2vsptprompt is None:
                print(f'No prompt v2 vs pt found for trial {j}: {trails[j].trail_path}')
                continue
            
            skip_trial = False

            for iFile in range(len(trails[j].ryfiles)):
                significance = trails[j].hsignificances[iFile].GetBinContent(ipt)
                chi2 = trails[j].hchi2s[iFile].GetBinContent(ipt)

                if (chi2 > 3 and chi2 != 0) or (significance < 6 or significance > 600):
                    print(f'Skipping trial {j}: {trails[j].ryfiles[iFile]} for pt bin {ipt} with chi2 = {chi2} and significance = {significance}')
                    skip_trial = True
                    break
            
            if skip_trial:
                continue

            if trails[j].hv2vsptprompt.GetBinContent(ipt) == 0:
                print(f'v2 in pt bin {ipt} for trial {j} is 0: {trails[j].trail_path}')
                continue
            if trails[j].hv2vsptprompt.GetBinContent(ipt) < 0:
                print(f'v2 in pt bin {ipt} for trial {j} is negative: {trails[j].trail_path}')
                print(f'v2 = {trails[j].hv2vsptprompt.GetBinContent(ipt)}')
                continue

            hvn[-1].SetBinContent(j, trails[j].hv2vsptprompt.GetBinContent(ipt))
            hvn[-1].SetBinError(j, trails[j].hv2vsptprompt.GetBinError(ipt))
            print(i)
            hsyst[-1].Fill(trails[j].hv2vsptprompt.GetBinContent(ipt) - default_trail.hv2vsptprompt.GetBinContent(ipt))
        #input(f'pt bin {ipt} done. Number of trials skipped = {counter}')

        # Compute systematic uncertainty
        rms = hsyst[-1].GetRMS()
        mean = hsyst[-1].GetMean()
        syst = np.sqrt(rms**2 + mean**2)
        print(f'pt bin {ipt}, mean = {mean}, rms = {rms}, syst = {syst}')
        gsyst.append(TGraphAsymmErrors())
        gsyst[-1].SetPoint(0, 0, hsyst[-1].GetMaximum()*0.5)
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
        gref[-1].SetPoint(0, 0, default_trail.hv2vsptprompt.GetBinContent(ipt))
        gref[-1].SetPointError(0, 0, 0, default_trail.hv2vsptprompt.GetBinError(ipt),
                           default_trail.hv2vsptprompt.GetBinError(ipt))
        gref[-1].SetPoint(1, len(trails), default_trail.hv2vsptprompt.GetBinContent(ipt))
        gref[-1].SetPointError(1, 0, 0, default_trail.hv2vsptprompt.GetBinError(ipt),
                           default_trail.hv2vsptprompt.GetBinError(ipt))
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
        leg[-1].SetHeader(f'{default_trail.gvnsimfits[0].GetX()[i] - default_trail.gvnsimfits[0].GetEXlow()[i]} < #it{{p}}_{{T}} < {default_trail.gvnsimfits[0].GetX()[i] + default_trail.gvnsimfits[0].GetEXhigh()[i]} GeV/#it{{c}}')
        leg[-1].AddEntry(hvn[-1], f'v2 vs trial ({hvn[-1].GetEntries()} trials)', 'p')
        leg[-1].AddEntry(gsyst[-1], f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst:.3f}', 'f')
        leg[-1].AddEntry(gref[-1], 'v2 stat. unc.', 'f')

        # Define reference vertical line at 1
        gref_two.append(gref[-1].Clone())
        gref_two[-1].SetPoint(0, 0, hsyst[-1].GetMaximum()*0.5)
        gref_two[-1].SetPointError(0, default_trail.hv2vsptprompt.GetBinError(ipt),
                                   default_trail.hv2vsptprompt.GetBinError(ipt),
                                   hsyst[-1].GetMaximum()*0.5, hsyst[-1].GetMaximum()*0.5)
        SetObjectStyle(gref_two[-1], markerstyle=20, markercolor=kAzure+2,
                       markersize=1, linecolor=kAzure+2,
                       linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
        hsyst[-1].GetYaxis().SetRangeUser(0, hsyst[-1].GetMaximum()*1.8)
        hsyst[-1].GetXaxis().SetRangeUser(hsyst[-1].GetXaxis().GetXmin(), hsyst[-1].GetXaxis().GetXmax())
        hsyst_final.SetBinContent(ipt, syst)
        hsyst[-1].Draw('same')
        gsyst[-1].Draw('2')
        gref_two[-1].Draw('2')
        gsyst[-1].Draw('2')
        leg[-1].Draw()
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
        if len(canvas) == 1:
            suffix_pdf = ''
        canv.SaveAs(f'{outputdir}/SystPromptv2.pdf{suffix_pdf}')
    canvsyst.SaveAs(f'{outputdir}/SystPromptv2_vs_pt.pdf)')
    input()

    outdir = os.path.join(outputdir, 'syst_multitrial')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfilename = os.path.join(outdir, f'syst_prompt.root')
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
    parser.add_argument("--prompt", "-p", action="store_true", help="compute systematics for prompt D mesons")
    args = parser.parse_args()

    if args.prompt:
        SystMultitrial.prompt = True
        

    compute_syst_multitrial(args.rypathsyst,
                            args.ry_default,
                            args.outputdir)