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
    
def get_ry_files(rypath, cutvar_path):
    files = []
    for raw_file in os.listdir(path=f'{rypath}/{cutvar_path}/ry'):
        if raw_file.endswith('.root'):
            files.append(raw_file)
        files.sort()
    return files

SetGlobalStyle(titleoffsety=1.1, maxdigits=3, topmargin=0.1, bottommargin=0.4, leftmargin=0.3, rightmargin=0.15,
               labelsizey=0.04, setoptstat=0, setopttitle=0, setdecimals=True,titleoffsetx=0.74)

def compute_syst_multitrial(rypathsyst, ry_default, outputdir):
    #______________________________________________________________________________________
    # Collect all .root files from rypathsyst
    if SystMultitrial.prompt:
        prompt_files = []
        root_files_dict = {}
        trails_path = os.listdir(rypathsyst)
        for trail_path in trails_path:
            try:
                for file in os.listdir(path=f'{rypathsyst}/{trail_path}/V2VsFrac'):
                    if file.endswith('.root'):
                        prompt_files.append(f'{rypathsyst}/{trail_path}/V2VsFrac/{file}')
            except FileNotFoundError:
                print(f'No V2VsFrac directory found in {rypathsyst}/{trail_path}')
                continue

            try:
                files = get_ry_files(rypath=rypathsyst, cutvar_path=trail_path)
            except FileNotFoundError:
                print(f'No ry directory found in {rypathsyst}/{trail_path}')
                continue
            
            for iFile, file in enumerate(files):
                if iFile in root_files_dict:
                    root_files_dict[iFile].append(f'{rypathsyst}/{trail_path}/ry/{file}')
                else:
                    root_files_dict[iFile] = [f'{rypathsyst}/{trail_path}/ry/{file}']
                
    else:
        root_files = []
        for file in os.listdir(rypathsyst):
            if file.endswith('.root'):
                root_files.append(file)
        if not root_files:
            sys.exit(f'ERROR: No .root files found in {rypathsyst}')

    #______________________________________________________________________________________
    # Get vn vs mass from default .root file
    if SystMultitrial.prompt:
        #python3 compute_syst_multitrial.py 
        # /home/wuct/ALICE/local/Results/test/flow/systematics/trails/trails_pt_10_15 
        # /home/wuct/ALICE/local/Results/test/flow/systematics/pre 
        # -o /home/wuct/ALICE/local/Results/test/flow/systematics/sys -p
        
        ptmin = float(rypathsyst.split('/')[-1].split('_')[-2]) / 10.0
        ptmax = float(rypathsyst.split('/')[-1].split('_')[-1]) / 10.0
        print(f'ptmin = {ptmin}, ptmax = {ptmax}')
        SystMultitrial.ptmin = ptmin
        SystMultitrial.ptmax = ptmax
        
        gvn_vs_mass_defaults = []
        tfiles = get_ry_files(rypath=ry_default, cutvar_path='cutvar_central')
        print(f'tfiles = {tfiles}')
        for iFile in range(len(root_files_dict)):
            tfile = TFile().Open(ry_default + '/cutvar_central/ry/' + tfiles[iFile])
            gvn_vs_mass_defaults.append(tfile.Get('gvnSimFit'))
            tfile.Close()
        
        # need to be updated to the combined file
        # prompt_tfile = TFile().Open(f'{ry_default}/cutvar_correlated/V2VsFrac/V2VsFrac_correlated.root')
        prompt_tfile = TFile().Open('/home/wuct/ALICE/local/Results/BDT/k3050/approve_2_11/combined/V2VsPt.root')
        # hvn_vs_frac_default = prompt_tfile.Get('hV2VsPtPrompt')
        hvn_vs_frac_default = prompt_tfile.Get('hv2Prompt')
        hvn_vs_frac_default.SetDirectory(0)
        prompt_tfile.Close()

    else:
        tfile = TFile().Open(ry_default)
        gvn_vs_mass_default = tfile.Get('gvnSimFit')

    #______________________________________________________________________________________
    # Get vn vs mass from all .root files
    if SystMultitrial.prompt:
        gvn_vs_mass_dict = {}
        hchi2s_dict = {}
        hsignificances_dict = {}
        hprompt_v2s = []
        for iFile, root_file in enumerate(root_files_dict.values()):
            gvn_vs_mass_dict[iFile] = []
            hchi2s_dict[iFile] = []
            hsignificances_dict[iFile] = []
            for file in root_file:
                tfile = TFile().Open(file)
                # check if the file has the right keys
                if not tfile.GetListOfKeys().Contains('gvnSimFit'):
                    print(f'ERROR: {file} does not contain gvnSimFit')
                    continue
                gvn_vs_mass_dict[iFile].append(tfile.Get('gvnSimFit'))
                hchi2s_dict[iFile].append(tfile.Get('hRedChi2SimFit'))
                hchi2s_dict[iFile][-1].SetDirectory(0)
                hsignificances_dict[iFile].append(tfile.Get('hRawYieldsSignificanceSimFit'))
                hsignificances_dict[iFile][-1].SetDirectory(0)
                tfile.Close()
        for file in prompt_files:
            tfile = TFile().Open(file)
            # check if the file has the right keys
            if not tfile.GetListOfKeys().Contains('hV2VsPtPrompt'):
                print(f'ERROR: {file} does not contain hV2VsPtPrompt')
                continue
            hprompt_v2s.append(tfile.Get('hV2VsPtPrompt'))
            hprompt_v2s[-1].SetDirectory(0)
            tfile.Close()
    else:
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
    if SystMultitrial.prompt:
        outputdir = outputdir + f'/pt_{int(ptmin*10)}_{int(ptmax*10)}'
        os.makedirs(outputdir, exist_ok=True)
        for iFile in range(len(root_files_dict)):
            compute_systematics(outputdir, gvn_vs_mass_defaults[iFile], gvn_vs_mass_dict[iFile], hchi2s_dict[iFile], hsignificances_dict[iFile], iFile)
            compute_systematics_prompt(outputdir, hvn_vs_frac_default, hprompt_v2s)
    else:
        compute_systematics(outputdir, gvn_vs_mass_default, gvn_vs_mass, hchi2, hsignificance)

def compute_systematics(outputdir, gvn_vs_mass_default, gvn_vs_mass, hchi2, hsignificance, iFile=0):
    hvn, hsyst, hchi2_vs_trial, hsignificance_vs_trial = [], [], [], []
    gref, gref_two, gsyst = [], [], []
    hsyst_final = hchi2[0].Clone('hsyst_final')
    hsyst_final.SetTitle('hsyst_final;#it{p}_{T} (GeV/#it{c});syst (%)')
    hsyst_final.Reset()
    canvas, leg = [], []
    for i in range(0, gvn_vs_mass[0].GetN()): # loop over pt bins
        if SystMultitrial.prompt:
            xPoints = gvn_vs_mass_default.GetX()
            for iPoint, xPt in enumerate(xPoints):
                if xPt > SystMultitrial.ptmin and xPt < SystMultitrial.ptmax:
                    break
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
            significance = hsignificance[j].GetBinContent(ipt)
            chi2 = hchi2[j].GetBinContent(ipt)

            # Skip chi2 higher than 5 and significance lower than 3
            if (chi2 > 10 and chi2 != 0) or (significance < 1 and significance != 0):
                print(f'Skipping trial {j} for pt bin {ipt} with chi2 = {chi2} and significance = {significance}')
                continue
            hchi2_vs_trial[-1].SetBinContent(j, chi2)
            hchi2_vs_trial[-1].SetBinError(j, hchi2[j].GetBinError(ipt))
            hsignificance_vs_trial[-1].SetBinContent(j, significance)
            hsignificance_vs_trial[-1].SetBinError(j, hsignificance[j].GetBinError(ipt))
            hvn[-1].SetBinContent(j, gvn_vs_mass[j].GetY()[i])
            hvn[-1].SetBinError(j, gvn_vs_mass[j].GetEYlow()[i])
            if SystMultitrial.prompt:
                hsyst[-1].Fill(gvn_vs_mass[j].GetY()[i] - gvn_vs_mass_default.GetY()[iPoint])
            else:
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
        if SystMultitrial.prompt:
            gref[-1].SetPoint(0, 0, gvn_vs_mass_default.GetY()[iPoint])
            gref[-1].SetPointError(0, 0, 0, gvn_vs_mass_default.GetEYlow()[iPoint],
                                   gvn_vs_mass_default.GetEYhigh()[iPoint])
            gref[-1].SetPoint(1, len(gvn_vs_mass), gvn_vs_mass_default.GetY()[iPoint])
            gref[-1].SetPointError(1, 0, 0, gvn_vs_mass_default.GetEYlow()[iPoint],
                                   gvn_vs_mass_default.GetEYhigh()[iPoint])
        else:
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
        if SystMultitrial.prompt:
            gref_two[-1].SetPointError(0, gvn_vs_mass_default.GetEYlow()[iPoint],
                                   gvn_vs_mass_default.GetEYhigh()[iPoint],
                                   hsyst[-1].GetMaximum()*0.5, hsyst[-1].GetMaximum()*0.5)
        else:
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
        if len(canvas) > 1:
            if icanv == 0:
                suffix_pdf = '('
            elif icanv == len(canvas) - 1:
                suffix_pdf = ')'
            else:
                suffix_pdf = ''
            canv.SaveAs(f'{outputdir}/SystMultitrial.pdf{suffix_pdf}')
    if SystMultitrial.prompt:
        canvas[-1].SaveAs(f'{outputdir}/SystMultitrial_ry_{iFile}.pdf')
        canvsyst.SaveAs(f'{outputdir}/SystMultitrial_ry_pt_{iFile}.pdf')
    else:
        canvas[-1].SaveAs(f'{outputdir}/SystMultitrial.pdf')
        canvsyst.SaveAs(f'{outputdir}/SystMultitrial_vs_pt.pdf)')
    input()

    outdir = os.path.join(outputdir, 'syst_multitrial')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    if SystMultitrial.prompt:
        outfilename = os.path.join(outdir, f'syst_multitrial_ry_{iFile}.root')
    else:
        outfilename = os.path.join(outdir, 'syst_multitrial.root')
    outfile = TFile(outfilename, 'recreate')
    for h in hvn:
        h.Write()
    for h in hsyst:
        h.Write()
    hsyst_final.Write()
    outfile.Close()

def compute_systematics_prompt(outputdir, hvn_vs_frac_default, hprompt_v2s):

    hvn, hsyst = [], []
    gref, gref_two, gsyst = [], [], []
    hsyst_final = hprompt_v2s[0].Clone('hsyst_final')
    hsyst_final.SetTitle('hsyst_final;#it{p}_{T} (GeV/#it{c});syst (%)')
    hsyst_final.Reset()
    canvas, leg = [], []

    canvas = TCanvas('sys_prompt_v2_multitrial', 'sys_prompt_v2_multitrial', 800, 800)
    canvas.cd().Divide(1, 2)
    hvn = TH1F('hvn', 'hvn;trial;#it{v}_{2}', len(hprompt_v2s), 0, len(hprompt_v2s)+1)
    hsyst = TH1F('hsyst', 'hsyst;#it{v}_{2}^{trial} - #it{v}_{2}^{default};"', 200, -0.2, 0.2)
    
    SetObjectStyle(hvn, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    SetObjectStyle(hsyst, markerstyle=20, markercolor=kBlack, markersize=1., linecolor=kBlack)
    
    iBin_default = hvn_vs_frac_default.FindBin((SystMultitrial.ptmin + SystMultitrial.ptmax) / 2)
    print(f'iBin_default = {iBin_default}')
    
    for i, prompt_v2 in enumerate(hprompt_v2s):    # loop over trials
        hvn.SetBinContent(i, prompt_v2.GetBinContent(1))
        hvn.SetBinError(i, prompt_v2.GetBinError(1))
        hsyst.Fill(prompt_v2.GetBinContent(1) - hvn_vs_frac_default.GetBinContent(iBin_default))
    #input(f'pt bin {ipt} done. Number of trials skipped = {counter}')
    
    rms = hsyst.GetRMS()
    mean = hsyst.GetMean()
    syst = np.sqrt(rms**2 + mean**2)
    print(f'mean = {mean}, rms = {rms}, syst = {syst}')
    gsyst = TGraphAsymmErrors()
    gsyst.SetPoint(0, 0, hsyst.GetMaximum()*0.5)
    gsyst.SetPointError(0, syst, syst,
                        hsyst.GetMaximum()*0.5, hsyst.GetMaximum()*0.5)
    SetObjectStyle(gsyst, markerstyle=20, markercolor=kOrange+2,
                    markersize=1, linecolor=kOrange+2,
                        linewidth=2, fillcolor=kOrange+2, fillstyle=3153,
                        fillalpha=0.5, linestyle=9)
    
    # Pad 1: prompt v2 vs trial
    canvas.cd(1).SetGrid()
    # Define reference line
    gref = TGraphAsymmErrors()
    gref.SetPoint(0, 0, hvn_vs_frac_default.GetBinContent(iBin_default))
    gref.SetPointError(0, 0, 0, hvn_vs_frac_default.GetBinError(iBin_default),
                        hvn_vs_frac_default.GetBinError(iBin_default))
    gref.SetPoint(1, len(hprompt_v2s), hvn_vs_frac_default.GetBinContent(iBin_default))
    gref.SetPointError(1, 0, 0, hvn_vs_frac_default.GetBinError(iBin_default),
                        hvn_vs_frac_default.GetBinError(iBin_default))
    SetObjectStyle(gref, markerstyle=20, markercolor=kAzure+2,
                    markersize=0, linecolor=kAzure+2,
                    linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    hvn.Draw('same')
    gref.Draw('c3 same')
    canvas.cd(2)
    
    # Legend
    leg = TLegend(0.6, 0.6, 0.9, 0.9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.SetHeader(f'{SystMultitrial.ptmin} < #it{{p}}_{{T}} < {SystMultitrial.ptmax} GeV/#it{{c}}')
    leg.AddEntry(hvn, f'v2 vs trial ({hvn.GetEntries()} trials)', 'p')
    leg.AddEntry(gsyst, f'#sqrt{{shift^{{2}} + rms^{{2}}}} = {syst:.3f}', 'f')
    leg.AddEntry(gref, 'v2 stat. unc.', 'f')
    
    # Define reference vertical line at 1
    gref_two = gref.Clone()
    gref_two.SetPoint(0, 0, hsyst.GetMaximum()*0.5)
    gref_two.SetPointError(0, hvn_vs_frac_default.GetBinError(iBin_default),
                        hvn_vs_frac_default.GetBinError(iBin_default),
                        hsyst.GetMaximum()*0.5, hsyst.GetMaximum()*0.5)
    SetObjectStyle(gref_two, markerstyle=20, markercolor=kAzure+2,
                    markersize=1, linecolor=kAzure+2,
                    linewidth=2, fillcolor=kAzure+2, fillstyle=3135, fillalpha=0.5, linestyle=9)
    hsyst.GetYaxis().SetRangeUser(0, hsyst.GetMaximum()*1.8)
    hsyst.GetXaxis().SetRangeUser(hsyst.GetXaxis().GetXmin(), hsyst.GetXaxis().GetXmax())
    hsyst_final.SetBinContent(1, syst)
    hsyst.Draw('same')
    gsyst.Draw('2')
    gref_two.Draw('2')
    gsyst.Draw('2')
    leg.Draw()
    canvas.Update()
    canvas.SaveAs(f'{outputdir}/SystMultitrial_Prompt_v2.pdf')
    input()
    
    canvsyst = TCanvas('sys_prompt_v2_multitrial', 'sys_prompt_v2_multitrial', 800, 800)
    canvsyst.cd()
    canvsyst.SetGrid()
    SetObjectStyle(hsyst_final, markerstyle=20, markercolor=kAzure+2,
                     markersize=1.,linecolor=kOrange+2,
                     linewidth=2, fillcolor=kOrange+2, fillstyle=3135, fillalpha=0.5)
    hsyst_final.Draw('same e')
    canvsyst.SaveAs(f'{outputdir}/SystMultitrial_Prompt_v2_pt.pdf')
    
    outdir = os.path.join(outputdir, 'syst_multitrial')
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfilename = os.path.join(outdir, 'syst_multitrial.root')
    outfile = TFile(outfilename, 'recreate')
    hvn.Write()
    hsyst.Write()
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
