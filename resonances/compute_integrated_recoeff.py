'''
python script to compute the integrated reco efficiency for resonances
tip: use after compute_reso_eff.py
run: python compute_integrated_recoeff.py
'''
import sys
import argparse
import yaml
import numpy as np
import pandas as pd
import uproot
from  multiprocessing import Process
from alive_progress import alive_bar
from particle import Particle
from ROOT import TFile, TH2F, TH1F, TCanvas, TLatex, kBlack, kRed, kAzure, TLegend, TGraphAsymmErrors, TEfficiency
sys.path.insert(0, '..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
SetGlobalStyle(padtopmargin=0.05,padbottommargin=0.15, palette=55, labelsize=0.03, titlesize=0.05,
               labeloffset=0.01, titleoffsety=1.2, titleoffsetx=1.2, titleoffsetz=1.,
               opttitle=0, optstat=0)

def main():
    ptmins = [2.]
    ptmaxs = [24.]
    resos = ['10433']
    mults = ['MB', 'HM'] # 'HM' or 'MB
    ptshapes = ['DsVeryHarderptshape', 'DsHarderptshape', 'DsVeryHarderptshape', 'DsHarderptshape'] # 'Dsptshape' or 'Lcptshape'

    for ieff, (ptmin, ptmax, reso, mult, ptshape) in enumerate(zip(ptmins, ptmaxs, resos)):
        for mult in mults:
            for ptshape in ptshapes:
                print(f'Processing {reso} - {mult} - {ptshape}')
                #______________________________________________________________________________
                # Load data
                if reso == '435':
                    path = f'./435_Ds2_Dplus_K0S/{mult}'
                else:
                    path = f'./10433_Ds1_Dstar_K0S/{mult}'
                print(f'Loading {path}/eff_reso{reso}_{mult}_{ptshape}.root')
                input_file = TFile(f'{path}/eff_reso{reso}_{mult}_{ptshape}.root')
                histo_name_prompt = 'hreco_pt'
                histo_name_nonprompt = 'hreco_pt_np'
                hprompt = input_file.Get(histo_name_prompt)
                hnonprompt = input_file.Get(histo_name_nonprompt)
                hden = input_file.Get('hgen_pt')
                SetObjectStyle(hprompt, color=kRed+1, markerstyle=20, markersize=1.2,
                               linecolor=kRed+1, markercolor=kRed+1, fillalpha=0.2, fillstyle=0,
                               linewidth=2)
                SetObjectStyle(hnonprompt, color=kAzure+1, markerstyle=20, markersize=1.2,
                               linecolor=kAzure+1, markercolor=kAzure+1, fillalpha=0.2, fillstyle=0,
                               linewidth=2)
                SetObjectStyle(hden, color=kBlack, markerstyle=20, markersize=1.2,
                               linecolor=kBlack, markercolor=kBlack, fillalpha=0.2, fillstyle=0,
                               linewidth=2)

                #______________________________________________________________________________
                # Compute efficiency
                outfile = TFile(f'{path}/integrated_recoeff_reso{reso}_{mult}_{ptshape}.root', 'recreate')
                geff_int_prompt = TGraphAsymmErrors()
                geff_int_prompt.SetName(f'geff_int_prompt_{reso}_{mult}_{ptshape}')
                geff_int_nonprompt = TGraphAsymmErrors()
                geff_int_nonprompt.SetName(f'geff_int_nonprompt_{reso}_{mult}_{ptshape}')
                # check ingredients
                canv = TCanvas(f'c_gen_reco_{reso}_{mult}_{ptshape}', '', 800, 600)
                canv.cd().SetLogy()
                hden.Draw('hist e')
                hprompt.Draw('same hist e')
                hnonprompt.Draw('same hist e')
                canv.Write()

                eff_prompt, eff_nonprompt, eff_unc_prompt, eff_unc_nonprompt,\
                gen, gen_unc, reco_prompt, reco_prompt_unc, reco_nonprompt, reco_nonprompt_unc = (0. for i in range(10))
                for ibin in range(1, hden.GetNbinsX()+1):
                    lowedge = hden.GetBinLowEdge(ibin)
                    upedge = hden.GetBinLowEdge(ibin) + hden.GetBinWidth(ibin)
                    if lowedge >= ptmin and upedge <= ptmax:
                        gen += hden.GetBinContent(ibin)*hden.GetBinWidth(ibin)
                        gen_unc += hden.GetBinError(ibin)**2

                for ibin in range(1, hprompt.GetNbinsX()+1):
                    lowedge = hprompt.GetBinLowEdge(ibin)
                    upedge = hprompt.GetBinLowEdge(ibin) + hprompt.GetBinWidth(ibin)
                    if lowedge >= ptmin and upedge <= ptmax:
                        reco_prompt += hprompt.GetBinContent(ibin)*hprompt.GetBinWidth(ibin)
                        reco_prompt_unc += hprompt.GetBinError(ibin)**2

                for ibin in range(1, hnonprompt.GetNbinsX()+1):
                    lowedge = hnonprompt.GetBinLowEdge(ibin)
                    upedge = hnonprompt.GetBinLowEdge(ibin) + hnonprompt.GetBinWidth(ibin)
                    if lowedge >= ptmin and upedge <= ptmax:
                        reco_nonprompt += hnonprompt.GetBinContent(ibin)*hnonprompt.GetBinWidth(ibin)
                        reco_nonprompt_unc += hnonprompt.GetBinError(ibin)**2

                reco_prompt_unc = np.sqrt(reco_prompt_unc)/reco_prompt if reco_prompt > 0 else 0
                reco_nonprompt_unc = np.sqrt(reco_nonprompt_unc)/reco_nonprompt if reco_nonprompt > 0 else 0
                gen_unc = np.sqrt(gen_unc)/gen if gen > 0 else 0
                print('gen: {} +- {}'.format(gen, gen_unc))
                print('reco_prompt: {} +- {}'.format(reco_prompt, np.sqrt(reco_prompt_unc)))
                print('reco_nonprompt: {} +- {}'.format(reco_nonprompt, np.sqrt(reco_nonprompt_unc)))

                eff_prompt += reco_prompt / gen if gen > 0 else 0
                eff_prompt /= (ptmax-ptmin)
                eff_nonprompt += reco_nonprompt / gen if gen > 0 else 0
                eff_nonprompt /= (ptmax-ptmin)

                eff_unc_prompt = np.sqrt(reco_prompt_unc**2 +gen_unc**2) * eff_prompt
                eff_unc_nonprompt = np.sqrt(reco_nonprompt_unc**2 +gen_unc**2) * eff_nonprompt
                print('eff_prompt: {} +- {}'.format(eff_prompt, eff_unc_prompt))
                print('eff_nonprompt: {} +- {}'.format(eff_nonprompt, eff_unc_nonprompt))
                #input('Press enter to continue')
        
                geff_int_prompt.SetPoint(0, 0.5*(ptmin+ptmax), eff_prompt)
                geff_int_prompt.SetPointError(0, 0.5*(ptmax-ptmin), 0.5*(ptmax-ptmin), eff_unc_prompt, eff_unc_prompt)
                geff_int_nonprompt.SetPoint(0, 0.5*(ptmin+ptmax), eff_nonprompt)
                geff_int_nonprompt.SetPointError(0, 0.5*(ptmax-ptmin),
                                                 0.5*(ptmax-ptmin),
                                                 eff_unc_nonprompt,
                                                 eff_unc_nonprompt)
                SetObjectStyle(geff_int_prompt, color=kRed+1,
                               markerstyle=20, markersize=1.2,
                               linecolor=kRed+1, markercolor=kRed+1,
                               fillalpha=0.2, fillstyle=0,
                               linewidth=2)
                SetObjectStyle(geff_int_nonprompt, color=kAzure+4,
                               markerstyle=20, markersize=1.2,
                               linecolor=kAzure+4, markercolor=kAzure+4,
                               fillalpha=0.2, fillstyle=0,
                               linewidth=2)
                canv = TCanvas('canv', 'canv', 800, 800)
                leg = TLegend(0.2, 0.7, 0.8, 0.9)
                leg.SetBorderSize(0)
                leg.SetTextSize(0.05)
                reso_label = Particle.from_pdgid(reso).name
                leg.SetHeader(f'{reso_label} - {mult} - {ptshape}')
                ymax = 1.#* max(effp, effnp)
                canv.cd().SetLogy()
                hFrame = canv.cd().DrawFrame(0, 1.e-5, 25, ymax, ';;#varepsilon^{int.}')
                geff_int_prompt.Draw('same pz')
                geff_int_nonprompt.Draw('same pz')
                leg.Draw()

                canv.Write()
                geff_int_prompt.Write()
                geff_int_nonprompt.Write()
                hprompt.Write()
                hnonprompt.Write()
                outfile.Close()
                canv.SaveAs(f'{path}/integrated_recoeff_reso{reso}_{mult}_{ptshape}.pdf')
    input('Press enter to exit')

main()
