'''
python script to compute the selection efficiency for resonances
run: python compute_reso_eff.py
'''
import sys
import os
import itertools
import numpy as np
import uproot
from alive_progress import alive_bar
from particle import Particle
from ROOT import TFile, TH2F, TH1F, TCanvas, kRed, kAzure, TLegend, gRandom # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

SetGlobalStyle(padtopmargin=0.05, padrightmargin=0.15, padleftmargin=0.15, padbottommargin=0.15,
               palette=55, labelsize=0.03, titlesize=0.03, labeloffset=0.01,
               titleoffsety=1.8, titleoffsetx=1.8, titleoffsetz=1., opttitle=0)

resonances = [413]  # 413, 435, 10433
mults = ['MB', 'HM']
ptshapes = ['Dsptshape']  # 'Dsptshape' or 'Lcptshape'
set_seed = 42
infile_path = '~/Desktop/Analyses/Ds_resonances'

def main():
    if set_seed:
        gRandom.SetSeed(set_seed)

    for reso in resonances:
        reso_name = Particle.from_pdgid(reso).name
        if reso == 435:
            meson = 'Dplus'
            decay_chain = f'{reso_name} #rightarrow D^{{+}} K^{{0}}_{{S}}'
        elif reso == 10433:
            meson = 'Dstar'
            decay_chain = f'{reso_name} #rightarrow D^{{*+}} K^{{0}}_{{S}}'
        elif reso == 413:
            meson = 'Dstar'
            decay_chain = f'{reso_name} #rightarrow D^{{0}} #pi^{{+}}'
        else:
            print(f'\033[1m\033[91mResonance {reso} not supported\033[0m')
            sys.exit()
        print(f'\033[1m\033[92mComputing efficiency for {reso_name} ({reso})\033[0m')

        for (mult, ptshape) in itertools.product(mults, ptshapes):
            # Load data
            data = os.path.join(infile_path, f'inputs/kinematics/Decays_{reso}_Pythia8_decayer_{mult}_{ptshape}.root')
            if reso == 413 and (mult == 'HM' or ptshape == 'Lcptshape'):
                print('\033[1m\033[91mERROR: D* efficiency for HM and LC not available\033[0m')
                sys.exit()

            print(f'\033[1mLoading kine data from {data}\033[0m')
            reso_df_kine = uproot.open(data)['tupleDreso'].arrays(library='pd')

            # phi is defined in [-pi, pi]
            reso_df_kine['phi_reso'] = reso_df_kine['phi_reso'].apply(lambda x: x + np.pi)
            reso_df_kine['phi_d'] = reso_df_kine['phi_d'].apply(lambda x: x + np.pi)
            reso_df_kine['phi_light'] = reso_df_kine['phi_light'].apply(lambda x: x + np.pi)

            if reso == 435:
                path = os.path.join(infile_path, 'eff/435_Ds2_Dplus_K0S')
            elif reso == 10433:
                path = os.path.join(infile_path, 'eff/10433_Ds1_Dstar_K0S')
            else:
                path = os.path.join(infile_path, 'check_efficiency_Dstar/compute_reso_eff')

            if reso != 413:
                dau_eff_file = TFile.Open(f'{path}/ResoEff_{reso}_{mult}_frac_{reso}.root')
                print(f'\033[1mLoading dau eff data from {dau_eff_file.GetName()}\033[0m')
                heff_v0 = dau_eff_file.Get(f'K0S_{meson}/hEff_{meson}_K0S')
                heff_D_prompt = dau_eff_file.Get(f'{meson}_Prompt/h3Eff{meson}Prompt')
                heff_D_nonprompt = dau_eff_file.Get(f'{meson}_NonPrompt/h3Eff{meson}NonPrompt')
            else:
                pi_eff_file = TFile.Open(f'{path}/EfftimesAcc_pTEtaPhiMap_SoftPion_GenKine.root')
                print(f'\033[1mLoading pi eff data from {pi_eff_file.GetName()}\033[0m')
                d_eff_file = TFile.Open(f'{path}/EfftimesAcc_pTYPhiMap_Dzero.root')
                print(f'\033[1mLoading D eff data from {d_eff_file.GetName()}\033[0m')
                heff_v0 = pi_eff_file.Get('Pion_eff_pTEtaPhiMap')
                heff_D_prompt = d_eff_file.Get('Dzero_eff_pTYPhiMap')
                heff_D_nonprompt = d_eff_file.Get('Dzero_eff_pTYPhiMap')

            # ______________________________________________________________________________
            # Compute efficiency
            nbins = 30
            xbins = np.array(list(range(nbins+1)), dtype=np.float32)
            hgen_pt = TH1F('hgen_pt', '', nbins, xbins)
            hreco_pt = hgen_pt.Clone('hreco_pt')
            hreco_pt_low = hgen_pt.Clone('hreco_pt_low')
            hreco_pt_high = hgen_pt.Clone('hreco_pt_high')
            hreco_pt_np = hgen_pt.Clone('hreco_pt_np')
            hreco_pt_np_low = hgen_pt.Clone('hreco_pt_np_low')
            hreco_pt_np_high = hgen_pt.Clone('hreco_pt_np_high')

            SetObjectStyle(hreco_pt, color=kRed+1, markerstyle=20, markersize=1.,
                           linecolor=kRed+1, fillalpha=0.2, fillstyle=3345)
            SetObjectStyle(hreco_pt_low, color=kRed+1, markerstyle=20, markersize=1.,
                           linecolor=kRed+1, fillalpha=0.2, linestyle=2)
            SetObjectStyle(hreco_pt_high, color=kRed+1, markerstyle=20, markersize=1.,
                           linecolor=kRed+1, fillalpha=0.2, linestyle=8)
            SetObjectStyle(hreco_pt_np, color=kAzure+4, markerstyle=20, markersize=1.,
                           fillalpha=0.2, fillcolor=kAzure+4, fillstyle=3354)
            SetObjectStyle(hreco_pt_np_low, color=kAzure+4, markerstyle=20, markersize=1.,
                           linecolor=kAzure+4, fillalpha=0.2, linestyle=2)
            SetObjectStyle(hreco_pt_np_high, color=kAzure+4, markerstyle=20, markersize=1.,
                           linecolor=kAzure+4, fillalpha=0.2, linestyle=8)

            hkine_mot_v0 = TH2F('hkine_mot_v0',
                                f';#it{{p}}_{{T}}(V0) (GeV/#it{{c}});#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}})',
                                nbins*10, xbins[0], xbins[-1], nbins*10, xbins[0], xbins[-1])
            hkine_mot_D = TH2F('hkine_mot_D',
                               f';#it{{p}}_{{T}}(D) (GeV/#it{{c}});#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}})',
                                nbins*10, xbins[0], xbins[-1], nbins*10, xbins[0], xbins[-1])
            if reso == 413:
                hkine_mot_v0.GetXaxis().SetTitle('#it{p}_{T}(#pi) (GeV/#it{c})')
                hkine_mot_v0.GetYaxis().SetTitle('#it{p}_{T}(D*^{+}) (GeV/#it{c})')
                hkine_mot_D.GetXaxis().SetTitle('#it{p}_{T}(D) (GeV/#it{c})')
                hkine_mot_D.GetYaxis().SetTitle('#it{p}_{T}(D*^{+}) (GeV/#it{c})')

            with alive_bar(len(reso_df_kine), bar='smooth', spinner='classic', title='Computing reso efficiency') as time_bar:
                for pt_reso, y_reso, pt_D, y_D, phi_D, pt_v0, y_v0, phi_v0, eta_v0, angle_pis in zip(reso_df_kine['pt_reso'].to_numpy(),
                                                                                                     reso_df_kine['y_reso'].to_numpy(),
                                                                                                     reso_df_kine['pt_d'].to_numpy(),
                                                                                                     reso_df_kine['y_d'].to_numpy(),
                                                                                                     reso_df_kine['phi_d'].to_numpy(),
                                                                                                     reso_df_kine['pt_light'].to_numpy(),
                                                                                                     reso_df_kine['y_light'].to_numpy(),
                                                                                                     reso_df_kine['phi_light'].to_numpy(),
                                                                                                     reso_df_kine['eta_light'].to_numpy(),
                                                                                                     reso_df_kine['theta_pis'].to_numpy()):
                    hkine_mot_D.Fill(pt_D, pt_reso)
                    hkine_mot_v0.Fill(pt_v0, pt_reso)

                    bin_D = heff_D_prompt.FindBin(pt_D, y_D, phi_D)
                    eff_D = heff_D_prompt.GetBinContent(bin_D)
                    eff_D_err = heff_D_prompt.GetBinError(bin_D)
                    eff_Dnp = heff_D_nonprompt.GetBinContent(bin_D)
                    eff_Dnp_err = heff_D_nonprompt.GetBinError(bin_D)

                    if reso == 413:
                        # for D* we use eta instead of y
                        bin_v0 = heff_v0.FindBin(pt_v0, eta_v0, phi_v0)
                    else:
                        bin_v0 = heff_v0.FindBin(pt_v0, y_v0, phi_v0)
                    eff_v0 = heff_v0.GetBinContent(bin_v0)

                    rnd_D = gRandom.Rndm()
                    rnd_v0 = gRandom.Rndm()

                    # cut on generated mother |y| < 0.5
                    if abs(y_reso) < 0.5:
                        hgen_pt.Fill(pt_reso)

                    # for D* cut on the angle between pion and D0 to reproduce AOD filtering
                    if angle_pis > 0.5 and reso == 413:
                        time_bar()
                        continue
                    if rnd_v0 > eff_v0:
                        time_bar()
                        continue

                    if rnd_D < eff_D:
                        hreco_pt.Fill(pt_reso)
                    if rnd_D < eff_D - eff_D_err:
                        hreco_pt_low.Fill(pt_reso)
                    if rnd_D < eff_D + eff_D_err:
                        hreco_pt_high.Fill(pt_reso)

                    if rnd_D < eff_Dnp:
                        hreco_pt_np.Fill(pt_reso)
                    if rnd_D < eff_Dnp - eff_Dnp_err:
                        hreco_pt_np_low.Fill(pt_reso)
                    if rnd_D < eff_Dnp + eff_Dnp_err:
                        hreco_pt_np_high.Fill(pt_reso)

                    time_bar()

            heff = hreco_pt.Clone('heff_prompt')
            heff_low = hreco_pt_low.Clone('heff_prompt_low')
            heff_high = hreco_pt_high.Clone('heff_prompt_high')
            heff_np = hreco_pt_np.Clone('heff_nonprompt')
            heff_np_low = hreco_pt_np_low.Clone('heff_nonprompt_low')
            heff_np_high = hreco_pt_np_high.Clone('heff_nonprompt_high')
            heff.Divide(hreco_pt, hgen_pt, 1., 1., 'B')
            heff_low.Divide(hreco_pt_low, hgen_pt, 1., 1., 'B')
            heff_high.Divide(hreco_pt_high, hgen_pt, 1., 1., 'B')
            heff_np.Divide(hreco_pt_np, hgen_pt, 1., 1., 'B')
            heff_np_low.Divide(hreco_pt_np_low, hgen_pt, 1., 1., 'B')
            heff_np_high.Divide(hreco_pt_np_high, hgen_pt, 1., 1., 'B')

            canvas = TCanvas(f'canvas_eff_reso{reso}_{mult}_{ptshape}', '', 800, 800)
            canvas.cd().SetLogy()
            eff_title = f';#it{{p}}_{{T}} (GeV/#it{{c}}); Acc #times #epsilon ({reso_name}) {mult}'
            canvas.cd().DrawFrame(0., 1.e-4, 25., 10, eff_title)
            heff.Draw('e hist same')
            heff_np.Draw('e hist same')
            legend = TLegend(0.2, 0.83, 0.7, 0.93)
            legend.SetHeader(f'{decay_chain}')
            legend.SetNColumns(2)
            legend.AddEntry(heff, 'Prompt', 'pl')
            legend.AddEntry(heff_np, 'Non-prompt', 'pl')
            legend.Draw()

            canvas_kine = TCanvas(f'canvas_kine_reso{reso}_{mult}_{ptshape}', '', 1600, 800)
            canvas_kine.Divide(2, 1)
            for i, histo in enumerate([hkine_mot_D, hkine_mot_v0]):
                canvas_kine.cd(i+1).SetLogz()
                corr_title = f';#it{{p}}_{{T}} (GeV/#it{{c}}); #it{{p}}_{{T}} ({reso_name}) (GeV/#it{{c}})'
                canvas_kine.cd(i+1).DrawFrame(0., 0., 25., 25., corr_title)
                histo.Draw('colz')

            # ______________________________________________________________________________
            # Save output
            outfile = TFile.Open(f'{path}/eff_reso{reso}_{mult}_{ptshape}.root', 'recreate')
            hgen_pt.Write()
            hreco_pt.Write()
            hreco_pt_low.Write()
            hreco_pt_high.Write()
            hreco_pt_np.Write()
            hreco_pt_np_low.Write()
            hreco_pt_np_high.Write()
            heff.Write('heff_prompt')
            heff_low.Write('heff_prompt_low')
            heff_high.Write('heff_prompt_high')
            heff_np.Write('heff_nonprompt')
            heff_np_low.Write('heff_nonprompt_low')
            heff_np_high.Write('heff_nonprompt_high')
            hkine_mot_v0.Write()
            hkine_mot_D.Write()
            canvas.Write()
            canvas.SaveAs(f'{path}/eff_reso{reso}_{mult}_{ptshape}.pdf')
            canvas_kine.Write()
            canvas_kine.SaveAs(f'{path}/kine_reso{reso}_{mult}_{ptshape}.pdf')
            outfile.Close()
            input('Press enter to continue...')

    input('Press enter to exit...')

main()
