'''
python script to compute the selection efficiency for resonances
'''

import argparse
import sys
import os
import numpy as np
import uproot
from alive_progress import alive_bar
from particle import Particle
from ROOT import TFile, TNtuple, TH2F, TH1F, TCanvas, kRed, kAzure, TLegend, gRandom, gROOT # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

SetGlobalStyle(padtopmargin=0.05, padrightmargin=0.15, padleftmargin=0.15, padbottommargin=0.15,
               palette=55, labelsize=0.03, titlesize=0.03, labeloffset=0.01,
               titleoffsety=1.8, titleoffsetx=1.8, titleoffsetz=1., opttitle=0)

def compute_efficiency(
    input_file,
    kinefile,
    pt_min_reso,
    pt_max_reso,
    input_file_pi,
    output_dir,
    suffix,
    seed
):
    """
    function for the propagation of the efficiency
    from the daughters to the D resonance mother

    Parameters
    ----------

    - input_file (str): input ROOT file with the 3D efficiency maps of the daugthers
    - kinefile (str): input ROOT file with D* -> D kinematics
    - pt_min_reso (float): min pT for pT integrated efficiency
    - pt_max_reso (float): max pT for pT integrated efficiency
    - input_file_pi (str): additional input ROOT file for pion in case of D*(2010)
    - output_dir (str): output directory
    - suffix (str): suffix for output file
    - seed (str): seed for random sampling
    """

    gROOT.SetBatch(True)
    gRandom.SetSeed(seed)

    if '435' in kinefile:
        reso = 435
        reso_label = Particle.from_pdgid(reso).name
        reso_name = 'Ds2starplus'
        meson = 'Dplus'
        decay_chain = f'{reso_label} #rightarrow D^{{+}} K^{{0}}_{{S}}'
    elif '10433' in kinefile:
        reso = 10433
        reso_label = Particle.from_pdgid(reso).name
        reso_name = 'Ds1plus'
        meson = 'Dstar'
        decay_chain = f'{reso_label} #rightarrow D^{{*+}} K^{{0}}_{{S}}'
    elif '413' in kinefile:
        reso = 413
        reso_label = Particle.from_pdgid(reso).name
        reso_name = 'Dstar'
        meson = 'Dzero'
        decay_chain = f'{reso_label} #rightarrow D^{{0}} #pi^{{+}}'
    else:
        print(f'\033[1m\033[91mResonance not supported\033[0m')
        sys.exit()
    print(f'\033[1m\033[92mComputing efficiency for {reso_label} ({reso})\033[0m')

    ptshape = "Dsptshape"
    if "Lcptshape" in kinefile:
        ptshape = "Lcptshape"
    elif "DsHarderptshape" in kinefile:
        ptshape = "DsHarderptshape"
    elif "DsVeryHarderptshape" in kinefile:
        ptshape = "DsVeryHarderptshape"
    elif "DsSofterptshape" in kinefile:
        ptshape = "DsSofterptshape"
    elif "DsVerySofterptshape" in kinefile:
        ptshape = "DsVerySofterptshape"

    trigger = ""
    if "MB" in kinefile:
        trigger = "MB"
    elif "HM" in kinefile:
        trigger = "HM"

    mult_weights_suffix = ""
    if "_multweights_all" in input_file:
        mult_weights_suffix = "_multweights_all"
    elif "_multweights_cand" in input_file:
        mult_weights_suffix = "_multweights_cand"
    elif "_multweights_candinmass" in input_file:
        mult_weights_suffix = "_multweights_candinmass"

    reso_label = Particle.from_pdgid(reso).name

    # Load data
    print(f'\033[1mLoading kine data from {kinefile}\033[0m')
    reso_df_kine = uproot.open(kinefile)['tupleDreso'].arrays(library='pd')

    # phi is defined in [-pi, pi]
    reso_df_kine['phi_reso'] = reso_df_kine['phi_reso'].apply(lambda x: x + np.pi)
    reso_df_kine['phi_d'] = reso_df_kine['phi_d'].apply(lambda x: x + np.pi)
    reso_df_kine['phi_light'] = reso_df_kine['phi_light'].apply(lambda x: x + np.pi)

    if reso != 413:
        dau_eff_file = TFile.Open(input_file)
        print(f'\033[1mLoading dau eff data from {dau_eff_file.GetName()}\033[0m')
        heff_v0 = dau_eff_file.Get(f'K0S_{meson}/hEff_{meson}_K0S')
        heff_D_prompt = dau_eff_file.Get(f'{meson}_Prompt/h3Eff{meson}Prompt')
        heff_D_nonprompt = dau_eff_file.Get(f'{meson}_NonPrompt/h3Eff{meson}NonPrompt')
        heff_v0.SetDirectory(0)
        heff_D_prompt.SetDirectory(0)
        heff_D_nonprompt.SetDirectory(0)
        dau_eff_file.Close()
    else:
        pi_eff_file = TFile.Open(input_file_pi)
        print(f'\033[1mLoading pi eff data from {pi_eff_file.GetName()}\033[0m')
        d_eff_file = TFile.Open(input_file)
        print(f'\033[1mLoading D eff data from {d_eff_file.GetName()}\033[0m')
        heff_v0 = pi_eff_file.Get('Pion_eff_pTEtaPhiMap')
        heff_D_prompt = d_eff_file.Get('Dzero_eff_pTYPhiMap')
        heff_D_nonprompt = d_eff_file.Get('Dzero_eff_pTYPhiMap')
        heff_v0.SetDirectory(0)
        heff_D_prompt.SetDirectory(0)
        heff_D_nonprompt.SetDirectory(0)
        pi_eff_file.Close()
        d_eff_file.Close()

    # ______________________________________________________________________________
    # Create output file
    outfile = TFile.Open(
        os.path.join(
            output_dir, f"eff_times_acc_{reso_name}_{trigger}_{ptshape}{mult_weights_suffix}_propagated{suffix}.root"),
        'recreate'
    )

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
                        f';#it{{p}}_{{T}}(V0) (GeV/#it{{c}});#it{{p}}_{{T}}({reso_label}) (GeV/#it{{c}})',
                        nbins*10, xbins[0], xbins[-1], nbins*10, xbins[0], xbins[-1])
    hkine_mot_D = TH2F('hkine_mot_D',
                       f';#it{{p}}_{{T}}(D) (GeV/#it{{c}});#it{{p}}_{{T}}({reso_label}) (GeV/#it{{c}})',
                        nbins*10, xbins[0], xbins[-1], nbins*10, xbins[0], xbins[-1])
    if reso == 413:
        hkine_mot_v0.GetXaxis().SetTitle('#it{p}_{T}(#pi) (GeV/#it{c})')
        hkine_mot_v0.GetYaxis().SetTitle('#it{p}_{T}(D*^{+}) (GeV/#it{c})')
        hkine_mot_D.GetXaxis().SetTitle('#it{p}_{T}(D) (GeV/#it{c})')
        hkine_mot_D.GetYaxis().SetTitle('#it{p}_{T}(D*^{+}) (GeV/#it{c})')

    # create TNtuple to store information of "selected" candidates
    reco_tree = TNtuple('recoTreePrompt', '', 'pt_reso:pt_d:pt_light')
    reco_tree_np = TNtuple('recoTreeNonPrompt', '', 'pt_reso:pt_d:pt_light')

    with alive_bar(len(reso_df_kine), bar='smooth', spinner='classic', title='Computing reso efficiency') as time_bar:
        for pt_reso, y_reso, pt_D, y_D, phi_D, pt_v0, y_v0, phi_v0, eta_v0, angle_pis in zip(
            reso_df_kine['pt_reso'].to_numpy(),
            reso_df_kine['y_reso'].to_numpy(),
            reso_df_kine['pt_d'].to_numpy(),
            reso_df_kine['y_d'].to_numpy(),
            reso_df_kine['phi_d'].to_numpy(),
            reso_df_kine['pt_light'].to_numpy(),
            reso_df_kine['y_light'].to_numpy(),
            reso_df_kine['phi_light'].to_numpy(),
            reso_df_kine['eta_light'].to_numpy(),
            reso_df_kine['theta_pis'].to_numpy()
        ):

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
                reco_tree.Fill(pt_reso, pt_D, pt_v0)
            if rnd_D < eff_D - eff_D_err:
                hreco_pt_low.Fill(pt_reso)
            if rnd_D < eff_D + eff_D_err:
                hreco_pt_high.Fill(pt_reso)

            if rnd_D < eff_Dnp:
                hreco_pt_np.Fill(pt_reso)
                reco_tree_np.Fill(pt_reso, pt_D, pt_v0)
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

    canvas = TCanvas(f'cEff_reso{reso}', '', 800, 800)
    canvas.cd().SetLogy()
    eff_title = f';#it{{p}}_{{T}} (GeV/#it{{c}}); Acc #times #epsilon ({reso_label}) {trigger}'
    canvas.cd().DrawFrame(0., 1.e-4, 25., 10, eff_title)
    heff.Draw('esame')
    heff_np.Draw('esame')
    legend = TLegend(0.2, 0.83, 0.7, 0.93)
    legend.SetHeader(f'{decay_chain}')
    legend.SetNColumns(2)
    legend.AddEntry(heff, 'Prompt', 'pl')
    legend.AddEntry(heff_np, 'Non-prompt', 'pl')
    legend.Draw()

    canvas_kine = TCanvas(f'canvas_kine_reso{reso}', '', 1600, 800)
    canvas_kine.Divide(2, 1)
    for i, histo in enumerate([hkine_mot_D, hkine_mot_v0]):
        canvas_kine.cd(i+1).SetLogz()
        corr_title = f';#it{{p}}_{{T}} (GeV/#it{{c}}); #it{{p}}_{{T}} ({reso_label}) (GeV/#it{{c}})'
        canvas_kine.cd(i+1).DrawFrame(0., 0., 25., 25., corr_title)
        histo.Draw('colz')

    # ______________________________________________________________________________
    # Compute integrated efficiency
    hgen_pt_int = TH1F('hgen_ptint', '', 1, np.array([pt_min_reso, pt_max_reso], dtype=np.float32))
    hreco_pt_int = hgen_pt_int.Clone('hreco_ptint')
    hreco_pt_int_np = hgen_pt_int.Clone('hreco_ptint')
    heff_pt_int = hgen_pt_int.Clone('heff_prompt_ptint')
    heff_pt_int_np = hgen_pt_int.Clone('heff_nonprompt_ptint')

    gen_int = hgen_pt.Integral(hgen_pt.FindBin(pt_min_reso*1.001),
                               hgen_pt.FindBin(pt_max_reso*0.999))
    reco_int = hreco_pt.Integral(hreco_pt.FindBin(pt_min_reso*1.001),
                                 hreco_pt.FindBin(pt_max_reso*0.999))
    reco_int_np = hreco_pt_np.Integral(hreco_pt_np.FindBin(pt_min_reso*1.001),
                                      hreco_pt_np.FindBin(pt_max_reso*0.999))
    hgen_pt_int.SetBinContent(1, gen_int)
    hreco_pt_int.SetBinContent(1, reco_int)
    hreco_pt_int_np.SetBinContent(1, reco_int_np)
    hgen_pt_int.Sumw2(False)
    hreco_pt_int.Sumw2(False)
    hreco_pt_int_np.Sumw2(False)

    heff_pt_int.Divide(hreco_pt_int, hgen_pt_int, 1., 1., 'B')
    heff_pt_int_np.Divide(hreco_pt_int_np, hgen_pt_int, 1., 1., 'B')
    SetObjectStyle(heff_pt_int, color=kRed+1, markerstyle=20, markersize=1.2,
                   linecolor=kRed+1, markercolor=kRed+1, linewidth=2)
    SetObjectStyle(heff_pt_int_np, color=kAzure+4, markerstyle=20, markersize=1.2,
                   linecolor=kAzure+4, markercolor=kAzure+4, linewidth=2)

    print((f'Prompt eff in {pt_min_reso:.0f}-{pt_max_reso:.0f}: '
           f'{heff_pt_int.GetBinContent(1):.3e} +/- {heff_pt_int.GetBinError(1):.3e}'))
    print((f'NonPrompt eff in {pt_min_reso:.0f}-{pt_max_reso:.0f}: '
          f'{heff_pt_int_np.GetBinContent(1):.3e} +/- {heff_pt_int_np.GetBinError(1):.3e}'))

    cEffTot = TCanvas(f'cEff_reso{reso}_Tot', '', 800, 800)
    cEffTot.cd().SetLogy()
    cEffTot.cd().DrawFrame(pt_min_reso, 1.e-4, pt_max_reso, 1, eff_title)
    heff_pt_int.Draw('same')
    heff_pt_int_np.Draw('same')
    legend.Draw()

    # ______________________________________________________________________________
    # Save output
    hgen_pt.Write()
    hreco_pt.Write()
    hreco_pt_low.Write()
    hreco_pt_high.Write()
    hreco_pt_np.Write()
    hreco_pt_np_low.Write()
    hreco_pt_np_high.Write()
    heff.Write()
    heff_low.Write()
    heff_high.Write()
    heff_np.Write()
    heff_np_low.Write()
    heff_np_high.Write()
    hkine_mot_v0.Write()
    hkine_mot_D.Write()
    hgen_pt_int.Write()
    hreco_pt_int.Write()
    hreco_pt_int_np.Write()
    heff_pt_int.Write()
    heff_pt_int_np.Write()
    canvas.Write()
    canvas_kine.Write()
    cEffTot.Write()
    reco_tree.Write()
    reco_tree_np.Write()

    canvas.SaveAs(
        os.path.join(
            output_dir,
            f'eff_times_acc_{reso_name}_{trigger}_{ptshape}{mult_weights_suffix}_propagated{suffix}.pdf'
        )
    )
    canvas_kine.SaveAs(
        os.path.join(
            output_dir,
            f'kine_reso_{reso_name}_{trigger}_{ptshape}{mult_weights_suffix}_propagated{suffix}.pdf'
        )
    )
    cEffTot.SaveAs(
        os.path.join(
            output_dir,
            f'eff_times_acc_ptint_{reso_name}_{trigger}_{ptshape}{mult_weights_suffix}_propagated{suffix}.pdf'
        )
    )

    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("input_file", metavar="text",
                        default="efficiencies.root", help="input ROOT file with 3D maps")
    parser.add_argument("kinefile", metavar="text",
                        default="kine.root", help="ROOT file with kinematics")
    parser.add_argument("--pt_min_reso", type=float,
                        default=2., help="min pT for pT integrated fraction")
    parser.add_argument("--pt_max_reso", type=float,
                        default=24., help="max pT for pT integrated fraction")
    parser.add_argument("--input_file_pi", "-i", metavar="text",
                        default="effpi.root", help="input file with 3D maps for pion")
    parser.add_argument("--output_dir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--seed", type=int,
                        default=42, help="seed for random sampling")
    args = parser.parse_args()

    compute_efficiency(
        args.input_file,
        args.kinefile,
        args.pt_min_reso,
        args.pt_max_reso,
        args.input_file_pi,
        args.output_dir,
        args.suffix,
        args.seed
    )
