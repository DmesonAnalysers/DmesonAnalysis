'''
python script for the propagation of the daughter tracking uncertainties to Ds resonances
run: python propagate_tracksyst.py
'''

import sys
import os
import itertools
import yaml
import numpy as np
import uproot
from particle import Particle
from ROOT import TFile, TCanvas, TH1F, TH2F, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import kRed, kOrange, kAzure, kCividis, kFullCircle, kFullSquare # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.AnalysisUtils import ApplyHistoEntriesToColumn #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle #pylint: disable=wrong-import-position,import-error

SetGlobalStyle(palette=kCividis, padleftmargin=0.14, padrightmargin=0.14, padbottommargin=0.14, titleoffsety=1.3)

resonances = [435, 10433]  # 413, 435, 10433
mults = ['MB', 'HM']
ptshapes = ['Dsptshape']  # 'Dsptshape' or 'Lcptshape'
infile_path = '~/Desktop/Analyses/Ds_resonances'

for reso in resonances:
    reso_name = Particle.from_pdgid(reso).name
    if reso == 435:
        reso_tag = 'Ds2starplus'
        meson = 'Dplus'
        decay_chain = f'{reso_name} #rightarrow D^{{+}} K^{{0}}_{{S}}'
        path = os.path.join(infile_path, 'eff')
    elif reso == 10433:
        reso_tag = 'Ds1plus'
        meson = 'Dstar'
        decay_chain = f'{reso_name} #rightarrow D^{{*+}} K^{{0}}_{{S}}'
        path = os.path.join(infile_path, 'eff')
    elif reso == 413:
        reso_tag = 'Dstar'
        meson = 'Dzero'
        decay_chain = f'{reso_name} #rightarrow D^{{0}} #pi^{{+}}'
        path = os.path.join(infile_path, 'check_efficiency_Dstar/compute_reso_eff')
    else:
        print(f'\033[1m\033[91mResonance {reso} not supported\033[0m')
        sys.exit()
    print(f'\033[1m\033[92mComputing efficiency for {reso_name} ({reso})\033[0m')

    # load input systematic uncertainties
    with open('config_tracksyst.yml', 'r', encoding='utf-8') as ymlCfgFile:
        inputTrackSyst = yaml.load(ymlCfgFile, yaml.FullLoader)

    ptMins = inputTrackSyst[reso]['D']['ptmins']
    ptMaxs = inputTrackSyst[reso]['D']['ptmaxs']
    systs = inputTrackSyst[reso]['D']['syst']
    ptLims = ptMins.copy()
    ptLims.append(ptMaxs[-1])
    hTrackSystD = TH1F('hTrSyst_D', ';#it{p}_{T} (GeV/#it{c});Rel. track. syst.', len(ptMins), np.array(ptLims, 'f'))
    for iPt, syst in enumerate(systs):
        hTrackSystD.SetBinContent(iPt+1, syst)
        hTrackSystD.SetBinError(iPt+1, 1.e-20)

    ptMins = inputTrackSyst[reso]['V0']['ptmins']
    ptMaxs = inputTrackSyst[reso]['V0']['ptmaxs']
    systs = inputTrackSyst[reso]['V0']['syst']
    systMats = inputTrackSyst[reso]['V0']['systMat']
    ptLims = ptMins.copy()
    ptLims.append(ptMaxs[-1])
    hTrackSystV0 = TH1F('hTrSyst_V0', ';#it{p}_{T} (GeV/#it{c});Rel. track. syst.', len(ptMins), np.array(ptLims, 'f'))
    for iPt, (syst, systMat) in enumerate(zip(systs, systMats)):
        hTrackSystV0.SetBinContent(iPt+1, np.sqrt(syst**2 + systMat**2))
        hTrackSystV0.SetBinError(iPt+1, 1.e-20)

    SetObjectStyle(hTrackSystD, color=kOrange+2, fillstyle=0, markersize=0.8)
    SetObjectStyle(hTrackSystV0, color=kOrange+7, fillstyle=0, markersize=0.8)

    for (mult, ptshape) in itertools.product(mults, ptshapes):
        # Load reco kine trees
        if reso == 413 and (mult == 'HM' or ptshape == 'Lcptshape'):
            print('\033[1m\033[91mERROR: D* efficiency for HM and LC not available\033[0m')
            sys.exit()

        data = f'{path}/eff_times_acc_{reso_tag}_{mult}_{ptshape}_multweights_candinmass_propagated.root'
        print(f'\033[1mLoading reco kine data from {data}\033[0m')
        df_kine = uproot.open(data)['recoTreePrompt'].arrays(library='pd')
        df_kine_np = uproot.open(data)['recoTreeNonPrompt'].arrays(library='pd')

        df_kine['track_syst_D'] = ApplyHistoEntriesToColumn(df_kine, 'pt_d', hTrackSystD)
        df_kine['track_syst_V0'] = ApplyHistoEntriesToColumn(df_kine, 'pt_light', hTrackSystV0)
        df_kine_np['track_syst_D'] = ApplyHistoEntriesToColumn(df_kine_np, 'pt_d', hTrackSystD)
        df_kine_np['track_syst_V0'] = ApplyHistoEntriesToColumn(df_kine_np, 'pt_light', hTrackSystV0)

        df_kine['track_syst_reso'] = df_kine.apply(lambda row: np.sqrt(row['track_syst_D']**2 + row['track_syst_V0']**2), axis=1)
        df_kine_np['track_syst_reso'] = df_kine_np.apply(lambda row: np.sqrt(row['track_syst_D']**2 + row['track_syst_V0']**2), axis=1)

        nPtBins = 24
        ptLims = np.array(list(range(nPtBins+1)), dtype=np.float32)

        hSystVsPtPrompt = TH2F('hSystVsPtPrompt', f';#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}});Rel. syst. unc.',
                               nPtBins, np.array(ptLims, 'd'), 18, 0., 0.18)
        hSystVsPtNonPrompt = TH2F('hSystVsPtNonPrompt', f';#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}});Rel. syst. unc.',
                                  nPtBins, np.array(ptLims, 'd'), 18, 0., 0.18)
        hPtDVsPtPrompt = TH2F('hPtDVsPtPrompt',
                              f';#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}});#it{{p}}_{{T}}(D) (GeV/#it{{c}})',
                              nPtBins, ptLims[0], ptLims[-1], nPtBins, ptLims[0], ptLims[-1])
        hPtV0VsPtPrompt = TH2F('hPtV0VsPtPrompt',
                               f';#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}});#it{{p}}_{{T}}(V0) (GeV/#it{{c}})',
                               nPtBins, ptLims[0], ptLims[-1], nPtBins, ptLims[0], ptLims[-1])
        hPtDVsPtNonPrompt = TH2F('hPtDVsPtNonPrompt',
                                 f';#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}});#it{{p}}_{{T}}(D) (GeV/#it{{c}})',
                                 nPtBins, ptLims[0], ptLims[-1], nPtBins, ptLims[0], ptLims[-1])
        hPtV0VsPtNonPrompt = TH2F('hPtV0VsPtNonPrompt',
                                  f';#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}});#it{{p}}_{{T}}(V0) (GeV/#it{{c}})',
                                  nPtBins, ptLims[0], ptLims[-1], nPtBins, ptLims[0], ptLims[-1])

        for pt_reso, pt_D, pt_V0, syst in zip(df_kine['pt_reso'].to_numpy(), df_kine['pt_d'].to_numpy(),
                                              df_kine['pt_light'].to_numpy(), df_kine['track_syst_reso'].to_numpy()):
            hSystVsPtPrompt.Fill(pt_reso, syst)
            hPtDVsPtPrompt.Fill(pt_reso, pt_D)
            hPtV0VsPtPrompt.Fill(pt_reso, pt_V0)
        for pt_reso, pt_D, pt_V0, syst in zip(df_kine_np['pt_reso'].to_numpy(), df_kine_np['pt_d'].to_numpy(),
                                              df_kine_np['pt_light'].to_numpy(),
                                              df_kine_np['track_syst_reso'].to_numpy()):
            hSystVsPtNonPrompt.Fill(pt_reso, syst)
            hPtDVsPtNonPrompt.Fill(pt_reso, pt_D)
            hPtV0VsPtNonPrompt.Fill(pt_reso, pt_V0)

        hSystMeanPrompt = hSystVsPtPrompt.ProfileX()
        hSystMeanNonPrompt = hSystVsPtNonPrompt.ProfileX()
        SetObjectStyle(hSystMeanPrompt, linewidht=3, fillstyle=0, markersize=0.8)
        SetObjectStyle(hSystMeanNonPrompt, linewidht=3, fillstyle=0, markersize=0.8)

        hPtDMeanPrompt = hPtDVsPtPrompt.ProfileX()
        hPtV0MeanPrompt = hPtV0VsPtPrompt.ProfileX()
        hPtDMeanNonPrompt = hPtDVsPtNonPrompt.ProfileX()
        hPtV0MeanNonPrompt = hPtV0VsPtNonPrompt.ProfileX()
        SetObjectStyle(hPtDMeanPrompt, fillstyle=0, markersize=0.5)
        SetObjectStyle(hPtV0MeanPrompt, fillstyle=0, markersize=0.5)
        SetObjectStyle(hPtDMeanNonPrompt, fillstyle=0, markersize=0.5)
        SetObjectStyle(hPtV0MeanNonPrompt, fillstyle=0, markersize=0.5)


        hTotSystPrompt = TH1F('hTotSystPrompt', '', 1, 2., 24.)
        hTotSystNonPrompt = TH1F('hTotSystNonPrompt', '', 1, 2., 24.)
        hTotSystPrompt.SetBinContent(1, hSystVsPtPrompt.ProfileY("_pfy", 3, 24).GetMean())
        hTotSystNonPrompt.SetBinContent(1, hSystVsPtNonPrompt.ProfileY("_pfy", 3, 24).GetMean())
        hTotSystPrompt.SetBinError(1, 1.e-20)
        hTotSystNonPrompt.SetBinError(1, 1.e-20)
        SetObjectStyle(hTotSystPrompt, color=kRed+1, linewidht=3, fillstyle=0, markerstyle=kFullCircle, markersize=1)
        SetObjectStyle(hTotSystNonPrompt, color=kAzure+4, linewidht=3, fillstyle=0, markerstyle=kFullSquare,
                       markersize=1)

        legAverage = TLegend(0.15, 0.75, 0.7, 0.9)
        legAverage.SetTextSize(0.045)
        legAverage.SetFillStyle(0)
        legAverage.AddEntry(hSystMeanPrompt, 'Average', 'pl')
        legAverage.AddEntry(hTotSystPrompt, 'p_{T} integrated', 'pl')

        legAverageNonPrompt = TLegend(0.15, 0.75, 0.7, 0.9)
        legAverageNonPrompt.SetTextSize(0.045)
        legAverageNonPrompt.SetFillStyle(0)
        legAverageNonPrompt.AddEntry(hSystMeanNonPrompt, 'Average', 'pl')
        legAverageNonPrompt.AddEntry(hTotSystNonPrompt, 'p_{T} integrated', 'pl')

        leg = TLegend(0.15, 0.75, 0.7, 0.9)
        leg.SetTextSize(0.035)
        leg.SetFillStyle(0)
        leg.AddEntry(hSystMeanPrompt, 'Resonance track. syst.', 'pl')
        leg.AddEntry(hTotSystPrompt, 'p_{T} integrated', 'pl')
        leg.AddEntry(hTrackSystD, 'D-meson track. syst.', 'pl')
        leg.AddEntry(hTrackSystV0, 'V0 track. syst. (w. mat. budget)', 'pl')

        cPrompt = TCanvas('cPrompt', 'Prompt', 1800, 600)
        cPrompt.Divide(3, 1)
        cPrompt.cd(1).SetLogz()
        hSystVsPtPrompt.Draw('colz')
        hSystMeanPrompt.DrawCopy('PH][ E0 same')
        hTotSystPrompt.DrawCopy('PH][ E0 same')
        legAverage.Draw()
        cPrompt.cd(2).SetLogz()
        hPtDVsPtPrompt.Draw('colz')
        hPtDMeanPrompt.DrawCopy('PH][ E0 same')
        cPrompt.cd(3).SetLogz()
        hPtV0VsPtPrompt.Draw('colz')
        hPtV0MeanPrompt.DrawCopy('PH][ E0 same')
        cPrompt.Modified()
        cPrompt.Update()

        cNonPrompt = TCanvas('cNonPrompt', 'NonPrompt', 1800, 600)
        cNonPrompt.Divide(3, 1)
        cNonPrompt.cd(1).SetLogz()
        hSystVsPtNonPrompt.Draw('colz')
        hSystMeanNonPrompt.DrawCopy('PH][ E0 same')
        hTotSystNonPrompt.DrawCopy('PH][ E0 same')
        legAverageNonPrompt.Draw()
        cNonPrompt.cd(2).SetLogz()
        hPtDVsPtNonPrompt.Draw('colz')
        hPtDMeanNonPrompt.DrawCopy('PH][ E0 same')
        cNonPrompt.cd(3).SetLogz()
        hPtV0VsPtNonPrompt.Draw('colz')
        hPtV0MeanNonPrompt.DrawCopy('PH][ E0 same')
        cNonPrompt.Modified()
        cNonPrompt.Update()

        cFinalSyst = TCanvas('cFinalSyst', '', 600, 600)
        cFinalSyst.DrawFrame(0., 0., 24., 0.18, f';#it{{p}}_{{T}}({reso_name}) (GeV/#it{{c}}); Rel. syst. unc.')
        hTotSystPrompt.DrawCopy('PH][ E0 same')
        hSystMeanPrompt.DrawCopy('PH][ E0 same')
        hTrackSystD.DrawCopy('PH][ E0 same')
        hTrackSystV0.DrawCopy('PH][ E0 same')
        leg.Draw()
        cFinalSyst.Modified()
        cFinalSyst.Update()

        outFile = TFile.Open(f'{infile_path}/systematics/tracking/tracking_syst_{reso}_{mult}_{ptshape}.root',
                             'recreate')
        hTrackSystD.Write()
        hTrackSystV0.Write()
        hSystVsPtPrompt.Write()
        hSystMeanPrompt.Write()
        hTotSystPrompt.Write()
        hPtDVsPtPrompt.Write()
        hPtV0VsPtPrompt.Write()
        hSystVsPtNonPrompt.Write()
        hSystMeanNonPrompt.Write()
        hTotSystNonPrompt.Write()
        hPtDVsPtNonPrompt.Write()
        hPtV0VsPtNonPrompt.Write()
        cPrompt.Write()
        cNonPrompt.Write()
        cFinalSyst.Write()
        outFile.Close()

        outFileNamePDF = f'{infile_path}/systematics/tracking/tracking_syst_{reso}_{mult}_{ptshape}.pdf'
        cPrompt.SaveAs(f'{outFileNamePDF}_Prompt.pdf')
        cNonPrompt.SaveAs(f'{outFileNamePDF}_NonPrompt.pdf')
        cFinalSyst.SaveAs(f'{outFileNamePDF}_FinalSyst.pdf')

        print('\nTotal tracking systematic uncertainty:')
        print('\nPrompt')
        print(f'\t\t pT (GeV/c): [2 - 24] syst: {hTotSystPrompt.GetBinContent(1):.3f}')
        print('\nFD')
        print(f'\t\t pT (GeV/c): [2 - 24] syst: {hTotSystNonPrompt.GetBinContent(1):.3f}')

        input('\nPress enter to exit')
