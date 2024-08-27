'''
python script to compute the prompt fraction in v2 analyses
'''
import sys
import argparse
import yaml
import numpy as np
from ROOT import TFile, TGraphAsymmErrors, TH1, TH1F, kOrange, kAzure, TCanvas, TLegend
sys.path.append('../../')
from utils.StyleFormatter import SetObjectStyle
from utils.ReadModel import ReadFONLL, ReadTAMU
from utils.AnalysisUtils import GetFractionNb, GetPromptFDFractionFc

def get_nbfraction(ry, eff_fd, xsc_pp, delta_pt, br, nev, raa, taa):
    '''
    Compute the prompt fraction in v2 analyses with Nb method

    Args:
    ry: float
        raw yield
    eff_fd: float
        feed-down efficiency
    xsc_pp: list
        cross section in pp interpolated from FONLL (Central, Min, Max)
    delta_pt: float
        pT bin width
    br: float
        branching ratio
    nev: float
        number of events
    raa: list
        nuclear modification factor (Central, Min, Max)
    taa: float
        nuclear overlap function
    '''

    # Compute prompt fraction
    #_____________________________________________________________________________________
    # prompt fraction
    fd_frac_central = taa * xsc_pp['Central'] * raa['Central'] * \
                      2 * delta_pt * eff_fd * br * nev / ry['Central']
    fd_frac = []
    for Xsec  in [xsc_pp['Min'], xsc_pp['Max']]:
        for Raa  in [raa['Min'], raa['Max']]:
            fd_frac.append(taa * Xsec * Raa * \
                           2 * delta_pt * eff_fd * \
                           br * nev / ry['Central'])
    fd_frac.sort()

    return [fd_frac[0], fd_frac_central, fd_frac[-1]]

def compute_nbfraction(configfile_name, efffile, outputdir, suffix):

    # Read configuration file
    with open(configfile_name, 'r', encoding='utf8') as ymlconfig:
        config = yaml.load(ymlconfig, yaml.FullLoader)

    # read global configuration
    # pt
    ptMins = config['ptmins']
    ptMaxs = config['ptmaxs']
    ptLims = list(ptMins)
    nPtBins = len(ptMins)
    ptLims.append(ptMaxs[-1])
    ptLims = np.asarray(ptLims, 'f')
    # fonll
    Dmeson = config['Dmeson']
    if Dmeson == 'D0': Dmeson = 'Dzero'
    fonllpred = {}
    fonllpred_prompt = {}
    fonllfile = TFile.Open(config['fonllfile'])
    for key in ['Central', 'Min', 'Max']:
        fonllpred[key] = fonllfile.Get(f'NonPrompt/{Dmeson}/hFonllNonPrompt{Dmeson}{key}')
        fonllpred_prompt[key] = fonllfile.Get(f'Prompt/Dplus/hFonllPromptDplus{key}')
        fonllpred[key].SetName(f"fonll_{key}")
        fonllpred[key].SetDirectory(0)
        fonllpred_prompt[key].SetName(f"fonll_prompt_{key}")
        fonllpred_prompt[key].SetDirectory(0)
    fonllfile.Close()
    # raa
    raa = {}
    if 'raa_tamu' in config:
        spline_raa, _, _, _ = ReadTAMU(config['raa_tamu'])
    else:
        print("\033[91mERROR: raa model not recognized. Exit.\033[0m")
        sys.exit(1)

    # Normalization
    # taa, Taa: 0.000000003917 #pb-1
    taa = 0.000000003917 #pb-1 (from David, https://alice-notes.web.cern.ch/system/files/notes/analysis/1541/2024-04-30-Centrality_Studies_2023%20%281%29.pdf)
    br = config['br']
    nevents = config['nevents']
    hist_norm = TH1F("hist_norm", "hist_norm", 3, 0, 3)
    hist_norm.GetXaxis().SetBinLabel(1, "taa")
    hist_norm.GetXaxis().SetBinLabel(2, "br")
    hist_norm.GetXaxis().SetBinLabel(3, "nevents")
    hist_norm.SetBinContent(1, taa)
    hist_norm.SetBinContent(2, br)
    hist_norm.SetBinContent(3, nevents)

    # Efficiency
    eff_prompt_hist = TFile.Open(efffile).Get('effPrompt')
    eff_fd_hist = TFile.Open(efffile).Get('effFD')

    # Raw yield
    ryield_filename = config['ryield_filename']
    ry_file = TFile.Open(ryield_filename)
    ry = ry_file.Get('hRawYieldsSimFit')
    ry.SetDirectory(0)
    ry_file.Close()

    gfd_frac, gprompt_frac = TGraphAsymmErrors(nPtBins), TGraphAsymmErrors(nPtBins)
    gprompt_frac_min = TGraphAsymmErrors(nPtBins)
    gprompt_frac_max = TGraphAsymmErrors(nPtBins)
    hprompt_frac = TH1F("hprompt_frac", "hprompt_frac", nPtBins, np.array(ptLims))
    hprompt_frac_min = TH1F("hprompt_frac_min", "hprompt_frac_min", nPtBins, np.array(ptLims))
    hprompt_frac_max = TH1F("hprompt_frac_max", "hprompt_frac_max", nPtBins, np.asarray(ptLims))
    hprompt_frac.SetDirectory(0)
    hprompt_frac_min.SetDirectory(0)
    hprompt_frac_max.SetDirectory(0)
    for ipt, (ptmin, ptmax) in enumerate(zip(ptMins, ptMaxs)):
        deltaPt = ptmax - ptmin
        ptcent = (ptmin + ptmax) / 2.
        print(f"ptmin: {ptmin}, ptmax: {ptmax}")

        # Collect inputs
        #_____________________________________________________________________________________
        # cross section in pp
        # fonll
        xsec = {}
        for pred in fonllpred:
            xsec[pred] = fonllpred[pred].GetBinContent(fonllpred[pred].FindBin(ptcent + 1.e-9))
        print("____________________________________")
        print(f"Cross section in pp: {xsec['Central']} - {xsec['Min']} + {xsec['Max']}")

        # efficiency
        eff_prompt = {'Central': eff_prompt_hist.GetBinContent(eff_prompt_hist.FindBin(ptcent + 1.e-9)),
                      'Min': eff_prompt_hist.GetBinError(eff_prompt_hist.FindBin(ptcent + 1.e-9)),
                      'Max': eff_prompt_hist.GetBinError(eff_prompt_hist.FindBin(ptcent + 1.e-9))}
        eff_fd = {'Central': eff_fd_hist.GetBinContent(eff_fd_hist.FindBin(ptcent + 1.e-9)),
                    'Min': eff_fd_hist.GetBinError(eff_fd_hist.FindBin(ptcent + 1.e-9)),
                    'Max': eff_fd_hist.GetBinError(eff_fd_hist.FindBin(ptcent + 1.e-9))}
        print("Efficiency")
        print(f"Prompt: {eff_prompt['Central']} - {eff_prompt['Min']} + {eff_prompt['Max']}")
        print(f"FD: {eff_fd['Central']} - {eff_fd['Min']} + {eff_fd['Max']}")

        # raw yield
        ryield = {'Central': abs(ry.GetBinContent(ry.FindBin(ptcent + 1.e-9))),
                  'Min': abs(ry.GetBinContent(ry.FindBin(ptcent + 1.e-9)) - ry.GetBinError(ry.FindBin(ptcent + 1.e-9))),
                  'Max': abs(ry.GetBinContent(ry.FindBin(ptcent + 1.e-9)) + ry.GetBinError(ry.FindBin(ptcent + 1.e-9)))}
        print(f"Raw yield: {ryield['Central']} - {ryield['Min']} + {ryield['Max']}")

        # raa
        raa['Central'] = spline_raa['yCent'](ptcent)
        if 'yMin' in spline_raa and 'yMax' in spline_raa:
            raa['Min'] = spline_raa['yMin'](ptcent)
            raa['Max'] = spline_raa['yMax'](ptcent)
        else:
            # raa +-50%
            raa['Min'] = raa['Central'] - 0.5*raa['Central']
            raa['Max'] = raa['Central'] + 0.5*raa['Central']
            raa['Min'] = raa['Central']
            raa['Max'] = raa['Central']
        print("____________________________________")
        print(f"Raa: {raa['Central']} - {raa['Min']} + {raa['Max']}")

        #_____________________________________________________________________________________
        # normalization
        print(f"Normalization: taa = {taa}, br = {br}, nevents = {nevents}")

        # Compute prompt fraction
        #_____________________________________________________________________________________
        # prompt fraction (0: Min, 1: Central, 2: Max)
        fd_frac = get_nbfraction(ryield,
                                 eff_fd['Central'],
                                 xsec,
                                 deltaPt, br,
                                 nevents, raa, taa)
        prompt_frac = [1 - fd_frac[2], 1 - fd_frac[1], 1 - fd_frac[0]]
        print(f"FD fraction: {fd_frac[1]:.3f} - {fd_frac[0]:.3f} + {fd_frac[2]:.3f}")
        print(f"Prompt fraction: {prompt_frac[1]:.3f} - {prompt_frac[0]:.3f} + {prompt_frac[2]:.3f}")

        # Fill TGraph
        gfd_frac.SetPoint(ipt, ptcent, fd_frac[1])
        gfd_frac.SetPointError(ipt, ptcent - ptmin,
                               ptmax - ptcent,
                               fd_frac[1] - fd_frac[0],
                               fd_frac[2] - fd_frac[1])
        gprompt_frac.SetPoint(ipt, ptcent, 1 - fd_frac[1])
        gprompt_frac.SetPointError(ipt, ptcent - ptmin,
                                   ptmax - ptcent,
                                   fd_frac[2] - fd_frac[1],
                                   fd_frac[1] - fd_frac[0])
        SetObjectStyle(gprompt_frac, markerstyle=20, markercolor=kOrange+1,
                       markersize=1.5, linecolor=kOrange+1)
        SetObjectStyle(gfd_frac, markerstyle=21, markercolor=kAzure+4,
                       markersize=1.5, linecolor=kAzure+4)

        hprompt_frac.SetBinContent(ipt + 1, prompt_frac[1])
        hprompt_frac.SetBinError(ipt + 1, max(prompt_frac[1] - prompt_frac[0], prompt_frac[2] - prompt_frac[1]))
        hprompt_frac_min.SetBinContent(ipt + 1, prompt_frac[0])
        hprompt_frac_min.SetBinError(ipt + 1, 0)
        hprompt_frac_max.SetBinContent(ipt + 1, prompt_frac[2])
        hprompt_frac_max.SetBinError(ipt + 1, 0)

    # Plot
    #_____________________________________________________________________________________
    gfd_frac.SetName("gfd_frac")
    gprompt_frac.SetName("gprompt_frac")
    gprompt_frac_min.SetName("gprompt_frac_min")
    gprompt_frac_max.SetName("gprompt_frac_max")
    canv = TCanvas("canv", "canv", 500, 500)
    hframe = canv.DrawFrame(0, 0, 50, 1.1, ";#it{p}_{T} (GeV/#it{c}); Fraction")
    leg = TLegend(0.15, 0.3, 0.5, 0.5)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(gprompt_frac, 'Prompt', 'p')
    leg.AddEntry(gfd_frac, 'Feed-down', 'p')
    gfd_frac.Draw("PZ same")
    gprompt_frac.Draw("PZ same")
    leg.Draw()
    canv.Update()

    input("Press Enter to continue...")

    # Save to file
    #_____________________________________________________________________________________
    outputfile_name = f"{outputdir}/promptfrac{suffix}.root"
    outfile = TFile(outputfile_name, "recreate")
    ry.Write()
    eff_prompt_hist.Write()
    eff_fd_hist.Write()
    hist_norm.Write()
    gfd_frac.Write()
    gprompt_frac.Write()
    gprompt_frac_min.Write()
    gprompt_frac_max.Write()
    hprompt_frac.Write()
    hprompt_frac_min.Write()
    hprompt_frac_max.Write()
    canv.Write()
    outfile.Close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument('configFileName', metavar='text', default='config_Ds_Fit.yml', help='configuration file')
    parser.add_argument('efffile', metavar='text', default='eff.root', help='efficiency file')
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    compute_nbfraction(
        args.configFileName,
        args.efffile,
        args.outputdir,
        args.suffix
    )
