import ROOT
import yaml
import argparse
import sys
import ctypes
from flow_analysis_utils import get_centrality_bins
sys.path.append('../../')
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle
from utils.AnalysisUtils import ComputeEfficiency

SetGlobalStyle(titleoffsety=1.1, maxdigits=3, topmargin=0.1,
               bottommargin=0.4, leftmargin=0.3, rightmargin=0.15,
               labelsizey=0.04, setoptstat=0, setopttitle=0,
               setdecimals=True, titleoffsetx=0.91, titlesizex=0.05)

def compute_eff(config_file, centclass, outputdir, suffix):
    '''
    Compute efficiency for prompt and feed-down D mesons

    Args:
    config_file: str
        configuration file
    centclass: str
        centrality class
    outputdir: str
        output directory
    suffix: str
        suffix for output files

    Returns:
    output ROOT file with efficiency histograms
    '''

    #_____________________________________________________________________________________
    # Read configuration file
    with open(config_file, 'r', encoding='utf8') as ymlconfig:
        config = yaml.load(ymlconfig, yaml.FullLoader)

    #_____________________________________________________________________________________
    # Load input files
    eff_filename = config['task_filename']
    charm_hadron = config['Dmeson']
    infile = ROOT.TFile.Open(eff_filename)

    #_____________________________________________________________________________________
    # Check centrality selection
    if charm_hadron in ['Dplus', 'Ds']:
        cent_hist = infile.Get('hf-candidate-creator-3prong/hSelCollisionsCent')
        _, centMinMax = get_centrality_bins(centclass)
        for ipt in range(cent_hist.GetNbinsX()):
            if not (cent_hist.GetBinContent(centMinMax[0] + 1) > 0 \
                and cent_hist.GetBinContent(centMinMax[0]) == 0 \
                and cent_hist.GetBinContent(centMinMax[1]) > 0 \
                and cent_hist.GetBinContent(centMinMax[1] + 1) == 0):
                print(f'\033[91mFATAL: Invalid centrality class: {centclass}. Exit!\033[0m')
                sys.exit(1)
    elif charm_hadron == 'Dzero':
        cent_hist = infile.Get('hf-candidate-creator-2prong/hSelCollisionsCent')
        _, centMinMax = get_centrality_bins(centclass)
        for ipt in range(cent_hist.GetNbinsX()):
            if not (cent_hist.GetBinContent(centMinMax[0] + 1) > 0 \
                    and cent_hist.GetBinContent(centMinMax[0]) == 0 \
                    and cent_hist.GetBinContent(centMinMax[1]) > 0 \
                    and cent_hist.GetBinContent(centMinMax[1] + 1) == 0):
                print(f'\033[91mFATAL: Invalid centrality class: {centclass}. Exit!\033[0m')
                sys.exit(1)
    else:
        print('\033[93mWARNING: Invalid charm hadron for centrality check.\033[0m')

    if charm_hadron == 'Dplus':
        hgen_prompt = infile.Get('hf-task-dplus/hPtGenPrompt')
        hgen_fd = infile.Get('hf-task-dplus/hPtGenNonPrompt')
        hrec_prompt = infile.Get('hf-task-dplus/hPtRecSigPrompt')
        hrec_fd = infile.Get('hf-task-dplus/hPtRecSigNonPrompt')
    elif charm_hadron == 'Ds':
        hgen_prompt = infile.Get('hf-task-ds/MC/Ds/Prompt/hPtGen')
        hgen_fd = infile.Get('hf-task-ds/MC/Ds/NonPrompt/hPtGen')
        hrec_prompt = infile.Get('hf-task-ds/MC/Ds/Prompt/hPtRecSig')
        hrec_fd = infile.Get('hf-task-ds/MC/Ds/NonPrompt/hPtRecSig')
    #TODO: Add Dzero and Lc
    else:
        sys.exit(f'\033[91mFATAL: Invalid charm_hadron: {charm_hadron}. Exit!\033[0m')
    infile.Close()

    #_____________________________________________________________________________________
    # Compute efficiency
    heff_prompt = hgen_prompt.Clone('effPrompt')
    heff_prompt.GetXaxis().SetTitle('#it{p}_{T} (GeV/#it{c})')
    heff_prompt.GetYaxis().SetTitle('Acc.#times#varepsilon') # unlike run 2, in run 3 we directly have the Eff.#times Acc. from the HF task
    heff_prompt.Reset()
    heff_prompt.SetDirectory(0)
    heff_fd = hgen_fd.Clone('effFD')
    heff_fd.Reset()
    heff_fd.SetDirectory(0)
    for ipt in range(hgen_prompt.GetNbinsX()):
        # get unweighted yields (for uncertainty)
        nRecoPromptUnc, nGenPromptUnc, nRecoFDUnc, nGenFDUnc = (ctypes.c_double() for _ in range(4))
        nRecoPrompt = hrec_prompt.IntegralAndError(ipt, ipt+1, nRecoPromptUnc)
        nGenPrompt = hgen_prompt.IntegralAndError(ipt, ipt+1, nGenPromptUnc)
        nRecoFD = hrec_fd.IntegralAndError(ipt, ipt+1, nRecoFDUnc)
        nGenFD = hgen_fd.IntegralAndError(ipt, ipt+1, nGenFDUnc)

        effPrompt, effPromptUnc = ComputeEfficiency(nRecoPrompt, nGenPrompt,
                                                    nRecoPromptUnc.value, nGenPromptUnc.value)
        effFD, effFDUnc = ComputeEfficiency(nRecoFD, nGenFD, nRecoFDUnc.value, nGenFDUnc.value)
        heff_prompt.SetBinContent(ipt+1, effPrompt)
        heff_prompt.SetBinError(ipt+1, effPromptUnc)
        heff_fd.SetBinContent(ipt+1, effFD)
        heff_fd.SetBinError(ipt+1, effFDUnc)

    #_____________________________________________________________________________________
    # Draw histograms
    c = ROOT.TCanvas('c', 'c', 600, 600)
    c.cd().SetLogy()
    SetObjectStyle(heff_prompt, color=ROOT.kOrange+1, markerstyle=20, fillalpha=0.)
    SetObjectStyle(heff_fd, color=ROOT.kAzure+4, markerstyle=22, markersize=1.5,
                   linewidh=2, linestyle=7, fillalpha=0.)
    heff_prompt.GetYaxis().SetRangeUser(1.e-4, 1.1)
    heff_prompt.Draw('same')
    heff_fd.Draw('same')
    leg = ROOT.TLegend(0.15, 0.7, 0.4, 0.85)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextSize(0.03)
    leg.AddEntry(heff_prompt, 'Prompt', 'lp')
    leg.AddEntry(heff_fd, 'Feed-down', 'lp')
    leg.Draw()

    #_____________________________________________________________________________________
    # Save output
    outfile_name = f'{outputdir}/eff{suffix}.root'
    outfile = ROOT.TFile(outfile_name, 'RECREATE')
    outfile.cd()
    c.Write()
    heff_prompt.Write()
    heff_fd.Write()
    hgen_prompt.Write()
    hgen_fd.Write()
    hrec_prompt.Write()
    hrec_fd.Write()
    outfile.Close()
    input(f'Saving efficiency histograms to {outfile_name}. Press any key to exit.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument('--centclass', '-c', metavar='text', default='')
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    compute_eff(
        config_file=args.config,
        centclass=args.centclass,
        outputdir=args.outputdir,
        suffix=args.suffix)
