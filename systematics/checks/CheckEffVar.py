'''
python script to check the efficiency of a variable cut between two different "configuration"
(e.g. Pythia6 vs Pythia8, Improver vs No Improver)
run: python CheckEffVar.py cfgFileName.yml outFileDir
'''

import sys
import argparse
import yaml
import numpy as np
from root_numpy import fill_hist
from ROOT import TH1F, kRed, kAzure, kFullCircle, TCanvas, TLegend # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
#pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.DfUtils import LoadDfFromRootOrParquet
from utils.StyleFormatter import SetGlobalStyle, SetObjectStyle

def main():
    """
    Main function of the script
    """
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml', help='config file name')
    parser.add_argument('outFileDir', metavar='text', default='./', help='output file directory')
    args = parser.parse_args()

    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    #Load data
    dfPromptP6 = LoadDfFromRootOrParquet(inputCfg['input']['prompt_files'][0])
    dfPromptP8 = LoadDfFromRootOrParquet(inputCfg['input']['prompt_files'][1])
    dfFDP6 = LoadDfFromRootOrParquet(inputCfg['input']['fd_files'][0])
    dfFDP8 = LoadDfFromRootOrParquet(inputCfg['input']['fd_files'][1])

    #Select pt bin
    ptMin = inputCfg['pt_bin'][0]
    ptMax = inputCfg['pt_bin'][1]
    dfPromptP6 = dfPromptP6.query(f'{ptMin} < pt_cand < {ptMax}')
    dfPromptP8 = dfPromptP8.query(f'{ptMin} < pt_cand < {ptMax}')
    dfFDP6 = dfFDP6.query(f'{ptMin} < pt_cand < {ptMax}')
    dfFDP8 = dfFDP8.query(f'{ptMin} < pt_cand < {ptMax}')

    SetGlobalStyle(padbottommargin=0.14, padleftmargin=0.18, padrightmargin=0.06, titleoffsety=1.6)
    varTitle = inputCfg['scan_variable']['title']
    nBins = inputCfg['scan_variable']['histo_bins']
    binLims = inputCfg['scan_variable']['histo_lims']
    varName = inputCfg['scan_variable']['name']
    hPromptP6 = TH1F('hPromptP6', f';{varTitle};Counts', nBins, binLims[0], binLims[1])
    hPromptP8 = TH1F('hPromptP8', f';{varTitle};Counts', nBins, binLims[0], binLims[1])
    hFDP6 = TH1F('hFDP6', f';{varTitle};Counts', nBins, binLims[0], binLims[1])
    hFDP8 = TH1F('hFDP8', f';{varTitle};Counts', nBins, binLims[0], binLims[1])
    scaleFactor = inputCfg['scan_variable']['rescale_factor']
    dfPromptP6[varName] = dfPromptP6[varName] * scaleFactor
    dfPromptP8[varName] = dfPromptP8[varName] * scaleFactor
    dfFDP6[varName] = dfFDP6[varName] * scaleFactor
    dfFDP8[varName] = dfFDP8[varName] * scaleFactor
    fill_hist(hPromptP6, dfPromptP6[varName].values)
    fill_hist(hPromptP8, dfPromptP8[varName].values)
    fill_hist(hFDP6, dfFDP6[varName].values)
    fill_hist(hFDP8, dfFDP8[varName].values)
    SetObjectStyle(hPromptP6, color=kAzure+4, marker=kFullCircle)
    SetObjectStyle(hPromptP8, color=kRed+1, marker=kFullCircle)
    SetObjectStyle(hFDP6, color=kAzure+4, marker=kFullCircle)
    SetObjectStyle(hFDP8, color=kRed+1, marker=kFullCircle)
    hPromptP6.GetXaxis().SetNdivisions(505)
    hFDP6.GetXaxis().SetNdivisions(505)
    hPromptP8.GetXaxis().SetNdivisions(505)
    hFDP8.GetXaxis().SetNdivisions(505)

    scanRange = inputCfg['scan_variable']['scan_range']
    scanStep = inputCfg['scan_variable']['scan_step']
    nEffBins = round((scanRange[1] - scanRange[0]) / scanStep)
    hEffPromptP6 = TH1F('hEffPromptP6', f';{varTitle} >;Efficiency', nEffBins, scanRange[0], scanRange[1])
    hEffPromptP8 = TH1F('hEffPromptP8', f';{varTitle} >;Efficiency', nEffBins, scanRange[0], scanRange[1])
    hEffFDP6 = TH1F('hEffFDP6', f';{varTitle} >;Efficiency', nEffBins, scanRange[0], scanRange[1])
    hEffFDP8 = TH1F('hEffFDP8', f';{varTitle} >;Efficiency', nEffBins, scanRange[0], scanRange[1])
    SetObjectStyle(hEffPromptP6, color=kAzure+4, marker=kFullCircle)
    SetObjectStyle(hEffPromptP8, color=kRed+1, marker=kFullCircle)
    SetObjectStyle(hEffFDP6, color=kAzure+4, marker=kFullCircle)
    SetObjectStyle(hEffFDP8, color=kRed+1, marker=kFullCircle)

    effPromptP6, effPromptP8, effFDP6, effFDP8 = ([] for _ in range(4))
    labelsConf = inputCfg['legend']['conf_labels']
    legPrompt = TLegend(0.25, 0.2, 0.6, 0.4)
    legPrompt.SetBorderSize(0)
    legPrompt.SetFillStyle(0)
    legPrompt.SetHeader('Prompt')
    legPrompt.AddEntry(hEffPromptP6, labelsConf[0], 'p')
    legPrompt.AddEntry(hEffPromptP8, labelsConf[1], 'p')
    legFD = TLegend(0.25, 0.2, 0.65, 0.4)
    legFD.SetBorderSize(0)
    legFD.SetFillStyle(0)
    legFD.SetHeader('Non-prompt')
    legFD.AddEntry(hEffPromptP6, labelsConf[0], 'p')
    legFD.AddEntry(hEffPromptP8, labelsConf[1], 'p')

    for iBin, cut in enumerate(np.arange(scanRange[0], scanRange[1], scanStep)):
        dfPromptP6Sel = dfPromptP6.query(f'{varName} > {cut}')
        dfPromptP8Sel = dfPromptP8.query(f'{varName} > {cut}')
        dfFDP6Sel = dfFDP6.query(f'{varName} > {cut}')
        dfFDP8Sel = dfFDP8.query(f'{varName} > {cut}')

        effPromptP6.append(float(len(dfPromptP6Sel) / len(dfPromptP6)))
        effPromptP8.append(float(len(dfPromptP8Sel) / len(dfPromptP8)))
        effFDP6.append(float(len(dfFDP6Sel) / len(dfFDP6)))
        effFDP8.append(float(len(dfFDP8Sel) / len(dfFDP8)))

        hEffPromptP6.SetBinContent(iBin+1, effPromptP6[-1])
        hEffPromptP8.SetBinContent(iBin+1, effPromptP8[-1])
        hEffFDP6.SetBinContent(iBin+1, effFDP6[-1])
        hEffFDP8.SetBinContent(iBin+1, effFDP8[-1])

        hEffPromptP6.SetBinError(iBin+1, 1.e-20)
        hEffPromptP8.SetBinError(iBin+1, 1.e-20)
        hEffFDP6.SetBinError(iBin+1, 1.e-20)
        hEffFDP8.SetBinError(iBin+1, 1.e-20)

    hEffPromptP6.GetXaxis().SetNdivisions(505)
    hEffFDP6.GetXaxis().SetNdivisions(505)
    hEffPromptP8.GetXaxis().SetNdivisions(505)
    hEffFDP8.GetXaxis().SetNdivisions(505)

    hEffPromptRatio = hEffPromptP8.Clone('hEffPromptRatio')
    hEffPromptRatio.Divide(hEffPromptP6)
    hEffPromptRatio.GetYaxis().SetTitle(f'Prompt eff ratio {labelsConf[1]} / {labelsConf[0]}')
    hEffFDRatio = hEffFDP8.Clone('hEffFDRatio')
    hEffFDRatio.Divide(hEffFDP6)
    hEffFDRatio.GetYaxis().SetTitle(f'Non-prompt eff ratio {labelsConf[1]} / {labelsConf[0]}')

    hEffPromptRatio.GetXaxis().SetNdivisions(505)
    hEffFDRatio.GetXaxis().SetNdivisions(505)

    cDistributions = TCanvas('cDistributions', '', 1920, 1080)
    cDistributions.Divide(2, 1)
    cDistributions.cd(1).SetLogy()
    hPromptP8.Draw('e')
    hPromptP6.Draw('esame')
    legPrompt.Draw()
    cDistributions.cd(2).SetLogy()
    hFDP8.Draw('e')
    hFDP6.Draw('esame')
    legFD.Draw()

    cEfficiency = TCanvas('cEfficiency', '', 1920, 1080)
    cEfficiency.Divide(2, 1)
    cEfficiency.cd(1).SetLogy()
    hEffPromptP6.Draw('e')
    hEffPromptP8.Draw('esame')
    legPrompt.Draw()
    cEfficiency.cd(2).SetLogy()
    hEffFDP6.Draw('e')
    hEffFDP8.Draw('esame')
    legFD.Draw()

    cEfficiencyRatio = TCanvas('cEfficiencyRatio', '', 1920, 1080)
    cEfficiencyRatio.Divide(2, 1)
    cEfficiencyRatio.cd(1)
    hEffPromptRatio.Draw('e')
    cEfficiencyRatio.cd(2)
    hEffFDRatio.Draw('e')

    tag = f'{labelsConf[0]}Vs{labelsConf[1]}_pT{ptMin}_{ptMax}'
    cDistributions.SaveAs(f'{args.outFileDir}/{varName}_Distr_{tag}.pdf')
    cEfficiency.SaveAs(f'{args.outFileDir}/{varName}_CutEff_{tag}.pdf')
    cEfficiencyRatio.SaveAs(f'{args.outFileDir}/{varName}_CutEffRatio_{tag}.pdf')

    print('Press any key to exit!')
    input()

main()
