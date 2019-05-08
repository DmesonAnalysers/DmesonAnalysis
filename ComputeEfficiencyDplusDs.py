#*************************************************************************************************************#
# python script for the computation of the D+ and Ds+ efficiency mesons from ProjectDplusDsSparse.py output   #
# run: python ComputeEfficiencyDplusDs.py fitConfigFileName.yml inputFileName.root outFileName.root              #
# author: Fabrizio Grosa, fabrizio.grosa@to.infn.it ,INFN Torino                                              #
#*************************************************************************************************************#

from ROOT import TFile, TCanvas, TH1F, TLegend # pylint: disable=import-error,no-name-in-module
from ROOT import gROOT, gStyle, kRed, kBlue, kFullCircle, kOpenSquare # pylint: disable=import-error,no-name-in-module
import yaml, array, math, string, argparse, six

def ComputeEfficiency(recoCounts, genCounts, recoCountsError, genCountsError):
    hTmpNum = TH1F('hTmpNum', '', 1,0,1)
    hTmpDen = TH1F('hTmpDen', '', 1,0,1)
    hTmpNum.SetBinContent(1, recoCounts)
    hTmpDen.SetBinContent(1, genCounts)
    hTmpNum.SetBinError(1, recoCountsError)
    hTmpDen.SetBinError(1, genCountsError)
    hTmpNum.Divide(hTmpNum,hTmpDen,1.,1,'B')
    return hTmpNum.GetBinContent(1), hTmpNum.GetBinError(1)

def SetHistoStyle(histo, color, marker, markersize=1.5, linewidth=2, linestyle=1):
    histo.SetMarkerColor(color)
    histo.SetLineColor(color)
    histo.SetLineStyle(linestyle)
    histo.SetLineWidth(linewidth)
    histo.SetMarkerStyle(marker)
    histo.SetMarkerSize(markersize)

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('fitConfigFileName', metavar='text', default='config_Ds_Fit.yml')
parser.add_argument('centClass', metavar='text', default='')
parser.add_argument('inFileName', metavar='text', default='')
parser.add_argument('outFileName', metavar='text', default='')
parser.add_argument('--ptweights', metavar=('text','text'), nargs=2, required=False, help='First path of the pT weights file, second name of the pT weights histogram')
parser.add_argument("--batch", help="suppress video output", action="store_true")
args = parser.parse_args()

with open(args.fitConfigFileName, 'r') as ymlfitConfigFile:
    fitConfig = yaml.load(ymlfitConfigFile)

cent = ''
if args.centClass == 'k010':
    cent = 'Cent010'
elif args.centClass == 'k3050':
    cent = 'Cent3050'

gROOT.SetBatch(args.batch)

PtMin = fitConfig[cent]['PtMin']
PtMax = fitConfig[cent]['PtMax']
PtLims = list(PtMin)
nPtBins = len(PtMin)
PtLims.append(PtMax[nPtBins - 1])

hEffPrompt = TH1F('hEffPrompt',';#it{p}_{T} (GeV/#it{c});Efficiency',nPtBins,array.array('f',PtLims))
hEffFD = TH1F('hEffFD',';#it{p}_{T} (GeV/#it{c});Efficiency',nPtBins,array.array('f',PtLims))
hYieldPromptGen = TH1F('hYieldPromptGen',';#it{p}_{T} (GeV/#it{c}); # Generated MC',nPtBins,array.array('f',PtLims))
hYieldFDGen = TH1F('hYieldFDGen',';#it{p}_{T} (GeV/#it{c}); # Generated MC',nPtBins,array.array('f',PtLims))
hYieldPromptReco = TH1F('hYieldPromptReco',';#it{p}_{T} (GeV/#it{c}); # Reco MC',nPtBins,array.array('f',PtLims))
hYieldFDReco = TH1F('hYieldFDReco',';#it{p}_{T} (GeV/#it{c}); # Reco MC',nPtBins,array.array('f',PtLims))
SetHistoStyle(hEffPrompt,kRed,kFullCircle)
SetHistoStyle(hEffFD,kBlue,kOpenSquare,1.5,2,7)
SetHistoStyle(hYieldPromptGen,kRed,kFullCircle)
SetHistoStyle(hYieldFDGen,kBlue,kOpenSquare,1.5,2,7)
SetHistoStyle(hYieldPromptReco,kRed,kFullCircle)
SetHistoStyle(hYieldFDReco,kBlue,kOpenSquare,1.5,2,7)

if args.ptweights:
    infileWeigts = TFile.Open(args.ptweights[0])
    hPtWeights = infileWeigts.Get(args.ptweights[1])

hRecoPrompt, hRecoFD, hGenPrompt, hGenFD = ([] for iHisto in range(4))

infile = TFile(args.inFileName)
for iPt in range(len(PtMin)):
    hRecoPrompt.append(infile.Get('hPromptPt_%0.f_%0.f' % (PtMin[iPt], PtMax[iPt])))
    hRecoFD.append(infile.Get('hFDPt_%0.f_%0.f' % (PtMin[iPt], PtMax[iPt])))
    hGenPrompt.append(infile.Get('hPromptGenPt_%0.f_%0.f' % (PtMin[iPt], PtMax[iPt])))
    hGenFD.append(infile.Get('hFDGenPt_%0.f_%0.f' % (PtMin[iPt], PtMax[iPt])))

    #get unweighted yields (for uncertainty)
    nRecoPrompt = hRecoPrompt[iPt].Integral()
    nGenPrompt = hGenPrompt[iPt].Integral()
    nRecoFD = hRecoFD[iPt].Integral()
    nGenFD = hGenFD[iPt].Integral()

    #get weighted yields
    nRecoPromptWeighted, nGenPromptWeighted, nRecoFDWeighted, nGenFDWeighted = (0 for i in range(4))
    for iBin in range(hRecoPrompt[iPt].GetNbinsX()):
        if args.ptweights:
            binweigths = hPtWeights.GetXaxis().FindBin(hRecoPrompt[iPt].GetBinCenter(iBin+1))
            weight = hPtWeights.GetBinContent(binweigths)
        else:
            weight = 1
        nRecoPromptWeighted += weight*hRecoPrompt[iPt].GetBinContent(iBin+1)
        nGenPromptWeighted += weight*hGenPrompt[iPt].GetBinContent(iBin+1)
        nRecoFDWeighted += weight*hRecoFD[iPt].GetBinContent(iBin+1)
        nGenFDWeighted += weight*hGenFD[iPt].GetBinContent(iBin+1)

    effPrompt, effPromptUnc = ComputeEfficiency(nRecoPromptWeighted,nGenPromptWeighted,nRecoPromptWeighted/math.sqrt(nRecoPrompt),nGenPromptWeighted/math.sqrt(nGenPrompt))
    effFD, effFDUnc = ComputeEfficiency(nRecoFDWeighted,nGenFDWeighted,nRecoFDWeighted/math.sqrt(nRecoFD),nGenFDWeighted/math.sqrt(nGenFD))
    hEffPrompt.SetBinContent(iPt+1,effPrompt)
    hEffPrompt.SetBinError(iPt+1,effPromptUnc)
    hEffFD.SetBinContent(iPt+1,effFD)
    hEffFD.SetBinError(iPt+1,effFDUnc)

    hYieldPromptGen.SetBinContent(iPt+1, nGenPromptWeighted)
    hYieldPromptGen.SetBinError(iPt+1, nGenPromptWeighted/math.sqrt(nGenPrompt))
    hYieldFDGen.SetBinContent(iPt+1, nGenFDWeighted)
    hYieldFDGen.SetBinError(iPt+1, nGenFDWeighted/math.sqrt(nGenFD))
    hYieldPromptReco.SetBinContent(iPt+1, nRecoPromptWeighted)
    hYieldPromptReco.SetBinError(iPt+1, nRecoPromptWeighted/math.sqrt(nRecoPrompt))
    hYieldFDReco.SetBinContent(iPt+1, nRecoFDWeighted)
    hYieldFDReco.SetBinError(iPt+1, nRecoFDWeighted/math.sqrt(nRecoFD))

gStyle.SetPadRightMargin(0.035)
gStyle.SetPadLeftMargin(0.14)
gStyle.SetPadTopMargin(0.035)
gStyle.SetTitleSize(0.045,'xy')
gStyle.SetLabelSize(0.040,'xy')
gStyle.SetPadTickX(1)
gStyle.SetPadTickY(1)
gStyle.SetLegendBorderSize(0)
gStyle.SetOptStat(0)

leg = TLegend(0.6,0.2,0.8,0.4)
leg.SetTextSize(0.045)
leg.SetFillStyle(0)
leg.AddEntry(hEffPrompt,"Prompt","p")
leg.AddEntry(hEffFD,"Feed-down","p")

hEffPrompt.SetDirectory(0)
hEffFD.SetDirectory(0)
cEff = TCanvas('cEff','',800,800)
cEff.DrawFrame(PtMin[0],1.e-5,PtMax[nPtBins-1],1.,';#it{p}_{T} (GeV/#it{c});Efficiency;')
cEff.SetLogy()
hEffPrompt.Draw('same')
hEffFD.Draw('same')
leg.Draw()

outFile = TFile(args.outFileName,'recreate')
hEffPrompt.Write()
hEffFD.Write()
hYieldPromptGen.Write()
hYieldFDGen.Write()
hYieldPromptReco.Write()
hYieldFDReco.Write()
outFile.Close()

if not args.batch:
    if six.PY2:
        outFileNamePDF = string.replace(args.outFileName,'.root','.pdf')
        cEff.SaveAs(outFileNamePDF)
        raw_input('Press enter to exit')
    elif six.PY3:
        outFileNamePDF = args.outFileName.replace('.root','.pdf')
        cEff.SaveAs(outFileNamePDF)
        input('Press enter to exit')
