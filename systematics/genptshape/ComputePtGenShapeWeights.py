'''
Script for the computation of pT shape weights
run: python3 ComputePtGenWeights.py inFileMC.root outFile.root [--Dspecie Dname] [--Bspecie Bname] [--PbPb]
'''

import sys
import argparse
from ROOT import TFile  # pylint: disable=import-error,no-name-in-module
sys.path.insert(0, '../..')
#pylint: disable=wrong-import-position,import-error,no-name-in-module
from utils.ReadModel import ReadFONLL, ReadTAMU, ReadPHSD, ReadCatania, ReadMCatsHQ

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('fileNameMC', metavar='text', default='inFileMC.root',
                    help='input MC root file name')
parser.add_argument('outFileName', metavar='text', default='outFile.root',
                    help='output root file name')
parser.add_argument('--Dspecie', metavar='text', default='Ds',
                    help='D specie for which compute the weights (Ds, Dplus, Dzero, Dstar, Lc)')
parser.add_argument('--Bspecie', metavar='text', required=False,
                    help='B specie for which compute the weights (Bs, Bplus, Bzero, Lb)')
parser.add_argument('--PbPb', action='store_true', default=False,
                    help='flag to activate pT shapes with Raa')
args = parser.parse_args()

if not args.Dspecie:
    print('ERROR: Dspecie argument is mandatory! Exit')

if args.Dspecie not in ['Ds', 'Dplus', 'Dzero', 'Lc']:
    print(f'ERROR: D specie {args.Dspecie} not supported! Only Ds, Dplus, Dzero, Dstar, Lc are supported! Exit')

# load MC input
infileGenPtShape = TFile.Open(args.fileNameMC)
listGenPtShape = infileGenPtShape.Get('HFMCCheck/clistHFMCCheck')
if args.Dspecie == 'Ds':
    hPtYGenD = listGenPtShape.FindObject('hyptDsprompt')
elif args.Dspecie == 'Dplus':
    hPtYGenD = listGenPtShape.FindObject('hyptDplusprompt')    
elif args.Dspecie == 'Dzero':
    hPtYGenD = listGenPtShape.FindObject('hyptD0prompt')
elif args.Dspecie == 'Dstar':
    hPtYGenD = listGenPtShape.FindObject('hyptDstarprompt')
elif args.Dspecie == 'Lc':
    hPtYGenD = listGenPtShape.FindObject('hyptLcprompt')
hPtYGenD.SetDirectory(0)
hPtGenD = hPtYGenD.ProjectionX('hPtGenD')
hPtGenD.SetDirectory(0)
hPtGenD.Sumw2()
hPtGenD.Scale(1./hPtGenD.Integral())

if args.Bspecie:
    if args.Bspecie not in ['Bs', 'Bplus', 'Bzero', 'Lb']:
        print(f'ERROR: B specie {args.Bspecie} not supported! Only Bs, Bplus, Bzero, Lb are supported! Exit')
    if args.Bspecie == 'Bs':
        hYPtGenB = listGenPtShape.FindObject('hyptBsAllDecay')
    elif args.Bspecie == 'Bplus':
        hYPtGenB = listGenPtShape.FindObject('hyptBplusAllDecay')
    elif args.Bspecie == 'Bzero':
        hYPtGenB = listGenPtShape.FindObject('hyptB0AllDecay')
    elif args.Bspecie == 'Lb':
        hYPtGenB = listGenPtShape.FindObject('hyptLbAllDecay')
    hYPtGenB.SetDirectory(0)
    hPtGenB = hYPtGenB.ProjectionX('hPtGenB')
    hPtGenB.SetDirectory(0)
    hPtGenB.Sumw2()
    hPtGenB.Scale(1./hPtGenB.Integral())

infileGenPtShape.Close()

# default models
# TODO: consider to make them arguments
sFONLLD, _ = ReadFONLL('../../models/fonll/FONLL_D_pp5_y05.txt', True)
if args.Bspecie:
    sFONLLB, _ = ReadFONLL('../../models/fonll/FONLL_B_pp5_y05.txt', True)

if args.PbPb:
    if args.Dspecie == 'Ds':
        sTAMU, _ = ReadTAMU('../../models/tamu/Ds_TAMU_RAA_5TeV_010.txt')
        sPHSD, _ = ReadPHSD('../../models/phsd/Ds_PHSD_RAA_5TeV_010.txt')
        sCatania, _ = ReadCatania('../../models/catania/Ds_Catania_RAA_5TeV_010.dat')
    else:
        sTAMU, _ = ReadTAMU('../../models/tamu/D_TAMU_RAA_5TeV_010.txt')
        sPHSD, _ = ReadPHSD('../../models/phsd/D_PHSD_RAA_5TeV_010.txt')
        sCatania, _ = ReadCatania('../../models/catania/D_Catania_RAA_5TeV_010.dat')
    sGossiaux, _ = ReadMCatsHQ('../../models/mcatshq/D_Gossiaux_RAA_5TeV_010.txt')

# TODO: add FONLLxRaa weights for B
histoDNames = ['hPtFONLLDcent', 'hPtFONLLDmin', 'hPtFONLLDmax']
histoBNames = ['hPtFONLLBcent', 'hPtFONLLBmin', 'hPtFONLLBmax']
modelPred = ['yCent', 'yMin', 'yMax']

hPtFONLLD, hPtFONLLB, hPtFONLLtimesTAMUD, hPtFONLLtimesPHSDD, \
    hPtFONLLtimesGossiauxD, hPtFONLLtimesCataniaD = ([] for _ in range(6))

hPtWeightsFONLLD, hPtWeightsFONLLB, hPtWeightsFONLLtimesTAMUD, hPtWeightsFONLLtimesPHSDD, \
    hPtWeightsFONLLtimesGossiauxD, hPtWeightsFONLLtimesCataniaD = ([] for _ in range(6))

# D meson weights
for histoName, pred in zip(histoDNames, modelPred):
    hPtFONLLD.append(hPtGenD.Clone(histoName))
    if args.PbPb:
        hPtFONLLtimesTAMUD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesTAMU')))
        hPtFONLLtimesPHSDD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesPHSD')))
        hPtFONLLtimesGossiauxD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesGossiaux')))
        hPtFONLLtimesCataniaD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesCatania')))

    for iPt in range(1, hPtFONLLD[-1].GetNbinsX()+1):
        ptCent = hPtFONLLD[-1].GetBinCenter(iPt)
        hPtFONLLD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent))
        if args.PbPb:
            hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptCent))
            hPtFONLLtimesPHSDD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sPHSD['yCent'](ptCent))
            hPtFONLLtimesGossiauxD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sGossiaux['yCent'](ptCent))
            hPtFONLLtimesCataniaD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sCatania['yCent'](ptCent))

    hPtFONLLD[-1].Sumw2()
    hPtFONLLD[-1].Scale(1./hPtFONLLD[-1].Integral())
    hPtWeightsFONLLD.append(hPtFONLLD[-1].Clone(histoName.replace('Pt', 'PtWeights')))
    hPtWeightsFONLLD[-1].Divide(hPtFONLLD[-1], hPtGenD)
    if args.PbPb:
        hPtFONLLtimesTAMUD[-1].Scale(1./hPtFONLLtimesTAMUD[-1].Integral())
        hPtFONLLtimesPHSDD[-1].Scale(1./hPtFONLLtimesPHSDD[-1].Integral())
        hPtFONLLtimesGossiauxD[-1].Scale(1./hPtFONLLtimesGossiauxD[-1].Integral())
        hPtFONLLtimesCataniaD[-1].Scale(1./hPtFONLLtimesCataniaD[-1].Integral())

        hPtWeightsFONLLtimesTAMUD[-1].append(hPtFONLLtimesTAMUD[-1].Clone(histoName.replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesTAMUD[-1].Divide(hPtFONLLtimesTAMUD[-1], hPtGenD)
        hPtWeightsFONLLtimesPHSDD.append(hPtFONLLtimesPHSDD[-1].Clone(histoName.replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesPHSDD[-1].Divide(hPtFONLLtimesPHSDD[-1], hPtGenD)
        hPtWeightsFONLLtimesGossiauxD.append(hPtFONLLtimesGossiauxD[-1].Clone(histoName.replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesGossiauxD[-1].Divide(hPtFONLLtimesGossiauxD[-1], hPtGenD)
        hPtWeightsFONLLtimesCataniaD.append(hPtFONLLtimesCataniaD[-1].Clone(histoName.replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesCataniaD[-1].Divide(hPtFONLLtimesCataniaD[-1], hPtGenD)

# B meson weights
if args.Bspecie:
    for histoName, pred in zip(histoBNames, modelPred):
        hPtFONLLB.append(hPtGenB.Clone(histoName))

        for iPt in range(1, hPtFONLLB[-1].GetNbinsX()+1):
            ptCent = hPtFONLLB[-1].GetBinCenter(iPt)
            hPtFONLLB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent))

        hPtFONLLB[-1].Sumw2()
        hPtFONLLB[-1].Scale(1./hPtFONLLB[-1].Integral())
        hPtWeightsFONLLB.append(hPtFONLLB[-1].Clone(histoName.replace('Pt', 'PtWeights')))
        hPtWeightsFONLLB[-1].Divide(hPtFONLLB[-1], hPtGenB)

outfile = TFile(args.outFileName, 'recreate')
hPtGenD.Write()
if args.Bspecie:
    hPtGenB.Write()
for iHisto, _ in enumerate(hPtFONLLD):
    hPtFONLLD[iHisto].Write()
    hPtWeightsFONLLD[iHisto].Write()
    if args.PbPb:
        hPtFONLLtimesTAMUD[iHisto].Write()
        hPtFONLLtimesPHSDD[iHisto].Write()
        hPtFONLLtimesGossiauxD[iHisto].Write()
        hPtFONLLtimesCataniaD[iHisto].Write()
        hPtWeightsFONLLtimesTAMUD[iHisto].Write()
        hPtWeightsFONLLtimesPHSDD[iHisto].Write()
        hPtWeightsFONLLtimesGossiauxD[iHisto].Write()
        hPtWeightsFONLLtimesCataniaD[iHisto].Write()
    if args.Bspecie:
        hPtFONLLB[iHisto].Write()
        hPtWeightsFONLLB[iHisto].Write()
outfile.Close()
