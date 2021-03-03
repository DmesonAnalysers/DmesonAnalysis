'''
Script for the computation of pT shape weights
run: python3 ComputePtGenWeights.py inFileMC.root outFile.root
                                    [--Dspecie Dname] [--Bspecie Bname]
                                    [--PbPb] [--rebin] [--smooth]
'''

import sys
import argparse
import yaml
from ROOT import TFile  # pylint: disable=import-error,no-name-in-module
sys.path.append('../..')
from utils.ReadModel import ReadFONLL, ReadTAMU, ReadPHSD, ReadCatania, ReadMCatsHQ  #pylint: disable=wrong-import-position,import-error

parser = argparse.ArgumentParser(description='Arguments to pass')
parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='input yaml cibfugfile name')
args = parser.parse_args()

with open(args.cfgFileName, 'r') as ymlfitConfigFile:
    inputCfg = yaml.load(ymlfitConfigFile, yaml.FullLoader)

Dspecie = inputCfg['meson']['Dspecie']
Bspecie = inputCfg['meson']['Bspecie']

fileNameMC = inputCfg['inputMC']['filename']
suffixCF = inputCfg['inputMC']['suffixCF']

rebin = inputCfg['options']['rebin']
smooth = inputCfg['options']['smooth']

shapesD = inputCfg['shapes']['D']
shapesB = inputCfg['shapes']['B']

outFileName = inputCfg['outputfile']

if not Dspecie:
    print('ERROR: Dspecie argument is mandatory! Exit')

if Dspecie not in ['Ds', 'Dplus', 'Dzero', 'Lc']:
    print(f'ERROR: D specie {Dspecie} not supported! Only Ds, Dplus, Dzero, Dstar, Lc are supported! Exit')

# load MC input
infileGenPtShape = TFile.Open(fileNameMC)
listGenPtShape = infileGenPtShape.Get('HFMCCheck/clistHFMCCheck')
if listGenPtShape:
    if Dspecie == 'Ds':
        hPtYGenD = listGenPtShape.FindObject('hyptDsprompt')
    elif Dspecie == 'Dplus':
        hPtYGenD = listGenPtShape.FindObject('hyptDplusprompt')
    elif Dspecie == 'Dzero':
        hPtYGenD = listGenPtShape.FindObject('hyptD0prompt')
    elif Dspecie == 'Dstar':
        hPtYGenD = listGenPtShape.FindObject('hyptDstarprompt')
    elif Dspecie == 'Lc':
        hPtYGenD = listGenPtShape.FindObject('hyptLcprompt')
    hPtYGenD.SetDirectory(0)
    hPtGenD = hPtYGenD.ProjectionX('hPtGenD')
elif suffixCF != '':
    if Dspecie == 'Ds':
        dirPtShape = infileGenPtShape.Get(f'PWG3_D2H_CFtaskDstoKKpi_CommonFramework_Phi{suffixCF}')
        contPtShape = dirPtShape.Get(f'CFHFccontainer0_3ProngDstoKKpi_CommonFramework_Phi{suffixCF}')
        hPtGenD = contPtShape.Project(0, 0)
        hPtGenD.SetName('hPtGenD')
    else:
        print(f'ERROR: D specie {Dspecie} not implemented for CF outputs! Exit')
hPtGenD.SetDirectory(0)
hPtGenD.Sumw2()
hPtGenD.Rebin(rebin)
hPtGenD.Scale(1./hPtGenD.Integral())

if Bspecie:
    if Bspecie not in ['Bs', 'Bplus', 'Bzero', 'Lb', 'BsBmix']:
        print(f'ERROR: B specie {Bspecie} not supported! Only Bs, Bplus, Bzero, Lb, BsBmix are supported! Exit')
        sys.exit()
    if listGenPtShape:
        if Bspecie in ['Bs', 'BsBmix']:
            hYPtGenB = listGenPtShape.FindObject('hyptBsAllDecay')
        elif Bspecie == 'Bplus':
            hYPtGenB = listGenPtShape.FindObject('hyptBplusAllDecay')
        elif Bspecie == 'Bzero':
            hYPtGenB = listGenPtShape.FindObject('hyptB0AllDecay')
        elif Bspecie == 'Lb':
            hYPtGenB = listGenPtShape.FindObject('hyptLbAllDecay')
        hYPtGenB.SetDirectory(0)
        hPtGenB = hYPtGenB.ProjectionX('hPtGenB')
        hPtGenB.SetName('hPtGenB')
        if Bspecie == 'BsBmix':
            hYPtGenB2 = listGenPtShape.FindObject('hyptB0AllDecay')
            hPtGenB2 = hYPtGenB.ProjectionX('hPtGenB2')
            hPtGenB2.SetName('hPtGenB2')
    elif suffixCF != '':
        print(f'ERROR: B specie {Bspecie} not implemented for CF outputs! Exit')
    hPtGenB.SetDirectory(0)
    hPtGenB.Sumw2()
    hPtGenB.Rebin(rebin)
    hPtGenB.Scale(1./hPtGenB.Integral())
    if Bspecie == 'BsBmix':
        hPtGenB2.SetDirectory(0)
        hPtGenB2.Sumw2()
        hPtGenB2.Rebin(rebin)
        hPtGenB2.Scale(1./hPtGenB2.Integral())
        hPtGenB.Add(hPtGenB2) # assuming 50% Bs and 50% B, reasonable for non-prompt Ds
        hPtGenB.Scale(1./2)

infileGenPtShape.Close()

# default models
isPbPb = False
sFONLLD, _, ptMinFONLL, ptMaxFONLL = ReadFONLL(shapesD['fonll']['file'], True)
if Bspecie:
    sFONLLB, _, ptMinFONLLB, ptMaxFONLLB = ReadFONLL(shapesB['fonll']['file'], True)

if 'tamu' in shapesD and shapesD['tamu']['enabled']:
    sTAMU, _, ptMinTAMU, ptMaxTAMU = ReadTAMU(shapesD['tamu']['file'])
if 'phsd' in shapesD and shapesD['phsd']['enabled']:
    sPHSD, _, ptMinPHSD, ptMaxPHSD = ReadPHSD(shapesD['phsd']['file'])
if 'catania' in shapesD and shapesD['catania']['enabled']:
    sCatania, _, ptMinCatania, ptMaxCatania = ReadCatania(shapesD['catania']['file'])
if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
    sGossiaux, _, ptMinGoss, ptMaxGoss = ReadMCatsHQ(shapesD['mc@shq']['file'])

if 'tamu' in shapesB and shapesB['tamu']['enabled']:
    if Bspecie not in ['Bplus', 'Bzero', 'Bs', 'BsBmix']:
        print(f'WARNING: no models available for {Bspecie} RAA, using B RAA as defult')
        sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(shapesB['tamu']['file']['B'])
    elif Bspecie in ['Bplus', 'Bzero']:
        sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(shapesB['tamu']['file']['B'])
    elif Bspecie == 'Bs':
        sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(shapesB['tamu']['file']['Bs'])
    elif Bspecie == 'BsBmix':
        sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(shapesB['tamu']['file']['B'])
        sTAMUBs, _, ptMinTAMUBs, ptMaxTAMUBs = ReadTAMU(shapesB['tamu']['file']['Bs'])

histoDNames = ['hPtFONLLDcent', 'hPtFONLLDmin', 'hPtFONLLDmax']
histoBNames = ['hPtFONLLBcent', 'hPtFONLLBmin', 'hPtFONLLBmax']
modelPred = ['yCent', 'yMin', 'yMax']

hPtFONLLD, hPtFONLLB, hPtFONLLtimesTAMUD, hPtFONLLtimesTAMUB, hPtFONLLtimesPHSDD, \
    hPtFONLLtimesGossiauxD, hPtFONLLtimesCataniaD = ([] for _ in range(7))

hPtWeightsFONLLD, hPtWeightsFONLLB, hPtWeightsFONLLtimesTAMUD, hPtWeightsFONLLtimesTAMUB, \
    hPtWeightsFONLLtimesPHSDD, hPtWeightsFONLLtimesGossiauxD, hPtWeightsFONLLtimesCataniaD = ([] for _ in range(7))

# D meson weights
for histoName, pred in zip(histoDNames, modelPred):
    hPtFONLLD.append(hPtGenD.Clone(histoName))
    if 'tamu' in shapesD and shapesD['tamu']['enabled']:
        hPtFONLLtimesTAMUD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesTAMU')))
    if 'phsd' in shapesD and shapesD['phsd']['enabled']:
        hPtFONLLtimesPHSDD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesPHSD')))
    if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
        hPtFONLLtimesGossiauxD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesGossiaux')))
    if 'catania' in shapesD and shapesD['catania']['enabled']:
        hPtFONLLtimesCataniaD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesCatania')))

    for iPt in range(1, hPtFONLLD[-1].GetNbinsX()+1):
        ptCent = hPtFONLLD[-1].GetBinCenter(iPt)
        if ptMinFONLL < ptCent < ptMaxFONLL:
            hPtFONLLD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent))
        elif ptCent > ptMaxFONLL:
            print(f'WARNING: Results for pT > {ptMaxFONLL} not reliable!')
            continue
        else:
            print(f'WARNING: Results for pT < {ptMinFONLL} not reliable!')
            continue
        if 'tamu' in shapesD and shapesD['tamu']['enabled']:
            if ptMinTAMU < ptCent < ptMaxTAMU:
                hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptCent))
            elif ptCent > ptMaxTAMU:
                hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptMaxTAMU))
            else:
                hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptMinTAMU))

        if 'phsd' in shapesD and shapesD['phsd']['enabled']:
            if ptMinPHSD < ptCent < ptMaxPHSD:
                hPtFONLLtimesPHSDD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sPHSD['yCent'](ptCent))
            elif ptCent > ptMaxPHSD:
                hPtFONLLtimesPHSDD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sPHSD['yCent'](ptMaxPHSD))
            else:
                hPtFONLLtimesPHSDD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sPHSD['yCent'](ptMinPHSD))

        if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
            if ptMinGoss < ptCent < ptMaxGoss:
                hPtFONLLtimesGossiauxD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sGossiaux['yCent'](ptCent))
            elif ptCent > ptMaxGoss:
                hPtFONLLtimesGossiauxD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sGossiaux['yCent'](ptMaxGoss))
            else:
                hPtFONLLtimesGossiauxD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sGossiaux['yCent'](ptMinGoss))

        if 'catania' in shapesD and shapesD['catania']['enabled']:
            if ptMinCatania < ptCent < ptMaxCatania:
                hPtFONLLtimesCataniaD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sCatania['yCent'](ptCent))
            elif ptCent > ptMaxCatania:
                hPtFONLLtimesCataniaD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sCatania['yCent'](ptMaxCatania))
            else:
                hPtFONLLtimesCataniaD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sCatania['yCent'](ptMinCatania))

    hPtFONLLD[-1].Sumw2()
    hPtFONLLD[-1].Scale(1./hPtFONLLD[-1].Integral())
    hPtWeightsFONLLD.append(hPtFONLLD[-1].Clone(histoName.replace('Pt', 'PtWeights')))
    hPtWeightsFONLLD[-1].Divide(hPtFONLLD[-1], hPtGenD)
    hPtWeightsFONLLD[-1].Smooth(smooth)
    if 'tamu' in shapesD and shapesD['tamu']['enabled']:
        hPtFONLLtimesTAMUD[-1].Scale(1./hPtFONLLtimesTAMUD[-1].Integral())
        hPtWeightsFONLLtimesTAMUD.append(
            hPtFONLLtimesTAMUD[-1].Clone(hPtFONLLtimesTAMUD[-1].GetName().replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesTAMUD[-1].Divide(hPtFONLLtimesTAMUD[-1], hPtGenD)
        hPtWeightsFONLLtimesTAMUD[-1].Smooth(smooth)
    if 'phsd' in shapesD and shapesD['phsd']['enabled']:
        hPtFONLLtimesPHSDD[-1].Scale(1./hPtFONLLtimesPHSDD[-1].Integral())
        hPtWeightsFONLLtimesPHSDD.append(
            hPtFONLLtimesPHSDD[-1].Clone(hPtFONLLtimesPHSDD[-1].GetName().replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesPHSDD[-1].Divide(hPtFONLLtimesPHSDD[-1], hPtGenD)
        hPtWeightsFONLLtimesPHSDD[-1].Smooth(smooth)
    if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
        hPtFONLLtimesGossiauxD[-1].Scale(1./hPtFONLLtimesGossiauxD[-1].Integral())
        hPtWeightsFONLLtimesGossiauxD.append(
            hPtFONLLtimesGossiauxD[-1].Clone(hPtFONLLtimesGossiauxD[-1].GetName().replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesGossiauxD[-1].Divide(hPtFONLLtimesGossiauxD[-1], hPtGenD)
        hPtWeightsFONLLtimesGossiauxD[-1].Smooth(smooth)
    if 'catania' in shapesD and shapesD['catania']['enabled']:
        hPtFONLLtimesCataniaD[-1].Scale(1./hPtFONLLtimesCataniaD[-1].Integral())
        hPtWeightsFONLLtimesCataniaD.append(
            hPtFONLLtimesCataniaD[-1].Clone(hPtFONLLtimesCataniaD[-1].GetName().replace('Pt', 'PtWeights')))
        hPtWeightsFONLLtimesCataniaD[-1].Divide(hPtFONLLtimesCataniaD[-1], hPtGenD)
        hPtWeightsFONLLtimesCataniaD[-1].Smooth(smooth)

# B meson weights
if Bspecie:
    for histoName, pred in zip(histoBNames, modelPred):
        hPtFONLLB.append(hPtGenB.Clone(histoName))
        if 'tamu' in shapesB and shapesB['tamu']['enabled']:
            hPtFONLLtimesTAMUB.append(hPtGenB.Clone(histoName.replace('FONLL', 'FONLLtimesTAMU')))

        for iPt in range(1, hPtFONLLB[-1].GetNbinsX()+1):
            ptCent = hPtFONLLB[-1].GetBinCenter(iPt)
            if ptMinFONLLB < ptCent < ptMaxFONLLB:
                hPtFONLLB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent))
            elif ptCent > ptMaxFONLLB:
                print(f'WARNING: Results for pT > {ptMaxFONLLB} not reliable! Set weights to 0')
                continue
            else:
                print(f'WARNING: Results for pT < {ptMinFONLLB} not reliable! Set weights to 0')
                continue
            if 'tamu' in shapesB and shapesB['tamu']['enabled']:
                if Bspecie != 'BsBmix':
                    if ptMinTAMUB < ptCent < ptMaxTAMUB:
                        hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptCent))
                    elif ptCent > ptMaxTAMUB:
                        hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptMaxTAMUB))
                    else:
                        hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptMinTAMUB))
                else:
                    ptMaxMix = min([ptMaxTAMUB, ptMaxTAMUBs])
                    ptMinMix = max([ptMinTAMUB, ptMinTAMUBs])
                    if ptMinMix < ptCent < ptMaxMix:
                        rAAMix = (sTAMUB['yCent'](ptCent) + sTAMUBs['yCent'](ptCent)) / 2
                    elif ptCent > ptMaxMix:
                        rAAMix = (sTAMUB['yCent'](ptMaxMix) + sTAMUBs['yCent'](ptMaxMix)) / 2
                    else:
                        rAAMix = (sTAMUB['yCent'](ptMinMix) + sTAMUBs['yCent'](ptMinMix)) / 2
                    hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * rAAMix)

        hPtFONLLB[-1].Sumw2()
        hPtFONLLB[-1].Scale(1./hPtFONLLB[-1].Integral())
        hPtWeightsFONLLB.append(hPtFONLLB[-1].Clone(histoName.replace('Pt', 'PtWeights')))
        hPtWeightsFONLLB[-1].Divide(hPtFONLLB[-1], hPtGenB)
        hPtWeightsFONLLB[-1].Smooth(smooth)
        if 'tamu' in shapesB and shapesB['tamu']['enabled']:
            hPtFONLLtimesTAMUB[-1].Scale(1./hPtFONLLtimesTAMUB[-1].Integral())
            hPtWeightsFONLLtimesTAMUB.append(
                hPtFONLLtimesTAMUB[-1].Clone(hPtFONLLtimesTAMUB[-1].GetName().replace('Pt', 'PtWeights')))
            hPtWeightsFONLLtimesTAMUB[-1].Divide(hPtFONLLtimesTAMUB[-1], hPtGenD)
            hPtWeightsFONLLtimesTAMUB[-1].Smooth(smooth)

outfile = TFile(outFileName, 'recreate')
hPtGenD.Write()
if Bspecie:
    hPtGenB.Write()
for iHisto, _ in enumerate(hPtFONLLD):
    hPtFONLLD[iHisto].Write()
    hPtWeightsFONLLD[iHisto].Write()
    if 'tamu' in shapesD and shapesD['tamu']['enabled']:
        hPtFONLLtimesTAMUD[iHisto].Write()
        hPtWeightsFONLLtimesTAMUD[iHisto].Write()
    if 'phsd' in shapesD and shapesD['phsd']['enabled']:
        hPtFONLLtimesPHSDD[iHisto].Write()
        hPtWeightsFONLLtimesPHSDD[iHisto].Write()
    if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
        hPtFONLLtimesGossiauxD[iHisto].Write()
        hPtWeightsFONLLtimesGossiauxD[iHisto].Write()
    if 'catania' in shapesD and shapesD['catania']['enabled']:
        hPtFONLLtimesCataniaD[iHisto].Write()
        hPtWeightsFONLLtimesCataniaD[iHisto].Write()
    if Bspecie:
        hPtFONLLB[iHisto].Write()
        hPtWeightsFONLLB[iHisto].Write()
        if 'tamu' in shapesB and shapesB['tamu']['enabled']:
            hPtFONLLtimesTAMUB[iHisto].Write()
            hPtWeightsFONLLtimesTAMUB[iHisto].Write()
outfile.Close()
