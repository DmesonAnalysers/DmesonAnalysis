'''
Script for the computation of pT shape weights
run: python3 ComputePtGenWeights.py cfgFileName.yml
'''
import os
import sys
import argparse
import yaml
from ROOT import TFile, TCanvas, TLegend  # pylint: disable=import-error,no-name-in-module
from ROOT import kBlack, kRed, kAzure, kGray, kOrange, kFullCircle, kFullSquare, kOpenCircle
sys.path.append('../../..')
from utils.ReadModel import ReadFONLL, ReadTAMU, ReadPHSD, ReadCatania, ReadMCatsHQ, ReadLIDO, ReadLGR  #pylint: disable=wrong-import-position,import-error
from utils.StyleFormatter import SetObjectStyle #pylint: disable=wrong-import-position,import-error
from sparse_dicts import get_sparses


def computePtWeights(config, outputDir, suffix):

    # load input configuration
    with open(config, 'r') as ymlCfgFile:
        config = yaml.load(ymlCfgFile, yaml.FullLoader)

    Dspecie = config['Dmeson']
    Bspecie = config['Bspecie']
    # 'Ds', 'Dplus', 'Dzero', 'Lc'
    if Dspecie not in ['Ds', 'Dplus', 'Dzero', 'Lc']:
        print(f'ERROR: D specie {Dspecie} not supported! Only Ds, Dplus, Dzero, Lc is supported! Exit')
        sys.exit()

    rebin = config['rebin']
    smooth = config['smooth']

    shapesD = config['shapes']['D']
    shapesB = config['shapes']['B']

    _, _, sparsesGen, axes = get_sparses(config, False, False, True)
    # load thnSparse
    sparseGenD = sparsesGen['GenPrompt']
    sparseGenD.SetName('sparseGenD')
    sparseGenB = sparsesGen['GenFD']
    sparseGenB.SetName('sparseGenB')

    hPtGenD = sparseGenD.Projection(axes['GenPrompt']['Pt'])
    hPtGenD.SetDirectory(0)
    hPtGenD.SetName('hPtGenD')
    hPtGenD.Sumw2()
    hPtGenD.Rebin(rebin)
    hPtGenD.Scale(1./hPtGenD.Integral())

    if Bspecie:
        hPtGenB = sparseGenB.Projection(axes['GenFD']['Pt'])
        if Dspecie == 'Ds':
            sparseGenBPlusBZero = sparseGenB.Clone('sparseGenBPlusBZero')
            sparseGenBPlusBZero.GetAxis(axes['GenFD']['origin']).SetRange(1, 2)
            sparseGenLambdaBZero = sparseGenB.Clone('sparseGenLambdaBZero')
            sparseGenLambdaBZero.GetAxis(axes['GenFD']['origin']).SetRange(4, 4)
            hPtGenB = sparseGenLambdaBZero.Projection(axes['GenFD']['Pt'])
            hPtGenB.Add(sparseGenBPlusBZero.Projection(axes['GenFD']['Pt']))
        
        #TODO: modifications for other B mesons
        hPtGenB.SetDirectory(0)
        hPtGenB.SetName('hPtGenB')
        hPtGenB.Sumw2()
        hPtGenB.Rebin(rebin)
        hPtGenB.Scale(1./hPtGenB.Integral())
    
    if Bspecie == 'BsBmix':
        sparseGenB.GetAxis(axes['GenFD']['origin']).SetRange('the bin number of Bs flag, the bin number of Bs flag')
        hPtGenBs = sparseGenB.Projection(axes['GenFD']['Pt'])
        hPtGenBs.SetDirectory(0)
        hPtGenBs.SetName('hPtGenBs')
        hPtGenBs.Rebin(rebin)
        hPtGenBs.Scale(1./2 * hPtGenBs.Integral())
        hPtGenB.Scale(1./2)
        hPtGenB.Add(hPtGenBs) # assuming 50% Bs and 50% B, reasonable for non-prompt Ds

    print('INFO: MC input loaded')
    # load models predictions
    #___________________________________________________________________________________________________________________________
    sFONLLD, _, ptMinFONLL, ptMaxFONLL = ReadFONLL(shapesD['fonll'], True, Dspecie)
    if Bspecie:
        sFONLLB, _, ptMinFONLLB, ptMaxFONLLB = ReadFONLL(shapesB['fonll'], True, 'B')
        # this ReadFONLL method should be modified to accept B meson species, but depends on the FONLL we will use
        # now it will load all B meson

    if 'tamu' in shapesD and shapesD['tamu']['enabled']:
        sTAMU, _, ptMinTAMU, ptMaxTAMU = ReadTAMU(shapesD['tamu']['file'])
    # if 'phsd' in shapesD and shapesD['phsd']['enabled']:
    #     sPHSD, _, ptMinPHSD, ptMaxPHSD = ReadPHSD(shapesD['phsd']['file'])
    # if 'catania' in shapesD and shapesD['catania']['enabled']:
    #     sCatania, _, ptMinCatania, ptMaxCatania = ReadCatania(shapesD['catania']['file'])
    # if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
    #     sGossiaux, _, ptMinGoss, ptMaxGoss = ReadMCatsHQ(shapesD['mc@shq']['file'])
    # if 'lido' in shapesD and shapesD['lido']['enabled']:
    #     sLIDO, _, ptMinLIDO, ptMaxLIDO = ReadLIDO(shapesD['lido']['file'])
    # if 'lgr' in shapesD and shapesD['lgr']['enabled']:
    #     sLGR, _, ptMinLGR, ptMaxLGR = ReadLGR(shapesD['lgr']['file'])
    # if 'datashape' in shapesD and shapesD['datashape']['enabled']:
    #     fData, _, _ = FitSpectra(shapesD['datashape']['file'], shapesD['datashape']['hSpectraName'],
    #                             shapesD['datashape']['systErrName'], shapesD['datashape']['parConfig'],
    #                             shapesD['datashape']['controlFile'])

    if 'tamu' in shapesB and shapesB['tamu']['enabled']:
        if Bspecie in ['Ball', 'BsBmix']:
            print(f'WARNING: no TAMU model available for {Bspecie} RAA, using B RAA as defult')
            sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(shapesB['tamu']['file']['B'])
        # elif Bspecie in ['Bplus', 'Bzero']:
        #     sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(shapesB['tamu']['file']['B'])
        # elif Bspecie == 'Bs':
        #     sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(shapesB['tamu']['file']['Bs'])
        elif Bspecie == 'BsBmix':
            sTAMUB, _, ptMinTAMUB, ptMaxTAMUB = ReadTAMU(shapesB['tamu']['file']['B'])
            sTAMUBs, _, ptMinTAMUBs, ptMaxTAMUBs = ReadTAMU(shapesB['tamu']['file']['Bs'])
    # if 'lido' in shapesB and shapesB['lido']['enabled']:
    #     if Bspecie not in ['Bplus', 'Bzero']:
    #         print(f'WARNING: no LIDO model available for {Bspecie} RAA, using B RAA as defult')
    #     sLIDOB, _, ptMinLIDOB, ptMaxLIDOB = ReadLIDO(shapesB['lido']['file'])

    histoDNames = ['hPtFONLLDcent', 'hPtFONLLDmin', 'hPtFONLLDmax']
    histoBNames = ['hPtFONLLBcent', 'hPtFONLLBmin', 'hPtFONLLBmax']
    modelPred = ['yCent', 'yMin', 'yMax']

    hPtFONLLD, hPtFONLLB, hPtFONLLtimesTAMUD, hPtFONLLtimesTAMUB, hPtFONLLtimesPHSDD, \
        hPtFONLLtimesGossiauxD, hPtFONLLtimesCataniaD, hPtFONLLtimesLIDOD, hPtFONLLtimesLIDOB, \
            hPtFONLLtimesLGRD, hPtDataShape = ([] for _ in range(11))

    hPtWeightsFONLLD, hPtWeightsFONLLB, hPtWeightsFONLLtimesTAMUD, hPtWeightsFONLLtimesTAMUB, \
        hPtWeightsFONLLtimesPHSDD, hPtWeightsFONLLtimesGossiauxD, hPtWeightsFONLLtimesCataniaD, \
            hPtWeightsFONLLtimesLIDOD, hPtWeightsFONLLtimesLIDOB, hPtWeightsFONLLtimesLGRD, \
                hPtWeightsDataShape = ([] for _ in range(11))

    print('INFO: Start computing pT weights')
    # D meson weights
    #___________________________________________________________________________________________________________________________
    for histoName, pred in zip(histoDNames, modelPred):
        hPtFONLLD.append(hPtGenD.Clone(histoName))
        if 'tamu' in shapesD and shapesD['tamu']['enabled']:
            hPtFONLLtimesTAMUD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesTAMU')))
        # if 'phsd' in shapesD and shapesD['phsd']['enabled']:
        #     hPtFONLLtimesPHSDD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesPHSD')))
        # if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
        #     hPtFONLLtimesGossiauxD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesGossiaux')))
        # if 'catania' in shapesD and shapesD['catania']['enabled']:
        #     hPtFONLLtimesCataniaD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesCatania')))
        # if 'lido' in shapesD and shapesD['lido']['enabled']:
        #     hPtFONLLtimesLIDOD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesLIDO')))
        # if 'lgr' in shapesD and shapesD['lgr']['enabled']:
        #     hPtFONLLtimesLGRD.append(hPtGenD.Clone(histoName.replace('FONLL', 'FONLLtimesLGR')))
        # if 'datashape' in shapesD and shapesD['datashape']['enabled']:
        #     hPtDataShape.append(hPtGenD.Clone(histoName.replace('FONLL', 'DataShape')))

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
                if ptMinTAMU <= ptCent <= ptMaxTAMU:
                    hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptCent))
                elif ptCent > ptMaxTAMU:
                    hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptMaxTAMU))
                else:
                    hPtFONLLtimesTAMUD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sTAMU['yCent'](ptMinTAMU))

            # if 'phsd' in shapesD and shapesD['phsd']['enabled']:
            #     if ptMinPHSD <= ptCent <= ptMaxPHSD:
            #         hPtFONLLtimesPHSDD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sPHSD['yCent'](ptCent))
            #     elif ptCent > ptMaxPHSD:
            #         hPtFONLLtimesPHSDD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sPHSD['yCent'](ptMaxPHSD))
            #     else:
            #         hPtFONLLtimesPHSDD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sPHSD['yCent'](ptMinPHSD))

            # if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
            #     if ptMinGoss <= ptCent <= ptMaxGoss:
            #         hPtFONLLtimesGossiauxD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sGossiaux['yCent'](ptCent))
            #     elif ptCent > ptMaxGoss:
            #         hPtFONLLtimesGossiauxD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sGossiaux['yCent'](ptMaxGoss))
            #     else:
            #         hPtFONLLtimesGossiauxD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sGossiaux['yCent'](ptMinGoss))

            # if 'catania' in shapesD and shapesD['catania']['enabled']:
            #     if ptMinCatania <= ptCent <= ptMaxCatania:
            #         hPtFONLLtimesCataniaD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sCatania['yCent'](ptCent))
            #     elif ptCent > ptMaxCatania:
            #         hPtFONLLtimesCataniaD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sCatania['yCent'](ptMaxCatania))
            #     else:
            #         hPtFONLLtimesCataniaD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sCatania['yCent'](ptMinCatania))

            # if 'lido' in shapesD and shapesD['lido']['enabled']:
            #     if ptMinLIDO <= ptCent <= ptMaxLIDO:
            #         hPtFONLLtimesLIDOD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sLIDO['yCent'](ptCent))
            #     elif ptCent > ptMaxLIDO:
            #         hPtFONLLtimesLIDOD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sLIDO['yCent'](ptMaxLIDO))
            #     else:
            #         hPtFONLLtimesLIDOD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sLIDO['yCent'](ptMinLIDO))

            # if 'lgr' in shapesD and shapesD['lgr']['enabled']:
            #     if ptMinLGR <= ptCent <= ptMaxLGR:
            #         hPtFONLLtimesLGRD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sLGR['yCent'](ptCent))
            #     elif ptCent > ptMaxLGR:
            #         hPtFONLLtimesLGRD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sLGR['yCent'](ptMaxLGR))
            #     else:
            #         hPtFONLLtimesLGRD[-1].SetBinContent(iPt, sFONLLD[pred](ptCent) * sLGR['yCent'](ptMinLGR))

            # if 'datashape' in shapesD and shapesD['datashape']['enabled']:
            #     # In this case we have to trust the data shape outside the fit range
            #     hPtDataShape[-1].SetBinContent(iPt, fData(ptCent))

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
        # if 'phsd' in shapesD and shapesD['phsd']['enabled']:
        #     hPtFONLLtimesPHSDD[-1].Scale(1./hPtFONLLtimesPHSDD[-1].Integral())
        #     hPtWeightsFONLLtimesPHSDD.append(
        #         hPtFONLLtimesPHSDD[-1].Clone(hPtFONLLtimesPHSDD[-1].GetName().replace('Pt', 'PtWeights')))
        #     hPtWeightsFONLLtimesPHSDD[-1].Divide(hPtFONLLtimesPHSDD[-1], hPtGenD)
        #     hPtWeightsFONLLtimesPHSDD[-1].Smooth(smooth)
        # if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
        #     hPtFONLLtimesGossiauxD[-1].Scale(1./hPtFONLLtimesGossiauxD[-1].Integral())
        #     hPtWeightsFONLLtimesGossiauxD.append(
        #         hPtFONLLtimesGossiauxD[-1].Clone(hPtFONLLtimesGossiauxD[-1].GetName().replace('Pt', 'PtWeights')))
        #     hPtWeightsFONLLtimesGossiauxD[-1].Divide(hPtFONLLtimesGossiauxD[-1], hPtGenD)
        #     hPtWeightsFONLLtimesGossiauxD[-1].Smooth(smooth)
        # if 'catania' in shapesD and shapesD['catania']['enabled']:
        #     hPtFONLLtimesCataniaD[-1].Scale(1./hPtFONLLtimesCataniaD[-1].Integral())
        #     hPtWeightsFONLLtimesCataniaD.append(
        #         hPtFONLLtimesCataniaD[-1].Clone(hPtFONLLtimesCataniaD[-1].GetName().replace('Pt', 'PtWeights')))
        #     hPtWeightsFONLLtimesCataniaD[-1].Divide(hPtFONLLtimesCataniaD[-1], hPtGenD)
        #     hPtWeightsFONLLtimesCataniaD[-1].Smooth(smooth)
        # if 'lido' in shapesD and shapesD['lido']['enabled']:
        #     hPtFONLLtimesLIDOD[-1].Scale(1./hPtFONLLtimesLIDOD[-1].Integral())
        #     hPtWeightsFONLLtimesLIDOD.append(
        #         hPtFONLLtimesLIDOD[-1].Clone(hPtFONLLtimesLIDOD[-1].GetName().replace('Pt', 'PtWeights')))
        #     hPtWeightsFONLLtimesLIDOD[-1].Divide(hPtFONLLtimesLIDOD[-1], hPtGenD)
        #     hPtWeightsFONLLtimesLIDOD[-1].Smooth(smooth)
        # if 'lgr' in shapesD and shapesD['lgr']['enabled']:
        #     hPtFONLLtimesLGRD[-1].Scale(1./hPtFONLLtimesLGRD[-1].Integral())
        #     hPtWeightsFONLLtimesLGRD.append(
        #         hPtFONLLtimesLGRD[-1].Clone(hPtFONLLtimesLGRD[-1].GetName().replace('Pt', 'PtWeights')))
        #     hPtWeightsFONLLtimesLGRD[-1].Divide(hPtFONLLtimesLGRD[-1], hPtGenD)
        #     hPtWeightsFONLLtimesLGRD[-1].Smooth(smooth)
        # if 'datashape' in shapesD and shapesD['datashape']['enabled']:
        #     hPtDataShape[-1].Scale(1./hPtDataShape[-1].Integral())
        #     hPtWeightsDataShape.append(
        #         hPtDataShape[-1].Clone(hPtDataShape[-1].GetName().replace('Pt', 'PtWeights')))
        #     hPtWeightsDataShape[-1].Divide(hPtDataShape[-1], hPtGenD)
        #     hPtWeightsDataShape[-1].Smooth(smooth)

    print('INFO: pT weights calculated')
    # B meson weights
    #___________________________________________________________________________________________________________________________
    if Bspecie:
        for histoName, pred in zip(histoBNames, modelPred):
            hPtFONLLB.append(hPtGenB.Clone(histoName))
            if 'tamu' in shapesB and shapesB['tamu']['enabled']:
                hPtFONLLtimesTAMUB.append(hPtGenB.Clone(histoName.replace('FONLL', 'FONLLtimesTAMU')))
            # if 'lido' in shapesB and shapesB['lido']['enabled']:
            #     hPtFONLLtimesLIDOB.append(hPtGenB.Clone(histoName.replace('FONLL', 'FONLLtimesLIDO')))

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
                        if ptMinTAMUB <= ptCent <= ptMaxTAMUB:
                            hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptCent))
                        elif ptCent > ptMaxTAMUB:
                            hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptMaxTAMUB))
                        else:
                            hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sTAMUB['yCent'](ptMinTAMUB))
                    else:
                        ptMaxMix = min([ptMaxTAMUB, ptMaxTAMUBs])
                        ptMinMix = max([ptMinTAMUB, ptMinTAMUBs])
                        if ptMinMix <= ptCent <= ptMaxMix:
                            rAAMix = (sTAMUB['yCent'](ptCent) + sTAMUBs['yCent'](ptCent)) / 2
                        elif ptCent > ptMaxMix:
                            rAAMix = (sTAMUB['yCent'](ptMaxMix) + sTAMUBs['yCent'](ptMaxMix)) / 2
                        else:
                            rAAMix = (sTAMUB['yCent'](ptMinMix) + sTAMUBs['yCent'](ptMinMix)) / 2
                        hPtFONLLtimesTAMUB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * rAAMix)
                # if 'lido' in shapesB and shapesB['lido']['enabled']:
                #     if ptMinLIDOB <= ptCent <= ptMaxLIDOB:
                #         hPtFONLLtimesLIDOB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sLIDOB['yCent'](ptCent))
                #     elif ptCent > ptMaxLIDOB:
                #         hPtFONLLtimesLIDOB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sLIDOB['yCent'](ptMaxLIDOB))
                #     else:
                #         hPtFONLLtimesLIDOB[-1].SetBinContent(iPt, sFONLLB[pred](ptCent) * sLIDOB['yCent'](ptMinLIDOB))

            hPtFONLLB[-1].Sumw2()
            hPtFONLLB[-1].Scale(1./hPtFONLLB[-1].Integral())
            hPtWeightsFONLLB.append(hPtFONLLB[-1].Clone(histoName.replace('Pt', 'PtWeights')))
            hPtWeightsFONLLB[-1].Divide(hPtFONLLB[-1], hPtGenB)
            hPtWeightsFONLLB[-1].Smooth(smooth)
            if 'tamu' in shapesB and shapesB['tamu']['enabled']:
                hPtFONLLtimesTAMUB[-1].Scale(1./hPtFONLLtimesTAMUB[-1].Integral())
                hPtWeightsFONLLtimesTAMUB.append(
                    hPtFONLLtimesTAMUB[-1].Clone(hPtFONLLtimesTAMUB[-1].GetName().replace('Pt', 'PtWeights')))
                hPtWeightsFONLLtimesTAMUB[-1].Divide(hPtFONLLtimesTAMUB[-1], hPtGenB)
                hPtWeightsFONLLtimesTAMUB[-1].Smooth(smooth)
            # if 'lido' in shapesB and shapesB['lido']['enabled']:
            #     hPtFONLLtimesLIDOB[-1].Scale(1./hPtFONLLtimesLIDOB[-1].Integral())
            #     hPtWeightsFONLLtimesLIDOB.append(
            #         hPtFONLLtimesLIDOB[-1].Clone(hPtFONLLtimesLIDOB[-1].GetName().replace('Pt', 'PtWeights')))
            #     hPtWeightsFONLLtimesLIDOB[-1].Divide(hPtFONLLtimesLIDOB[-1], hPtGenB)
            #     hPtWeightsFONLLtimesLIDOB[-1].Smooth(smooth)

    print('INFO: pT weights calculated')
    # save output
    #___________________________________________________________________________________________________________________________
    os.makedirs(f'{outputDir}/ptweights', exist_ok=True)
    outfile = TFile(f'{outputDir}/ptweights/pTweight_{suffix}.root', 'recreate')
    hPtGenD.Write()
    if Bspecie:
        hPtGenB.Write()
    for iHisto, _ in enumerate(hPtFONLLD):
        hPtFONLLD[iHisto].Write()
        hPtWeightsFONLLD[iHisto].Write()
        if 'tamu' in shapesD and shapesD['tamu']['enabled']:
            hPtFONLLtimesTAMUD[iHisto].Write()
            hPtWeightsFONLLtimesTAMUD[iHisto].Write()
        # if 'phsd' in shapesD and shapesD['phsd']['enabled']:
        #     hPtFONLLtimesPHSDD[iHisto].Write()
        #     hPtWeightsFONLLtimesPHSDD[iHisto].Write()
        # if 'mc@shq' in shapesD and shapesD['mc@shq']['enabled']:
        #     hPtFONLLtimesGossiauxD[iHisto].Write()
        #     hPtWeightsFONLLtimesGossiauxD[iHisto].Write()
        # if 'catania' in shapesD and shapesD['catania']['enabled']:
        #     hPtFONLLtimesCataniaD[iHisto].Write()
        #     hPtWeightsFONLLtimesCataniaD[iHisto].Write()
        # if 'lido' in shapesD and shapesD['lido']['enabled']:
        #     hPtFONLLtimesLIDOD[iHisto].Write()
        #     hPtWeightsFONLLtimesLIDOD[iHisto].Write()
        # if 'lgr' in shapesD and shapesD['lgr']['enabled']:
        #     hPtFONLLtimesLGRD[iHisto].Write()
        #     hPtWeightsFONLLtimesLGRD[iHisto].Write()
        # if 'datashape' in shapesD and shapesD['datashape']['enabled']:
        #     hPtDataShape[iHisto].Write()
        #     hPtWeightsDataShape[iHisto].Write()
        if Bspecie:
            hPtFONLLB[iHisto].Write()
            hPtWeightsFONLLB[iHisto].Write()
            if 'tamu' in shapesB and shapesB['tamu']['enabled']:
                hPtFONLLtimesTAMUB[iHisto].Write()
                hPtWeightsFONLLtimesTAMUB[iHisto].Write()
            # if 'lido' in shapesB and shapesB['lido']['enabled']:
            #     hPtFONLLtimesLIDOB[iHisto].Write()
            #     hPtWeightsFONLLtimesLIDOB[iHisto].Write()

    # pT shape D
    #___________________________________________________________________________________________________________________________
    canvPtshape = TCanvas('pTshape', 'pTshape', 2000, 900)
    canvPtshape.Divide(2, 1)
    ptD = [0, 36]
    canvPtshape.cd(1).DrawFrame(ptD[0], 0.0000001, ptD[1], 1,
                        f';#it{{p}}_{{T}} (GeV/c);Prompt {config["Dmeson"]}')
    canvPtshape.cd(1)
    canvPtshape.cd(1).SetLogy()

    leg = TLegend(0.5, 0.63, 0.7, 0.83)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.04)

    SetObjectStyle(hPtGenD, color=kRed, markersize=0.5)
    leg.AddEntry(hPtGenD, 'Gen Prompt', 'lp')
    SetObjectStyle(hPtFONLLtimesTAMUD[0], color=kBlack, markersize=0.5)
    leg.AddEntry(hPtFONLLtimesTAMUD[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtFONLLD[0], color=kAzure, markersize=0.5)
    #leg.AddEntry(hPtFONLLD[0], 'FONLL', 'lp')

    hPtFONLLtimesTAMUD[0].Draw('same')
    #hPtFONLLD[0].Draw('same')
    hPtGenD.Draw('same')
    leg.Draw()

    canvPtshape.cd(2).DrawFrame(0, 0., ptD[1], 5,
                        f';#it{{p}}_{{T}} (GeV/c);Ratio')
    canvPtshape.cd(2)

    legR = TLegend(0.5, 0.63, 0.7, 0.83)
    legR.SetFillStyle(0)
    legR.SetBorderSize(0)
    legR.SetTextSize(0.04)

    SetObjectStyle(hPtWeightsFONLLtimesTAMUD[0], color=kBlack, markersize=1)
    legR.AddEntry(hPtWeightsFONLLtimesTAMUD[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtWeightsFONLLD[0], color=kAzure, markersize=1)
    #legR.AddEntry(hPtWeightsFONLLD[0], 'FONLL', 'lp')

    hPtWeightsFONLLtimesTAMUD[0].Draw('same')
    #hPtWeightsFONLLD[0].Draw('same')
    legR.Draw()

    canvPtshape.Write()
    canvPtshape.SaveAs(f'{outputDir}/ptweights/pTweight.png')
    canvPtshape.SaveAs(f'{outputDir}/ptweights/pTweight.pdf')

    # pT shape B
    #___________________________________________________________________________________________________________________________
    canvPtshapeB = TCanvas('pTshapeB', 'pTshapeB', 2000, 900)
    canvPtshapeB.Divide(2, 1)
    ptB = [0, 72]
    canvPtshapeB.cd(1).DrawFrame(0, 0.0000001, ptB[1], 1,
                        f';#it{{p}}_{{T}}^{{B}} (GeV/c); B hadron')
    canvPtshapeB.cd(1)
    canvPtshapeB.cd(1).SetLogy()

    legB = TLegend(0.5, 0.63, 0.7, 0.83)
    legB.SetFillStyle(0)
    legB.SetBorderSize(0)
    legB.SetTextSize(0.04)

    SetObjectStyle(hPtGenB, color=kRed, markersize=0.5)
    legB.AddEntry(hPtGenB, 'Gen B', 'lp')
    SetObjectStyle(hPtFONLLtimesTAMUB[0], color=kBlack, markersize=0.5)
    legB.AddEntry(hPtFONLLtimesTAMUB[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtFONLLB[0], color=kAzure, markersize=0.5)
    #legB.AddEntry(hPtFONLLB[0], 'FONLL', 'lp')

    hPtFONLLtimesTAMUB[0].Draw('same')
    #hPtFONLLB[0].Draw('same')
    hPtGenB.Draw('same')
    legB.Draw()

    canvPtshapeB.cd(2).DrawFrame(0, 0., ptB[1], 10,
                        f';#it{{p}}_{{T}}^{{B}} (GeV/c);Ratio')
    canvPtshapeB.cd(2)

    legBR = TLegend(0.5, 0.63, 0.7, 0.83)
    legBR.SetFillStyle(0)
    legBR.SetBorderSize(0)
    legBR.SetTextSize(0.04)

    SetObjectStyle(hPtWeightsFONLLtimesTAMUB[0], color=kBlack, markersize=1)
    legBR.AddEntry(hPtWeightsFONLLtimesTAMUB[0], 'FONLL #times TAMU (R_{AA})', 'lp')
    SetObjectStyle(hPtWeightsFONLLB[0], color=kAzure, markersize=1)
    #legBR.AddEntry(hPtWeightsFONLLB[0], 'FONLL', 'lp')

    hPtWeightsFONLLtimesTAMUB[0].Draw('same')
    #hPtWeightsFONLLB[0].Draw('same')
    legBR.Draw()

    canvPtshapeB.Write()
    canvPtshapeB.SaveAs(f'{outputDir}/ptweights/pTweightB.png')
    canvPtshapeB.SaveAs(f'{outputDir}/ptweights/pTweightB.pdf')

    outfile.Close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Arguments')
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="flow configuration file")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    args = parser.parse_args()

    computePtWeights(args.config, args.outputdir, args.suffix)