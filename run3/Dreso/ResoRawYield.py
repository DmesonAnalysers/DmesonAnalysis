'''
python script to upload tree from resonance task output and extract raw yield
'''
import argparse
import sys
import numpy as np
import yaml
sys.path.append('../../utils')
from DfUtils import LoadDfFromRootOrParquet #pylint: disable=wrong-import-position,import-error
import os
import pandas as pd
import uproot
import ROOT
import matplotlib.pyplot as plt
from ResoFitUtils import fit_invariant_mass_with_background, roofit_plot_with_matplotlib
from particle import Particle


#Utility to merge trees from AO2D and save it as pd.dataframe
def LoadTreeFromAO2D(inFileNames, inDirName, inTreeName):
    dataframe = pd.DataFrame()
    if not isinstance(inFileNames, list):
            inFileNames = [inFileNames]
    if not inDirName:
        for file in inFileNames:
            print(f'Loading trees form {file}')
            f = uproot.open(file)
            dirnames = f.keys(recursive = False)
            print(f"keys found in file: {dirnames}")
            for dirname in dirnames:
                if 'parent' in dirname:
                    continue
                else:
                    dataframe = pd.concat([dataframe, LoadDfFromRootOrParquet(file, dirname, inTreeName)], ignore_index=True)
                    print(f"processed key: {dirname}")
    else:  
        dataframe = LoadDfFromRootOrParquet(inFileNames, inDirName, inTreeName)
    return dataframe

# utilities for mass cut
def shiftedMass(mass, pT):
    return mass -0.0025 + 0.0001*pT
def sigma_par(pT):
    return 0.00796 + 0.00176*pT

#Parametrized pt-Differential cut on D daughter invariant mass (parameters from HFFilterHelpers.h)
def ApplyDMassCut(df, cfg):
    ptLabel = cfg['cutvarsDDaughter']['pt']['namevar']
    massLabel = cfg['cutvarsDDaughter']['invMass']['namevar']
    nsigma = cfg['cutvarsDDaughter']['invMass']['nsigma']
    mass = 1.86966
    condition = np.where(df[ptLabel]< 10, 
                         np.abs(shiftedMass(mass, df[ptLabel]) - df[massLabel]) < nsigma * sigma_par(df[ptLabel]),
                         np.abs(mass - df[massLabel]) < nsigma * sigma_par(df[ptLabel]))
    filteredDf = df[condition]
    return filteredDf

#Parametrized pt-Differential cut on D daughter invariant mass (parameters from HFFilterHelpers.h)
def ApplyV0MassCut(df, cfg):
    
    massLabel = cfg['cutvarsV0Daughter']['invMass']['namevar']
    massmin = cfg['cutvarsV0Daughter']['invMass']['min']
    massmax = cfg['cutvarsV0Daughter']['invMass']['max']
    condition = (df[massLabel] > massmin) &  (df[massLabel] < massmax)
    filteredDf = df[condition]
    return filteredDf

def CheckDaughtersDist(df, cfg, outdir):
    ptBinsMin = cfg['cutvarsDDaughter']['pt']['min']
    ptBinsMax = cfg['cutvarsDDaughter']['pt']['max']
    massD = df['fInvMassProng0']
    ptD = df['fPtProng0']
    massV0 = df['fInvMassProng1']
    ptV0 = df['fPtProng1']
    plotLimsD = [cfg['cutvarsDDaughter']['invMass']['min'],cfg['cutvarsDDaughter']['invMass']['max']]
    plotLimsV0 = [cfg['cutvarsV0Daughter']['invMass']['min'],cfg['cutvarsV0Daughter']['invMass']['max']]
    nBins = cfg['cutvarsDDaughter']['invMass']['nbins']
    nsigma = cfg['cutvarsDDaughter']['invMass']['nsigma']
    # mass = 1.86966
    mass = Particle.from_pdgid(411).mass*1e-3
    for ptMin, ptMax in zip(ptBinsMin,ptBinsMax):
        selectedDMass = massD[(ptD>ptMin) & (ptD<ptMax)]
        selectedV0Mass = massV0[(ptV0>ptMin) & (ptV0<ptMax)]
        histD, bin_edgesD = np.histogram(selectedDMass, bins=nBins, range=plotLimsD)
        binctrsD = [(bin_edgesD[i] + bin_edgesD[i+1]) / 2 for i in range(len(bin_edgesD) -1)]
        histV0, bin_edgesV0 = np.histogram(selectedV0Mass, bins=nBins, range=plotLimsV0)
        binctrsV0 = [(bin_edgesV0[i] + bin_edgesV0[i+1]) / 2 for i in range(len(bin_edgesV0) -1)]
        fig, ax = plt.subplots(ncols=2, figsize=(16,10))
        ax[0].errorbar(binctrsD, histD, yerr=np.sqrt(histD), fmt='o', label=f'D Invariant mass Spectrum {ptMin}< pT < {ptMax}', markersize=5, color='black')
        ax[0].set_title(f'D Invariant mass Spectrum {ptMin}< pT < {ptMax}', fontsize=15)
        ax[0].set_xlabel('invariant mass GeV/c²', fontsize=15)
        ax[0].set_ylabel(f'counts per {np.round((plotLimsD[1] - plotLimsD[0])/nBins * 1000)} MeV/c²', fontsize=15)
        ax[0].axvline(x=shiftedMass(mass, (ptMax + ptMin)/2.) - nsigma * sigma_par((ptMax + ptMin)/2.) , color='r', linestyle='--')
        ax[0].axvline(x=shiftedMass(mass, (ptMax + ptMin)/2.) + nsigma * sigma_par((ptMax + ptMin)/2.) , color='r', linestyle='--')
        ax[0].axvline(x=shiftedMass(mass, (ptMax + ptMin)/2.) - 4 * sigma_par((ptMax + ptMin)/2.) , color='g', linestyle='--')
        ax[0].axvline(x=shiftedMass(mass, (ptMax + ptMin)/2.) + 4 * sigma_par((ptMax + ptMin)/2.) , color='g', linestyle='--')
        ax[0].axvline(x=shiftedMass(mass, (ptMax + ptMin)/2.) - 6 * sigma_par((ptMax + ptMin)/2.) , color='g', linestyle='--')
        ax[0].axvline(x=shiftedMass(mass, (ptMax + ptMin)/2.) + 6 * sigma_par((ptMax + ptMin)/2.) , color='g', linestyle='--')
        ax[1].errorbar(binctrsV0, histV0, yerr=np.sqrt(histV0), fmt='o', label=f'V0 Invariant mass Spectrum {ptMin}< pT < {ptMax}', markersize=5, color='black')
        ax[1].set_title(f'V0 Invariant mass Spectrum {ptMin}< pT < {ptMax}', fontsize=15)
        ax[1].set_xlabel('invariant mass GeV/c²', fontsize=15)
        ax[1].set_ylabel(f'counts per {np.round((plotLimsV0[1] - plotLimsV0[0])/nBins * 1000)} MeV/c²', fontsize=15)
        fig.savefig(f'{outdir}/DaughtersInvMassSpectrum{ptMin}<pT<{ptMax}.png')
        ax[0].cla()
        ax[1].cla()

# Function to create a template for the background by mixing V0 and D+ from different events
def BkgFormSidebands(df, cfg):
    nsig_min = cfg['bkg']['nsig_min']
    nsig_max = cfg['bkg']['nsig_max']
    mass = Particle.from_pdgid(411).mass*1e-3
    ptLabel = cfg['cutvarsDDaughter']['pt']['namevar']
    massLabel = cfg['cutvarsDDaughter']['invMass']['namevar']
    condition = np.where(df[ptLabel]< 10, 
                         (np.abs(shiftedMass(mass, df[ptLabel]) - df[massLabel]) > nsig_min * sigma_par(df[ptLabel])) & (np.abs(shiftedMass(mass, df[ptLabel]) - df[massLabel]) < nsig_max * sigma_par(df[ptLabel])),
                         (np.abs(mass - df[massLabel]) > nsig_min * sigma_par(df[ptLabel])) & (np.abs(mass - df[massLabel]) < nsig_max * sigma_par(df[ptLabel])))
    filteredDf = df[condition]
    return filteredDf



# main
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('configfile', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
    args = parser.parse_args()
    # Loading configurations for data handling
    with open(args.configfile, 'r') as ymlCfgFile:
        cfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    inFileNames = cfg['infile']['filename']
    inDirNames = cfg['infile']['dirname']
    inTreeName = cfg['infile']['treename']
    # Loading and Filtering Dataframe
    df = LoadTreeFromAO2D(inFileNames, inDirNames, inTreeName)
     # Creating output directory for plots
    outdir = cfg['outdir']
    if not os.path.exists(outdir):
            os.makedirs(outdir)
    # Check and apply invariant mass cuts on Daughters
    CheckDaughtersDist(df, cfg, outdir)
    dfMassCut = ApplyDMassCut(df, cfg)
    dfMassCut = ApplyV0MassCut(dfMassCut, cfg)
    # D pT cut
    dfMassCut = dfMassCut[dfMassCut[cfg['cutvarsDDaughter']['pt']['namevar']] > 1]
   
    # Creating Resonance invariant mass histogram and fitting
    ptReso = dfMassCut[cfg['resoHist']['pt']['namevar']]
    invMassReso = dfMassCut[cfg['resoHist']['invMass']['namevar']]
    minPtLims = cfg['resoHist']['pt']['min']
    maxPtLims = cfg['resoHist']['pt']['max']
    nBins = cfg['resoHist']['invMass']['nbins']
    plotLims = [cfg['resoHist']['invMass']['min'],cfg['resoHist']['invMass']['max']]
    fitLims = [cfg['resoFit']['min'],cfg['resoFit']['max']]
    pdgId = cfg['pdgReso']
    if pdgId == 666:
        nameRes = 'Xc3055'
    elif pdgId == 10433:
        nameRes = 'Ds1'
    elif pdgId == 435:
        nameRes = 'Ds2Star'
    else: 
        raise KeyError("Only pdgId 666, 10433, 435 are valid!")
    isDeltaM = cfg['deltaMass']
    if len(minPtLims) != len(maxPtLims):
        raise ValueError("pT bins bin minima and maxima have different length")
    for ptMin, ptMax in zip(minPtLims, maxPtLims):
        selectedMasses = invMassReso[(ptReso>ptMin) & (ptReso<ptMax)]
        fig, ax = plt.subplots(figsize=(16,10))
        # Don't Fit Xc2
        if pdgId != 666:
            th1 = ROOT.TH1F(f"{nameRes}InvMassSpectrum{ptMin}<pT<{ptMax}", "invariant mass histogram", nBins, plotLims[0], plotLims[1])
            for m in selectedMasses:
                th1.Fill(m)
            workspace, mass_frame = fit_invariant_mass_with_background(th1, fitLims[0], fitLims[1],pdgId, isDeltaM)
            results = roofit_plot_with_matplotlib(ax, th1, workspace, plotLims[0], plotLims[1], fitLims[0], fitLims[1], plot = True)
            ax.axvline(x=2.536 - 0.135 - 0.038*135/1870 , color='r', linestyle='--', label='Ds1 2536 - pi0')
            ax.axvline(x=2.569 , color='b', linestyle='--', label='Ds2* 2573')
            # ax.axvline(x=2.714 , color='g', linestyle='--', label='Ds1 2700')
            # ax.axvline(x=2.859 , color='orange', linestyle='--', label='Ds1 2860')
            ax.legend(loc='upper right', fontsize=15)
        else: 
            histReso, bin_edgesReso = np.histogram(selectedMasses, bins=nBins, range=plotLims)
            binctrsReso = [(bin_edgesReso[i] + bin_edgesReso[i+1]) / 2 for i in range(len(bin_edgesReso) -1)]
            ax.errorbar(binctrsReso, histReso, yerr=np.sqrt(histReso), fmt='o', label=f'Reso Invariant mass Spectrum {ptMin}< pT < {ptMax}', markersize=5, color='blue')
            ax.set_title(f'{nameRes} Reso Invariant mass Spectrum {ptMin}< pT < {ptMax}', fontsize=15)
            ax.set_xlabel('invariant mass GeV/c²', fontsize=15)
            ax.set_ylabel(f'counts per {np.round((plotLims[1] - plotLims[0])/nBins * 1000)} MeV/c²', fontsize=15)
            ax.axvline(x=3.055 , color='r', linestyle='--', label='Xc 3055')
            ax.axvline(x=3.080 , color='g', linestyle='--', label='Xc 3080')
            ax.legend(loc='upper right', fontsize=15)

        fig.savefig(f'{outdir}/{nameRes}InvMassSpectrum{ptMin}<pT<{ptMax}.png')
         # Temporary plot
        if cfg['bkg']['sidebands'] :
            SBdf = BkgFormSidebands(df, cfg)
            fig2, ax2 = plt.subplots(figsize=(16,10))
            ptLab = cfg['resoHist']['pt']['namevar']
            filteredDfPtCut = SBdf[(SBdf[cfg['cutvarsDDaughter']['pt']['namevar']] > 1) & (SBdf[ptLab] > 5) & (SBdf[ptLab] < 50)]#[SBdf[ptLab] > 5 & SBdf[ptLab] < 50] # &  
            invMassSB = filteredDfPtCut[cfg['resoHist']['invMass']['namevar']]
            histSideBands, bin_edgesSideBands = np.histogram(invMassSB, bins=nBins, range=plotLims)
            binctrsSideBands = [(bin_edgesSideBands[i] + bin_edgesSideBands[i+1]) / 2 for i in range(len(bin_edgesSideBands) -1)]
            ax2.errorbar(binctrsSideBands, histSideBands, yerr=np.sqrt(histSideBands), fmt='o', label=f'SideBands Invariant mass Spectrum {ptMin}< pT < {ptMax}', markersize=5, color='blue')
            ax2.set_title(f'{nameRes} Side Bands BKG Invariant mass Spectrum {ptMin}< pT < {ptMax}', fontsize=15)
            ax2.set_xlabel('invariant mass GeV/c²', fontsize=15)
            ax2.set_ylabel(f'counts per {np.round((plotLims[1] - plotLims[0])/nBins * 1000)} MeV/c²', fontsize=15)
            fig2.savefig(f'{outdir}/{nameRes}SidebandsBKGSpectrum{5}<pT<{50}.png')

        print(dfMassCut)
    