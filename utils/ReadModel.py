'''
Script with helper functions to load model predictions from txt files
'''

import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline

def InterpolateModel(ptCent, yCent, yMin=None, yMax=None):
    '''
    Helper function to interpolate model predictions.
    The returned splines will raise an error if applied out of the data boundary.

    Parameters
    -----------
    ptCent: list of pT centres to interpolate
    yCent: list of central values to interpolate
    yMin: list of min values to interpolate
    yMax: list of max values to interpolate

    Returns:
    -----------
    splinesAll: dictionary with splines {yCent, yMin, yMax}
    ptMin: minimum pt for which the interpolation is valid
    ptMax: maximum pt for which the interpolation is valid
    '''

    splinesAll = {}
    splinesAll['yCent'] = InterpolatedUnivariateSpline(ptCent, yCent, ext='raise', check_finite=True)

    if yMin is not None and yMin.any():
        splinesAll['yMin'] = InterpolatedUnivariateSpline(ptCent, yMin, ext='raise', check_finite=True)
    if yMax is not None and yMax.any():
        splinesAll['yMax'] = InterpolatedUnivariateSpline(ptCent, yMax, ext='raise', check_finite=True)

    return splinesAll, min(ptCent), max(ptCent)


def ReadFONLL(fileNameFONLL, isPtDiff=False):
    '''
    Helper function to read FONLL txt files

    Parameters
    -----------
    fileNameFONLL: FONLL file name

    Returns:
    -----------
    splineFONLL: dictionary with splines {yCent, yMin, yMax}
    dfFONLL: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfFONLL = pd.read_csv(fileNameFONLL, sep=" ", header=13).astype('float64')
    if not isPtDiff:
        dfFONLL['ptcent'] = (dfFONLL['ptmin']+dfFONLL['ptmax']) / 2
        dfFONLL['central_ptdiff'] = dfFONLL['central'] / (dfFONLL['ptmax']-dfFONLL['ptmin'])
        dfFONLL['min_ptdiff'] = dfFONLL['min'] / (dfFONLL['ptmax']-dfFONLL['ptmin'])
        dfFONLL['max_ptdiff'] = dfFONLL['max'] / (dfFONLL['ptmax']-dfFONLL['ptmin'])
    else:
        dfFONLL.rename(columns={'pt': 'ptcent'}, inplace=True)
        dfFONLL.rename(columns={'central': 'central_ptdiff'}, inplace=True)
        dfFONLL.rename(columns={'min': 'min_ptdiff'}, inplace=True)
        dfFONLL.rename(columns={'max': 'max_ptdiff'}, inplace=True)

    splineFONLL, ptMin, ptMax = InterpolateModel(dfFONLL['ptcent'], dfFONLL['central_ptdiff'],
                                                 dfFONLL['min_ptdiff'], dfFONLL['max_ptdiff'])

    return splineFONLL, dfFONLL, ptMin, ptMax


def ReadGMVFNS(fileNameGMVFNS, isSACOT=False):
    '''
    Helper function to read GVMFS txt files

    Parameters
    -----------
    fileNameGMVFNS: GVMFS file name
    isSACOT: flag to tell the function if it is the SACOT-mT version of GVMFS

    Returns:
    -----------
    splineGMVFNS: dictionary with splines {yCent, yMin, yMax}
    dfGMVFNS: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfGMVFNS = pd.read_csv(fileNameGMVFNS, sep=' ', comment='#')
    if not isSACOT:
        splineGMVFNS, ptMin, ptMax = InterpolateModel(dfGMVFNS['pT'], dfGMVFNS['cen'],
                                                      dfGMVFNS['min'], dfGMVFNS['max'])
    else:
        dfGMVFNS['xsec_min[mb]'] = dfGMVFNS['xsec[mb]'] - \
            np.sqrt(dfGMVFNS['PDFerr[mb]']**2 + dfGMVFNS['down.scale.err[mb]']**2)
        dfGMVFNS['xsec_max[mb]'] = dfGMVFNS['xsec[mb]'] + \
            np.sqrt(dfGMVFNS['PDFerr[mb]']**2 + dfGMVFNS['up.scale.err[mb]']**2)

        splineGMVFNS, ptMin, ptMax = InterpolateModel(dfGMVFNS['pT'], dfGMVFNS['xsec[mb]'],
                                                      dfGMVFNS['xsec_min[mb]'], dfGMVFNS['xsec_max[mb]'])

    return splineGMVFNS, dfGMVFNS, ptMin, ptMax


def ReadKtFact(fileNameKtFact):
    '''
    Helper function to read kT-factorisation txt files

    Parameters
    -----------
    fileNameKtFact: kT-factorisation file name

    Returns:
    -----------
    splineKtFact: dictionary with splines {yCent, yMin, yMax}
    dfKtFact: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfKtFact = pd.read_csv(fileNameKtFact, sep=' ', comment='#')
    dfKtFact['ptcent'] = (dfKtFact['ptmin']+dfKtFact['ptmax']) / 2
    splineKtFact, ptMin, ptMax = InterpolateModel(dfKtFact['ptcent'], dfKtFact['central'],
                                                  dfKtFact['lower'], dfKtFact['upper'])

    return splineKtFact, dfKtFact, ptMin, ptMax


def ReadTAMU(fileNameTAMU):
    '''
    Helper function to read TAMU txt files

    Parameters
    -----------
    fileNameTAMU: TAMU file name

    Returns:
    -----------
    splineTAMU: dictionary with splines {yCent, yMin, yMax}
    dfTAMU: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfTAMU = pd.read_csv(fileNameTAMU, sep=' ', comment='#')
    if 'R_AA_max' in dfTAMU and 'R_AA_min' in dfTAMU:
        dfTAMU['R_AA'] = (dfTAMU['R_AA_min'] + dfTAMU['R_AA_max']) / 2 #central value taken as average of min and max
        splineTAMU, ptMin, ptMax = InterpolateModel(dfTAMU['PtCent'], dfTAMU['R_AA'],
                                                    dfTAMU['R_AA_min'], dfTAMU['R_AA_max'])
    else:
        splineTAMU, ptMin, ptMax = InterpolateModel(dfTAMU['PtCent'], dfTAMU['R_AA'])

    return splineTAMU, dfTAMU, ptMin, ptMax


def ReadTAMUv2(fileNameTAMUv2):
    '''
    Helper function to read TAMU v2 txt files

    Parameters
    -----------
    fileNameTAMU: TAMU file name

    Returns:
    -----------
    splineTAMU: dictionary with splines {yCent, yMin, yMax}
    dfTAMU: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfTAMU = pd.read_csv(fileNameTAMUv2, sep=' ', comment='#')
    if 'v2min' in dfTAMU and 'v2max' in dfTAMU:
        dfTAMU['v2'] = (dfTAMU['v2min'] + dfTAMU['v2max']) / 2 #central value taken as average of min and max
        splineTAMU, ptMin, ptMax = InterpolateModel(dfTAMU['pT'], dfTAMU['v2'],
                                                    dfTAMU['v2min'], dfTAMU['v2max'])
    else:
        splineTAMU, ptMin, ptMax = InterpolateModel(dfTAMU['pT'], dfTAMU['v2'])

    return splineTAMU, dfTAMU, ptMin, ptMax


def ReadPHSD(fileNamePHSD):
    '''
    Helper function to read PHSD txt files

    Parameters
    -----------
    fileNamePHSD: PHSD file name

    Returns:
    -----------
    splinePHSD: dictionary with splines {yCent}
    dfPHSD: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfPHSD = pd.read_csv(fileNamePHSD, sep=' ', comment='#')
    splinePHSD, ptMin, ptMax = InterpolateModel(dfPHSD['pt'], dfPHSD['Raa'])

    return splinePHSD, dfPHSD, ptMin, ptMax


def ReadMCatsHQ(fileNameMCatsHQ):
    '''
    Helper function to read MCatsHQ txt files

    Parameters
    -----------
    fileNameMCatsHQ: MCatsHQ file name

    Returns:
    -----------
    splineMCatsHQ: dictionary with splines {yCent, yMin, yMax}
    dfMCatsHQ: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfMCatsHQ = pd.read_csv(fileNameMCatsHQ, sep=' ', comment='#')
    dfMCatsHQ['Raa_min'] = dfMCatsHQ[['RAAcolK1.5', 'RAAcolradLPMK0.8', 'RAAcolradLPMgludampK0.8']].min(axis=1)
    dfMCatsHQ['Raa_max'] = dfMCatsHQ[['RAAcolK1.5', 'RAAcolradLPMK0.8', 'RAAcolradLPMgludampK0.8']].max(axis=1)
    dfMCatsHQ['Raa_cent'] = (dfMCatsHQ['Raa_min'] + dfMCatsHQ['Raa_max']) / 2
    splineMCatsHQ, ptMin, ptMax = InterpolateModel(dfMCatsHQ['pt'], dfMCatsHQ['Raa_cent'],
                                                   dfMCatsHQ['Raa_min'], dfMCatsHQ['Raa_max'])

    return splineMCatsHQ, dfMCatsHQ, ptMin, ptMax


def ReadCatania(fileNameCatania):
    '''
    Helper function to read Catania txt files

    Parameters
    -----------
    fileNameCatania: Catania file name

    Returns:
    -----------
    splineCatania: dictionary with splines {yCent}
    dfCatania: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''
    dfCatania = pd.read_csv(fileNameCatania, sep=' ', comment='#')
    splineCatania, ptMin, ptMax = InterpolateModel(dfCatania['pt'], dfCatania['Raa'])

    return splineCatania, dfCatania, ptMin, ptMax


def ReadLIDO(fileName):
    '''
    Method to read LIDO Raa files

    Inputs
    ----------
    - fileName: file name
    - obs: observable (Raa or v2)

    Returns
    ----------
    splineLIDO: dictionary of splines with LIDO predictions {yCent, yMin, yMax}
    dfLIDO: pandas dataframe with original values
    ptMin: minimum pt for which the model is valid
    ptMax: maximum pt for which the model is valid
    '''

    dfLIDO = pd.read_csv(fileName, sep=' ')
    dfLIDO['Raa_min'] = dfLIDO['Raa'] - dfLIDO['Raa-error']
    dfLIDO['Raa_max'] = dfLIDO['Raa'] + dfLIDO['Raa-error']

    splineLIDO, ptMin, ptMax = InterpolateModel(dfLIDO['pT'], dfLIDO['Raa'], dfLIDO['Raa_min'], dfLIDO['Raa_max'])

    return splineLIDO, dfLIDO, ptMin, ptMax
