'''
Script with helper functions to load models from txt files
'''

import numpy as np
import pandas as pd
from scipy.interpolate import InterpolatedUnivariateSpline

def InterpolateModel(ptCent, yCent, yMin=None, yMax=None):
    '''
    Helper function to interpolate model

    Parameters
    -----------
    ptRange: list with min and max pt values for interpolation
    ptStep: pt step for interpolation
    ptCent: list of pT centres to interpolate
    yCent: list of central values to interpolate
    yMin: list of min values to interpolate
    yMax: list of max values to interpolate

    Returns:
    -----------
    splinesAll: dictionary with values of splines {yCent, yMin, yMax}
    '''

    splinesAll = {}
    splinesAll['yCent'] = InterpolatedUnivariateSpline(ptCent, yCent)

    if yMin is not None and yMin.any():
        splinesAll['yMin'] = InterpolatedUnivariateSpline(ptCent, yMin)
    if yMax is not None and yMax.any():
        splinesAll['yMax'] = InterpolatedUnivariateSpline(ptCent, yMax)

    return splinesAll


def ReadFONLL(fileNameFONLL):
    '''
    Helper function to read FONLL txt files

    Parameters
    -----------
    fileNameFONLL: FONLL file name

    Returns:
    -----------
    splineFONLL: dictionary with values of splines {yCent, yMin, yMax}
    dfFONLL: pandas dataframe with original values
    '''
    dfFONLL = pd.read_csv(fileNameFONLL, sep=' ', comment='#')
    dfFONLL['ptcent'] = (dfFONLL['ptmin']+dfFONLL['ptmax']) / 2
    dfFONLL['central_ptdiff'] = dfFONLL['central'] / (dfFONLL['ptmax']-dfFONLL['ptmin'])
    dfFONLL['min_ptdiff'] = dfFONLL['min'] / (dfFONLL['ptmax']-dfFONLL['ptmin'])
    dfFONLL['max_ptdiff'] = dfFONLL['max'] / (dfFONLL['ptmax']-dfFONLL['ptmin'])

    splineFONLL = InterpolateModel(dfFONLL['ptcent'], dfFONLL['central_ptdiff'],
                                   dfFONLL['min_ptdiff'], dfFONLL['max_ptdiff'])

    return splineFONLL, dfFONLL


def ReadGMVFNS(fileNameGMVFNS, isSACOT=False):
    '''
    Helper function to read GVMFS txt files

    Parameters
    -----------
    fileNameGMVFNS: GVMFS file name
    isSACOT: flag to tell the function if it is the SACOT-mT version of GVMFS

    Returns:
    -----------
    splineGMVFNS: dictionary with values of splines {yCent, yMin, yMax}
    dfGMVFNS: pandas dataframe with original values
    '''
    dfGMVFNS = pd.read_csv(fileNameGMVFNS, sep=' ', comment='#')
    if not isSACOT:
        splineGMVFNS = InterpolateModel(dfGMVFNS['pT'], dfGMVFNS['cen'],
                                        dfGMVFNS['min'], dfGMVFNS['max'])
    else:
        dfGMVFNS['xsec_min[mb]'] = dfGMVFNS['xsec[mb]'] - \
            np.sqrt(dfGMVFNS['PDFerr[mb]']**2 + dfGMVFNS['down.scale.err[mb]']**2)
        dfGMVFNS['xsec_max[mb]'] = dfGMVFNS['xsec[mb]'] + \
            np.sqrt(dfGMVFNS['PDFerr[mb]']**2 + dfGMVFNS['up.scale.err[mb]']**2)

        splineGMVFNS = InterpolateModel(dfGMVFNS['pT'], dfGMVFNS['xsec[mb]'],
                                        dfGMVFNS['xsec_min[mb]'], dfGMVFNS['xsec_max[mb]'])

    return splineGMVFNS, dfGMVFNS


def ReadKtFact(fileNameKtFact):
    '''
    Helper function to read kT-factorisation txt files

    Parameters
    -----------
    fileNameKtFact: kT-factorisation file name

    Returns:
    -----------
    splineKtFact: dictionary with values of splines {yCent, yMin, yMax}
    dfKtFact: pandas dataframe with original values
    '''
    dfKtFact = pd.read_csv(fileNameKtFact, sep=' ', comment='#')
    dfKtFact['ptcent'] = (dfKtFact['ptmin']+dfKtFact['ptmax']) / 2
    splineKtFact = InterpolateModel(dfKtFact['ptcent'], dfKtFact['central'],
                                    dfKtFact['lower'], dfKtFact['upper'])

    return splineKtFact, dfKtFact


def ReadTAMU(fileNameTAMU):
    '''
    Helper function to read TAMU txt files

    Parameters
    -----------
    fileNameTAMU: TAMU file name

    Returns:
    -----------
    splineTAMU: dictionary with values of splines {yCent, yMin, yMax}
    dfTAMU: pandas dataframe with original values
    '''
    dfTAMU = pd.read_csv(fileNameTAMU, sep=' ', comment='#')
    dfTAMU['Cent'] = (dfTAMU['Min'] + dfTAMU['Max']) / 2 #central value taken as average of min and max
    splineTAMU = InterpolateModel(dfTAMU['PtCent'], dfTAMU['Cent'],
                                  dfTAMU['Min'], dfTAMU['Max'])

    return splineTAMU, dfTAMU


def ReadPHSD(fileNamePHSD):
    '''
    Helper function to read PHSD txt files

    Parameters
    -----------
    fileNamePHSD: PHSD file name

    Returns:
    -----------
    splinePHSD: dictionary with values of splines {yCent}
    dfPHSD: pandas dataframe with original values
    '''
    dfPHSD = pd.read_csv(fileNamePHSD, sep=' ', comment='#')
    splinePHSD = InterpolateModel(dfPHSD['pt'], dfPHSD['Raa'])

    return splinePHSD, dfPHSD


def ReadMCatsHQ(fileNameMCatsHQ):
    '''
    Helper function to read MCatsHQ txt files

    Parameters
    -----------
    fileNameMCatsHQ: MCatsHQ file name

    Returns:
    -----------
    splineMCatsHQ: dictionary with values of splines {yCent, yMin, yMax}
    dfMCatsHQ: pandas dataframe with original values
    '''
    dfMCatsHQ = pd.read_csv(fileNameMCatsHQ, sep=' ', comment='#')
    dfMCatsHQ['Raa_min'] = dfMCatsHQ[['RAAcolK1.5', 'RAAcolradLPMK0.8', 'RAAcolradLPMgludampK0.8']].min(axis=1)
    dfMCatsHQ['Raa_max'] = dfMCatsHQ[['RAAcolK1.5', 'RAAcolradLPMK0.8', 'RAAcolradLPMgludampK0.8']].max(axis=1)
    dfMCatsHQ['Raa_cent'] = (dfMCatsHQ['Raa_min'] + dfMCatsHQ['Raa_max']) / 2
    splineMCatsHQ = InterpolateModel(dfMCatsHQ['pt'], dfMCatsHQ['Raa_cent'],
                                     dfMCatsHQ['Raa_min'], dfMCatsHQ['Raa_max'])

    return splineMCatsHQ, dfMCatsHQ


def ReadCatania(fileNameCatania):
    '''
    Helper function to read Catania txt files

    Parameters
    -----------
    fileNameCatania: Catania file name

    Returns:
    -----------
    splineCatania: dictionary with values of splines {yCent}
    dfCatania: pandas dataframe with original values
    '''
    dfCatania = pd.read_csv(fileNameCatania, sep=' ', comment='#')
    splineCatania = InterpolateModel(dfCatania['pt'], dfCatania['Raa'])

    return splineCatania, dfCatania
