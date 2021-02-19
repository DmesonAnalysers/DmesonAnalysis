'''
Script with utils methods for managment and operations on pandas dataframes
'''
import pandas as pd
import uproot

def GetMaskOfBits(bits):
    '''
    Helper method to get bit mask from bits

    Arguments
    ----------
    - list of bits

    Returns
    ----------
    - mask corresponding to the input bits
    '''
    mask = 0
    for bit in bits:
        mask += 2**bit

    return mask


def FilterBitDf(dfToFilter, column, bitsToTest, logic='or'):
    '''
    Method to apply selection testing one or more bits

    Arguments
    ----------
    - pandas dataframe to filter
    - colum with bitmap
    - list of bits to test
    - logic to combine the bits (and, or)

    Returns
    ----------
    - filtered pandas dataframe
    '''
    maskOfBits = GetMaskOfBits(bitsToTest)
    flags = dfToFilter[column].astype(int) & maskOfBits
    if logic == 'or':
        flags = flags.astype('bool')
    elif logic == 'and':
        flags -= maskOfBits
        flags = ~flags.astype('bool')
    elif logic == 'not':
        flags = ~flags.astype('bool')
    else:
        print('Error: only and, or, and not logics are supported for bitwise operations')
        return None

    dfFilt = dfToFilter[flags.to_numpy()]

    return dfFilt


def LoadDfFromRootOrParquet(inFileNames, inDirNames=None, inTreeNames=None):
    '''
    Helper method to load a pandas dataframe from either root or parquet files

    Arguments
    ----------
    - input file name of list of input file names
    - input dir name of list of input dir names (needed only in case of root files)

    Returns
    ----------
    - loaded pandas dataframe
    '''

    if not isinstance(inFileNames, list):
        inFileNames = [inFileNames]
    if not isinstance(inDirNames, list):
        inDirName = inDirNames
        inDirNames = [inDirName] * len(inFileNames)
    if not isinstance(inTreeNames, list):
        inTreeName = inTreeNames
        inTreeNames = [inTreeName] * len(inFileNames)
    dfOut = pd.DataFrame()

    for inFile, inDir, inTree in zip(inFileNames, inDirNames, inTreeNames):
        if '.root' in inFile:
            path = f'{inFile}:{inDir}/{inTree}' if inDir else f'{inFile}:{inTree}'
            dfOut = dfOut.append(uproot.open(path).arrays(library='pd'), ignore_index=True)
        elif '.parquet' in inFile:
            dfOut = dfOut.append(pd.read_parquet(inFile), ignore_index=True)
        else:
            print('ERROR: only root or parquet files are supported! Returning empty dataframe')
            return pd.DataFrame()

    return dfOut


def GetMind0(ptList, d0List, ptThrs):
    '''
    Helper method to get minimum impact parameter for given pt threshold as in AOD filtering

    Arguments
    ----------
    - list of pt of daughter tracks
    - list of impact parameters of daughter tracks
    - pt threshold (selected less than threshold)

    Returns
    ----------
    - minimum impact parameter for given pt threshold
    '''
    d0SelList = []
    for d0, pt in zip(d0List, ptList):
        if pt < ptThrs:
            d0SelList.append(abs(d0))

    if len(d0SelList) == 0: # no daughters with pt < pt threshold
        return 999.

    return min(d0SelList)
