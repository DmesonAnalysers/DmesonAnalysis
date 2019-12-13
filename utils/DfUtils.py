import numpy as np
import uproot


def WriteTree(df, cols, treeName, fileName):
    '''
    Helper method to write a tree with uproot

    Arguments
    ----------

    - pandas data frame to be written as tree in a root file
    - name of the columns
    - name of the output tree
    - name of the output file
    '''
    outBranches = {}
    for colName in cols:
        outBranches[colName] = np.float32 #define all branches as float for the moment
    with uproot.recreate(fileName, compression=uproot.LZ4(4)) as outFile:
        outFile[treeName] = uproot.newtree(outBranches, compression=uproot.LZ4(4))
        outFile[treeName].extend(dict(df[cols]))


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
    flags = dfToFilter[column] & maskOfBits
    if logic == 'or':
        flags = flags.astype('bool')
    elif logic == 'and':
        flags -= maskOfBits
        flags = ~flags.astype('bool')
    else:
        print('Error: only and and or logics are supported for bitwise operations')
        return None

    dfFilt = dfToFilter[flags.values]

    return dfFilt
