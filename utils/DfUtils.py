'''
Script with utils methods for managment and operations on pandas dataframes
'''
import pandas as pd
import os
import sys
import uproot
from ROOT import TFile, TTree, TList, TDirectoryFile
import numpy as np
from alive_progress import alive_bar

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
    - logic to combine the bits (and, or, not)

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

def GetObjectFromFile(inFile, pathToObj):
    '''
    Function to extract an object inside a root file.
    Supports nested containers with the following Data Types:
     - TFile
     - TList

    Parameters
    -----------
    inFile: TFile or name of the input file
    pathToObj: path of the object inside the root file

    Returns:
    -----------
    outObj: target root object
    '''

    pathToObj = os.path.normpath(pathToObj)
    pathElements = pathToObj.split(os.sep)

    outObjOld = None
    if isinstance(inFile, str):
        outObjOld = TFile.Open(inFile, 'read')

    # if isinstance(inFile, str):
    #     print('lallero', outObjOld, type(outObjOld))

    #     if outObjOld == None:
    #         print('\033[31mError\033[0m: File does not exist. Exit!')
    #         sys.exit()
    elif isinstance(inFile, TFile):
        outObjOld = inFile
        # pass
    else:
        print('\033[31mError\033[0m: input file must be TFile or str. Exit!')
        sys.exit()

    # print(containerName, outObjOld.Get('HM_CharmFemto_Dpion_TrackCuts0'))
    for iContainer, containerName in enumerate(pathElements):
        print(f'\n\nStart of loop. container: {outObjOld}. name of target obj: {containerName}')
        if isinstance(outObjOld, TFile):
            print("---------------- 1")
            outObjOld.ls()
            outObjOld = outObjOld.Get(containerName)
            print("---------------- 2")
            outObjOld.ls()
            print(outObjOld)
        elif isinstance(outObjOld, TList):
            outObjOld = outObjOld.FindObject(containerName)
        elif isinstance(outObjOld, TDirectoryFile):
            outObjOld = outObjOld.Get(containerName)
        else:
            print(f'\033[31mError\033[0m: instance of {type(outObjOld)} not implemented. Exit!')
            sys.exit()
        print(f'End of loop. container: {outObjOld}. name of target obj: {containerName}')
        
    return outObjOld

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

def ConvertParquet2Root(parquetFile, outputName, treename='tree'):
    '''
    Converts a pandas dataframe to a TTree and saves it in a root file.
    Warning: this function is very slow.
    Arguments
    ----------
    - parquetFile: the absolute path of the file that you wanto to convert
    - outputName: the name of the root file

    '''

    oFile = TFile(outputName, 'recreate')
    tree = TTree(treename, treename)

    df = pd.read_parquet(parquetFile, engine='pyarrow')

    print(f'\n\033[94mConversion of {parquetFile}\033[0m')

    data = []
    for iCol, colname in enumerate(list(df)):
        data.append(np.ones((1), dtype="float32"))
        tree.Branch(colname, data[iCol], colname+"/F")

    with alive_bar(len(df)) as bar:
        for i in range(len(df)):
            for iCol, colname in enumerate(list(df)):
                data[iCol][0] = df.values[i, iCol]
            tree.Fill()
            bar()

    oFile.Write()
