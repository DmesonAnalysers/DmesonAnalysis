# !/usr/bin/python3
'''
Script to merge unmerged outputs of tasks run on the grid
run: python MergeTrainOutputs.py files_to_merge.yml
'''

import os
import argparse
import yaml
from ROOT import TGrid, TFileMerger #pylint: disable=import-error, no-name-in-module
from ROOT import gROOT #pylint: disable=import-error, no-name-in-module


def Merge(merger, outfilename, objtomerge=None, mode=None):
    '''
    general function for merging
    '''

    merger.OutputFile(outfilename)
    if objtomerge is not None and objtomerge:
        merger.AddObjectNames(objtomerge)
        if mode is not None:
            merger.PartialMerge(mode)
    else:
        merger.Merge()
    merger.Reset()
    print(f'\33[32mMerged files in {outfilename}\33[0m')

parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgfile', metavar='text',
                    default='files_to_merge.yml', help='input yml config file')
args = parser.parse_args()

with open(args.cfgfile, 'r') as ymlinputCfg:
    inputCfg = yaml.load(ymlinputCfg, yaml.FullLoader)

# download files from alien
runs = inputCfg['MergeOptions']['RunNumbers']
runsToRemove = []
nPerRun = inputCfg['MergeOptions']['NfilesPerRun']
inDirs = []
grid = TGrid.Connect('alien://')
for iRun, run in enumerate(runs):
    dirName = inputCfg['DataPath']
    inDirs.append([])
    if run is not None:
        if isinstance(run, int):
            if inputCfg['MergeOptions']['IsMC']:
                dirName = os.path.join(dirName, f'{run:d}')
            else:
                dirName = os.path.join(dirName, f'{run:09d}')
        else:
            dirName = os.path.join(dirName, f'{run}')

    if inputCfg['RecoPass'] is not None:
        dirName = os.path.join(dirName, inputCfg['RecoPass'])
    if inputCfg['TrainName'] is not None:
        dirName = os.path.join(dirName, inputCfg['TrainName'])

    if not grid.Cd(dirName.replace('alien://', '')):
        print(f'\33[31mNo outputs for run {run}, continue\33[0m')
        runsToRemove.append(run)
        continue

    listOfFiles = grid.Ls(dirName.replace('alien://', ''))
    for iFile in range(listOfFiles.GetEntries()):
        if listOfFiles.GetFileName(iFile).isdigit():
            inDirs[iRun].append(listOfFiles.GetFileName(iFile))

    print(f'\33[32mStart download of files for run {run}\33[0m')
    if not inputCfg['MergeOptions']['MergeByRun']:
        for inDir in inDirs[iRun]:
            outName = os.path.join(inputCfg["OutputPath"], f'{run}')
            os.system(f'alien_cp -select {inputCfg["InputFileName"]} -y 2 -T 32 {dirName} file://{outName}/')
    else:
        outName = os.path.join(inputCfg["OutputPath"], f'{run}')
        os.system(f'alien_cp -T 32 {dirName}/{inputCfg["InputFileName"]} file://{outName}/')

# remove runs with no output
for run in runsToRemove:
    runs.remove(run)

# merge options
fileMerger = TFileMerger()
fileMerger.SetPrintLevel(0)

if not inputCfg['MergeOptions']['MergeTrees']:
    fileMerger.SetNotrees(True)

defMode = (TFileMerger.kAll | TFileMerger.kIncremental)
Mode = (defMode | TFileMerger.kOnlyListed)
objToMerge = ''
for obj in inputCfg['MergeOptions']['ObjectsToMerge']:
    if obj is not None:
        objToMerge += f'{obj} '

# do partial merge within runs
if not inputCfg['MergeOptions']['MergeByRun']:
    for iRun, run in enumerate(runs):
        print(f'\33[32mPartially merging files for run {run}\33[0m')
        nBunch = 0
        for dirNum, inDir in enumerate(inDirs[iRun]):
            inName = os.path.join(inputCfg["OutputPath"], f'{run}', f'{inDir}')
            fileMerger.AddFile(f'{inName}/{inputCfg["InputFileName"]}')
            if (dirNum%inputCfg['MergeOptions']['NfilesPerChunk'] == 0 and dirNum != 0) or inDir == max(inDirs[iRun]):
                outBunchName = os.path.join(inputCfg['OutputPath'], f'{run}',
                                            inputCfg['OutputFileName'].replace('.root', f'_{nBunch:04d}.root'))
                Merge(fileMerger, outBunchName, objToMerge, Mode)
                nBunch += 1
        # remove partial files
        for dirNum, inDir in enumerate(inDirs[iRun]):
            inName = os.path.join(inputCfg["OutputPath"], f'{run}', f'{inDir}')
            os.remove(f'{inName}/{inputCfg["InputFileName"]}')
            os.rmdir(inName)

fileMerger.Reset()

# merge partially merged files (in steps, if needed)
if inputCfg['MergeOptions']['DoTotalMerge']:
    print('\33[32mMerging partial files\33[0m')
    filesToMerge = []
    for run in runs:
        inRunDir = os.path.join(inputCfg['OutputPath'], f'{run}')
        for fileName in os.listdir(inRunDir):
            if os.path.isdir(os.path.join(inRunDir, fileName)):
                for fileNameSubdir in os.listdir(os.path.join(inRunDir, fileName)):
                    filesToMerge.append(os.path.join(inRunDir, fileName, fileNameSubdir))
            else:
                filesToMerge.append(os.path.join(inRunDir, fileName))
    if len(filesToMerge) > 10: # do final merge in bunches
        nBunch = 0
        for iFile, fileToMerge in enumerate(filesToMerge):
            fileMerger.AddFile(fileToMerge)
            if (iFile > 0 and iFile % 10 == 0) or iFile == len(filesToMerge)-1:
                outName = os.path.join(inputCfg['OutputPath'],
                                       inputCfg['OutputFileName'].replace('.root', f'_{nBunch:02d}.root'))
                Merge(fileMerger, outName, objToMerge, Mode)
                nBunch += 1
        fileMerger.Reset()
        for bunch in range(nBunch):
            inName = os.path.join(inputCfg['OutputPath'],
                                  inputCfg['OutputFileName'].replace('.root', f'_{bunch:02d}.root'))
            fileMerger.AddFile(inName)
        Merge(fileMerger, os.path.join(inputCfg['OutputPath'], inputCfg['OutputFileName']))
        # clenup of intermediate steps
        for bunch in range(nBunch):
            inName = os.path.join(inputCfg['OutputPath'],
                                  inputCfg['OutputFileName'].replace('.root', f'_{bunch:02d}.root'))
            os.remove(inName)
    elif 1 < len(filesToMerge) <= 10:
        for iFile, fileToMerge in enumerate(filesToMerge):
            fileMerger.AddFile(fileToMerge)
        Merge(fileMerger, os.path.join(inputCfg['OutputPath'], inputCfg['OutputFileName']))
    else:
        if os.path.isfile(filesToMerge[0]):
            os.rename(filesToMerge[0], os.path.join(inputCfg['OutputPath'], inputCfg['OutputFileName']))

    print(f'\33[32mTotal merged output: {os.path.join(inputCfg["OutputPath"], inputCfg["OutputFileName"])}\33[0m')

    # clenup of intermediate steps
    for fileToMerge in filesToMerge:
        if os.path.isfile(fileToMerge):
            os.remove(fileToMerge)
    for run in runs:
        inRunDir = os.path.join(inputCfg['OutputPath'], f'{run}')
        for subDir in os.listdir(inRunDir):
            if os.path.isdir(os.path.join(inRunDir, subDir)):
                os.rmdir(os.path.join(inRunDir, subDir))
        os.rmdir(os.path.join(inRunDir))

gROOT.Reset()
