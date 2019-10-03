# !/usr/bin/python3
import os
import argparse
import six
import yaml
from ROOT import TGrid, TFileMerger #pylint: disable=import-error, no-name-in-module
from ROOT import gROOT #pylint: disable=import-error, no-name-in-module


def merge(merger, outfilename, objtomerge=None, mode=None):
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
    print('Merged files in ', outfilename)


#pylint: disable=invalid-name
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgfile', metavar='text', \
    default='files_to_merge.yml', help='input yml config file')
args = parser.parse_args()

with open(args.cfgfile, 'r') as ymlinputCfg:
    if six.PY2:
        inputCfg = yaml.load(ymlinputCfg)
    elif six.PY3:
        inputCfg = yaml.safe_load(ymlinputCfg)

grid = TGrid.Connect('alien://')
fileMerger = TFileMerger()
fileMerger.SetPrintLevel(0)

if not inputCfg['MergeOptions']['MergeTrees']:
    fileMerger.SetNotrees(True)

defMode = (TFileMerger.kAll | TFileMerger.kIncremental)
Mode = (defMode | TFileMerger.kOnlyListed)
objToMerge = ''
for obj in inputCfg['MergeOptions']['ObjectsToMerge']:
    if obj is not None:
        objToMerge += '{0} '.format(obj)

#partially merged files from alien
nBunch = 0
runs = inputCfg['MergeOptions']['RunNumbers']
nPerRun = inputCfg['MergeOptions']['NfilesPerRun']
for iRun, run in enumerate(runs):
    dirname = inputCfg['DataPath']
    indirs = []
    if run is not None:
        if inputCfg['MergeOptions']['IsMC']:
            dirname = os.path.join(dirname, '{:d}'.format(run))
        else:
            dirname = os.path.join(dirname, '{:09d}'.format(run))
    if (inputCfg['RecoPass'] is not None) and (inputCfg['TrainName'] is not None):
        dirname = os.path.join(dirname, inputCfg['RecoPass'], inputCfg['TrainName'])

    if not grid.Cd(dirname.replace('alien://', '')):
        print('No outputs for this run, continue')
        continue
    else:
        listoffiles = grid.Ls(dirname.replace('alien://', ''))
        for iFile in range(listoffiles.GetEntries()):
            if listoffiles.GetFileName(iFile).isdigit():
                indirs.append(listoffiles.GetFileName(iFile))

    dirnum = 1
    if not inputCfg['MergeOptions']['MergeByRun']:
        for indir in indirs:
            inname = os.path.join(dirname, '{0}'.format(indir), inputCfg['InputFileName'])
            fileMerger.AddFile(inname)
            if dirnum%inputCfg['MergeOptions']['NfilesPerChunk'] == 0 or indir == max(indir):
                print('Merging up to ', inname)
                outname = inputCfg['OutputFileName'].replace('.root', '_{:04d}.root'.format(nBunch))
                outname = os.path.join(inputCfg['OutputPath'], outname)
                merge(fileMerger, outname, objToMerge, Mode)
                nBunch += 1
            dirnum += 1
        print('Merging up to ', inname)
        outname = inputCfg['OutputFileName'].replace('.root', '_{:04d}.root'.format(nBunch))
        outname = os.path.join(inputCfg['OutputPath'], outname)
        merge(fileMerger, outname, objToMerge, Mode)
        nBunch += 1
    else:
        inname = os.path.join(dirname, inputCfg['InputFileName'])
        fileMerger.AddFile(inname)
        if iRun%inputCfg['MergeOptions']['NfilesPerChunk'] == 0 or iRun == len(iRun)-1:
            print('Merging up to ', inname)
            outname = inputCfg['OutputFileName'].replace('.root', '_{:04d}.root'.format(nBunch))
            outname = os.path.join(inputCfg['OutputPath'], outname)
            merge(fileMerger, outname, objToMerge, Mode)
            nBunch += 1

fileMerger.Reset()

# merge partially merged files
if inputCfg['MergeOptions']['DoTotalMerge']:
    print('Merging partial files')
    outname = os.path.join(inputCfg['OutputPath'], inputCfg['OutputFileName'])
    if nBunch > 1:
        for iMergedFile in range(nBunch):
            inname = inputCfg['OutputFileName'].replace('.root', '_{:04d}.root'.format(iMergedFile))
            inname = os.path.join(inputCfg['OutputPath'], inname)
            fileMerger.AddFile(inname)
        merge(fileMerger, outname)
        for iMergedFile in range(nBunch):
            inname = inputCfg['OutputFileName'].replace('.root', '_{:04d}.root'.format(iMergedFile))
            inname = os.path.join(inputCfg['OutputPath'], inname)
            if os.path.isfile(inname):
                os.remove(inname)
    else:
        inname = inputCfg['OutputFileName'].replace('.root', '_0000.root')
        inname = os.path.join(inputCfg['OutputPath'], inname)
        if os.path.isfile(inname):
            os.rename(inname, outname)
    print('Total merged output: ', outname)

gROOT.Reset()
