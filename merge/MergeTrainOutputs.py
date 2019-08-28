# !/usr/bin/python3
import os
import argparse
import six
import yaml
from ROOT import TGrid, TFileMerger, gROOT #pylint: disable=import-error, no-name-in-module

#pylint: disable=invalid-name
parser = argparse.ArgumentParser(description='Arguments')
parser.add_argument('cfgfile', metavar='text', \
    default='files_to_merge.yml', help='input yml config file')
args = parser.parse_args()

with open(args.cfgfile, 'r') as ymlinputCfg:
    if six.PY2:
        inputCfg = yaml.load(ymlinputCfg)
    else:
        inputCfg = yaml.safe_load(ymlinputCfg)

TGrid.Connect('alien://')
fileMerger = TFileMerger()

if not inputCfg['MergeOptions']['MergeTrees']:
    fileMerger.SetNotrees(True)

if inputCfg['MergeOptions']['MergeByRun']:
    runs = inputCfg['MergeOptions']['RunNumbers']
else:
    runs = [('%03d' % iNum) for iNum in range(inputCfg['MergeOptions']['NfilesTot'])]

defMode = (TFileMerger.kAll | TFileMerger.kIncremental)
Mode = (defMode | TFileMerger.kOnlyListed)
objToMerge = ''
for obj in inputCfg['MergeOptions']['ObjectsToMerge']:
    if obj is not None:
        objToMerge += '%s ' % obj

#partially merged files from alien
nBunch = 0
for iRun, run in enumerate(runs):
    if inputCfg['MergeOptions']['IsMC']:
        inname = os.path.join(inputCfg['DataPath'], '%s' % run)
    else:
        inname = os.path.join(inputCfg['DataPath'], '000%s' % run)
    if inputCfg['RecoPass'] is not None and inputCfg['TrainName'] is not None:
        inname = os.path.join(inname, inputCfg['RecoPass'], inputCfg['TrainName'])
    inname = os.path.join(inname, inputCfg['InputFileName'])
    fileMerger.AddFile(inname)
    if iRun%inputCfg['MergeOptions']['NfilesPerChunk'] == 0 or iRun == len(runs)-1:
        print('Merging up to %s' % inname)
        outname = inputCfg['OutputFileName'].replace('.root', '_%04d.root' % nBunch)
        outname = os.path.join(inputCfg['OutputPath'], outname)
        fileMerger.OutputFile(outname)
        if objToMerge:
            fileMerger.AddObjectNames(objToMerge)
            fileMerger.PartialMerge(Mode)
        else:
            fileMerger.Merge()
        fileMerger.Reset()
        nBunch += 1
fileMerger.Reset()

# merge partially merged files
if inputCfg['MergeOptions']['DoTotalMerge']:
    print('Merging partial files')
    outname = os.path.join(inputCfg['OutputPath'], inputCfg['OutputFileName'])
    if nBunch > 1:
        fileMerger.OutputFile(outname)
        for iMergedFile in range(nBunch):
            inname = inputCfg['OutputFileName'].replace('.root', '_%04d.root' % iMergedFile)
            inname = os.path.join(inputCfg['OutputPath'], inname)
            fileMerger.AddFile(inname)
        fileMerger.Merge()
        fileMerger.Reset()
        for iMergedFile in range(nBunch):
            inname = inputCfg['OutputFileName'].replace('.root', '_%04d.root' % iMergedFile)
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
