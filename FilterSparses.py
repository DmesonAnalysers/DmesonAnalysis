'''
python script to filter too large sparses (sparses in input files are overwritten with filtered ones)
run: python FilterSparse cfgFileName.yml filtFileName.yml
with option --suffix a suffix is added to the output file name, otherwise the input file is overwritten
'''

import argparse
import array
import six
import yaml
from ROOT import TFile, TDirectoryFile, TCanvas  # pylint: disable=import-error,no-name-in-module
from TaskFileLoader import LoadSparseFromTask, LoadListFromTask, LoadNormObjFromTask, LoadCutObjFromTask


def FilterSparses(sparsesOrig, cutvars, axestokeep):
    '''
    function that gets a dictionary of sparses and a dictionary of selections
    and returns a dictionary of filtered sparses
    '''
    sparseFiltPart, sparsesFilt = {}, {}
    for sparsetype in sparsesOrig:
        # apply selections to sparses
        for iPt in range(0, len(cutvars['Pt']['min'])):
            for iVar in cutvars:
                if iVar == 'InvMass':
                    continue
                binMin = sparsesOrig[sparsetype].GetAxis(
                    cutvars[iVar]['axisnum']).FindBin(cutvars[iVar]['min'][iPt]*1.0001)
                binMax = sparsesOrig[sparsetype].GetAxis(
                    cutvars[iVar]['axisnum']).FindBin(cutvars[iVar]['max'][iPt]*0.9999)
                sparsesOrig[sparsetype].GetAxis(
                    cutvars[iVar]['axisnum']).SetRange(binMin, binMax)

            # project on a new sparse for that pT bin
            naxes = len(axestokeep)
            sparseFiltPart[sparsetype] = sparsesOrig[sparsetype].Projection(
                naxes, array.array('i', axestokeep), 'O')
            if iPt == 0:
                sparsesFilt[sparsetype] = sparseFiltPart[sparsetype]
            else:
                sparsesFilt[sparsetype].Add(sparseFiltPart[sparsetype])

        sparsesFilt[sparsetype].SetName(sparseReco[sparsetype].GetName())

    return sparsesFilt


def PlotFiltVarsVsPt(sparses, canvasname):
    '''
    function that gets a dictionary of sparses and returns a dictionary of canvases with scatter plots
    '''
    cVars, histos = {}, {}
    for sparsetype in sparses:
        histos[sparsetype] = []
        cVars[sparsetype] = TCanvas(
            '{0}{1}'.format(canvasname, sparsetype), '', 1920, 1080)
        cVars[sparsetype].Divide(
            int(sparses[sparsetype].GetNdimensions()/3)+1, 3)
        for iVar in range(sparses[sparsetype].GetNdimensions()):
            histos[sparsetype].append(sparses[sparsetype].Projection(iVar, 1, 'O'))
            histos[sparsetype][iVar].SetName('{0}{1}{2}'.format(canvasname, sparsetype, iVar))
            cVars[sparsetype].cd(iVar+1).SetLogz()
            histos[sparsetype][iVar].Draw('colz')

        cVars[sparsetype].SaveAs('{0}{1}.pdf'.format(canvasname, sparsetype))

# main function
PARSER = argparse.ArgumentParser(description='Arguments to pass')
PARSER.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                    help='config file name with root input files')
PARSER.add_argument('cutSetFileName', metavar='text', default='filtFileName.yml',
                    help='input file with filtering selections')
PARSER.add_argument('--suffix', metavar='text',
                    help='suffix to be added to output file')
PARSER.add_argument('--plot', action='store_true',
                    help='flag to enable plots')
ARGS = PARSER.parse_args()

with open(ARGS.cfgFileName, 'r') as ymlCfgFile:
    inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

with open(ARGS.cutSetFileName, 'r') as ymlCutSetFile:
    cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
cutVars = cutSetCfg['cutvars']
axesToKeep = cutSetCfg['axestokeep']

infilenames = inputCfg['filename']
if not isinstance(infilenames, list):
    infilenames = [infilenames]
isMC = inputCfg['isMC']
enableSecPeak = inputCfg['enableSecPeak']

for iFile, infilename in enumerate(infilenames):
    sparseReco, sparseGen = LoadSparseFromTask(infilename, inputCfg)

    # filter sparses
    sparseFiltReco = FilterSparses(sparseReco, cutVars, axesToKeep)
    if inputCfg['isMC']:
        sparseFiltGen = FilterSparses(sparseGen, cutVars, axesToKeep)

    # plot filtered sparses (each variable vs pt)
    if ARGS.plot:
        PlotFiltVarsVsPt(sparseFiltReco, 'cReco{0}'.format(iFile))
        if inputCfg['isMC']:
            PlotFiltVarsVsPt(sparseFiltGen, 'cGen{0}'.format(iFile))

    # get other objects from original file
    inlist = LoadListFromTask(infilename, inputCfg)
    hEv, normCounter = LoadNormObjFromTask(infilename, inputCfg)
    cutObj, cutObjName = LoadCutObjFromTask(infilename, inputCfg)
    for sparse in sparseFiltReco.values():
        sparsetodel = inlist.FindObject(sparse.GetName())
        inlist.Remove(sparsetodel)
        inlist.Add(sparse)
    if inputCfg['isMC']:
        for sparse in sparseFiltGen.values():
            sparsetodel = inlist.FindObject(sparse.GetName())
            inlist.Remove(sparsetodel)
            inlist.Add(sparse)

    # save new file with filtered sparses
    outfilename = infilename
    if ARGS.suffix:
        outfilename = outfilename.replace(
            '.root', '{0}.root'.format(ARGS.suffix))
    print('Saving filtered ThnSparses in file', outfilename)
    outfile = TFile(outfilename, 'recreate')
    outdir = TDirectoryFile(inputCfg['dirname'], inputCfg['dirname'])
    outdir.Write(inputCfg['dirname'])
    outdir.cd()
    inlist.Write(inputCfg['listname'], 1)
    cutObj.Write(cutObjName)
    normCounter.Write()
    outfile.Close()
