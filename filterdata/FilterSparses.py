'''
python script to filter too large sparses (sparses in input files are overwritten with filtered ones)
and manage sparses for v2 analyses with ML
run: python FilterSparse cfgFileName.yml [filtFileName.yml, ...]
with option --suffix a suffix is added to the output file name, otherwise the input file is overwritten
'''

import sys
import argparse
import array
import yaml
from ROOT import TFile, TDirectoryFile, TCanvas  # pylint: disable=import-error,no-name-in-module
sys.path.append('..')
from utils.TaskFileLoader import LoadSparseFromTask, LoadListFromTask, LoadNormObjFromTask, LoadCutObjFromTask # pylint: disable=wrong-import-position,import-error
from utils.TaskFileLoader import LoadSparseFromTaskV2, LoadListFromTaskV2 # pylint: disable=wrong-import-position,import-error

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
                sparsesOrig[sparsetype].GetAxis(cutvars[iVar]['axisnum']).SetRange(binMin, binMax)

            # project on a new sparse for that pT bin
            naxes = len(axestokeep)
            sparseFiltPart[sparsetype] = sparsesOrig[sparsetype].Projection(
                naxes, array.array('i', axestokeep), 'O')
            if iPt == 0:
                sparsesFilt[sparsetype] = sparseFiltPart[sparsetype]
            else:
                sparsesFilt[sparsetype].Add(sparseFiltPart[sparsetype])

        sparsesFilt[sparsetype].SetName(sparsesOrig[sparsetype].GetName())

    return sparsesFilt

def FilterSparsesV2(sparse, cutvars, axestokeep):
    '''
    function that gets a dictionary of sparses and a dictionary of selections
    and returns a dictionary of filtered sparses
    '''
    sparseFiltPart = []
    # apply selections to sparses
    for iPt in range(0, len(cutvars['Pt']['min'])):
        for iVar in cutvars:
            if iVar == 'InvMass':
                continue
            binMin = sparse.GetAxis(
                cutvars[iVar]['axisnum']).FindBin(cutvars[iVar]['min'][iPt]*1.0001)
            binMax = sparse.GetAxis(
                cutvars[iVar]['axisnum']).FindBin(cutvars[iVar]['max'][iPt]*0.9999)
            sparse.GetAxis(cutvars[iVar]['axisnum']).SetRange(binMin, binMax)

        # project on a new sparse for that pT bin
        naxes = len(axestokeep)
        sparseFiltPart = sparse.Projection(naxes, array.array('i', axestokeep), 'O')
        if iPt == 0:
            sparsesFilt = sparseFiltPart
        else:
            sparsesFilt.Add(sparseFiltPart)
    sparsesFilt.SetName(sparse.GetName())

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

def FiltFunc(args):
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    with open(args.cutSetFileNames[0], 'r') as ymlCutSetFile:
        cutSetCfg = yaml.load(ymlCutSetFile, yaml.FullLoader)
    cutVars = cutSetCfg['cutvars']
    axesToKeep = cutSetCfg['axestokeep']

    infilenames = inputCfg['filename']
    if not isinstance(infilenames, list):
        infilenames = [infilenames]

    for iFile, infilename in enumerate(infilenames):
        sparseReco, sparseGen = LoadSparseFromTask(infilename, inputCfg)

        # filter sparses
        sparseFiltReco = FilterSparses(sparseReco, cutVars, axesToKeep)
        if inputCfg['isMC']:
            sparseFiltGen = FilterSparses(sparseGen, cutVars, axesToKeep)

        # plot filtered sparses (each variable vs pt)
        if args.plot:
            PlotFiltVarsVsPt(sparseFiltReco, 'cReco{0}'.format(iFile))
            if inputCfg['isMC']:
                PlotFiltVarsVsPt(sparseFiltGen, 'cGen{0}'.format(iFile))

        # get other objects from original file
        inlist = LoadListFromTask(infilename, inputCfg)
        _, normCounter = LoadNormObjFromTask(infilename, inputCfg)
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
        if args.suffix:
            outfilename = outfilename.replace(
                '.root', '{0}.root'.format(args.suffix))
        print('Saving filtered ThnSparses in file', outfilename)
        outfile = TFile(outfilename, 'recreate')
        outdir = TDirectoryFile(inputCfg['dirname'], inputCfg['dirname'])
        outdir.Write(inputCfg['dirname'])
        outdir.cd()
        inlist.Write(inputCfg['listname'], 1)
        cutObj.Write(cutObjName)
        normCounter.Write()
        outfile.Close()

def FiltFuncV2(args):
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)

    cutSetsCfg = []
    for cutFile in args.cutSetFileNames:
        with open(cutFile, 'r') as file:
            cutCfg = yaml.load(file, yaml.FullLoader)
            cutSetsCfg.append(cutCfg)

    if not len(cutSetsCfg) == len(inputCfg['dirname']) == len(inputCfg['listname']):
        print('Wrong number of cut sets, dir names or list names')
        return

    sparseList = LoadSparseFromTaskV2(inputCfg)
    sparseListFilt = []

    for sparse, cutSet in zip(sparseList, cutSetsCfg):
        sparseFilt = FilterSparsesV2(sparse, cutSet['cutvars'], cutSet['axestokeep'])
        sparseListFilt.append(sparseFilt)

    for sparse in sparseListFilt[1:]:
        sparseListFilt[0].Add(sparse)

    # get other objects from original file
    inlist = LoadListFromTaskV2(inputCfg['filename'], inputCfg['dirname'][0], inputCfg['listname'][0])
    sparsetodel = inlist.FindObject(sparseList[0].GetName())
    inlist.Remove(sparsetodel)
    inlist.Add(sparseListFilt[0])

    # save new file with filtered sparses
    outfilename = inputCfg['outfilename']
    print('Saving filtered ThnSparses in file', outfilename)
    outfile = TFile(outfilename, 'recreate')
    outdir = TDirectoryFile(inputCfg['outdirname'], inputCfg['outdirname'])
    outdir.Write(inputCfg['outdirname'])
    outdir.cd()
    inlist.Write(inputCfg['outlistname'], 1)
    outfile.Close()

def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='cfgFileName.yml',
                        help='config file name with root input files')
    parser.add_argument('cutSetFileNames', nargs='*', type=str, metavar='list',
                        default=['filtFileName.yml'], help='list of files with filtering selections')
    parser.add_argument('--v2', action='store_true', help='flag to filter for v2 analysis')
    parser.add_argument('--suffix', metavar='text', help='suffix to be added to output file')
    parser.add_argument('--plot', action='store_true', help='flag to enable plots')
    args = parser.parse_args()

    if not args.v2:
        FiltFunc(args)
    else:
        FiltFuncV2(args)

main()
