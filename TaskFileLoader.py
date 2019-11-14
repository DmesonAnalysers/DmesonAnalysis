'''
python script with helper functions to load objects from task
'''

from ROOT import TFile  # pylint: disable=import-error,no-name-in-module

def LoadSparseFromTask(infilename, inputCfg):
    print('Loading THnSparses from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])
    sparses, sparsesGen = {}, {}
    sparses['RecoAll'] = inlistData.FindObject(inputCfg['sparsenameAll'])
    if inputCfg['isMC']:
        sparses['RecoPrompt'] = inlistData.FindObject(
            inputCfg['sparsenamePrompt'])
        sparses['RecoFD'] = inlistData.FindObject(inputCfg['sparsenameFD'])
        sparsesGen['GenPrompt'] = inlistData.FindObject(
            inputCfg['sparsenameGenPrompt'])
        sparsesGen['GenFD'] = inlistData.FindObject(
            inputCfg['sparsenameGenFD'])
        if inputCfg['enableSecPeak']:
            sparses['RecoSecPeakPrompt'] = inlistData.FindObject(
                inputCfg['sparsenamePromptSecPeak'])
            sparses['RecoSecPeakFD'] = inlistData.FindObject(
                inputCfg['sparsenameFDSecPeak'])
            sparsesGen['GenSecPeakPrompt'] = inlistData.FindObject(
                inputCfg['sparsenameGenPromptSecPeak'])
            sparsesGen['GenSecPeakFD'] = inlistData.FindObject(
                inputCfg['sparsenameGenFDSecPeak'])

    return sparses, sparsesGen

def LoadSingleSparseFromTask(infilename, inputCfg, sparsetype='sparsenameBkg'):
    print('Loading THnSparse from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])
    sparse = inlistData.FindObject(inputCfg[sparsetype])

    return sparse


def LoadNormObjFromTask(infilename, inputCfg):
    print('Loading norm objects from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])
    normCounter = indirData.Get(inputCfg['normname'])
    hEv = inlistData.FindObject(inputCfg['histoevname'])

    return hEv, normCounter


def LoadListFromTask(infilename, inputCfg):
    print('Loading TList from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])

    return inlistData


def LoadCutObjFromTask(infilename, inputCfg):
    print('Loading cut object from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    cutobjname = inputCfg['listname'].replace('coutputDs', 'coutputDsCuts')
    cutobjname = cutobjname.replace('coutputDplus', 'coutputDplusCuts')
    cutobj = indirData.Get(cutobjname)

    return cutobj, cutobjname


def LoadPIDTH3(infilename, inputCfg):
    print('Loading Nsigma distributions from file', infilename)
    hNsigmaVsPtVsML = {}
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])
    species = ['K', 'Pi']
    detectors = ['TPC', 'TOF', 'Comb']
    for det in detectors:
        hNsigmaVsPtVsML[det] = {}
        for spe in species:
            hNsigmaVsPtVsML[det][spe] = {}
            for iprong in range(3):
                hNsigmaVsPtVsML[det][spe]['{0}'.format(iprong)] = inlistData.FindObject(
                    'fHistNsigma{0}{1}Prong{2}'.format(det, spe, iprong))

    return hNsigmaVsPtVsML


def LoadPIDSparses(infilename, inputCfg):
    print('Loading PID THnSparses from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])
    sparsePIDNsigma = inlistData.FindObject('fnSparsePID')
    sparsePIDNsigmaComb = inlistData.FindObject('fnSparsePIDcomb')

    #dictionary of sparse axes with detectors, mass hypothesis, daughter number
    axes = {'TPC': {'Pi':{'0':2, '1':6, '2':10}, 'K':{'0':3, '1':7, '2':11}},
            'TOF': {'Pi':{'0':4, '1':8, '2':12}, 'K':{'0':5, '1':9, '2':13}},
            'Comb': {'Pi':{'0':2, '1':4, '2':6}, 'K':{'0':3, '1':5, '2':7}}}

    return sparsePIDNsigma, sparsePIDNsigmaComb, axes

def LoadSparseFromTaskV2(inputCfg):
    infilename = inputCfg['filename']
    print(f'Loading THnSparses from file {infilename}')
    sparses = []
    infileData = TFile(infilename)

    for dirname, listname in zip(inputCfg['dirname'], inputCfg['listname']):
        indirData = infileData.Get(dirname)
        if not indirData:
            print(f'Directory {dirname} not found!')
            return []
        inlistData = indirData.Get(listname)
        if not inlistData:
            print(f'List {listname} not found!')
            return []
        sparse = inlistData.FindObject(inputCfg['sparsename'])
        if not sparse:
            print('Sparse not found!')
            return []
        sparses.append(sparse)

    return sparses

def LoadListFromTaskV2(infilename, dirname, listname):
    print('Loading TList from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(dirname)
    inlistData = indirData.Get(listname)

    return inlistData
