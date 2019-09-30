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
