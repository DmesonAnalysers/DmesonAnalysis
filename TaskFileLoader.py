'''
python script with helper functions to load objects from task
'''

from ROOT import TFile  # pylint: disable=import-error,no-name-in-module


def LoadSparseFromTask(infilename, inputCfg):
    '''
    Method to retrieve sparses from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - list of sparses with reconstructed quantities for all candidates, and prompt, FD D mesons and bkg candidates (only if MC)
    - list of sparses with generated quantities for prompt and FD D mesons (only if MC)
    '''
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
    '''
    Method to retrieve single sparse from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file
    - sparse that should be returned

    Returns
    ----------
    - selected sparse from file
    '''
    print('Loading THnSparse from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])
    sparse = inlistData.FindObject(inputCfg[sparsetype])

    return sparse


def LoadNormObjFromTask(infilename, inputCfg):
    '''
    Method to retrieve normalisation objects from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - histo with event info and normalisation counter
    '''
    print('Loading norm objects from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])
    normCounter = indirData.Get(inputCfg['normname'])
    hEv = inlistData.FindObject(inputCfg['histoevname'])

    return hEv, normCounter


def LoadListFromTask(infilename, inputCfg):
    '''
    Method to retrieve list of objects from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - TList of objects from input file
    '''
    print('Loading TList from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])

    return inlistData


def LoadCutObjFromTask(infilename, inputCfg):
    '''
    Method to retrieve D-meson cut object from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - D-meson cut object
    '''
    print('Loading cut object from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    cutobjname = inputCfg['listname'].replace('coutputDs', 'coutputDsCuts')
    cutobjname = cutobjname.replace('coutputDplus', 'coutputDplusCuts')
    cutobj = indirData.Get(cutobjname)

    return cutobj, cutobjname


def LoadPIDSparses(infilename, inputCfg):
    '''
    Method to retrieve PID sparses from output task file

    Inputs
    ----------
    - input root file name
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - sparse of NsigmaTPC and NsigmaTOF variables
    - sparse of NsigmaComb variables
    - dictionary with variable : sparse-axis number
    '''
    print('Loading PID THnSparses from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(inputCfg['dirname'])
    inlistData = indirData.Get(inputCfg['listname'])
    sparsePIDNsigma = inlistData.FindObject('fnSparsePID')
    sparsePIDNsigmaComb = inlistData.FindObject('fnSparsePIDcomb')

    # dictionary of sparse axes with detectors, mass hypothesis, daughter number
    axes = {'TPC': {'Pi': {'0': 2, '1': 6, '2': 10}, 'K': {'0': 3, '1': 7, '2': 11}},
            'TOF': {'Pi': {'0': 4, '1': 8, '2': 12}, 'K': {'0': 5, '1': 9, '2': 13}},
            'Comb': {'Pi': {'0': 2, '1': 4, '2': 6}, 'K': {'0': 3, '1': 5, '2': 7}}}

    return sparsePIDNsigma, sparsePIDNsigmaComb, axes


def LoadSparseFromTaskV2(inputCfg):
    '''
    Method to retrieve sparses from D-meson vn output task file

    Inputs
    ----------
    - config dictionary from yaml file with name of objects in root file

    Returns
    ----------
    - list of sparses with reconstructed quantities for all candidates
    '''
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
    '''
    Method to retrieve list of objects from D-meson vn output task file

    Inputs
    ----------
    - input root file name
    - name of directory in input root file
    - name of list in directory in input root file

    Returns
    ----------
    - TList of objects from input file
    '''
    print('Loading TList from file', infilename)
    infileData = TFile(infilename)
    indirData = infileData.Get(dirname)
    inlistData = indirData.Get(listname)

    return inlistData
