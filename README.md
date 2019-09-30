# DmesonStdAnalysisPbPb
Code for the measurement of Ds and D+ meson pT-differential yields starting from the outputs of the [AliPhysics](https://github.com/alisw/AliPhysics) tasks [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/AliAnalysisTaskSEDs.cxx) and [AliAnalysisTaskSEDplus.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/AliAnalysisTaskSEDplus.cxx), using rectangular or ML selections

## Significance optimisation

* Compute expected significance for all combinations of different selection criteria:
```
python ScanSignificanceSparse.py configfile.yml output.root
```
where ```configfile.yml``` is a configuration file such as [config_Ds_010_SignOpt.yml](configfiles/config_Ds_010_SignOpt.yml) 

* Project ntuple with expected significance as a function of relevant variables:
```
python ProjectSignifNtuple.py configfile.yml input.root PtMin PtMax minSignificance maxSignificance minEffPrompt maxEffPrompt
```
where the input file ```input.root``` is the one produced in the previous step

* Produce plots for expected significance as a function of pT:
```
python GetExpectedSignificance.py configfile.yml cutset.yml outputname
```
where the configuration file ```cutset.yml``` contains the selection that you want to apply, such as [cutset_010_central_2018.yml](configfiles/cutset_010_central_2018.yml) if you want to apply rectangular selections or [cutset_010_ML_test.yml](configfiles/cutset_010_ML_test.yml) if you want to apply a selection on the ML output.

## Main analysis 

### Projection of invariant-mass distributions from THnSparse
* Project the THnSparse with the desired selections into invariant-mass distributions (TH1F):
```
python ProjectDplusDsSparse.py configfile.yml cutset.yml output.root
```

### Raw yield extraction
* Perform raw-yield extraction (root)
```
GetRawYieldsDplusDs.C+ (int cent, bool isMC = false, TString infilename = "distributions.root", TString cfgfilename = "config_Fit.yml", TString outFileName = "output.root")
```
where ```distributions.root``` is the file obtained projecting the data or MC THnSparse and ```config_Fit.yml``` is a configuration file with the inputs needed to perform the invariant-mass fits such as [config_Ds_Fit.yml](configfiles/config_Ds_Fit.yml)

### Efficiency-times-acceptance computation
The efficiency-times-acceptance computation is done in two steps:
* Efficiency computation:
```
python ComputeEfficiencyDplusDs.py config_Fit.yml distributionsMC.root output.root
```
where ```distributions.root``` is the file obtained projecting the MC THnSparse and ```config_Fit.yml``` is the same config file used for the raw-yield extraction needed to have the same pT binning

To apply pT weights the ```--ptweights``` argument followed by the name of the input file with the pT weights and the name of the pT-weights histogram

* Acceptance and efficiency combination:
```
python CombineAccTimesEff.py effFileName.root accFileName.root outFileName.root
```
where ```effFileName.root``` is the file with the efficiencies computed in the previous step and ```accFileName.root``` is the file with the acceptance computed using the [ComputeAcceptance.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/ComputeAcceptance.C) macro from [AliPhysics](https://github.com/alisw/AliPhysics)

both can be run with the ```--batch``` argument to avoid the canvas window

### Cross section

* For the computation of the cross section, a modified version of [HFPtSpectrum.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/HFPtSpectrum.C) present in this repository, is used:

### Nuclear modification factor

* For the computation of the nuclear modification factor, a modified version of [HFPtSpectrumRaa.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/HFPtSpectrumRaa.C) present in this repository, is used

### Corrected yield

* For the computation of the pT-differential corrected yields, a modified version of [ComputeDmesonYield.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/ComputeDmesonYield.C) present in this repository, is used

### Run full analysis
* To run the full analysis, from the raw-yield extraction to the nuclear modification factor and the corrected yields the script
```
sh RunFullAnalysis.sh
```
can be used by setting some hard-coded parameters 

## Systematic uncertainties
### Selection efficiency
* For the cut-variation studies, the configuration files for each set of selection criteria can be created using:
```
MakeCutsFilesForSyst.py
```
the variables and the ranges should be set hard-coded in the script
Once the configuration files are created they can be used to repeat the main analysis with the different selection criteria

### Raw-yield extraction
* For the raw-yield extraction uncertainty a multi-trial study can be run with:
```
RawYieldSystematics(TString outfilerawname = "output.root")
```

## Creation of files with selections to be applied on the tasks
* In the [cutobjects](https://github.com/fgrosa/DmesonStdAnalysisPbPb/tree/master/cutobjects) folder all the macros needed to produce the cut-object files used in the tasks are stored

## Train output merge
* The by-hand merge of unmerged outputs of a [ALICE analysis train](http://alimonitor.cern.ch/map.jsp) can be performed with the script in the ```merge``` folder:
```
python MergeTrainOutputs.py files_to_merge.yml
```
where ```files_to_merge.yml``` is the configuration file containing the information about the outputs that has to be merged such as [files_to_merge_LHC18q.yml](merge/files_to_merge_LHC18q.yml)

## Pre-filter ThnSparses
* If the train outputs are too large and cannot be merged, they can be pre-filtered with:
```
python FilterSparse configfile.yml cutset.yml
```
where ```configfile.yml``` is the config file with the info of the input files and ```cutset.yml``` is the set of selections to be applied in the filtering. It creates output files as the input ones, with the ThnSparses filtered. With the option --suffix ```suffix```, a suffix is added to the output file names, otherwise the input files are overwritten

## Test and validation of code for production of trees used in ML studies
The validation of the code for production of trees used in ML studies can be done using the scripts in the ```treecreator``` folder

* To run the [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/AliAnalysisTaskSEDs.cxx) and [AliAnalysisTaskSEHFTreeCreator.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/treeHF/AliAnalysisTaskSEHFTreeCreator.cxx) on the same files:
```
RunAnalysisTreeCreator.C (TString configfilename = configfile.yml, TString runMode = "full", bool mergeviajdl = true)
```
where ```configfile.yml``` is a configuration file with the information about the dataset and the AliPhysics version to be used

* To run the validation of the output:
```
python ValidateTreeCreator.py inputfile inputdir inputlist
```
where ```inputfile```, ```inputdir```, and ```inputlist``` are the root file produced by ```RunAnalysisTreeCreator.C```, the name of the TDirectoryFile and the TList inside the root file
