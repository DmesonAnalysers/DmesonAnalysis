# DmesonStdAnalysisPbPb
Code for the measurement of Ds and D+ meson pT-differential yields starting from the outputs of the [AliPhysics](https://github.com/alisw/AliPhysics) tasks [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx) and [AliAnalysisTaskSEDplus.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDplus.cxx), using rectangular or ML selections

## Run analysis tasks

### Creation of files with selections to be applied on the tasks
* In the [cutobjects](https://github.com/fgrosa/DmesonStdAnalysisPbPb/tree/master/cutobjects) folder all the macros needed to produce the cut-object files used in the tasks are stored

### Run D+ and Ds tasks with private jobs
The [AliAnalysisTaskSEDplus.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDplus.cxx) and [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx) tasks can be run with private jobs using the ```RunAnalysisDsDplusTask.C``` script in the ```runanalysistask``` folder:

```
root -l RunAnalysisDsDplusTask.cc+(TString configfilename = configfile.yml, TString runMode = "full", bool mergeviajdl = true)
```
where ```configfile.yml``` is a configuration file (such as [runAnalysis_config_LHC17p_cent.yml](runanalysistask/runAnalysis_config_LHC17p_cent.yml)) with the information about the dataset, the AliPhysics version, and the task options to be used. The tasks options include the possibility to create a tree for the ML studies or apply a ML model trained with [xgboost](https://xgboost.readthedocs.io/en/latest/) or [scikit learn](https://scikit-learn.org/stable/).

### Train output merge
* The by-hand merge of unmerged outputs of a [ALICE analysis train](http://alimonitor.cern.ch/map.jsp) or private jobs can be performed with the script in the ```merge``` folder:
```
python MergeTrainOutputs.py files_to_merge.yml
```
where ```files_to_merge.yml``` is the configuration file containing the information about the outputs that has to be merged such as [files_to_merge_LHC18q.yml](merge/files_to_merge_LHC18q.yml)

## Main analysis with THnSparses

### Pre-filter ThnSparses
* The THnSparse in the task outputs can be pre-filtered to reduce the file size (useful if the train outputs are too large and cannot be merged) with the ```FilterSparse.py``` script in the ```filterdata``` folder:
```
python FilterSparse.py configfile.yml cutset.yml
```
where ```configfile.yml``` is the config file with the info of the input files and ```cutset.yml``` is the set of selections to be applied in the filtering. It creates output files as the input ones, with the ThnSparses filtered. 
With the option ```--suffix suffixname```, a suffix is added to the output file names, otherwise the input files are overwritten.
With the option ```--plot``` it creates control plots that are saved in .pdf files 

### Projection of invariant-mass distributions from THnSparse
* Project the THnSparse with the desired selections into invariant-mass distributions (TH1F):
```
python ProjectDplusDsSparse.py configfile.yml cutset.yml output.root
```
## Main analysis with TTrees or dataframes

### Filter trees to prepare data sets for ML studies
To filter trees produced with the Ds and D+ tasks and divide each category (data, MC prompt D, MC feed-down D, MC background) in a separated file (tree or dataframe) to prepare the datasets for the ML analyses, the ```FilterTrees4ML.cc``` and ```FilterTrees4ML.py``` scripts in the ```filterdata``` folder can be used:
```
root -l FilterTrees4ML.cc+(TString configfilename = configfile.yml)
```
or 
```
python3 FilterTrees4ML.py configfile.yml
```
where ```configfile.yml``` is a configuration file (such as [config_Dplus_data_skim_pp5TeV.yml](filterdata/config_Dplus_data_skim_pp5TeV.yml)) that contains the information about the input files, the preselections to apply, the features to keep and the output files. The output files are by default ```root``` files. In the case of the python script, if the ```--parquet``` option is used, the output data are saved into ```parquet``` files instead of ```root``` files. 

### Projection of invariant-mass distributions from THnSparse
* Project the TTree or dataframe with the desired selections into invariant-mass distributions (TH1F):
```
python ProjectDplusDsTree.py configfile.yml cutset.yml output.root
```
It autodetects whether the input files are ```root``` files containing TTrees or ```parquet``` files containing pandas dataframes.

## Common analysis
The following steps can be performed after having projected THnSparse or TTree (dataframe) objects

### Raw yield extraction
* Perform raw-yield extraction (root)
```
root -l GetRawYieldsDplusDs.C+(int cent, bool isMC = false, TString infilename = "distributions.root", TString cfgfilename = "config_Fit.yml", TString outFileName = "output.root")
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
where ```effFileName.root``` is the file with the efficiencies computed in the previous step and ```accFileName.root``` is the file with the acceptance computed using the [ComputeAcceptance.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/ComputeAcceptance.C)

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

## Significance optimisation

### Optimisation with THnSparse

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
where the configuration file ```cutset.yml``` contains the selection that you want to apply, such as [cutset_010_central_2018.yml](configfiles/cutsets/cutset_010_central_2018.yml) if you want to apply rectangular selections or [cutset_010_ML_test.yml](configfiles/cutsets/cutset_010_ML_test.yml) if you want to apply a selection on the ML output.

## Systematic uncertainties
### Selection efficiency
* For the cut-variation studies with THnSparses, the configuration files for each set of selection criteria can be created using:
```
MakeCutsFilesForSyst.py
```
the variables and the ranges should be set hard-coded in the script
Once the configuration files are created they can be used to repeat the main analysis with the different selection criteria

### Raw-yield extraction
* For the raw-yield extraction uncertainty a multi-trial study can be run with:
```
root -l RawYieldSystematics.C+(TString outfilerawname = "output.root")
```

## Test and validation of code for production of trees ([AliAnalysisTaskSEHFTreeCreator.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/treeHF/AliAnalysisTaskSEHFTreeCreator.cxx))
The validation of the code for production of trees used in ML studies can be done using the scripts in the ```runanalysistask``` folder

* To run the [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx) and [AliAnalysisTaskSEHFTreeCreator.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/treeHF/AliAnalysisTaskSEHFTreeCreator.cxx) on the same files:
```
root -l RunAnalysisTreeCreator.cc+(TString configfilename = configfile.yml, TString runMode = "full", bool mergeviajdl = true)
```
where ```configfile.yml``` is a configuration file with the information about the dataset and the AliPhysics version to be used

* To run the validation of the output:
```
python ValidateTreeCreator.py inputfile inputdir inputlist
```
where ```inputfile```, ```inputdir```, and ```inputlist``` are the root file produced by ```RunAnalysisTreeCreator.C```, the name of the TDirectoryFile and the TList inside the root file
