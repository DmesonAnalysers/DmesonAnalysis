# D<sub>s</sub><sup>+</sup> and D<sup>+</sup>-meson Analysis code

Code for the measurement of D<sub>s</sub><sup>+</sup> and D<sup>+</sup>-meson p<sub>T</sub>-differential yields starting from the outputs of the [AliPhysics](https://github.com/alisw/AliPhysics) tasks [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx) and [AliAnalysisTaskSEDplus.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDplus.cxx), using rectangular or ML selections

## Run analysis tasks

### Creation of files with selections to be applied on the tasks
* In the [cutobjects](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/cutobjects) folder all the macros needed to produce the cut-object files used in the tasks are stored

### Run D<sup>+</sup> and D<sub>s</sub><sup>+</sup>tasks with private jobs
The [AliAnalysisTaskSEDplus.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDplus.cxx) and [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx) tasks can be run with private jobs using the ```RunAnalysisDsDplusTask.C``` script in the ```runanalysistask``` folder:

```cpp
root -l RunAnalysisDsDplusTask.cc+(TString configfilename = configfile.yml, TString runMode = "full", bool mergeviajdl = true)
```
where ```configfile.yml``` is a configuration file (such as [runAnalysis_config_LHC17p_cent.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/runanalysistask/runAnalysis_config_LHC17p_cent.yml)) with the information about the dataset, the AliPhysics version, and the task options to be used. The tasks options include the possibility to create a tree for the ML studies or apply a ML model trained with [xgboost](https://xgboost.readthedocs.io/en/latest/) or [scikit learn](https://scikit-learn.org/stable/).

### Train output merge
* The by-hand merge of unmerged outputs of a [ALICE analysis train](http://alimonitor.cern.ch/map.jsp) or private jobs can be performed with the script in the ```merge``` folder:
```python
python3 MergeTrainOutputs.py files_to_merge.yml
```
where ```files_to_merge.yml``` is the configuration file containing the information about the outputs that has to be merged such as [files_to_merge_LHC18q.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/merge/files_to_merge_LHC18q.yml)

## Machine Learning analsyis for D-meson candidate selections
*To be added*

## Main analysis with THnSparses

### Pre-filter ThnSparses
* The THnSparse in the task outputs can be pre-filtered to reduce the file size (useful if the train outputs are too large and cannot be merged) with the ```FilterSparse.py``` script in the ```filterdata``` folder:
```python
python3 FilterSparse.py configfile.yml cutset.yml
```
where ```configfile.yml``` is the config file with the info of the input files and ```cutset.yml``` is the set of selections to be applied in the filtering. It creates output files as the input ones, with the ThnSparses filtered. 
With the option ```--suffix suffixname```, a suffix is added to the output file names, otherwise the input files are overwritten.
With the option ```--plot``` it creates control plots that are saved in .pdf files 

### Projection of invariant-mass distributions from THnSparses
* Project the THnSparse with the desired selections into invariant-mass distributions (TH1F):
```python
python3 ProjectDplusDsSparse.py configfile.yml cutset.yml output.root
```
where ```configfile.yml``` is a configuration file with the info of the input files (such as [config_Dplus_pp_data_tree.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/config_Ds_MC_3050.yml)), while ```cutset.yml``` is the set of selections to be applied.

## Main analysis with TTrees or dataframes

### Filter trees to prepare data sets for ML studies
To filter trees produced with the D<sub>s</sub><sup>+</sup>and D<sup>+</sup> tasks and divide each category (data, MC prompt D, MC feed-down D, MC background) in a separated file (tree or dataframe) to prepare the datasets for the ML analyses, the ```FilterTrees4ML.cc``` and ```FilterTrees4ML.py``` scripts in the ```filterdata``` folder can be used:
```cpp
root -l FilterTrees4ML.cc+(TString configfilename = configfile.yml)
```
or 
```python
python3 FilterTrees4ML.py configfile.yml
```
where ```configfile.yml``` is a configuration file (such as [config_Dplus_data_skim_pp5TeV.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/filterdata/config_Dplus_data_skim_pp5TeV.yml)) that contains the information about the input files, the preselections to apply, the features to keep and the output files. The output files are by default ```root``` files. In the case of the python script, if the ```--parquet``` option is used, the output data are saved into ```parquet``` files instead of ```root``` files. 

### Projection of invariant-mass distributions from TTrees
* Project the TTree or dataframe with the desired selections into invariant-mass distributions (TH1F):
```python
python3 ProjectDplusDsTree.py configfile.yml cutset.yml output.root
```
where ```configfile.yml``` is a configuration file with the info of the input files, including the original task output (such as [config_Dplus_pp_data_tree.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/config_Dplus_pp_data_tree.yml)), while ```cutset.yml``` is the set of selections to be applied.
It autodetects whether the input files are ```root``` files containing TTrees or ```parquet``` files containing pandas dataframes.

## Common analysis
The following steps can be performed after having projected THnSparse or TTree (dataframe) objects

### Raw yield extraction
To perform raw-yield extraction either a ROOT or a python script can be used.
* ROOT:
```cpp
root -l GetRawYieldsDplusDs.C+(int cent, bool isMC = false, TString infilename = "distributions.root", TString cfgfilename = "config_Fit.yml", TString outFileName = "output.root")
```
* python:
```python
python3 GetRawYieldsDplusDs.py config_Fit.yml centName distributions.root output.root
```
where ```distributions.root``` is the file obtained projecting the data or MC THnSparse and ```config_Fit.yml``` is a configuration file with the inputs needed to perform the invariant-mass fits such as [config_Ds_Fit.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/fit/config_Ds_Fit.yml) and ```output.root``` is the name of the output ```.root``` file name. In case of the python script, the ```--isMC``` option can be used to specify if the input distributions are from MC simulations and the ```--batch``` option can be used to execute the script in batch mode.

### Efficiency-times-acceptance computation
The efficiency-times-acceptance computation is done in two steps:
* Efficiency computation:
```python
python3 ComputeEfficiencyDplusDs.py config_Fit.yml centName distributionsMC.root output.root
```
where ```distributions.root``` is the file obtained projecting the MC THnSparse and ```config_Fit.yml``` is the same config file used for the raw-yield extraction needed to have the same p<sub>T</sub> binning. The ```--batch``` option can be used to execute the script in batch mode.

To apply p<sub>T</sub> weights the ```--ptweights``` argument followed by the name of the input file with the p<sub>T</sub> weights and the name of the p<sub>T</sub>-weights histogram

* Acceptance and efficiency combination:
```python
python3 CombineAccTimesEff.py effFileName.root accFileName.root outFileName.root
```
where ```effFileName.root``` is the file with the efficiencies computed in the previous step and ```accFileName.root``` is the file with the acceptance computed using the [ComputeAcceptance.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/ComputeAcceptance.C)

both can be run with the ```--batch``` argument to avoid the canvas window

## Standard analysis with theory-driven prompt fraction evaluation
### Cross section

* For the computation of the cross section, a modified version of [HFPtSpectrum.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/HFPtSpectrum.C) present in this repository, is used:

### Nuclear modification factor

* For the computation of the nuclear modification factor, a modified version of [HFPtSpectrumRaa.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/HFPtSpectrumRaa.C) present in this repository, is used

### Corrected yield

* For the computation of the p<sub>T</sub>-differential corrected yields, a modified version of [ComputeDmesonYield.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/ComputeDmesonYield.C) present in this repository, is used

### Run full analysis
* To run the full analysis, from the raw-yield extraction to the nuclear modification factor and the corrected yields the script
```sh
sh RunFullAnalysis.sh
```
can be used by setting some hard-coded parameters

## Analysis with data-driven evaluation of prompt / feed-down fraction
### Prompt / feed-down fraction
* The evaluation of the prompt / feed-down fractions can be performed with the *cut-variation* method with the script:

```python
python3 ComputeCutVarPromptFrac.py cfgFileName.yml outFileName.root
```
where ```cfgFileName.yml``` is a configuration file such as [config_Dplus_PromptFrac_pp5TeV.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/datadrivenfprompt/config_Dplus_PromptFrac_pp5TeV.yml)). The method requires several raw yields and efficiency files obtained with different topological selections applied to enrich/reduce the prompt or the feed-down contribution.

* Cross section
*To be added*

## Significance optimisation
### Optimisation with TTrees
* The script [ScanSelectionsTree.py](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/optimisation/ScanSelectionsTree.py) can be used to compute expected quantities (i.e. signal, background, significance, S/B, prompt fraction) for all combinations of different selection criteria:
```python
python3 ScanSelectionsTree.py cfgFileName.yml outFileName.root
```
where ```cfgFileName.yml``` is a yaml config file containing all the information about the input data to be used and the selections to be tested, such as [config_Dplus_pp5TeV_Optimisation.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/optimisation/config_Dplus_pp5TeV_Optimisation.yml). If the number of variables tested are less or equal 2 (i.e. ML outputs), the script produces plots with expected quantities as a function of the applied selections. In any case, a ntuple with all the expected quantities and the values of applied selections is produced and stored in the output file. 

### Optimisation with THnSparses
*To be updated*

* Compute expected significance for all combinations of different selection criteria:
```python
python3 ScanSignificanceSparse.py configfile.yml output.root
```
where ```configfile.yml``` is a configuration file such as [config_Ds_010_SignOpt.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/config_Ds_010_SignOpt.yml) 

* Project ntuple with expected significance as a function of relevant variables:
```python
python3 ProjectSignifNtuple.py configfile.yml input.root PtMin PtMax minSignificance maxSignificance minEffPrompt maxEffPrompt
```
where the input file ```input.root``` is the one produced in the previous step

## Systematic uncertainties
### Selection efficiency
* For the cut-variation studies with THnSparses, the configuration files for each set of selection criteria can be created using:
```python
python3 MakeCutsFilesForSyst.py
```
the variables and the ranges should be set hard-coded in the script
Once the configuration files are created they can be used to repeat the main analysis with the different selection criteria

### Raw-yield extraction
* For the raw-yield extraction uncertainty a multi-trial study can be run with:
```cpp
root -l RawYieldSystematics.C+(TString outfilerawname = "output.root")
```

## Test and validation of alternative code for production of TTrees ([AliAnalysisTaskSEHFTreeCreator.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/treeHF/AliAnalysisTaskSEHFTreeCreator.cxx))
The validation of the code for production of trees used in ML studies can be done using the scripts in the ```runanalysistask``` folder

* To run the [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx) and [AliAnalysisTaskSEHFTreeCreator.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/treeHF/AliAnalysisTaskSEHFTreeCreator.cxx) on the same files:
```cpp
root -l RunAnalysisTreeCreator.cc+(TString configfilename = configfile.yml, TString runMode = "full", bool mergeviajdl = true)
```
where ```configfile.yml``` is a configuration file with the information about the dataset and the AliPhysics version to be used

* To run the validation of the output:
```python
python3 ValidateTreeCreator.py inputfile inputdir inputlist
```
where ```inputfile```, ```inputdir```, and ```inputlist``` are the root file produced by ```RunAnalysisTreeCreator.C```, the name of the TDirectoryFile and the TList inside the root file
