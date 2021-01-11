# D<sub>s</sub><sup>+</sup>, D<sup>+</sup>-meson, and Λ<sub>c</sub><sup>+</sup> Analysis code

Code for the measurement of D<sub>s</sub><sup>+</sup>, D<sup>+</sup>-meson p<sub>T</sub>, and Λ<sub>c</sub><sup>+</sup>-hadron p<sub>T</sub>-differential yields starting from the outputs of the [AliPhysics](https://github.com/alisw/AliPhysics) tasks [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx), [AliAnalysisTaskSEDplus.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDplus.cxx), and [AliAnalysisTaskSENonPromptLc.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSENonPromptLc.cxx), using rectangular or ML selections

## Run analysis tasks

### Creation of files with selections to be applied on the tasks
* In the [cutobjects](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/cutobjects) folder all the macros needed to produce the cut-object files used in the tasks are stored

### Run D<sup>+</sup>, D<sub>s</sub><sup>+</sup>, and Λ<sub>c</sub><sup>+</sup> tasks with private jobs
The [AliAnalysisTaskSEDplus.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDplus.cxx), [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx) and [AliAnalysisTaskSENonPromptLc.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSENonPromptLc.cxx) tasks can be run with private jobs using the ```RunAnalysisDplusDsLcTask.cc``` script in the ```runanalysistask``` folder:

```cpp
root -l RunAnalysisDplusDsLcTask.cc+(TString configfilename = configfile.yml, TString runMode = "full", bool mergeviajdl = true)
```
where ```configfile.yml``` is a configuration file (such as [runAnalysis_config_LHC17p_cent.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/runanalysistask/runAnalysis_config_LHC17p_cent.yml)) with the information about the dataset, the AliPhysics version, and the task options to be used. The tasks options include the possibility to create a tree for the ML studies or apply a ML model trained with [xgboost](https://xgboost.readthedocs.io/en/latest/) or [scikit learn](https://scikit-learn.org/stable/). The ML model application is not supported by the Λ<sub>c</sub><sup>+</sup> task.

### Train output merge
* The by-hand merge of unmerged outputs of a [ALICE analysis train](http://alimonitor.cern.ch/map.jsp) or private jobs can be performed with the script in the ```merge``` folder:
```python3
python3 MergeTrainOutputs.py files_to_merge.yml
```
where ```files_to_merge.yml``` is the configuration file containing the information about the outputs that has to be merged such as [files_to_merge_LHC18q.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/merge/files_to_merge_LHC18q.yml)

## Main analysis with THnSparses

### Pre-filter ThnSparses
* The THnSparse in the task outputs can be pre-filtered to reduce the file size (useful if the train outputs are too large and cannot be merged) with the ```FilterSparse.py``` script in the ```filterdata``` folder:
```python3
python3 FilterSparse.py configfile.yml cutset.yml
```
where ```configfile.yml``` is the config file with the info of the input files and ```cutset.yml``` is the set of selections to be applied in the filtering. It creates output files as the input ones, with the ThnSparses filtered. 
With the option ```--suffix suffixname```, a suffix is added to the output file names, otherwise the input files are overwritten.
With the option ```--plot``` it creates control plots that are saved in .pdf files 

### Projection of invariant-mass distributions from THnSparses
* Project the THnSparse with the desired selections into invariant-mass distributions (TH1F):
```python3
python3 ProjectDplusDsSparse.py configfile.yml cutset.yml output.root
```
where ```configfile.yml``` is a configuration file with the info of the input files (such as [config_Dplus_pp_data_tree.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/config_Ds_MC_3050.yml)), while ```cutset.yml``` is the set of selections to be applied.

To apply *p*<sub>T</sub> weights in case of MC the ```--ptweights``` argument followed by the name of the input file with the *p*<sub>T</sub> weights and the name of the *p*<sub>T</sub>-weights histogram should be parsed. In this case, the *p*<sub>T</sub> weights are applied to both the prompt and the FD distributions. If also the ```--ptweightsB``` argument followed by the name of the input file with the *p*<sub>T</sub><sup>B</sup> weights and the name of the *p*<sub>T</sub><sup>B</sup>-weights histogram is parsed, the *p*<sub>T</sub> weights for the FD are computed from the B-mother *p*<sub>T</sub>

## Main analysis with TTrees or dataframes

### Filter trees to prepare data sets for ML studies
To filter trees produced with the D<sub>s</sub><sup>+</sup>and D<sup>+</sup> tasks and divide each category (data, MC prompt D, MC feed-down D, MC background) in a separated file (tree or dataframe) to prepare the datasets for the ML analyses, the ```FilterTrees4ML.cc``` and ```FilterTrees4ML.py``` scripts in the ```filterdata``` folder can be used:
```cpp
root -l FilterTrees4ML.cc+(TString configfilename = configfile.yml)
```
or 
```python3
python3 FilterTrees4ML.py configfile.yml
```
where ```configfile.yml``` is a configuration file (such as [config_Dplus_data_skim_pp5TeV.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/filterdata/config_Dplus_data_skim_pp5TeV.yml)) that contains the information about the decay channel, the input files, the preselections to apply, the features to keep and the output files. The output files are by default ```root``` files in the c++ script and ```parquet``` in the python script. If the ```--root``` option is used, the output data are saved into ```root``` files instead of ```parquet``` files. 

## Machine Learning analsyis for D-meson candidate selections
*To be added*

### Projection of invariant-mass distributions from TTrees
* Project the TTree or dataframe with the desired selections into invariant-mass distributions (TH1F):
```python3
python3 ProjectDplusDsTree.py configfile.yml cutset.yml output.root
```
where ```configfile.yml``` is a configuration file with the info of the input files, including the original task output (such as [config_Dplus_pp_data_tree.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/config_Dplus_pp_data_tree.yml)), while ```cutset.yml``` is the set of selections to be applied.
It autodetects whether the input files are ```root``` files containing TTrees or ```parquet``` files containing pandas dataframes.

To apply *p*<sub>T</sub> weights in case of MC the ```--ptweights``` argument followed by the name of the input file with the *p*<sub>T</sub> weights and the name of the *p*<sub>T</sub>-weights histogram should be parsed. In this case, the *p*<sub>T</sub> weights are applied to both the prompt and the FD distributions. If also the ```--ptweightsB``` argument followed by the name of the input file with the *p*<sub>T</sub><sup>B</sup> weights and the name of the *p*<sub>T</sub><sup>B</sup>-weights histogram is parsed, the *p*<sub>T</sub> weights for the FD are computed from the B-mother *p*<sub>T</sub>

## Common analysis
The following steps can be performed after having projected THnSparse or TTree (dataframe) objects

### Raw yield extraction
To perform raw-yield extraction either a ROOT or a python script can be used.
* c++:
```cpp
root -l GetRawYieldsDplusDs.C+(int cent, bool isMC = false, TString infilename = "distributions.root", TString cfgfilename = "config_Fit.yml", TString outFileName = "output.root")
```
* python:
```python3
python3 GetRawYieldsDplusDs.py config_Fit.yml centName distributions.root output.root
```
where ```distributions.root``` is the file obtained projecting the data or MC THnSparse and ```config_Fit.yml``` is a configuration file with the inputs needed to perform the invariant-mass fits such as [config_Ds_Fit.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/fit/config_Ds_Fit.yml) and ```output.root``` is the name of the output ```.root``` file name. In case of the python script, the ```--isMC``` option can be used to specify if the input distributions are from MC simulations and the ```--batch``` option can be used to execute the script in batch mode.

### Efficiency-times-acceptance computation
The efficiency-times-acceptance computation is done in two steps:
* Efficiency computation:
```python3
python3 ComputeEfficiencyDplusDs.py config_Fit.yml centName distributionsMC.root output.root
```
where ```distributions.root``` is the file obtained projecting the MC THnSparse and ```config_Fit.yml``` is the same config file used for the raw-yield extraction needed to have the same *p*<sub>T</sub> binning. The ```--batch``` option can be used to execute the script in batch mode.

* Acceptance and efficiency combination:
```python3
python3 CombineAccTimesEff.py effFileName.root accFileName.root outFileName.root
```
where ```effFileName.root``` is the file with the efficiencies computed in the previous step and ```accFileName.root``` is the file with the acceptance computed using the [ComputeAcceptance.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/ComputeAcceptance.C)

both can be run with the ```--batch``` argument to avoid the canvas window

## Standard analysis with theory-driven prompt fraction evaluation
### Cross section

* For the computation of the cross section, a modified version of [HFPtSpectrum.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/HFPtSpectrum.C) present in this repository, is used

### Nuclear modification factor

* For the computation of the nuclear modification factor, a modified version of [HFPtSpectrumRaa.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/HFPtSpectrumRaa.C) present in this repository, is used

### Corrected yield

* For the computation of the *p*<sub>T</sub>-differential corrected yields, a modified version of [ComputeDmesonYield.C](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/macros/ComputeDmesonYield.C) present in this repository, is used

## Analysis with data-driven evaluation of prompt / feed-down fraction
### Prompt / feed-down fraction
* The evaluation of the prompt / feed-down fractions can be performed with the *cut-variation* method with the script:

```python3
python3 ComputeCutVarPromptFrac.py cfgFileName.yml outFileName.root
```
where ```cfgFileName.yml``` is a configuration file such as [config_Dplus_PromptFrac_pp5TeV.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/configfiles/datadrivenfprompt/config_Dplus_PromptFrac_pp5TeV.yml)). The method requires several raw yields and efficiency files obtained with different topological selections applied to enrich/reduce the prompt or the feed-down contribution.

### Cross section
* The computation of the prompt / feed-down *p*<sub>T</sub>-differential cross sections can be performed with the script:
```python3
python3 ComputeDataDrivenCrossSection.py rawYieldFile.root effAccFile.root fracFile.root outFile.root [--prompt] [--FD] [--Dplus] [--Ds] [--system] [--energy] [--batch]
```
where ```rawYieldFile.root```, ```effAccFile.root```, ```fracFile.root``` are the root files containing the raw yields, the acceptance-times-efficiency factors, and the fraction of prompt (feed-down) D mesons estimated with the cut-variation method (previous paragraph), while ```outFile.root``` is the ROOT output file. The optional parameters are needed to define wether the prompt or the feed-down cross section should be computed for the D<sub>s</sub><sup>+</sup> or D<sup>+</sup> meson, the system (```pp``` or ```Pb-Pb```) and the centre-of-mass energy.

## Run full analysis
* To run the full analysis escept for the ML part, from the raw-yield extraction to the nuclear modification factor and the corrected yields, the script
```sh
sh RunFullAnalysis.sh
```
can be used by setting some hard-coded parameters

## Significance optimisation
### Optimisation with TTrees
* The script [ScanSelectionsTree.py](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/optimisation/ScanSelectionsTree.py) can be used to compute expected quantities (i.e. signal, background, significance, S/B, prompt fraction) for all combinations of different selection criteria:
```python3
python3 ScanSelectionsTree.py cfgFileName.yml outFileName.root
```
where ```cfgFileName.yml``` is a yaml config file containing all the information about the input data to be used and the selections to be tested, such as [config_Dplus_pp5TeV_Optimisation.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/optimisation/config_Dplus_pp5TeV_Optimisation.yml). If the number of variables tested are less or equal 2 (i.e. ML outputs), the script produces plots with expected quantities as a function of the applied selections. In any case, a ntuple with all the expected quantities and the values of applied selections is produced and stored in the output file. 

## Systematic uncertainties
All the code for the evaluation of the systematic uncertainties is in the [systematics](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/systematics/) directory.
### Selection efficiency
* For the cut-variation studies with TTrees or THnSparses, the configuration files for each set of selection criteria can be created using:
```python3
python3 MakeCutsFilesForSyst.py
```
the variables and the ranges should be set hard-coded in the script.
Once the configuration files are created they can be used to repeat the main analysis with the different selection criteria. 

Once the analysis has been repeated for all the sets of selections, the systematic uncertainty can be evaluated using the script in the [systematics/seleff](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/systematics/seleff/) directory:
```cpp
root -l PlotCutVariationsOnePtBin.cc+(TString cfgFileName = "cfgFile.yml")
```
where the config file ```cfgFile.yml``` includes all the information of the sets of selections to be used and the quality criteria that has to be applied, such as [config_cutvar_DsFD_pp.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/systematics/seleff/config_cutvar_DsFD_pp.yml).

### Raw-yield extraction
* For the raw-yield extraction uncertainty a multi-trial study (with the usage of [AliHFInvMassMultiTrialFit.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/AliHFInvMassMultiTrialFit.cxx)) can be run with the script in the [systematics/rawyields](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/systematics/rawyields/) directory:
```cpp
root -l RawYieldSystematics.cc+(TString cfgFileName = "cfgFile.yml")
```
where the config file ```cfgFile.yml``` includes all the information of the variations that has to be applied, such as [config_multi_trial_DplusFD_pp.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/systematics/rawyields/config_multi_trial_DsFD_pp.yml)

### Generated MC *p*<sub>T</sub> shape
* The systematic uncertainty arising from the shape of the *p*<sub>T</sub> distributions in the MC simulation can be evaluated with the code in the [systematics/genptshape](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/systematics/genptshape/) directory.
    * The first step is the computation of the *p*<sub>T</sub> weights:
    
    ```python3
    python3 ComputePtGenShapeWeights.py inFileMC.root outFile.root [--Dspecie Dname] [--Bspecie Bname] [--PbPb] [--rebin] [--smooth]
    ```
    where the root file ```inFileMC.root``` can be the output of the [AliAnalysisTaskCheckHFMCProd.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/AliAnalysisTaskCheckHFMCProd.cxx) task or [AliCFTaskVertexingHF.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/AliCFTaskVertexingHF.cxx), ```outFile.root``` is the output file name, ```--Dspecie``` and ```--Bspecie``` is the argument to chose the D-meson and B-meson species to use, ```--PbPb``` is a flag to enable in case of Pb-Pb analysis while ```--rebin``` and ```--smooth``` are two flags to apply a rebin of the spectra and a smoothening of the weights.

    * The second step step is the computation of the efficiencies with and without *p*<sub>T</sub> weights, as described in the dedicated section.

    * The second step step is the evaluation of the systematic uncertainty:
    ```python3
    python3 GetPtWeightSyst.py cfgFileName.yml
    ```
    where ```cfgFileName.yml``` is a config file as [config_ptshape_syst.yml.yml](https://github.com/DmesonAnalysers/DmesonAnalysis/tree/master/systematics/genptshape/config_ptshape_syst.yml.yml)

## Test and validation of alternative code for production of TTrees ([AliAnalysisTaskSEHFTreeCreator.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/treeHF/AliAnalysisTaskSEHFTreeCreator.cxx))
The validation of the code for production of trees used in ML studies can be done using the scripts in the ```runanalysistask``` folder

* To run the [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx) and [AliAnalysisTaskSEHFTreeCreator.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/treeHF/AliAnalysisTaskSEHFTreeCreator.cxx) on the same files:
```cpp
root -l RunAnalysisTreeCreator.cc+(TString configfilename = configfile.yml, TString runMode = "full", bool mergeviajdl = true)
```
where ```configfile.yml``` is a configuration file with the information about the dataset and the AliPhysics version to be used

* To run the validation of the output:
```python3
python3 ValidateTreeCreator.py inputfile inputdir inputlist
```
where ```inputfile```, ```inputdir```, and ```inputlist``` are the root file produced by ```RunAnalysisTreeCreator.C```, the name of the TDirectoryFile and the TList inside the root file
