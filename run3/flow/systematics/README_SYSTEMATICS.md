# Multitrial systematic
Code for the multitrial systematic starting from the outputs of the [taskFlowCharmHadrons.cxx](https://github.com/AliceO2Group/O2Physics/blob/master/PWGHF/D2H/Tasks/taskFlowCharmHadrons.cxx) task and the raw yield .root file from ``run_full_flow_analysis.py``.


## Analysis steps
The analysis is organized as follows:
- `modifications_config.yml`: Config file with the variables to be changed in multitrial or resolution systematics. All possible combination of the variations are considered. For each varation, a dedicated config file is saved in the given output directory.

### Multitrial
- `multitrial.sh`: Starting from a set of varaibles the user MUST set, the analysis chain up to the raw yield extraction is performed for each configuration file. Finally, the systematic is computed with the `compute_syst_multitrial.py` script 
- `compute_syst_multitrial.py`: Collect each raw yield file produced in the previous step and computes the systematic. A default cut on significance (higher than 3) and chi2 (lower than 2) are applied on the inputs. The script can be run stand alone with the following command: 

`python3 compute_syst_multitrial.py /path/to/multitrial/ry/folder/ /path/to/default/ry/ -o outputdir`

where `[OPTIONS]` are:
- `--outputdir` or `-o` (`-o output_directory`): set the output directory


### Resolution
- `reso_syst.sh`: Starting from a set of varaibles the user MUST set, the analysis chain up to the raw yield extraction for each configuration file. Finally, the systematic is computed with the `compute_syst_reso.py` script 
- `compute_syst_reso.py`: Collect each raw yield file produced in the previous step and computes the systematic. The script can be run stand alone with the following command: 
Compute the `vn` with different methods (simultaneous fit to the projected inv. mass and EP/SP, fit to invariant mass in-plane vs. out-of-plane)

`python3 compute_syst_reso.py /config/with/variations /config/default /path/to/default/ry/ /path/to/resolution/file/reso.root -c centrality_bin -o outputdir`

where `[OPTIONS]` are:
- `--outputdir` or `-o` (`-o output_directory`): set the output directory
- `--centrality` or `-c` (`-o kllhh`): set the centrality bin

