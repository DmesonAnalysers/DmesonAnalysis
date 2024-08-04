# Multitrial systematic
Code for the multitrial systematic starting from the outputs of the [taskFlowCharmHadrons.cxx](https://github.com/AliceO2Group/O2Physics/blob/master/PWGHF/D2H/Tasks/taskFlowCharmHadrons.cxx) task and the raw yield .root file from ``run_full_flow_analysis.py``.


## Analysis steps
The analysis is organized as follows:
- `modifications_config.yml`: Config fiel with the variables to be changed in multitrial. All possible combination are considered. For each varation, a dedicated config file is saved
- `multitrial.sh`: Starting from a set of varaibles the user MUST set, the analysis chain up to the raw yield extraction. Finally, the systematic is computed with the `compute_syst_multitrial.py` script 
- `compute_syst_multitrial.py`: Collect each raw yield filed produced in the previous step and computes the systematic. A default cut on significance (higher than 3) and chi2 (lower than 2) are applied on the inputs. The script can be run stand alone with the following command: 
Compute the `vn` with different methods (simultaneous fit to the projected inv. mass and EP/SP, fit to invariant mass in-plane vs. out-of-plane)

The full analysis chain can be run with ``run_full_flow_analysis.py``.

The python scripts necessitate a configuration file like ``config_flow.yml``.

For instance, the code can be run with the following command:

`python3 compute_syst_multitrial.py /path/to/multitrial/ry/folder/ /path/to/default/ry/raw_yields.root -o outputdir`

where `[OPTIONS]` are:
- `--outputdir` or `-o` (`-o output_directory`): set the output directory