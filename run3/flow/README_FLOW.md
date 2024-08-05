# Flow Analysis in Run 3

Code for the measurement of D meson flow starting from the outputs of the [taskFlowCharmHadrons.cxx](https://github.com/AliceO2Group/O2Physics/blob/master/PWGHF/D2H/Tasks/taskFlowCharmHadrons.cxx) task using rectangular or ML selections performed with the [hipe4ml](https://github.com/hipe4ml/hipe4ml) package.

## Analysis steps
The analysis is organized as follows:
- `compute_reso.py`: compute SP/EP resolution term
- `project_thnsparse.py`: project ThnSparse as a function of *p*<sub>T</sub>/centrality from AnalysisResults.root
- `get_vn_vs_mass.py`: Compute the `vn` with different methods (simultaneous fit to the projected inv. mass and EP/SP, fit to invariant mass in-plane vs. out-of-plane)
- `compute_efficiency.py`: Compute the `skip_efficiency` with from the AnalysisResults.root from the D meson task (NOT the flow task)

The full analysis chain can be run with ``run_full_flow_analysis.py``.

The python scripts necessitate a configuration file like ``config_flow.yml``.

For instance, the code can be run with the following command:

`python3 run_full_flow_analysis.py config_flow.yml path_to_AnalysisResults [OPTIONS]`

where `[OPTIONS]` are:
- `--centrality` or `-c` (`-c k3050`): set the centrality class
- `--resolution` or `-r` (`-r 1.`): set the resolution value or the path to the resolution file computed with `compute_reso.py`
- `--outputdir` or `-o` (`-o output_directory`): set the output directory (subdirectories will be automatically made)
- `--suffix` or `-s` (`-s suffix`): set the suffix of the output (if the input AnalysisResults.root already has a suffix, this will be used as a default)
- `--vn_method` or `-v` (`-vn sp`): set the `vn` method. Only `sp` (default), `ep`, `deltaphi` are available.
- `--wagon_id` or `-w` (`-w 0000`): set the index of the wagon inside the AnalysisResults.root to be considered. As a default, this is not needed
- `--skip_resolution`: avoid resolution estimation
- `--skip_projection`: avoid ThnSparse projection
- `--skip_vn`: avoid `vn` estimation
- `--skip_efficiency`: avoid `skip_efficiency` estimation