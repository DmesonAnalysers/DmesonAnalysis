# Flow Analysis in Run 3

Code for the measurement of D meson flow starting from the outputs of the [taskFlowCharmHadrons.cxx](https://github.com/AliceO2Group/O2Physics/blob/master/PWGHF/D2H/Tasks/taskFlowCharmHadrons.cxx) task using rectangular or ML selections performed with the [hipe4ml](https://github.com/hipe4ml/hipe4ml) package.

## Analysis steps
The analysis is organized as follows:
- `compute_reso.py: compute SP/EP resolution term
- `project_thnsparse.py`: projcect ThnSparse as a function of pT/centrality from AnalysisResults.root
- `get_vn_vs_mass.py`: Simultaneous fit to the projected inv. mass and EP/SP

The full analysis chain can be run with ``run_full_flow_analysis.py``.

The python scripts necessitate a configuration file like ``config_flow.yml``.