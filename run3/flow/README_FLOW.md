# Flow Analysis in Run 3

Code for the measurement of D meson flow starting from the outputs of the [taskFlowCharmHadrons.cxx](https://github.com/AliceO2Group/O2Physics/blob/master/PWGHF/D2H/Tasks/taskFlowCharmHadrons.cxx) task using rectangular or ML selections performed with the [hipe4ml](https://github.com/hipe4ml/hipe4ml) package.

## Prerequisites
The code in this repository requires having installed: 
- [python3](https://www.python.org)(>=3.6)
- [AliPhysics](https://github.com/alisw/AliPhysics) (until simultaneous fit is ported to flarefly/O2Physics)
- [hipe4ml](https://github.com/hipe4ml/hipe4ml)(>=0.0.10 or installed from dev branch)
- [pyaml](https://pypi.org/project/pyaml)
- [alive_progress](https://github.com/rsalmei/alive-progress)

## Analysis flow
The analysis flow is organized as follows:
.                
└── flow                        # The folder containing flow analysis material
    ├── compute_reso.py         # Compute SP/EP resolution term
    ├── project_thnsparse.py    # Projcect ThnSparse as a function of pT/centrality from AnalysisResults.root
    ├── get_vn_vs_mass.py       # Simultanoues fit to the project inv. mass and EP/SP

The full analysis chain can be run with run_full_flow_analysis.py. 

The python scripts necessitate of a configuration file like config_flow.yml
