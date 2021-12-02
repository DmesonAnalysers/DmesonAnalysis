[![](https://img.shields.io/github/license/DmesonAnalysers/DmesonAnalysis?color=blue)](https://github.com/DmesonAnalysers/DmesonAnalysis/blob/master/LICENSE)
![](https://img.shields.io/github/languages/count/DmesonAnalysers/DmesonAnalysis?color=green)
![](https://img.shields.io/github/last-commit/DmesonAnalysers/DmesonAnalysis?color=red)

# DmesonAnalysis

Code for the measurement of D<sup>0</sup>, D<sup>+</sup>, D<sub>s</sub><sup>+</sup>, D<sup>*+</sup>, and Λ<sub>c</sub><sup>+</sup>-baryon p<sub>T</sub>-differential yields starting from the outputs of the [AliPhysics](https://github.com/alisw/AliPhysics) tasks [AliAnalysisTaskSEDs.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDs.cxx), [AliAnalysisTaskSEDplus.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDplus.cxx), [AliAnalysisTaskSENonPromptLc.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSENonPromptLc.cxx), and [AliAnalysisTaskSEDmesonTree.cxx](https://github.com/alisw/AliPhysics/blob/master/PWGHF/vertexingHF/vHFML/AliAnalysisTaskSEDmesonTree.cxx) using rectangular or ML selections performed with the [hipe4ml](https://github.com/hipe4ml/hipe4ml) package.

## Prerequisites
The code in this repository requires to have installed: 
- [python3](https://www.python.org)(>=3.6)
- [AliPhysics](https://github.com/alisw/AliPhysics)
- [hipe4ml](https://github.com/hipe4ml/hipe4ml)(>=0.0.10 or installed from dev branch)
- [pyaml](https://pypi.org/project/pyaml)
- [alive_progress](https://github.com/rsalmei/alive-progress)

## Documentation
https://dmesonanalysers.github.io/DmesonAnalysis
