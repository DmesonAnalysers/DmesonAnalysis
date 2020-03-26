#!/bin/bash
#steps to be performed
DoDataProjection=true
DoMCProjection=true
DoDataRawYields=true
DoMCRawYields=false
DoEfficiency=true
DoAccEff=true
DoHFPtSpec=false
DoHFPtSpecRaa=false
DoDmesonYield=false
DoDataDrivenCrossSection=false

#wheter you are projecting a tree or a sparse
ProjectTree=true

#PARAMETERS TO BE SET (use "" for parameters not needed)
################################################################################################
Meson="Ds" # whether it is Dplus or Ds analysis
Cent="kpp5TeVFD" # used also to asses prompt or non-prompt and system

cfgFileData="configfiles/config_Ds_pp_data_tree.yml"
cfgFileMC="configfiles/config_Ds_pp_MC_tree.yml"
cfgFileFit="configfiles/fit/config_Ds_Fit_pp5TeV.yml"

accFileName="accfiles/Acceptance_Toy_DsKKpi_yfidPtDep_etaDau09_ptDau100_FONLL5ptshape.root"
predFileName="models/D0DplusDstarPredictions_502TeV_y05_noYShift_all_191017_BDShapeCorrected.root"
pprefFileName="" #"ppreference/Ds_ppreference_pp5TeV_noyshift_pt_2_3_4_6_8_12_16_24_36_50.root"

PtWeightsFileName="" #"ptweights/PtWeigths_LHC19c3b.root"
PtWeightsHistoName="" #"hPtWeightsFONLLtimesTAMUcent"

DataDrivenFractionFileName=""

#assuming cutsets config files starting with "cutset" and are .yml
CutSetsDir="configfiles/cutsets/Ds/pp/data_driven_fprompt_FDen_cuts"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

OutDirRawyields="../../Analyses/pp5TeV/Ds_wML_mult/outputs/100320/data_driven_fprompt/raw_yield"
OutDirEfficiency="../../Analyses/pp5TeV/Ds_wML_mult/outputs/100320/data_driven_fprompt/eff_testw"
OutDirCrossSec=""
OutDirRaa=""
################################################################################################

if [ ${Meson} != "Dplus" ] && [ ${Meson} != "Ds" ]; then
  echo $(tput setaf 1) ERROR: only Ds and Dplus mesons are supported! $(tput sgr0)
  exit 2
fi

if [ ${Cent} != "k010" ] && [ ${Cent} != "k3050" ] && [ ${Cent} != "k6080" ] && [ ${Cent} != "kpp5TeVFD" ] && [ ${Cent} != "kpp5TeVPrompt" ]; then
  echo $(tput setaf 1) ERROR: system ${Cent} is not supported! $(tput sgr0)
  exit 2
fi

if [ ! -f "${cfgFileData}" ]; then
  echo $(tput setaf 1) ERROR: data config file "${cfgFileData}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ! -f "${cfgFileMC}" ]; then
  echo $(tput setaf 1) ERROR: MC config file "${cfgFileMC}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ! -f "${cfgFileFit}" ]; then
  echo $(tput setaf 1) ERROR: fit config file "${cfgFileFit}"does not exist! $(tput sgr0)
  exit 2
fi

if [ ! -f "${accFileName}" ]; then
  echo $(tput setaf 1) ERROR: acceptance file "${accFileName}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ! -f "${predFileName}" ]; then
  echo $(tput setaf 1) ERROR: FONLL file "${predFileName}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ! -f "${pprefFileName}" ] && [ "${pprefFileName}" != "" ]; then
  echo $(tput setaf 1) ERROR: pp reference file "${pprefFileName}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ! -f "${DataDrivenFractionFileName}" ] && [ ${DoDataDrivenCrossSection} ]; then
  echo $(tput setaf 1) ERROR: data-driven fraction file "${DataDrivenFractionFileName}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ! -f "${PtWeightsFileName}" ] && [ "${PtWeightsFileName}" != "" ]; then
  echo $(tput setaf 3) WARNING: pT-weights file "${PtWeightsFileName}" does not exist! $(tput sgr0)
fi

if [ ! -d "${OutDirRawyields}" ]; then
  mkdir ${OutDirRawyields}
fi

if [ ! -d "${OutDirEfficiency}" ]; then
  mkdir ${OutDirEfficiency}
fi

if [ ! -d "${OutDirCrossSec}" ] && [ "${OutDirCrossSec}" != "" ]; then
  mkdir ${OutDirCrossSec}
fi

if [ ! -d "${OutDirRaa}" ] && [ "${OutDirRaa}" != "" ]; then
  mkdir ${OutDirRaa}
fi

#get prompt or non-prompt, system
isPrompt=true
System=pp
if [ ${Cent} == "kpp5TeVFD" ]; then
  isPrompt=false
  System=pp
elif [ ${Cent} == "kpp5TeVPrompt" ]; then
  isPrompt=true
  System=pp
elif [ ${Cent} == "k010" ]; then
  isPrompt=true
  System=PbPb
elif [ ${Cent} == "k3050" ]; then
  isPrompt=true
  System=PbPb
elif [ ${Cent} == "k6080" ]; then
  isPrompt=true
  System=PbPb
fi

#project sparses or trees
ProjectScript="ProjectDplusDsSparse.py"
if $ProjectTree; then
  ProjectScript="ProjectDplusDsTree.py"
fi

if $DoDataProjection; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Projecting data distributions $(tput sgr0)
    python3 ${ProjectScript} ${cfgFileData} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml ${OutDirRawyields}/Distr_${Meson}_data${CutSets[$iCutSet]}.root
  done
fi

if $DoMCProjection; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Projecting MC distributions $(tput sgr0)
    python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Meson}_MC${CutSets[$iCutSet]}.root
  done
fi

#compute raw yields
if $DoMCRawYields; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Extract raw yields from ${OutDirEfficiency}/Distr_${Meson}_MC${CutSets[$iCutSet]}.root $(tput sgr0)
    echo '.x GetRawYieldsDplusDs.C+('${Cent}',true, "'${OutDirEfficiency}'/Distr_'${Meson}'_MC'${CutSets[$iCutSet]}'.root", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYields'${Meson}'_MC'${CutSets[$iCutSet]}'.root")' | root -l -b
    echo '.q'
  done
fi

if $DoDataRawYields; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Extract raw yields from ${OutDirRawyields}/Distr_${Meson}_data${CutSets[$iCutSet]}.root $(tput sgr0)
    echo '.x GetRawYieldsDplusDs.C+('${Cent}',false, "'${OutDirRawyields}'/Distr_'${Meson}'_data'${CutSets[$iCutSet]}'.root", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYields'${Meson}${CutSets[$iCutSet]}'.root")' | root -l -b
    echo '.q'
  done
fi

#compute efficiency
if $DoEfficiency; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute efficiency from ${OutDirEfficiency}/Distr_${Meson}_MC${CutSets[$iCutSet]}.root $(tput sgr0)
    if [ "${PtWeightsFileName}" == "" ] || [ "${PtWeightsHistoName}" == "" ]; then
      python3 ComputeEfficiencyDplusDs.py ${cfgFileFit} ${Cent} ${OutDirEfficiency}/Distr_${Meson}_MC${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Efficiency_${Meson}${CutSets[$iCutSet]}.root --batch
    else
      echo $(tput setaf 6) Using ${PtWeightsHistoName} pt weights from ${PtWeightsFileName}
      python3 ComputeEfficiencyDplusDs.py ${cfgFileFit} ${Cent} ${OutDirEfficiency}/Distr_${Meson}_MC${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Efficiency_${Meson}${CutSets[$iCutSet]}.root --ptweights ${PtWeightsFileName} ${PtWeightsHistoName} --batch
    fi
  done
fi

#compute efficiency times acceptance
if $DoAccEff; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute efficiency times acceptance $(tput sgr0)
    python3 CombineAccTimesEff.py ${OutDirEfficiency}/Efficiency_${Meson}${CutSets[$iCutSet]}.root ${accFileName} ${OutDirEfficiency}/Eff_times_Acc_${Meson}${CutSets[$iCutSet]}.root --batch
  done
fi

#compute cross section
if $DoDataDrivenCrossSection; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute cross section $(tput sgr0)
    if $isPrompt; then
      python3 ComputeDataDrivenCrossSection.py ${OutDirRawyields}/RawYields${Meson}${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Eff_times_Acc_${Meson}${CutSets[$iCutSet]}.root ${DataDrivenFractionFileName} ${OutDirCrossSec}/CrossSection${Meson}${CutSets[$iCutSet]}.root --prompt --${Meson} --system ${System} --energy 5.02 --batch
    else
      python3 ComputeDataDrivenCrossSection.py ${OutDirRawyields}/RawYields${Meson}${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Eff_times_Acc_${Meson}${CutSets[$iCutSet]}.root ${DataDrivenFractionFileName} ${OutDirCrossSec}/CrossSection${Meson}${CutSets[$iCutSet]}.root --FD --${Meson} --system ${System} --energy 5.02 --batch
    fi
  done
fi

#compute HFPtSpectrum
if $DoHFPtSpec; then
  Channel=""
  if $Meson == "Ds"; then
    Channel="kDsKKpi"
  else
    Channel="kDplusKpipi"
  fi
  
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute HFPtspectrum $(tput sgr0)
    echo '.x HFPtSpectrum.C+('${Channel}',"'${predFileName}'","'${OutDirEfficiency}'/Eff_times_Acc_'${Meson}${CutSets[$iCutSet]}'.root","'${OutDirRawyields}'/RawYields'${Meson}${CutSets[$iCutSet]}'.root","hRawYields","hAccEffPrompt","hAccEffFD","hEvForNorm","'${OutDirCrossSec}'/HFPtSpectrum'${Meson}${CutSets[$iCutSet]}'.root",kNb,1.,true,'${Cent}',k2018)' | root -l -b
    echo '.q'
  done
fi

#compute HFPtSpectrumRaa
if $DoHFPtSpecRaa; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute HFPtspectrumRaa $(tput sgr0)
    echo '.x HFPtSpectrumRaa.C+("'${pprefFileName}'","'${OutDirCrossSec}'/HFPtSpectrum'${Meson}${CutSets[$iCutSet]}'.root","'${OutDirRaa}'/HFPtSpectrumRaa'${Meson}${CutSets[$iCutSet]}'.root",4,1,kNb,'${Cent}',k2018,k5dot023,1./3,3,6,false,1)' | root -l -b
    echo '.q'
  done
fi

#compute corrected yield
if $DoDmesonYield; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute corrected yield $(tput sgr0)
    echo '.x ComputeDmesonYield.C+(k'${Meson}','${Cent}',2,1,"'${pprefFileName}'","'${OutDirCrossSec}'/HFPtSpectrum'${Meson}${CutSets[$iCutSet]}'.root","","'${OutDirRaa}'/HFPtSpectrumRaa'${Meson}${CutSets[$iCutSet]}'.root","","'${OutDirCrossSec}'","'${CutSets[$iCutSet]}'",1,1./3,3,false,1)' | root -l -b
    echo '.q'
  done
fi
