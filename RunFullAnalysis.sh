#!/bin/bash
#steps to be performed
DoDataProjection=true
DoMCProjection=true
DoDataRawYields=true
DoMCRawYields=false
DoEfficiency=true
DoAccEff=true
DoAccEffRw=false
DoHFPtSpec=false
DoHFPtSpecRaa=false
DoDmesonYield=false
DoDataDrivenCrossSection=false

#wheter you are projecting a tree or a sparse
ProjectTree=false

#PARAMETERS TO BE SET (use "" for parameters not needed)
################################################################################################
Particle="LctopKpi" # # whether it is Dplus, D0, Ds, LctopK0s or LctopKpi analysis
Cent="kpp13TeVFD" # used also to asses prompt or non-prompt and system

cfgFileData="configfiles/config_LctopKpi_data_tree_basic.yml"
cfgFileMC="configfiles/config_LctopKpi_MC_tree_basic.yml"
cfgFileFit="configfiles/fit/config_LctopKpi_Fit_pp13TeV_basic.yml"

#Config input MC files - specific for LctopKpi resonant channel
#WARNING: used only when Particle="LctopKpi"
cfgFileMCNonRes="configfiles/config_LctopKpi_NonRes_MC_tree_basic.yml"
cfgFileMCKStar="configfiles/config_LctopKpi_KStar_MC_tree_basic.yml"
cfgFileMCDelta="configfiles/config_LctopKpi_Delta_MC_tree_basic.yml"
cfgFileMCLambda1520="configfiles/config_LctopKpi_Lambda1520_MC_tree_basic.yml"
LctopKpi_reso_channel=('' '_NonRes' '_KStar' '_Delta' '_Lambda1520')

accFileName="accfiles/Acceptance_Toy_LcpKpi_yfidPtDep_etaDau09_ptDau100_promptD0FONLL13ptshape_FONLLy.root"
predFileName="models/fonll/feeddown/DmesonLcPredictions_13TeV_y05_FFptDepLHCb_BRpythia8_PDG2020_PromptLcMod.root"
pprefFileName="" #"ppreference/Ds_ppreference_pp5TeV_noyshift_pt_2_3_4_6_8_12_16_24_36_50.root"

PtWeightsDFileName=""
PtWeightsDHistoName=""
PtWeightsBFileName=""
PtWeightsBHistoName=""

MultWeightsFileName="systematics/genmultdistr/multweights/MultWeights_pp13TeV_MB_030_fnonprompt.root"
MultWeightsHistoName="hNtrklWeightsCandInMass"

DataDrivenFractionFileName=""

#assuming cutsets config files starting with "cutset" and are .yml

CutSetsDir="configfiles/cutsets/LctopKpi/pp/central/"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}



OutDirRawyields="rawyields/central"
OutDirEfficiency="efficiencies/central"
OutDirCrossSec=""
OutDirRaa=""
################################################################################################

if [ ${Particle} != "Dplus" ] && [ ${Particle} != "Ds" ]  && [ ${Particle} != "Dstar" ] && [ ${Particle} != "LctopK0s" ] && [ ${Particle} != "LctopKpi" ]; then
  echo $(tput setaf 1) ERROR: only Ds and Dplus mesons are supported! $(tput sgr0)
  exit 2
fi

if [ ${Cent} != "k010" ] && [ ${Cent} != "k3050" ] && [ ${Cent} != "k6080" ] && 
[ ${Cent} != "kpp5TeVFD" ] && [ ${Cent} != "kpp5TeVPrompt" ] &&
[ ${Cent} != "kpp13TeVPrompt" ] && [ ${Cent} != "kpp13TeVFD" ]; then
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

if [ ${Particle} == "LctopKpi" ] && [ ! -f "${cfgFileMCNonRes}" ]; then
  echo $(tput setaf 1) ERROR: MC config file "${cfgFileMCNonRes}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ${Particle} == "LctopKpi" ] && [ ! -f "${cfgFileMCKStar}" ]; then
  echo $(tput setaf 1) ERROR: MC config file "${cfgFileMCKStar}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ${Particle} == "LctopKpi" ] && [ ! -f "${cfgFileMCDelta}" ]; then
  echo $(tput setaf 1) ERROR: MC config file "${cfgFileMCDelta}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ${Particle} == "LctopKpi" ] && [ ! -f "${cfgFileMCLambda1520}" ]; then
  echo $(tput setaf 1) ERROR: MC config file "${cfgFileMCLambda1520}" does not exist! $(tput sgr0)
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

if [ ! -f "${DataDrivenFractionFileName}" ] && $DoDataDrivenCrossSection; then
  echo $(tput setaf 1) ERROR: data-driven fraction file "${DataDrivenFractionFileName}" does not exist! $(tput sgr0)
  exit 2
fi

if [ ! -f "${PtWeightsDFileName}" ] && [ "${PtWeightsDFileName}" != "" ]; then
  echo $(tput setaf 3) WARNING: pT-weights file "${PtWeightsDFileName}" does not exist! $(tput sgr0)
fi

if [ ! -f "${PtWeightsBFileName}" ] && [ "${PtWeightsBFileName}" != "" ]; then
  echo $(tput setaf 3) WARNING: pTB-weights file "${PtWeightsBFileName}" does not exist! $(tput sgr0)
fi

if [ ! -f "${MultWeightsFileName}" ] && [ "${MultWeightsFileName}" != "" ]; then
  echo $(tput setaf 3) WARNING: mult weights file "${MultWeightsFileName}" does not exist! $(tput sgr0)
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
elif [ ${Cent} == "kpp13TeVPrompt" ]; then
  isPrompt=true
  System=pp
elif [ ${Cent} == "kpp13TeVFD" ]; then
  isPrompt=false
  System=pp
else
  echo $(tput setaf 1) ERROR: system ${Cent} is not supported! $(tput sgr0)
  exit 2
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
    python3 ${ProjectScript} ${cfgFileData} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml ${OutDirRawyields}/Distr_${Particle}_data${CutSets[$iCutSet]}.root
  done
fi

if $DoMCProjection; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Projecting MC distributions $(tput sgr0)
    if [ ${Particle} == "LctopKpi" ]; then
      python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root
      python3 ${ProjectScript} ${cfgFileMCNonRes} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_NonRes_MC${CutSets[$iCutSet]}.root --LctopKpireso 1
      python3 ${ProjectScript} ${cfgFileMCKStar} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_KStar_MC${CutSets[$iCutSet]}.root --LctopKpireso 2
      python3 ${ProjectScript} ${cfgFileMCDelta} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_Delta_MC${CutSets[$iCutSet]}.root --LctopKpireso 3
      python3 ${ProjectScript} ${cfgFileMCLambda1520} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_Lambda1520_MC${CutSets[$iCutSet]}.root --LctopKpireso 4
    elif [ ${Particle} != "LctopKpi" ]; then    
      if [ "${PtWeightsDFileName}" == "" -o "${PtWeightsDHistoName}" == "" ] && [ "${PtWeightsBFileName}" == "" -o "${PtWeightsBHistoName}" == "" ] && [ "${MultWeightsFileName}" == "" -o "${MultWeightsHistoName}" == "" ]; then
        python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root
      elif [ "${PtWeightsDFileName}" != "" ] && [ "${PtWeightsDHistoName}" != "" ] && [ "${PtWeightsBFileName}" == "" -o "${PtWeightsBHistoName}" == "" ] && [ "${MultWeightsFileName}" == "" ] && [ "${MultWeightsHistoName}" == "" ]; then
          echo $(tput setaf 6) Using ${PtWeightsDHistoName} pt weights from ${PtWeightsDFileName} $(tput sgr0)
          python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root --ptweights ${PtWeightsDFileName} ${PtWeightsDHistoName}
      elif [ "${PtWeightsDFileName}" != "" ] && [ "${PtWeightsDHistoName}" != "" ] && [ "${PtWeightsBFileName}" != "" ] && [ "${PtWeightsBHistoName}" != "" ] && [ "${MultWeightsFileName}" == "" ] && [ "${MultWeightsHistoName}" == "" ]; then
          echo $(tput setaf 6) Using ${PtWeightsDHistoName} pt weights from ${PtWeightsDFileName} and ${PtWeightsBHistoName} ptB weights from ${PtWeightsBFileName} $(tput sgr0)
          python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root --ptweights ${PtWeightsDFileName} ${PtWeightsDHistoName} --ptweightsB ${PtWeightsBFileName} ${PtWeightsBHistoName}
      elif [ "${MultWeightsFileName}" != "" ] && [ "${MultWeightsHistoName}" != "" ] && [ "${PtWeightsDFileName}" == "" ] && [ "${PtWeightsDHistoName}" == "" ] && [ "${PtWeightsBFileName}" == "" ] && [ "${PtWeightsBHistoName}" == "" ]; then
          echo $(tput setaf 6) Using ${MultWeightsHistoName} mult weights from ${MultWeightsFileName} $(tput sgr0)
          python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root --multweights ${MultWeightsFileName} ${MultWeightsHistoName}
      elif [ "${MultWeightsFileName}" != "" ] && [ "${MultWeightsHistoName}" != "" ] && [ "${PtWeightsDFileName}" != "" ] && [ "${PtWeightsDHistoName}" != "" ] && [ "${PtWeightsBFileName}" == "" ] && [ "${PtWeightsBHistoName}" == "" ]; then
          echo $(tput setaf 6) Using ${MultWeightsHistoName} mult weights from ${MultWeightsFileName} and ${PtWeightsDHistoName} pt weights from ${PtWeightsDFileName} $(tput sgr0)
          python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root --multweights ${MultWeightsFileName} ${MultWeightsHistoName} --ptweights ${PtWeightsDFileName} ${PtWeightsDHistoName}
      elif [ "${MultWeightsFileName}" != "" ] && [ "${MultWeightsHistoName}" != "" ] && [ "${PtWeightsDFileName}" != "" ] && [ "${PtWeightsDHistoName}" != "" ] && [ "${PtWeightsBFileName}" != "" ] && [ "${PtWeightsBHistoName}" != "" ]; then
          echo $(tput setaf 6) Using ${MultWeightsHistoName} mult weights from ${MultWeightsFileName} and ${PtWeightsDHistoName} pt weights from ${PtWeightsDFileName} and ${PtWeightsBHistoName} ptB weights from ${PtWeightsBFileName} $(tput sgr0)
          python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root --multweights ${MultWeightsFileName} ${MultWeightsHistoName} --ptweights ${PtWeightsDFileName} ${PtWeightsDHistoName} --ptweightsB ${PtWeightsBFileName} ${PtWeightsBHistoName}
      fi
    fi
  done
fi

#compute raw yields
if $DoMCRawYields; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Extract raw yields from ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root $(tput sgr0)
    echo '.x GetRawYieldsDplusDs.C+('${Cent}',true, "'${OutDirEfficiency}'/Distr_'${Particle}'_MC'${CutSets[$iCutSet]}'.root", "", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYields'${Meson}'_MC'${CutSets[$iCutSet]}'.root")' | root -l -b
    echo '.q'
  done
fi

if $DoDataRawYields; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Extract raw yields from ${OutDirRawyields}/Distr_${Particle}_data${CutSets[$iCutSet]}.root $(tput sgr0)
    if [ ${Particle} != "D0" ]; then
      echo '.x GetRawYieldsDplusDs.C+('${Cent}',false, "'${OutDirRawyields}'/Distr_'${Particle}'_data'${CutSets[$iCutSet]}'.root", "", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYields'${Particle}${CutSets[$iCutSet]}'.root")' | root -l -b
      echo '.q'
    else
      echo '.x GetRawYieldsDplusDs.C+('${Cent}',false, "'${OutDirRawyields}'/Distr_'${Particle}'_data'${CutSets[$iCutSet]}'.root", "", "'${OutDirEfficiency}'/Distr_'${Particle}'_MC'${CutSets[$iCutSet]}'.root", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYields'${Particle}${CutSets[$iCutSet]}'.root")' | root -l -b
      echo '.q'
      fi
  done
fi

#compute efficiency
if $DoEfficiency; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    if [ ${Particle} == "LctopKpi" ]; then
      for Channel in "${LctopKpi_reso_channel[@]}";
      do
        echo $(tput setaf 4) Compute efficiency from ${OutDirEfficiency}/Distr_${Particle}${Channel}_MC${CutSets[$iCutSet]}.root $(tput sgr0)
        python3 ComputeEfficiencyDplusDs.py ${cfgFileFit} ${Cent} ${OutDirEfficiency}/Distr_${Particle}${Channel}_MC${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Efficiency_${Particle}${Channel}${CutSets[$iCutSet]}.root --batch
      done
    elif [ ${Particle} != "LctopKpi" ]; then
      echo $(tput setaf 4) Compute efficiency from ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root $(tput sgr0)
      python3 ComputeEfficiencyDplusDs.py ${cfgFileFit} ${Cent} ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Efficiency_${Particle}${CutSets[$iCutSet]}.root --batch
    fi
  done
fi

#compute efficiency times acceptance
if $DoAccEff; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    if [ ${Particle} == "LctopKpi" ]; then
    for Channel in "${LctopKpi_reso_channel[@]}";
    do
      echo $(tput setaf 4) Compute efficiency times acceptance ${Particle} ${Channel} $(tput sgr0)
      python3 CombineAccTimesEff.py ${OutDirEfficiency}/Efficiency_${Particle}${Channel}${CutSets[$iCutSet]}.root ${accFileName} ${OutDirEfficiency}/Eff_times_Acc_${Particle}${Channel}${CutSets[$iCutSet]}.root --batch
    done
    elif [ ${Particle} != "LctopKpi" ]; then
      echo $(tput setaf 4) Compute efficiency times acceptance $(tput sgr0)
      python3 CombineAccTimesEff.py ${OutDirEfficiency}/Efficiency_${Particle}${CutSets[$iCutSet]}.root ${accFileName} ${OutDirEfficiency}/Eff_times_Acc_${Particle}${CutSets[$iCutSet]}.root --batch
    fi
  done
fi

if $DoAccEffRw; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    if [ ${Particle} == "LctopKpi" ]; then
      echo $(tput setaf 4) Compute re-weighted efficiency times acceptance $(tput sgr0)
      python3 ComputeEffAccWeightedAvg.py ${OutDirEfficiency}/Eff_times_Acc_${Particle}_NonRes${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Eff_times_Acc_${Particle}_KStar${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Eff_times_Acc_${Particle}_Delta${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Eff_times_Acc_${Particle}_Lambda1520${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Eff_times_Acc_${Particle}${CutSets[$iCutSet]}_rw.root
    elif [ ${Particle} != "LctopKpi" ]; then
      echo $(tput setaf 4) This is not LctopKpi --> not computing average efficiency re-weighted with BR $(tput sgr0)
    fi
  done
fi

#compute cross section
if $DoDataDrivenCrossSection; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute cross section $(tput sgr0)
    if $isPrompt; then
      python3 ComputeDataDrivenCrossSection.py ${OutDirRawyields}/RawYields${Particle}${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Eff_times_Acc_${Particle}${CutSets[$iCutSet]}.root ${DataDrivenFractionFileName} ${OutDirCrossSec}/CrossSection${Particle}${CutSets[$iCutSet]}.root --prompt --${Particle} --system ${System} --energy 5.02 --batch
    else
      python3 ComputeDataDrivenCrossSection.py ${OutDirRawyields}/RawYields${Particle}${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Eff_times_Acc_${Particle}${CutSets[$iCutSet]}.root ${DataDrivenFractionFileName} ${OutDirCrossSec}/CrossSection${Particle}${CutSets[$iCutSet]}.root --FD --${Particle} --system ${System} --energy 5.02 --batch
    fi
  done
fi

#compute HFPtSpectrum
if $DoHFPtSpec; then
  Channel=""
  if [ ${Particle} == "Ds" ]; then 
    Channel="kDsKKpi"
  elif [ ${Particle} == "Dplus" ]; then 
    Channel="kDplusKpipi"
  elif [ ${Particle} == "D0" ]; then
    Channel="kD0Kpi"
  fi
  
  cc=""
  year=""
  sigma=1.
  if [ ${Cent} == "k010" -o "${Cent}" == "k3050" -o "${Cent}" == "k6080" ]; then
    cc=$Cent
    year="k2018"
    sigma=1.
  elif [ ${Cent} == "kpp5TeVPrompt" -o "${Cent}" == "kpp5TeVFD" ]; then
    cc="kpp5"
    year="k2017"
    sigma=50870000000.
  fi

  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute HFPtspectrum $(tput sgr0)
    echo '.x HFPtSpectrum.C+('${Channel}',"'${predFileName}'","'${OutDirEfficiency}'/Eff_times_Acc_'${Particle}${CutSets[$iCutSet]}'.root","'${OutDirRawyields}'/RawYields'${Particle}${CutSets[$iCutSet]}'.root","hRawYields","hAccEffPrompt","hAccEffFD","hEvForNorm","'${OutDirCrossSec}'/HFPtSpectrum'${Particle}${CutSets[$iCutSet]}'.root",kNb,'${sigma}',true,'${cc}','${year}')' | root -l -b
    echo '.q'
  done
fi

#compute HFPtSpectrumRaa
if $DoHFPtSpecRaa; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute HFPtspectrumRaa $(tput sgr0)
    echo '.x HFPtSpectrumRaa.C+("'${pprefFileName}'","'${OutDirCrossSec}'/HFPtSpectrum'${Particle}${CutSets[$iCutSet]}'.root","'${OutDirRaa}'/HFPtSpectrumRaa'${Particle}${CutSets[$iCutSet]}'.root",4,1,kNb,'${Cent}',k2018,k5dot023,1./3,3,6,false,1)' | root -l -b
    echo '.q'
  done
fi

#compute corrected yield
if $DoDmesonYield; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute corrected yield $(tput sgr0)
    echo '.x ComputeDmesonYield.C+(k'${Particle}','${Cent}',2,1,"'${pprefFileName}'","'${OutDirCrossSec}'/HFPtSpectrum'${Particle}${CutSets[$iCutSet]}'.root","","'${OutDirRaa}'/HFPtSpectrumRaa'${Particle}${CutSets[$iCutSet]}'.root","","'${OutDirCrossSec}'","'${CutSets[$iCutSet]}'",1,1./3,3,false,1)' | root -l -b
    echo '.q'
  done
fi
