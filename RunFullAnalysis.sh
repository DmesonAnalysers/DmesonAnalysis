#!/bin/bash
#steps to be performed
DoDataProjection=false
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
ProjectTree=false

#PARAMETERS TO BE SET (use "" for parameters not needed)
################################################################################################
Particle="Dplus" # whether it is Dplus,  Ds or Lc analysis
Cent="kpp13TeVFD" # used also to asses prompt or non-prompt and system

cfgFileData="../analyses/np_Dplus_vsmult_pp13TeV/configfiles/inputs/config_Dplus_pp13TeV_data_MB_multWeights.yml"
cfgFileMC="../analyses/np_Dplus_vsmult_pp13TeV/configfiles/inputs/config_Dplus_pp13TeV_MC_MB_multWeights.yml"
cfgFileFit="../analyses/np_Dplus_vsmult_pp13TeV/configfiles//fit/config_fit_Dplus_13TeV.yml"

accFileName="accfiles/Acceptance_Toy_DplusKpipi_yfidPtDep_etaDau09_ptDau100_promptDplusFONLL13ptshape_FONLLy.root"
predFileName="models/fonll/feeddown/DmesonLcPredictions_13TeV_y05_FFee_BRpythia8_SepContr_PDG2020.root"
pprefFileName="" #"ppreference/Ds_ppreference_pp5TeV_noyshift_pt_2_3_4_6_8_12_16_24_36_50.root"

PtWeightsDFileName=""
PtWeightsDHistoName=""
PtWeightsBFileName=""
PtWeightsBHistoName=""

MultWeightsFileName="systematics/genmultdistr/multweights/MultWeights_pp13TeV_MB_030_fnonprompt.root"
MultWeightsHistoName="hNtrklWeightsCandInMass"

DataDrivenFractionFileName=""

#assuming cutsets config files starting with "cutset" and are .yml
CutSetsDir="../analyses/np_Dplus_vsmult_pp13TeV/configfiles/cutsets/datadriven/030"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

OutDirRawyields="../analyses/np_Dplus_vsmult_pp13TeV/outputs/rawyields/030/"
OutDirEfficiency="../analyses/np_Dplus_vsmult_pp13TeV/outputs/efficiencies/030/"
OutDirCrossSec=""
OutDirRaa=""
################################################################################################

if [ ${Particle} != "Dplus" ] && [ ${Particle} != "Ds" ] && [ ${Particle} != "Lc" ]; then
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
    if [ "${PtWeightsDFileName}" == "" -o "${PtWeightsDHistoName}" == "" ] && [ "${PtWeightsBFileName}" == "" -o "${PtWeightsBHistoName}" == "" ] && [ "${MultWeightsFileName}" == "" -o "${MultWeightsHistoName}" == "" ]; then
      python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root
    elif [ "${PtWeightsDFileName}" != "" ] && [ "${PtWeightsDHistoName}" != "" ] && [ "${PtWeightsBFileName}" == "" -o "${PtWeightsBHistoName}" == "" ]; then
        echo $(tput setaf 6) Using ${PtWeightsDHistoName} pt weights from ${PtWeightsDFileName} $(tput sgr0)
        python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root --ptweights ${PtWeightsDFileName} ${PtWeightsDHistoName}
    elif [ "${PtWeightsDFileName}" != "" ] && [ "${PtWeightsDHistoName}" != "" ] && [ "${PtWeightsBFileName}" != "" ] && [ "${PtWeightsBHistoName}" != "" ]; then
        echo $(tput setaf 6) Using ${PtWeightsDHistoName} pt weights from ${PtWeightsDFileName} and ${PtWeightsBHistoName} ptB weights from ${PtWeightsBFileName} $(tput sgr0)
        python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root --ptweights ${PtWeightsDFileName} ${PtWeightsDHistoName} --ptweightsB ${PtWeightsBFileName} ${PtWeightsBHistoName}
    elif [ "${MultWeightsFileName}" != "" ] && [ "${MultWeightsHistoName}" != "" ]; then
        echo $(tput setaf 6) Using ${MultWeightsHistoName} mult weights from ${MultWeightsFileName} $(tput sgr0)
        python3 ${ProjectScript} ${cfgFileMC} ${CutSetsDir}/cutset${CutSets[$iCutSet]}.yml  ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root --multweights ${MultWeightsFileName} ${MultWeightsHistoName}
    fi
  done
fi

#compute raw yields
if $DoMCRawYields; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Extract raw yields from ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root $(tput sgr0)
    echo '.x GetRawYieldsDplusDs.C+('${Cent}',true, "'${OutDirEfficiency}'/Distr_'${Particle}'_MC'${CutSets[$iCutSet]}'.root", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYields'${Meson}'_MC'${CutSets[$iCutSet]}'.root")' | root -l -b
    echo '.q'
  done
fi

if $DoDataRawYields; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Extract raw yields from ${OutDirRawyields}/Distr_${Particle}_data${CutSets[$iCutSet]}.root $(tput sgr0)
    echo '.x GetRawYieldsDplusDs.C+('${Cent}',false, "'${OutDirRawyields}'/Distr_'${Particle}'_data'${CutSets[$iCutSet]}'.root", "'${cfgFileFit}'", "'${OutDirRawyields}'/RawYields'${Particle}${CutSets[$iCutSet]}'.root")' | root -l -b
    echo '.q'
  done
fi

#compute efficiency
if $DoEfficiency; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute efficiency from ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root $(tput sgr0)
    python3 ComputeEfficiencyDplusDs.py ${cfgFileFit} ${Cent} ${OutDirEfficiency}/Distr_${Particle}_MC${CutSets[$iCutSet]}.root ${OutDirEfficiency}/Efficiency_${Particle}${CutSets[$iCutSet]}.root --batch
  done
fi

#compute efficiency times acceptance
if $DoAccEff; then
  for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo $(tput setaf 4) Compute efficiency times acceptance $(tput sgr0)
    python3 CombineAccTimesEff.py ${OutDirEfficiency}/Efficiency_${Particle}${CutSets[$iCutSet]}.root ${accFileName} ${OutDirEfficiency}/Eff_times_Acc_${Particle}${CutSets[$iCutSet]}.root --batch
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
