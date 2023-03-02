#!/bin/bash

# Script to run the full analysis for a specific configuration
export nJobsPerStep=8
export nJobsPerConf=1

if (( nJobsPerConf * nJobsPerStep > $(nproc) )); then
  echo $(tput setaf 3) "WARNING: nJobsPerConf*nJobsPerStep = $((nJobsPerConf * nJobsPerStep)) is greater than the number of available cores" $(tput sgr0)
fi

export ConfFile="config_fit_Ds2starplus_rawyieldsyst_checkmultitrial.yml"
export InputFile="/home/stefano/Desktop/Analyses/DmesonAnalysis/resonances/multitrialcheck_Ds2starplus_HM.root"


echo "Performing analysis for ry ${ConfFile} --> ${InputFile}"
for (( iTrial=4; iTrial<30; iTrial++ )); # 100 is max number of trials
do
    Trials+=($iTrial)
done
echo "Trials: ${Trials[@]}"
export Trials

function multitrial () {
    echo "Prforming multitrial on ${InputFile} --> htrial_$1"
    python3 rawyield_syst_Dsreso.py ${ConfFile} ${InputFile} --histoname htrial_$1
    echo "Iteration $1 completed"
}
export -f multitrial
time parallel -j ${nJobsPerStep} --tag -k 'multitrial {1} 2>&1' ::: ${Trials[@]}