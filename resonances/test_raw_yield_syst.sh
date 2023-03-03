#!/bin/bash

# Script to run the full analysis for a specific configuration
export nJobsPerStep=8
export nJobsPerConf=1

#______________________________________________________________
# Configuration for multitrial generation and analysis
export RefFile="/home/stefano/Desktop/cernbox/Ds_resonances/raw_yield/mass_Ds2starplus_pt2.0-24.0_MB.root"
export doGenTrials=false # if true, generate multitrials, otherwise use the ones in InputTrialsFile
export TrialsOutName="multitrialcheck_Ds2starplus_MB.root" # multitrialcheck_{name_reso}_{trigger}.root
export ConfFile="config_fit_Ds2starplus_rawyieldsyst_checkmultitrial.yml"
export InputTrialsFile="/home/stefano/Desktop/Analyses/DmesonAnalysis/resonances/multitrialcheck_Ds2starplus_MB.root"  # in case generation is not needed
export nTrials=100 # number of trials to perform, max 100


if (( nJobsPerConf * nJobsPerStep > $(nproc) )); then
  echo $(tput setaf 3) "WARNING: nJobsPerConf*nJobsPerStep = $((nJobsPerConf * nJobsPerStep)) is greater than the number of available cores" $(tput sgr0)
fi

if [ "$doGenTrials" = true ]; then
    echo "Generating multitrials for ${RefFile} --> ${TrialsOutName}"
    python3 generate_multitrials.py ${RefFile} ${TrialsOutName}
    export InputTrialsFile=${TrialsOutName}
fi

echo "Performing analysis for ry ${ConfFile} --> ${InputFile}"
for (( iTrial=0; iTrial<$nTrials; iTrial++ )); # 100 is max number of trials
do
    Trials+=($iTrial)
done
echo "Trials: ${Trials[@]}"
export Trials

function multitrial () {
    echo "Prforming multitrial on ${InputFile} --> htrial_$1"
    python3 rawyield_syst_Dsreso.py ${ConfFile} ${InputTrialsFile} --histoname htrial_$1
    echo "Iteration $1 completed"
}
export -f multitrial
time parallel -j ${nJobsPerStep} --tag -k 'multitrial {1} 2>&1' ::: ${Trials[@]}