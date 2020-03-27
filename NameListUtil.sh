#!/bin/bash
CutSetsDir="configfiles/cutsets/Ds/pp/syst_cuts_FDen_1D"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo "\"${CutSets[$iCutSet]}\","
  done
echo " "