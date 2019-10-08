#!/bin/bash
CutSetsDir="configfiles/syst_cuts_Ds3050_MLlowpt"
declare -a CutSets=()
for filename in ${CutSetsDir}/*.yml; do
    tmp_name="$(basename -- ${filename} .yml)"
    tmp_name=${tmp_name:6}
    CutSets+=("${tmp_name}")
done
arraylength=${#CutSets[@]}

for (( iCutSet=0; iCutSet<${arraylength}; iCutSet++ ));
  do
    echo -n "\"${CutSets[$iCutSet]}\","
  done
echo " "