#!/bin/bash

declare -a arrayOfFiles=("outputs_valid" "event_stat.root" "*.xml" "*.jdl" "myAnalysis*" "Improver.root" "EventStat_temp.root" "stdout" "stderr" "*.d" "*.so" "*.pcm")

for FILE in "${arrayOfFiles[@]}"
do
    rm $FILE 2>/dev/null
done
