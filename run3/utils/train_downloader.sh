#!/bin/bash
TRAIN_RUN=XXXXXX
INPUT="files_to_download.txt"
IS_SLIM=true
MAX_FILES_TO_DOWNLOAD=2
START_FILE=0
SUFFIX="LHCXXYY"
SELECTIONS="a < X" # Can be null if no selection is needed
TRAIN_FRACTION=1
ISMC=true
MC_SELECTIONS=( "a:X > 1 and Y > 1" # Array of selections in the form "Name:selection", or "" if no selection is needed
                "b:X < 1 and Y > 1"
                "c:X > 1 and Y < 1"
                "d:X < 1 and Y < 1")

TREE_NAME=("tree1" "tree2" "tree3") # Can be an array of strings for multiple trees, or a single string for a single tree. "*" for all trees.
OUTPUT_DIRECTORY="Folder/For/Outputs" # Files will be downloaded in this/directory/MC(data)/TrainXXXXXX/
OPTIONAL_ARGS="" # --aod, --analysis, or --parquet


#_______________________________________________________________________________________________________________________

# Check if TREE_NAME is an array or a single string
if [ "$(declare -p TREE_NAME 2>/dev/null | grep 'declare -a')" ]; then
    # Convert array to JSON format
    tree_name_json=$(printf '%s\n' "${TREE_NAME[@]}" | jq -R . | jq -s .)
else
    # Single string, just wrap it in quotes
    tree_name_json=$(jq -n --arg tn "$TREE_NAME" '$tn')
fi

# Convert MC_SELECTIONS to JSON format using jq
if [ -z "$MC_SELECTIONS" ]; then
    selections_json="null"
else
    selections_json=$(printf '%s\n' "${MC_SELECTIONS[@]}" | awk -F: '{print ""$1": "$2""}' | jq -R -s 'split("\n") | map(select(length > 0)) | map(split(": ") | { (.[0]): .[1] }) | add')
fi

# Pass information to a json file using jq
# For strings, use --arg, for numbers, use --argjson
jq -n --argjson train_run "$TRAIN_RUN" \
      --arg input "$INPUT" \
      --argjson is_slim "$IS_SLIM" \
      --argjson max_files_to_download "$MAX_FILES_TO_DOWNLOAD" \
      --argjson start_file "$START_FILE" \
      --arg suffix "$SUFFIX" \
      --arg selections "$SELECTIONS" \
      --argjson train_fraction "$TRAIN_FRACTION" \
      --argjson isMC "$ISMC" \
      --argjson mc_selections "$selections_json" \
      --argjson tree_name "$tree_name_json" \
      --arg output_directory "$OUTPUT_DIRECTORY" \
      '{
          train_run: $train_run,
          input: $input,
          is_slim: $is_slim,
          max_files_to_download: $max_files_to_download,
          start_file: $start_file,
          suffix: $suffix,
          selections: $selections,
          train_fraction: $train_fraction,
          isMC: $isMC,
          mc_selections: $mc_selections,
          tree_name: $tree_name,
          output_directory: $output_directory
      }' > config_download.json

python3 train_output_downloader.py config_download.json $OPTIONAL_ARGS

rm -f config_download.json