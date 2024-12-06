#!/bin/bash

# Bash script to process input directories with a suffix dynamically picked from run_label.txt

# Define default parameters
TARGET_DIR="./inputs/"
OUTPUT_BASE_DIR="./output_reso"
INPUT_FILE="input_data_reso.txt"  # File with input directories (one per line)
RUN_LABEL="run_label.txt"         # File with corresponding suffixes for each input directory
FILES_TO_MERGE=("AnalysisResults")
CURRENT_DIR=$(pwd)
MAX_JOBS=4  # Limit the number of concurrent jobs

# Check if input files exist
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: Input file $INPUT_FILE not found."
    exit 1
fi
if [[ ! -f "$RUN_LABEL" ]]; then
    echo "Error: Run label file $RUN_LABEL not found."
    exit 1
fi

# Read input directories and suffixes into arrays
mapfile -t INPUT_DIRS < "$INPUT_FILE"
mapfile -t RUN_SUFFIXES < "$RUN_LABEL"

# Ensure both files have the same number of lines
if [[ ${#INPUT_DIRS[@]} -ne ${#RUN_SUFFIXES[@]} ]]; then
    echo "Error: Mismatch between number of input directories and run labels."
    exit 1
fi

# Create the output directory
mkdir -p "$OUTPUT_BASE_DIR"

# Function to process each input directory
process_input_dir() {
    local input_dir="$1"
    local suffix="$2"

    # Step 1: Download the file for the current input directory
    python3 download.py \
        --target_dir "$TARGET_DIR" \
        --input_dirs "$input_dir" \
        --suffix "$suffix" \
        --file_to_merge "${FILES_TO_MERGE[@]}"

    # Path to the downloaded AnalysisResults file
    local analysis_results_path="$CURRENT_DIR/inputs/AnalysisResults_${suffix}.root"

    # Step 2: Run the resolution
    python3 compute_reso.py "$analysis_results_path" -c "k0100" -o "$OUTPUT_BASE_DIR" -s "$suffix"

    # Step 3: (Eventually) Delete the original AnalysisResults file after use
    # rm "$analysis_results_path"
    echo "Processed and saved output for suffix $suffix."

    # Step 4: Process all reso together
    
}

# Loop through input directories and their corresponding suffixes
for i in "${!INPUT_DIRS[@]}"; do
    input_dir="${INPUT_DIRS[$i]}"
    suffix="${RUN_SUFFIXES[$i]}"

    process_input_dir "$input_dir" "$suffix" &  # Start the job in the background

    # Limit the number of concurrent jobs
    while (( $(jobs -r -p | wc -l) >= MAX_JOBS )); do
        sleep 1  # Wait for some jobs to finish
    done
done

# Wait for all background jobs to finish
wait

echo "All outputs saved in $OUTPUT_BASE_DIR."
