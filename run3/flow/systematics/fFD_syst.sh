#!/bin/bash
# Bash script to compute resolution systematics in flow analysis

########################################
# Define paths and configuration
########################################

# Default cut variation path (uncorrected)
export cutvar_default="path/to/default/V2VsFrac"

# Array of systematic variation paths (e.g., for "random" and "even" variations)
cutvar_systvariation_paths=(
  "path/to/Syst/random"
  "path/to/Syst/even"
  "path/to/Syst/odd"
  "path/to/Syst/1in4"
  "path/to/Syst/1in3"
  "path/to/Syst/minus3low"
  "path/to/Syst/minus3high"
  "path/to/Syst/minus3"
)

# Other important paths
export corrPath="path/to/correlated/folder"
export uncorrPath="path/to/uncorr/fodler/"
export output_dir="output/dir"  # Output directory for results
export suffix="suffix"
export config="../config_flow.yml"
export n_parallel=2  # Number of parallel jobs

########################################
# Parallel function to run analysis for each configuration
########################################

parallel_func() {
    local corrPath="$1"

    # Check if the efficiency folder already exists in corrPath.
    # If not, copy it from the uncorrected path.
    if [ -d "$corrPath/eff" ]; then
        echo "Efficiency folder already exists in $corrPath"
    else
        echo "Copying the efficiency files from $uncorrPath/eff to $corrPath/eff"
        cp -r "$uncorrPath/eff" "$corrPath/eff"
    fi
    
    # Compute a suffix from corrPath:
    # Example: If corrPath is "CutVarFrac_XXX.root", suffix becomes "CutVarFrac".
    local base
    base=$(basename "$corrPath")
    local suffix="${base%%_*}"   # Extract text before the first underscore

    # Determine the parent directory of corrPath's parent directory.
    # This could be used as the "cut variation" directory.
    local path_to_cutvar
    path_to_cutvar="$(dirname "$(dirname "$corrPath")")"
    echo "Path to cutvar: $path_to_cutvar"
    echo "Suffix: $suffix"

    ########################################
    # Build and execute Python commands
    ########################################

    # (Optional) Command to compute data-driven fractions (currently commented out)
    local cmd=(python3 "./ComputeDataDriFrac_flow.py" -i "$corrPath" -o "$corrPath" -s "$suffix" -b -comb -cc "$corrPath" -oc "$corrPath")
    echo "Executing: ${cmd[*]}"
    "${cmd[@]}"

    # Command to compute prompt FD v2
    local cmd=(python3 "./ComputeV2vsFDFrac.py" "$config" -i "$uncorrPath" -o "$corrPath" -s "$suffix" -comb -ic "$corrPath" -oc "$corrPath")
    echo "Executing: ${cmd[*]}"
    "${cmd[@]}"
}
export -f parallel_func

########################################
# Run the analysis for each configuration in parallel
########################################

# Change directory to the BDT folder (adjust as necessary)
cd ../BDT/ || { echo "Failed to change directory to ../BDT/"; exit 1; }

echo "Running parallel jobs with n_parallel: $n_parallel"

# Use GNU Parallel to run the function for each systematic variation path in the array
parallel -j "$n_parallel" parallel_func ::: "${cutvar_systvariation_paths[@]}"

########################################
# Further commands for computing systematic uncertainties can be added below
########################################
# Change directory to the BDT folder (adjust as necessary)
cd ../systematics/ || { echo "Failed to change directory to ../systematics/"; exit 1; }

echo "Computing systeamtic"
# Enable recursive globbing so that ** searches through all subdirectories
shopt -s globstar

# Declare an array to store all file paths
all_v2_files=()

# Add reference to cutvar_systvariation_paths
cutvar_systvariation_paths=(
  "$cutvar_default"
  "${cutvar_systvariation_paths[@]}"
)
# Iterate over each systematic variation path
for syst_path in "${cutvar_systvariation_paths[@]}"; do
    echo "Searching in: $syst_path"
    
    # Find all files matching the pattern V2VsFrac*.root in any subfolder
    mapfile -t files < <(find "$syst_path" -type f -name "V2VsFrac*.root")

    # Check if any files were found
    if [ ${#files[@]} -eq 0 ]; then
        echo "No V2VsFrac*.root files found in $syst_path"
    else
        # Append found files to the global list
        echo "V2VsFrac.root file found in $syst_path"
        all_v2_files+=("${files[@]}")
    fi
done

# Disable recursive globbing if no longer needed
shopt -u globstar

echo "Running python3 compute_syst_fFD.py -s $cutvar_systvariation_paths -o $output_dir"
python3 compute_syst_fFD.py ${all_v2_files[@]} -s $suffix -o $output_dir