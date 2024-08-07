#!/bin/bash
# bash script to compute resolution systematics in flow analysis

export config_default="/home/stefano/flowD/output/pass3/Ds/sp/config/config_flow_Ds_2060.yml_LHC23zzh_pass3_Ds_skimmed_FULL_2060_train33944.yml" # default configuration file
export config_modifies="/home/stefano/alice/DmesonAnalysis/run3/flow/systematics/modifications_config_syst_reso.yml" # configuration file with the modifications to be applied
export output_dir="/home/stefano/flowD/output/pass3/Ds/systematics/resolution" # output directory
export anres_dir="/home/stefano/flowD/input/LHC23/pass3/Ds/AnalysisResults_LHC23zzh_pass3_Ds_skimmed_FULL_2060_train33944.root"
export vn_method="sp"
export res_file="/home/stefano/flowD/output/pass3/Resolution/resosp_pass3_3050.root"
export wagon_id=""
export centmin=30
export centmax=49

export n_parallel=2 # number of parallel jobs

#______________________________________________________________________________________________________________________
# Make the list of configurations to be used with make_yaml_for_syst.py
echo "python3 make_yaml_for_syst.py $config_default -m $config_modifies -o $output_dir"
 python3 make_yaml_for_syst.py $config_default -m $config_modifies -o $output_dir

#______________________________________________________________________________________________________________________
# List of centrality bins 1% width (30-50)
export cent_bins=()
for i in $(seq $centmin 1 $centmax); do
    cent_bins+=("k${i}$(($i+1))")
done
# delta centrality k3050
export centup=$(($centmax+1))
export deltacent="k${centmin}${centup}"

#______________________________________________________________________________________________________________________
# Parallel loop to run the analysis for each configuration
parallel_func() {
    config_file=$1
    cent=$2
    export suffix=$(basename $config_file .yml)
    # add centrality to suffix
    suffix="${suffix}_${cent}"
    echo "suffix: $suffix"
    # enter if wagon_id is not empty
    if [ -n "$wagon_id" ]; then
        echo -e "\e[35m python3 run_full_flow_analysis.py $config_file $anres_dir -o $output_dir -c $cent -v $vn_method -r 1 $skip_resolution $skip_efficiency -w $wagon_id \e[0m"
        python3 run_full_flow_analysis.py $config_file $anres_dir -o $output_dir -c $cent -v $vn_method -r 1  -w $wagon_id -s $suffix --skip_resolution --skip_efficiency --batch
    else
        echo -e "\e[35m python3 run_full_flow_analysis.py $config_file $anres_dir -o $output_dir -c $cent $skip_resolution -v $vn_method -r $res_file $skip_efficiency \e[0m"
        python3 run_full_flow_analysis.py $config_file $anres_dir -o $output_dir -c $cent -v $vn_method -r 1 -s $suffix --skip_resolution --skip_efficiency --batch
    fi
}
export -f parallel_func

#______________________________________________________________________________________________________________________
# Run the analysis for each configuration file in parallel (wait)
# Get the list of configuration files
export outdir_config_files="$output_dir/config_syst"
config_files=$(ls $outdir_config_files/config_*.yml)
cd ../ # go to the parent directory
echo "config_files: $config_files"
echo "cent_bins: ${cent_bins[@]}"
echo "n_parallel: $n_parallel"
echo "parallel -j $n_parallel parallel_func ::: $config_files ::: ${cent_bins[@]}"
parallel -j $n_parallel parallel_func ::: $config_files ::: ${cent_bins[@]}

#______________________________________________________________________________________________________________________
# Compute the systematic uncertainties
export syst_ry_path="$output_dir/$vn_method/ry/"
echo "python3 compute_syst_reso.py $config_modifies $config_default $res_file -c $deltacent -o $output_dir"
python3 compute_syst_reso.py $config_modifies $config_default $syst_ry_path -c $deltacent $res_file -o $output_dir 
