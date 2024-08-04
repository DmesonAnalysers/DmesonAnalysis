#!/bin/bash
# bash script to compute raw yield and vn vs mass for multiple configurations in parallel

export config_default="/home/stefano/alice/DmesonAnalysis/run3/flow/config_flow_Dplus_2060_run2bins.yml" # default configuration file
export config_modifies="/home/stefano/alice/DmesonAnalysis/run3/flow/systematics/modifications_config.yml" # configuration file with the modifications to be applied
export ry_default="/home/stefano/flowD/output/pass3/Dplus/sp/ry/raw_yields_LHC23zzh_pass3_Dplus_skimmed_FULL_2060_train237912.root"
export output_dir="/home/stefano/alice/DmesonAnalysis/run3/flow/systematics" # output directory
export anres_dir="/home/stefano/flowD/input/LHC23/pass3/Dplus/AnalysisResults_LHC23zzh_pass3_Dplus_skimmed_FULL_2060_train237912.root"
export cent="k3050"
export skip_resolution="--skip_resolution"
export vn_method="sp"
export res_file="/home/stefano/flowD/output/pass3/Resolution/resosp_pass3_3050.root"
export skip_efficiency="--skip_efficiency"
export wagon_id=""

export n_parallel=3 # number of parallel jobs

#______________________________________________________________________________________________________________________
# Make the list of configurations to be used with make_yaml_for_syst.py
echo "python3 make_yaml_for_syst.py $config_default -m $config_modifies -o $output_dir"
 python3 make_yaml_for_syst.py $config_default -m $config_modifies -o $output_dir

#______________________________________________________________________________________________________________________
# Parallel loop to run the analysis for each configuration
parallel_func() {
    config_file=$1
    export suffix=$(basename $config_file .yml)
    echo "suffix: $suffix"
    # enter if wagon_id is not empty
    if [ -n "$wagon_id" ]; then
        # echo in magenta 
        echo -e "\e[35m python3 run_full_flow_analysis.py $config_file $anres_dir -o $output_dir -c $cent $skip_resolution -v $vn_method -r $res_file $skip_efficiency -w $wagon_id \e[0m"
        python3 run_full_flow_analysis.py $config_file $anres_dir -o $output_dir -c $cent $skip_resolution -v $vn_method -r $res_file $skip_efficiency -w $wagon_id -s $suffix 
    else
        # echo in magenta
        echo -e "\e[35m python3 run_full_flow_analysis.py $config_file $anres_dir -o $output_dir -c $cent $skip_resolution -v $vn_method -r $res_file $skip_efficiency \e[0m"
        python3 run_full_flow_analysis.py $config_file $anres_dir -o $output_dir -c $cent $skip_resolution -v $vn_method -r $res_file $skip_efficiency -s $suffix
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
echo "n_parallel: $n_parallel"
echo "parallel -j $n_parallel parallel_func ::: $config_files"
parallel -j $n_parallel parallel_func ::: $config_files

#______________________________________________________________________________________________________________________
# Compute the systematic uncertainties
export syst_ry_path = "$output_dir/$vn_method/ry"
echo "python3 compute_syst_multitrial.py $syst_ry_path $ry_default -o $output_dir"
python3 compute_syst_multitrial.py $syst_ry_path $ry_default -o $output_dir
