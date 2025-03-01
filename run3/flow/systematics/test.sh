#!/bin/bash
export config_modifies="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/systematics/modifications_config_fit.yml"
export config_default="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/config/config_flow.yml" # bdt cut of corr and uncorr, skip_cut, sysematic
export output_dir="/home/wuct/ALICE/local/Results/test/flow/systematics_1"

export sigle_pt=True
export n_parallel=1

########################################
# prepare the config files, pre-processed files, and run reference cut variation
########################################

mkdir -p $output_dir

# # # # produce the yaml files for the systematic: config_flow, config_pre, config_defalut, config_default_corelated
rm -rf $output_dir/config_sys
python3 ./make_yaml_for_syst.py $config_default -m $config_modifies -o $output_dir -mb

#* DONE
#____________________________________________________________________________________________________________________
# # produce the pre-processed files for the systematic
python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/tool/pre_process.py $output_dir/config_sys/config_pre.yml --out_dir $output_dir --pre_sys

# * DONE
#____________________________________________________________________________________________________________________
# # Uncorelated: provide the projection of MC
# # if you don't want the final results, delete: --do_vn --do_frac_cut_var --do_data_driven_frac --do_v2_vs_frac --do_merge_images
python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/reference/config_uncorrelated.yml \
        --use_preprocessed --do_calc_weights --do_make_yaml --do_projections --do_efficiency --do_vn \
        --do_frac_cut_var --do_data_driven_frac --do_v2_vs_frac --do_merge_images --do_sys_trail
# # python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/reference/config_uncorrelated.yml \
# #         --do_data_driven_frac --do_v2_vs_frac --do_merge_images --do_sys_trail

#* DONE
#____________________________________________________________________________________________________________________
# # Correlated: provide the fraction
# # if you don't want the final results, delete: --do_data_driven_frac --do_v2_vs_frac --do_merge_images
python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/reference/config_correlated.yml \
        --do_calc_weights --do_make_yaml --do_projections --do_efficiency --do_vn \
        --do_frac_cut_var --do_data_driven_frac --do_v2_vs_frac --do_merge_images --do_sys_trail

#* DONE
#____________________________________________________________________________________________________________________
# # Combined: provide the central value
# # if you rerun the combined instead of copy, add whole flag unless --do_sys_trail
cp -r ${output_dir}/pre_sys/cutvar_uncorr/ry $output_dir/pre_sys/cutvar_combined/ry
cp -r ${output_dir}/pre_sys/cutvar_uncorr/eff $output_dir/pre_sys/cutvar_combined/eff
python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/reference/config_combined.yml \
        --do_data_driven_frac --do_v2_vs_frac

########################################
# multi-trails
########################################

# --do_proj_mc actually means do not project MC
parallel_func() {
    config_file=$1
    suffix=$(basename $config_file .yml)
    echo "suffix: $suffix"
    python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py $config_file --use_preprocessed --do_projections --do_proj_mc --do_vn --do_data_driven_frac --do_v2_vs_frac --do_merge_images --do_sys_trail
}
export -f parallel_func


config_files=$(find $output_dir/config_sys/all_pt -name "config_*.yml")
echo "config_files: $config_files"
start_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "Start time: $start_time"
parallel -j $n_parallel parallel_func ::: $config_files
end_time=$(date +"%Y-%m-%d %H:%M:%S")
echo "End time: $end_time"

python3 compute_syst_multitrial_bdt.py $output_dir/trails/all_pt $output_dir/pre_sys/cutvar_combined -o $output_dir -p
