#!/bin/bash
export config_modifies="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/systematics/modifications_config_fit.yml"
export config_default="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/config/config_flow.yml" # bdt cut of corr and uncorr, skip_cut, sysematic
export res_file="/media/wuct/wulby/ALICE/AnRes/resolution/output_reso/resosp3050l_PASS4_full_PbPb_Reso.root"
export cent="k3050"
export output_dir="/home/wuct/ALICE/local/Results/test/flow/systematics_1"

export sigle_pt=True
export n_parallel=1

########################################
# prepare the config files, pre-processed files, and run reference cut variation
########################################

# mkdir -p $output_dir

# # # # produce the yaml files for the systematic: config_flow, config_pre, config_defalut, config_default_corelated
# rm -rf $output_dir/config_sys
# python3 ./make_yaml_for_syst.py $config_default -m $config_modifies -o $output_dir -mb

#* DONE
#____________________________________________________________________________________________________________________
# # produce the pre-processed files for the systematic
# python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/tool/pre_process.py $output_dir/config_sys/config_pre.yml --out_dir $output_dir --pre_sys

# * DONE
#____________________________________________________________________________________________________________________
# # Uncorelated: provide the projection of MC
# # if you don't want the final results, delete: --do_vn --do_frac_cut_var --do_data_driven_frac --do_v2_vs_frac --do_merge_images
# python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/reference/config_uncorrelated.yml \
#         --use_preprocessed --do_calc_weights --do_make_yaml --do_projections --do_efficiency --do_vn \
#         --do_frac_cut_var --do_data_driven_frac --do_v2_vs_frac --do_merge_images --do_sys_trail
# # python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/reference/config_uncorrelated.yml \
# #         --do_data_driven_frac --do_v2_vs_frac --do_merge_images --do_sys_trail

#* DONE
#____________________________________________________________________________________________________________________
# # Correlated: provide the fraction
# # if you don't want the final results, delete: --do_data_driven_frac --do_v2_vs_frac --do_merge_images
# python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/reference/config_correlated.yml \
#         --do_calc_weights --do_make_yaml --do_projections --do_efficiency --do_vn \
#         --do_frac_cut_var --do_data_driven_frac --do_v2_vs_frac --do_merge_images --do_sys_trail

#* DONE
#____________________________________________________________________________________________________________________
# # Combined: provide the central value
# # if you rerun the combined instead of copy, add whole flag unless --do_sys_trail
# cp -r ${output_dir}/pre_sys/cutvar_uncorr/ry $output_dir/pre_sys/cutvar_combined/ry
# cp -r ${output_dir}/pre_sys/cutvar_uncorr/eff $output_dir/pre_sys/cutvar_combined/eff
# python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/reference/config_combined.yml \
#         --do_data_driven_frac --do_v2_vs_frac

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

# if [ $sigle_pt = True ]; then
#     config_files_paths=$(find $output_dir/config_sys/pt_* -type d)
#     start_time=$(date +"%Y-%m-%d %H:%M:%S")
#     echo "Start time: $start_time"
#     for config_files_path in $config_files_paths; do
#         pt_bin=$(basename "$config_files_path" | grep -oE '[0-9]+_[0-9]+')
#         anres_dir=$(find "$output_dir/pre/AnRes" -name "AnalysisResults*${pt_bin}.root")
#         echo "anres_dir: $anres_dir"
#         echo "pt_bin: $pt_bin"

#         config_files=$(ls $config_files_path/config_*.yml)
#         parallel -j $n_parallel parallel_func ::: $config_files ::: $anres_dir ::: $pt_bin
#     done
#     end_time=$(date +"%Y-%m-%d %H:%M:%S")
#     echo "End time: $end_time"
# fi

# parallel_func() {
#     config_file=$1
#     anres_dir=$2
#     export pt_bin=$3
#     export suffix=$(basename $config_file .yml)
#     echo "suffix: $suffix"

#     start_time=$(date +"%Y-%m-%d %H:%M:%S")
#     echo "Start time: $start_time"
#     echo "python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py $config_file $anres_dir -c $cent -r $res_file -o $output_dir/trails/trails_pt_${pt_bin} -s $suffix -vn sp \
#     --preprocessed --skip_calc_weights --skip_make_yaml --skip_efficiency --skip_frac_cut_var \
#     --systematic --sys_pre $output_dir/pre"

#     python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py $config_file $anres_dir -c $cent -r $res_file -o $output_dir/trails/trails_pt_${pt_bin} -s $suffix -vn sp \
#     --preprocessed --skip_calc_weights --skip_make_yaml --skip_efficiency --skip_frac_cut_var \
#     --systematic --sys_pre $output_dir/pre

#     end_time=$(date +"%Y-%m-%d %H:%M:%S")
#     echo "End time: $end_time"
# }
# export -f parallel_func


# if [ $sigle_pt = True ]; then
#     config_files_paths=$(find $output_dir/config_sys/pt_* -type d)
#     start_time=$(date +"%Y-%m-%d %H:%M:%S")
#     echo "Start time: $start_time"
#     for config_files_path in $config_files_paths; do
#         pt_bin=$(basename "$config_files_path" | grep -oE '[0-9]+_[0-9]+')
#         anres_dir=$(find "$output_dir/pre/AnRes" -name "AnalysisResults*${pt_bin}.root")
#         echo "anres_dir: $anres_dir"
#         echo "pt_bin: $pt_bin"

#         config_files=$(ls $config_files_path/config_*.yml)
#         parallel -j $n_parallel parallel_func ::: $config_files ::: $anres_dir ::: $pt_bin
#     done
#     end_time=$(date +"%Y-%m-%d %H:%M:%S")
#     echo "End time: $end_time"
# fi



# # export outdir_config_files="$output_dir/config_syst"
# # config_files=$(ls $outdir_config_files/config_*.yml)
# # cd ../ # go to the parent directory
# # echo "config_files: $config_files"
# # echo "n_parallel: $n_parallel"
# # echo "parallel -j $n_parallel parallel_func ::: $config_files"
# # parallel -j $n_parallel parallel_func ::: $config_files