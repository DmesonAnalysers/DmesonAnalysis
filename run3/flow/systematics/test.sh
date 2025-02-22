#！bin/bash
export config_modifies="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/systematics/modifications_config_fit.yml"
export config_default="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/config/config_flow_D0_3050_mother.yml"
export res_file="/media/wuct/wulby/ALICE/AnRes/resolution/output_reso/resosp3050l_PASS4_full_PbPb_Reso.root"
export cent="k3050"
export output_dir="/home/wuct/ALICE/local/Results/test/flow/systematics"

export sigle_pt=True
export n_parallel=15


mkdir -p $output_dir

# produce the yaml files for the systematic: config_flow, config_pre, config_defalut, config_default_corelated
rm -rf $output_dir/config_sys
python3 ./make_yaml_for_syst.py $config_default -m $config_modifies -o $output_dir -mb

# produce the pre-processed files, 
# --preselected is optional, it means the pre-selection is applied, and flow/systematics/modifications_config_fit.yml should contain the AnRes file with single pt bin
# python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/tool/pre_process.py $output_dir/config_sys/config_pre.yml --out_dir $output_dir --pre --preselected

export anres_dir=$(ls $output_dir/pre/AnRes/AnalysisResults*.root)

# Uncorelated: provide the projection of MC, and the efficiency
# python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/config_default.yml \
#         $anres_dir -c $cent -r $res_file -o $output_dir/pre -s central -vn sp \
#         --preprocessed --skip_frac_cut_var --skip_data_driven_frac --skip_v2_vs_frac

# # Correlated: provide the fraction
# python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py ${output_dir}/config_sys/config_default_corelated.yml \
#         $anres_dir -c $cent -r $res_file -o $output_dir/pre -s correlated -vn sp \
#         --preprocessed --skip_calc_weights --skip_make_yaml --skip_proj_mc --skip_efficiency --skip_vn

parallel_func() {
    config_file=$1
    anres_dir=$2
    export pt_bin=$3
    export suffix=$(basename $config_file .yml)
    echo "suffix: $suffix"

    start_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "Start time: $start_time"
    echo "python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py $config_file $anres_dir -c $cent -r $res_file -o $output_dir/trails/trails_pt_${pt_bin} -s $suffix -vn sp \
    --preprocessed --skip_calc_weights --skip_make_yaml --skip_efficiency --skip_frac_cut_var \
    --systematic --sys_pre $output_dir/pre"

    python3 /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/BDT/run_cutvar.py $config_file $anres_dir -c $cent -r $res_file -o $output_dir/trails/trails_pt_${pt_bin} -s $suffix -vn sp \
    --preprocessed --skip_calc_weights --skip_make_yaml --skip_efficiency --skip_frac_cut_var \
    --systematic --sys_pre $output_dir/pre

    end_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "End time: $end_time"
}
export -f parallel_func


if [ $sigle_pt = True ]; then
    config_files_paths=$(find $output_dir/config_sys/pt_* -type d)
    start_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "Start time: $start_time"
    for config_files_path in $config_files_paths; do
        pt_bin=$(basename "$config_files_path" | grep -oE '[0-9]+_[0-9]+')
        anres_dir=$(find "$output_dir/pre/AnRes" -name "AnalysisResults*${pt_bin}.root")
        echo "anres_dir: $anres_dir"
        echo "pt_bin: $pt_bin"

        config_files=$(ls $config_files_path/config_*.yml)
        parallel -j $n_parallel parallel_func ::: $config_files ::: $anres_dir ::: $pt_bin
    done
    end_time=$(date +"%Y-%m-%d %H:%M:%S")
    echo "End time: $end_time"
fi



# export outdir_config_files="$output_dir/config_syst"
# config_files=$(ls $outdir_config_files/config_*.yml)
# cd ../ # go to the parent directory
# echo "config_files: $config_files"
# echo "n_parallel: $n_parallel"
# echo "parallel -j $n_parallel parallel_func ::: $config_files"
# parallel -j $n_parallel parallel_func ::: $config_files