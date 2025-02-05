#!/bin/bash

#Parameters
#----------
#- config (str): path of directory with config files
#- an_res_file (str): path of directory with analysis results
#- centrality (str): centrality class
#- resolution (str/int): resolution file or resolution value
#- outputdir (str): output directory
#- suffix (str): suffix for output files
#- vn_method (str): vn technique (sp, ep, deltaphi)
#- skip_calc_weights (bool): skip calculation of weights
#- skip_make_yaml (bool): skip make yaml
#- skip_cut_variation (bool): skip cut variation
#- skip_proj_mc (bool): skip projection for MC
#- skip_efficiency (bool): skip efficiency
#- skip_vn (bool): skip vn extraction
#- skip_frac_cut_var (bool): skip fraction by cut variation
#- skip_data_driven_frac (bool): skip fraction by data-driven method
#- skip_v2_vs_frac (bool): skip v2 vs FD fraction
#----------

# export config_flow="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/config_flow.yml"
export config_flow="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/config_flow.yml"
export anres_dir="/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s2_0.root \
    /media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s2_1.root \
    /media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s2_2.root \
    /media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s2_3.root \
    /media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/324130/temp_merged_s2_4.root"
export output_dir="/home/wuct/ALICE/local/Results/BDT/full/uncorrelated"
export cent="k3050"
export vn_method="sp"
export res_file="/media/wuct/wulby/ALICE/AnRes/resolution/output_reso/resospk3050_inte.root"
export suffix="pt2_3"

export use_prep=True # True or False (use pre-processed inputs for projections)
export sprep=True # True or False (perform pre-processing of grid output files)
export spw=True # True or False (skip calculation of weights)
export smy=False # True or False (skip make yaml)
export scv=True # True or False (skip cut variation), not used anymore
export spm=True # True or False (skip projection for MC)
export seff=True # True or False (skip efficiency)
export svn=True # True or False (skip vn extraction)
export sfcv=True # True or False (skip fraction by cut variation)
export sddf=True # True or False (skip fraction by data-driven method)
export sv2vf=True # True or False (skip v2 vs fraction)

if [ $use_prep = False ]; then
	export use_preprocessed=""
else
	export use_preprocessed="--preprocessed"
fi

if [ $sprep = False ]; then
    export skip_pre_process=""
else
    export skip_pre_process="--skip_pre_process"
fi

if [ $spw = False ]; then
	export skip_calc_weights=""
else
	export skip_calc_weights="--skip_calc_weights"
fi

if [ $smy = False ]; then
	export skip_make_yaml=""
else
	export skip_make_yaml="--skip_make_yaml"
fi

if [ $scv = False ]; then
	export skip_cut_variation=""
else
	export skip_cut_variation="--skip_cut_variation"
fi

if [ $spm = False ]; then
	export skip_proj_mc=""
else
	export skip_proj_mc="--skip_proj_mc"
fi

if [ $seff = False ]; then
	export skip_efficiency=""
else
	export skip_efficiency="--skip_efficiency"
fi

if [ $svn = False ]; then
	export skip_vn=""
else
	export skip_vn="--skip_vn"
fi

if [ $sfcv = False ]; then
	export skip_frac_cut_var=""
else
	export skip_frac_cut_var="--skip_frac_cut_var"
fi

if [ $sddf = False ]; then
	export skip_data_driven_frac=""
else
	export skip_data_driven_frac="--skip_data_driven_frac"
fi

if [ $sv2vf = False ]; then
	export skip_v2_vs_frac=""
else
	export skip_v2_vs_frac="--skip_v2_vs_frac"
fi

python3 run_cutvar.py $config_flow $anres_dir -c $cent -r $res_file -o $output_dir -s $suffix -vn $vn_method $use_preprocessed \
					  $skip_pre_process \
					  $skip_calc_weights \
					  $skip_make_yaml \
					  $skip_cut_variation \
					  $skip_proj_mc \
					  $skip_efficiency \
					  $skip_vn \
					  $skip_frac_cut_var \
					  $skip_data_driven_frac \
					  $skip_v2_vs_frac