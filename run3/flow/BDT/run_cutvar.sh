#!/bin/bash

#Parameters
#----------
# config_flow (str): path of directory with config files
#- usepreprocessed (bool): use pre-processed input
#- docw (bool): calculate pt weights
#- domy (bool): produce yaml config files
#- doproj (bool): project sparses (TODO: separate data and mc)
#- doeff (bool): perform efficiency calculation
#- dovn (bool): perform simultaneous fits
#- dofcv (bool): perform fraction by cut variation
#- doddf (bool): perform fraction by data-driven method
#- dov2vf (bool): perform v2 vs FD fraction
#- domergeimages (bool): perform cutvar images merging
#----------

export config_flow="/home/mdicosta/FlowDplus/FinalResults/templs_from_histo_parallel/config_3040_uncorrelated_templs_from_histo_parallel.yml"
export usepreprocessed=True
export docw=True
export domy=True
export doproj=True
export doeff=True
export dovn=True
export dofcv=False
export doddf=False
export dov2vf=False
export domergeimages=False

export usepreprocessed=$([ "$usepreprocessed" = "false" ] && echo "" || echo "--use_preprocessed")
export calc_weights=$([ "$docw" = "false" ] && echo "" || echo "--do_calc_weights")
export make_yaml=$([ "$domy" = "false" ] && echo "" || echo "--do_make_yaml")
export proj=$([ "$doproj" = "false" ] && echo "" || echo "--do_projections")
export efficiency=$([ "$doeff" = "false" ] && echo "" || echo "--do_efficiency")
export vn=$([ "$dovn" = "false" ] && echo "" || echo "--do_vn")
export frac_cut_var=$([ "$dofcv" = "false" ] && echo "" || echo "--do_frac_cut_var")
export data_driven_frac=$([ "$doddf" = "false" ] && echo "" || echo "--do_data_driven_frac")
export v2_vs_frac=$([ "$dov2vf" = "false" ] && echo "" || echo "--do_v2_vs_frac")
export merge_images=$([ "$domergeimages" = "false" ] && echo "" || echo "--do_merge_images")

python3 run_cutvar.py $config_flow \
					  $usepreprocessed \
					  $calc_weights \
					  $make_yaml \
					  $proj \
					  $efficiency \
					  $vn \
					  $frac_cut_var \
					  $data_driven_frac \
					  $v2_vs_frac \
					  $merge_images
