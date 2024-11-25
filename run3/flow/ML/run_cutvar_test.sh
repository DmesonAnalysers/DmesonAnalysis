#!/bin/bash

export config_flow="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/config_flow_d0_test_ml.yml"
export anres_dir="/media/wuct/wulby/ALICE/AnRes/D0_flow/pass4/ML/Results/AnalysisResults_BDT_DATA_test.root"
export output_dir="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/test"
export cent="k3050"
export vn_method="sp"
export res_file="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/Results/2060/k3050/small/sp/resolution/resosp3050s_291131_inte_gain_Reso.root"
export suffix="test"
export isMC="--isMC"

export MakeyamlPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/make_yaml_for_ml.py"
export SimFitPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/get_vn_vs_mass.py"
export ProjPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/proj_thn_mc.py"
export EffPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/compute_efficiency.py"
export CurVarFracPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/compute_frac_cut_var.py"
export DataDrivenFracPath="/home/wuct/ALICE/local/DmesonAnalysis/ComputeDataDrivenFraction.py"
export v2vsFDFracPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/ComputeV2vsFDFrac.py"

export output_dir="${output_dir}/cutvar_${suffix}"

# make yaml file
echo -e "\e[34mpython3 make_yaml.py ${config_flow} -o ${output_dir} -s ${suffix}\e[0m"
python3 ${MakeyamlPath} ${config_flow} -o ${output_dir} -s ${suffix}

# Get the number of projection files for different cutsets
CutVarNum=$(find "${output_dir}/config" -type f -name "cutset*.yml" | wc -l)
echo "Number of cutsets: ${CutVarNum}"

# # Cut variation (aply the cut and project)
# echo -e "\e[34mpython3 cut_variation.py ${config_flow} ${anres_dir} -c ${cent} -o ${output_dir} -s ${suffix}\e[0m"
# python3 cut_variation.py ${config_flow} ${anres_dir} -c ${cent} -r ${res_file} -o ${output_dir} -s ${suffix}

# # projection for MC and apply the ptweights
# for ((i=0; i<${CutVarNum}; i++)); do
# 	iCutSets=$(printf "%02d" $i)
# 	echo -e "\e[31mPyyhon3 proj_thn_mc.py ${config_flow} ${output_dir}/config/cutset_${suffix}_${iCutSets}.yml -o ${output_dir} -s ${suffix}_${iCutSets}\e[0m"
# 	echo -e "\e[31mProcessing cutset ${iCutSets}\e[0m"
# 	python3 ${ProjPath} ${config_flow} \
# 						${output_dir}/config/cutset_${suffix}_${iCutSets}.yml \
# 						-o ${output_dir} \
# 						-s ${suffix}_${iCutSets}
# done

# # compute the efficiency
# for ((i=0; i<${CutVarNum}; i++)); do
# 	iCutSets=$(printf "%02d" $i)
# 	echo -e "\e[31mProcessing cutset ${iCutSets}\e[0m"
# 	python3 ${EffPath} ${config_flow} \
# 			${output_dir}/proj_mc/proj_mc_${suffix}_${iCutSets}.root \
# 			-c ${cent} \
# 			-o ${output_dir} \
# 			-s ${suffix}_${iCutSets}
# done

# # do the simulation fit to get the raw yields
# if [ ! -d "${output_dir}/ry" ]; then mkdir -p ${output_dir}/ry; fi

# for ((i=0; i<${CutVarNum}; i++)); do
# 	iCutSets=$(printf "%02d" $i)
# 	echo -e "\e[31mProcessing cutset ${iCutSets}\e[0m"
# 	python3 ${SimFitPath} ${config_flow} ${cent} \
# 			${output_dir}/proj/cutset_${suffix}_${iCutSets}.root \
# 			-o ${output_dir}/ry \
# 			-s _${suffix}_${iCutSets} \
# 			-vn ${vn_method} \
# 			--batch
# done

# # compute the fraction by cut variation method
# python3 ${CurVarFracPath} ${config_flow} \
# 					${output_dir} \
# 					-o ${output_dir} \
# 					-s ${suffix}

# mkdir -p ${output_dir}/DataDrivenFrac
# # compute fraction by Data-driven method
# for ((i=0; i<${CutVarNum}; i++)); do
# 	iCutSets=$(printf "%02d" $i)
# 	echo -e "\e[31mProcessing cutset ${iCutSets}\e[0m"
# 	python3 ${DataDrivenFracPath} \
# 			${output_dir}/eff/eff_${suffix}.root \
# 			${output_dir}/CutVarFrac/CutVarFrac_${suffix}.root \
# 			${output_dir}/DataDrivenFrac/DataDrivenFrac_${suffix}_${iCutSets}.root
# done

mkdir -p ${output_dir}/v2vsFDFrac
# compute v2 vs FD fraction
python3 ${v2vsFDFracPath} /home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/config_D0_v2frac_PbPb3050_unc.yml \
					${output_dir}/v2vsFDFrac

#TODO: since the reflection needed to be add in simutaneouse fit, so the MC process should be done first
#_Process_MC__________________________________________________________________________________________________________

