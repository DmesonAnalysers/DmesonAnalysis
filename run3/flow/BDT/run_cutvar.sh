#!/bin/bash

export config_flow="path/to/config_flow.yml"
export anres_dir="path/to/analysis_results_data.root"
export output_dir="path/to/output/BDT"
export cent="k3050"
export vn_method="sp"
export res_file="path/to/resolution.root"
export suffix="suffix"

#___________________________________________________________________________________________________________________________
# Paths for the scripts
export MakeyamlPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/make_yaml_for_ml.py"
export SimFitPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/get_vn_vs_mass.py"
export ProjPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/proj_thn_mc.py"
export EffPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/compute_efficiency.py"
export CurVarFracPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/compute_frac_cut_var.py"
export DataDrivenFracPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ComputeDataDrivenFraction.py"
export v2vsFDFracPath="/home/wuct/ALICE/local/DmesonAnalysis/run3/flow/ML/ComputeV2vsFDFrac.py"

export output_dir="${output_dir}/cutvar_${suffix}"

#___________________________________________________________________________________________________________________________
# make yaml file
echo -e "\e[34mpython3 make_yaml.py ${config_flow} -o ${output_dir} -s ${suffix}\e[0m"
python3 ${MakeyamlPath} ${config_flow} -o ${output_dir} -s ${suffix}

#___________________________________________________________________________________________________________________________
# Get the number of projection files for different cutsets
CutVarNum=$(find "${output_dir}/config" -type f -name "cutset*.yml" | wc -l)
echo "Number of cutsets: ${CutVarNum}"

#___________________________________________________________________________________________________________________________
# Cut variation (aply the cut and project)
echo -e "\e[34mpython3 cut_variation.py ${config_flow} ${anres_dir} -c ${cent} -o ${output_dir} -s ${suffix}\e[0m"
python3 cut_variation.py ${config_flow} ${anres_dir} -c ${cent} -r ${res_file} -o ${output_dir} -s ${suffix}

#___________________________________________________________________________________________________________________________
# projection for MC and apply the ptweights
for ((i=0; i<${CutVarNum}; i++)); do
	iCutSets=$(printf "%02d" $i)
	echo -e "\e[34mpyhon3 proj_thn_mc.py ${config_flow} ${output_dir}/config/cutset_${suffix}_${iCutSets}.yml -o ${output_dir} -s ${suffix}_${iCutSets}\e[0m"
	python3 ${ProjPath} ${config_flow} \
						${output_dir}/config/cutset_${suffix}_${iCutSets}.yml \
						-o ${output_dir} \
						-s ${suffix}_${iCutSets}
done

#___________________________________________________________________________________________________________________________
# compute the efficiency
for ((i=0; i<${CutVarNum}; i++)); do
	iCutSets=$(printf "%02d" $i)
	echo -e "\e[34mpython3 compute_efficiency.py ${config_flow} ${output_dir}/proj_mc/proj_mc_${suffix}_${iCutSets}.root -c ${cent} -o ${output_dir} -s ${suffix}_${iCutSets}\e[0m"
	python3 ${EffPath} ${config_flow} \
			${output_dir}/proj_mc/proj_mc_${suffix}_${iCutSets}.root \
			-c ${cent} \
			-o ${output_dir} \
			-s ${suffix}_${iCutSets}
done

#___________________________________________________________________________________________________________________________
# do the simulation fit to get the raw yields
if [ ! -d "${output_dir}/ry" ]; then mkdir -p ${output_dir}/ry; fi

for ((i=0; i<${CutVarNum}; i++)); do
	iCutSets=$(printf "%02d" $i)
	echo -e "\e[34mpython3 get_vn_vs_mass.py ${config_flow} ${cent} ${output_dir}/proj/proj_${suffix}_${iCutSets}.root -o ${output_dir}/ry -s _${suffix}_${iCutSets} -vn ${vn_method} --batch\e[0m"
	python3 ${SimFitPath} ${config_flow} ${cent} \
			${output_dir}/proj/proj_${suffix}_${iCutSets}.root \
			-o ${output_dir}/ry \
			-s _${suffix}_${iCutSets} \
			-vn ${vn_method} \
			--batch
done

#___________________________________________________________________________________________________________________________
# compute the fraction by cut variation method
echo -e "\e[34mpython3 compute_frac_cut_var.py ${config_flow} ${output_dir} -o ${output_dir} -s ${suffix}\e[0m"
python3 ${CurVarFracPath} ${config_flow} \
					${output_dir} \
					-o ${output_dir} \
					-s ${suffix}

#___________________________________________________________________________________________________________________________
# compute fraction by Data-driven method
echo -e "\e[34mpython3 ComputeDataDrivenFraction.py -i ${output_dir} -o ${output_dir} -s ${suffix}\e[0m"
python3 ${DataDrivenFracPath} \
		-i ${output_dir} \
		-o ${output_dir} \
		-s ${suffix}

#___________________________________________________________________________________________________________________________
# compute v2 vs FD fraction
python3 ${v2vsFDFracPath} ${config_flow} \
		-i ${output_dir} \
		-o ${output_dir} \
		-s ${suffix}