import os
import sys
import numpy as np
import argparse
import yaml
import shutil
import concurrent.futures
import time
sys.path.append('..')
from flow_analysis_utils import get_cut_sets_config, cut_var_image_merger
from ComputeDataDriFrac_flow import main_data_driven_frac
from ComputeV2vsFDFrac import main_v2_vs_frac
from concurrent.futures import ProcessPoolExecutor

def check_dir(dir):

	if not os.path.exists(dir):
		print(f"\033[32m{dir} does not exist, it will be created\033[0m")
		os.makedirs(dir)
	else:
		print(f"\033[33m{dir} already exists, it will be removed and recreated\033[0m")
		shutil.rmtree(dir)
		os.makedirs(dir)

	return

def run_full_cut_variation(config_flow,
                           use_preprocessed, 
						   calc_weights=False,
						   make_yaml=False, 
						   proj=False,
						   efficiency=False,
						   vn = False,
						   frac_cut_var=False,
						   data_driven_frac=False,
						   v2_vs_frac=False,
         				   merge_images=False):    
#___________________________________________________________________________________________________________________________
	# Load and copy the configuration file
	with open(config_flow, 'r') as cfgFlow:
		config = yaml.safe_load(cfgFlow)

	anres_dir = config['anresdir'] 
	cent = config['centrality'] 
	res_file = config['res_file'] 
	output = config['out_dir'] 
	suffix = config['suffix'] 
	vn_method = config['vn_method']
	n_workers = config['nworkers']
	
	CutSets, _, _, _, _ = get_cut_sets_config(config_flow)
	# REVIEW: uniformize the max cutsets variable
	mCutSets = max(CutSets)

	print(f"\033[32mINFO: Number of cutsets: {mCutSets}\033[0m")

	output_dir = f"{output}/cutvar_{suffix}"
 
	os.system(f"mkdir -p {output_dir}")
  
	# copy the configuration file
	config_suffix = 0
	os.makedirs(f'{output_dir}/config_flow', exist_ok=True)
	while os.path.exists(f'{output_dir}/config_flow/config_flow_{suffix}_{config_suffix}.yml'):
		config_suffix = config_suffix + 1
	os.system(f'cp {config_flow} {output_dir}/{os.path.splitext(os.path.basename(config_flow))[0]}_{suffix}_{config_suffix}.yml')

	# backup the results into history
	file_to_check = f"{output_dir}/V2VsFrac/V2VsFrac_{suffix}.root"
	if os.path.exists(file_to_check):
		for sub_path in ['ry', 'CutVarFrac', 'V2VsFrac']:
			os.system(f"mkdir -p {output_dir}/history/{config_suffix}/{sub_path}")
			os.system(f"cp {output_dir}/{sub_path}/* {output_dir}/history/{config_suffix}/{sub_path}")
		os.system(f"cp {output_dir}/config_flow/{os.path.splitext(os.path.basename(config_flow))[0]}_{suffix}_{config_suffix-1}.yml {output_dir}/history/{config_suffix}")

#___________________________________________________________________________________________________________________________
	# calculate the pT weights
	if calc_weights:
		check_dir(f"{output_dir}/ptweights")
		CalcWeiPath = "./ComputePtWeights.py"

		print(f"\033[32mpython3 {CalcWeiPath} {config_flow} -o {output_dir} -s {suffix}\033[0m")
		os.system(f"python3 {CalcWeiPath} {config_flow} -o {output_dir} -s {suffix}")
	else:
		print("\033[33mWARNING: Calculation of weights will not be performed\033[0m")

	if 'ptWeights_path' in config and config['ptWeights_path'] != '':
		given_ptweights = True
		given_ptWeightsPath = config['ptWeights_path']
		print(f"\033[32mINFO: Given pt weights {given_ptWeightsPath} will be used\033[0m")
	else:
		given_ptweights = False

#___________________________________________________________________________________________________________________________
	# make yaml file
	if make_yaml:
		check_dir(f"{output_dir}/config")
		MakeyamlPath = './make_yaml_for_ml.py'
		pre_process = "--preprocessed" if use_preprocessed else ""

		print(f"\033[32mpython3 {MakeyamlPath} {config_flow} {pre_process} -o {output_dir} -s {suffix}\033[0m")
		os.system(f"python3 {MakeyamlPath} {config_flow} {pre_process} -o {output_dir} -s {suffix}")
	else:
		print("\033[33mWARNING: Make yaml will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# Projection for MC and apply the ptweights
	if proj:
		ProjPath = "./proj_thn.py"
		if use_preprocessed:
			pre_process = "--preprocessed"
		else:
			pre_process = ""
			anres_files = ' '.join(anres_dir)
     	check_dir(f"{output_dir}/proj")
		print('Projecting histograms')
    	
     	def run_projections(i):
			"""Run simulation fit for a given cutset index."""
			iCutSets = f"{i:02d}"
			print(f"\033[32mProcessing cutset {iCutSets}...\033[0m")
			if not os.path.exists(f'{output_dir}/ptweights/pTweight_{suffix}.root')  and not given_ptweights:
				# REVIEW: add the list of anres files
				print(f"\033[32mpython3 {ProjPath} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml {anres_files} {pre_process} -c {cent} -r {res_file} -o {output_dir} -s {suffix}_{iCutSets}\033[0m")
				os.system(f"python3 {ProjPath} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml {anres_files} {pre_process} -c {cent} -r {res_file} -o {output_dir} -s {suffix}_{iCutSets}")
			else:
				if given_ptweights:
					ptweightsPath = given_ptWeightsPath
				else:
					ptweightsPath = f'{output_dir}/ptweights/pTweight_{suffix}.root'

				print(
					f"\033[32mpython3 {ProjPath} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml {anres_files} {pre_process} "
					f"-w {ptweightsPath} hPtWeightsFONLLtimesTAMUDcent "
					f"-wb {ptweightsPath} hPtWeightsFONLLtimesTAMUBcent "
					f"-c {cent} -r {res_file} -o {output_dir} -s {suffix}_{iCutSets} \033[0m"
				)
    			# REVIEW: add the list of anres files
				os.system(f"python3 {ProjPath} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml {anres_files} {pre_process} "
						f"-w {ptweightsPath} hPtWeightsFONLLtimesTAMUDcent "
						f"-wb {ptweightsPath} hPtWeightsFONLLtimesTAMUBcent "
						f"-c {cent} -r {res_file} -o {output_dir} -s {suffix}_{iCutSets}")

		with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
			results_proj = list(executor.map(run_projections, range(mCutSets)))
	else:
		print("\033[33mWARNING: Projection for MC will not be performed\033[0m")							

#___________________________________________________________________________________________________________________________
	# Compute the efficiency
	if efficiency:
		check_dir(f"{output_dir}/eff")
		EffPath = "./../compute_efficiency.py"

		def run_efficiency(i):
			"""Run efficiency computation for a given cutset index."""
			iCutSets = f"{i:02d}"
			print(f"\033[32mpython3 {EffPath} {config_flow} {output_dir}/proj/proj_{suffix}_{iCutSets}.root -c {cent} -o {output_dir} -s {suffix}_{iCutSets}\033[0m")
			print(f"\033[32mProcessing cutset {iCutSets}\033[0m")
			os.system(f"python3 {EffPath} {config_flow} {output_dir}/proj/proj_{suffix}_{iCutSets}.root -c {cent} -o {output_dir} -s {suffix}_{iCutSets} --batch")
		
  		with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
			results_eff = list(executor.map(run_efficiency, range(mCutSets)))
	else:
		print("\033[33mWARNING: Efficiency will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# do the simulation fit to get the raw yields
	if vn:
		check_dir(f"{output_dir}/ry")
		SimFitPath = "./../get_vn_vs_mass.py"

		def run_simfit(i):
			"""Run simulation fit for a given cutset index."""
			iCutSets = f"{i:02d}"
			print(f"\033[32mpython3 {SimFitPath} {config_flow} {cent} {output_dir}/proj/proj_{suffix}.root -o {output_dir}/ry -s _{suffix}_{iCutSets} -vn {vn_method}\033[0m")
			print(f"\033[32mProcessing cutset {iCutSets}\033[0m")
			os.system(f"python3 {SimFitPath} {config_flow} {cent} {output_dir}/proj/proj_{suffix}_{iCutSets}.root -o {output_dir}/ry -s _{suffix}_{iCutSets} -vn {vn_method} --batch")
		
  		with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
			executor.map(run_simfit, range(mCutSets))
	else:
		print("\033[33mWARNING: vn extraction will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# Compute the fraction by cut variation method
	if frac_cut_var:
		check_dir(f"{output_dir}/CutVarFrac")
		CurVarFracPath = "./compute_frac_cut_var.py"

		print(f"\033[32mpython3 {CurVarFracPath} {config_flow} {output_dir} -o {output_dir} -s {suffix}\033[0m")
		os.system(f"python3 {CurVarFracPath} {config_flow} {output_dir} -o {output_dir} -s {suffix} --batch")
	else:
		print("\033[33mWARNING: Fraction by cut variation will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# Compute fraction by Data-driven method
	if data_driven_frac:
		check_dir(f"{output_dir}/DataDrivenFrac")
		DataDrivenFracPath = "./ComputeDataDriFrac_flow.py"

		# TODO: add the combined method with runing the correlated and uncorrelated at the same time
		combined = config['combined'] if 'combined' in config else False
		print(f"\033[32mCombined method: {combined}\033[0m")
		if combined:
			if config['minimisation']['correlated']:
				print("\033[31mOnly support combined method when the minimisation is uncorrelated\033[0m")
				print("\033[31mThe combined method will not be performed\033[0m")
				# run the data-driven method with the corelated results
				print(f"\033[32mpython3 {DataDrivenFracPath} -i {output_dir} -o {output_dir} -s {suffix} -b\033[0m")
				main_data_driven_frac(inputdir=output_dir, outputdir=output_dir, suffix=suffix, batch=True, combined=False)
			else:
				# the path of corresponding results with correlated cut method
				correlatedPath = config['correlatedPath'] if 'correlatedPath' in config else ''
				if correlatedPath == '':
					print("\033[31mPlease provide the path to the corresponding correlated cut method\033[0m")
					correlatedPath = input("Path(`output_dir` in run_cutvar.sh): ")
				correlatedCutVarPath = correlatedPath + '/cutvar_' + suffix
				
				# the path of the combined results
				if 'pass3' in correlatedPath:
					outputdir_corr = output_dir.replace(f'{output_dir.split("/")[-2]}', 'combined_pass3')
				else:
					outputdir_corr = output_dir.replace(f'{output_dir.split("/")[-2]}', 'combined')
				check_dir(f"{outputdir_corr}/DataDrivenFrac")
				
				# run the data-driven method with the combined results and the uncorrelated results
				print(f"\033[32mpython3 {DataDrivenFracPath} -i {output_dir} -o {output_dir} -s {suffix} -b -comb -cc {correlatedCutVarPath} -oc {outputdir_corr}\033[0m") 
				main_data_driven_frac(inputdir=output_dir, outputdir=output_dir, suffix=suffix, batch=True, \
										combined=True, correlatedCutVarPath=correlatedCutVarPath, outputdir_combined=outputdir_corr)

		# run the data-driven method with the uncorrelated or correlated results
		else:
			print(f"\033[32mpython3 {DataDrivenFracPath} -i {output_dir} -o {output_dir} -s {suffix} -b\033[0m")
			main_data_driven_frac(inputdir=output_dir, outputdir=output_dir, suffix=suffix, batch=True, combined=False)

	else:
		print("\033[33mWARNING: Fraction by Data-driven method will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# Compute v2 vs fraction
	if v2_vs_frac:
		check_dir(f"{output_dir}/V2VsFrac")
		v2vsFDFracPath = "./ComputeV2vsFDFrac.py"
		
		combined = config['combined'] if 'combined' in config else False
		
		if combined:
			if 'pass3' in output_dir:
				inputdir_combined = output_dir.replace(f'{output_dir.split("/")[-2]}', 'combined_pass3')
			else:
				inputdir_combined = output_dir.replace(f'{output_dir.split("/")[-2]}', 'combined')
			outputdir_combined = inputdir_combined
			check_dir(f"{outputdir_combined}/V2VsFrac")
			
			print(f"\033[32mthe combined method will be performed\033[0m")
			print(f"\033[32mpython3 {v2vsFDFracPath} {config_flow} -i {output_dir} -o {output_dir} -s {suffix} -comb -ic {inputdir_combined} -oc {outputdir_combined}\033[0m")
			main_v2_vs_frac(config=config_flow, inputdir=output_dir, outputdir=output_dir, suffix=suffix, \
								combined=True, inputdir_combined=inputdir_combined, outputdir_combined=outputdir_combined)
		else:
			print(f"\033[32mpython3 {v2vsFDFracPath} {config_flow} -i {output_dir} -o {output_dir} -s {suffix}\033[0m")
			main_v2_vs_frac(config=config_flow, inputdir=output_dir, outputdir=output_dir, suffix=suffix, combined=False)
	else:
		print("\033[33mWARNING: v2 vs fraction will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# Merge cut var figures in multipanel images
	if merge_images:
		print(f"\033[32m\nCut_var_image_merger({config_flow}, {output_dir}, {suffix})\033[0m")
		cut_var_image_merger(config, output_dir, suffix)
	return


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Arguments')
	parser.add_argument('flow_config', metavar='text', default='config_flow_d0.yml', help='configuration file')
	parser.add_argument("--preprocessed", "-prep", action="store_true", help="use preprocessed input")
	parser.add_argument("--do_calc_weights", "-cw", action="store_true", help="perform calculation of weights")
	parser.add_argument("--do_make_yaml", "-my", action="store_true", help="perform make yaml")
	parser.add_argument("--do_projections", "-pd", action="store_true", help="perform projections")
	parser.add_argument("--do_efficiency", "-e", action="store_true", help="perform efficiency")
	parser.add_argument("--do_vn", "-vn", action="store_true", help="perform vn extraction")
	parser.add_argument("--do_frac_cut_var", "-f", action="store_true", help="perform fraction by cut variation")
	parser.add_argument("--do_data_driven_frac", "-ddf", action="store_true", help="perform fraction by data-driven method")
	parser.add_argument("--do_v2_vs_frac", "-v2fd", action="store_true", help="perform v2 vs FD fraction")
	parser.add_argument("--do_merge_images", "-mergeim", action="store_true", help="perform v2 vs FD fraction")
	args = parser.parse_args()

	start_time = time.time()
	run_full_cut_variation(args.flow_config,
                           args.preprocessed,
						   args.do_calc_weights,
						   args.do_make_yaml, 
						   args.do_projections,
						   args.do_efficiency, 
						   args.do_vn,
						   args.do_frac_cut_var, 
						   args.do_data_driven_frac, 
						   args.do_v2_vs_frac,
						   args.do_merge_images)

	end_time = time.time()
	execution_time = end_time - start_time
	print(f"\033[34mTotal execution time: {execution_time:.2f} seconds\033[0m")
 