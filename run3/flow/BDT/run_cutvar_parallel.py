import os
import sys
import argparse
import yaml
import shutil
import concurrent.futures
import time
sys.path.append('..')
from flow_analysis_utils import get_cut_sets_config, cut_var_image_merger

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
						   proj_data=False,
						   proj_mc=False,
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

	cent = config['centrality']
	res_file = config['res_file']
	output = config['out_dir']
	suffix = config['suffix']
	vn_method = config['vn_method']
	n_workers = 6

	CutSets, _, _, _, _ = get_cut_sets_config(config_flow)
	nCutSets = max(CutSets)
	print(f"\033[32mNumber of cutsets: {nCutSets}\033[0m")

	output_dir = f"{output}/cutvar_{suffix}"
 
	if not os.path.exists(output_dir):
		print(f"Creating {output_dir}")
		os.makedirs(output_dir)
	else:
		print(f"Directory already exists: {output_dir}")

	# copy the configuration file
	config_suffix = 1
	while os.path.exists(f'{output_dir}/{os.path.splitext(os.path.basename(config_flow))[0]}_{suffix}_{config_suffix}.yml'):
		config_suffix = config_suffix + 1
	os.system(f'cp {config_flow} {output_dir}/{os.path.splitext(os.path.basename(config_flow))[0]}_{suffix}_{config_suffix}.yml')

#___________________________________________________________________________________________________________________________
	# calculate the pT weights
	if calc_weights:
		check_dir(f"{output_dir}/ptweights")
		CalcWeiPath = "./ComputePtWeights.py"

		print(f"\033[32mpython3 {CalcWeiPath} {config_flow} -o {output_dir} -s {suffix}\033[0m")
		os.system(f"python3 {CalcWeiPath} {config_flow} -o {output_dir} -s {suffix}")
	else:
		print("\033[33mWARNING: Calculation of weights will not be performed\033[0m")

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
	#TODO: 1.keep the yaml file for the user to check 2.modify the proj_thn 3.use make_combination in proj_thn.py


#___________________________________________________________________________________________________________________________
	ProjPath = "./proj_thn.py"
	pre_process = "--preprocessed" if use_preprocessed else ""
	proj_data = "--proj_data" if proj_data else ""
	proj_mc = "--proj_mc" if proj_mc else ""
	def run_projections(i):
		"""Run simulation fit for a given cutset index."""
		iCutSets = f"{i:02d}"
		print('CIAOOOOO')
		print(f"\033[32mProcessing cutset {iCutSets}...\033[0m")
		if not os.path.exists(f'{output_dir}/ptweights/pTweight_{suffix}.root'):
			print('USE PT WEIGHTS')
			print(f"\033[32mpython3 {ProjPath} {proj_data} {proj_mc} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml {pre_process} -c {cent} -r {res_file} -o {output_dir} -s {suffix}_{iCutSets}\033[0m")
			os.system(f"python3 {ProjPath} {proj_data} {proj_mc} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml {pre_process} -c {cent} -r {res_file} -o {output_dir} -s {suffix}_{iCutSets}")
		else:
			print('NOT USE PT WEIGHTS')
			print(
				f"\033[32mpython3 {ProjPath} {proj_data} {proj_mc} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml {pre_process} "
				f"-w {output_dir}/ptweights/pTweight_{suffix}.root hPtWeightsFONLLtimesTAMUDcent "
				f"-wb {output_dir}/ptweights/pTweight_{suffix}.root hPtWeightsFONLLtimesTAMUBcent "
				f"-c {cent} -r {res_file} -o {output_dir} -s {suffix}_{iCutSets} \033[0m"
			)
			os.system(f"python3 {ProjPath} {proj_data} {proj_mc} {config_flow} {output_dir}/config/cutset_{suffix}_{iCutSets}.yml {pre_process} \
					-w {output_dir}/ptweights/pTweight_{suffix}.root hPtWeightsFONLLtimesTAMUDcent \
					-wb {output_dir}/ptweights/pTweight_{suffix}.root hPtWeightsFONLLtimesTAMUBcent -c {cent} -r {res_file} -o {output_dir} -s {suffix}_{iCutSets}")
		print('CIAO END')
  
	if proj_mc or proj_data:
		print('Projecting histograms')
		with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
			results_proj = list(executor.map(run_projections, range(nCutSets)))
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
			print(f"\033[32mProcessing cutset {iCutSets} for efficiency...\033[0m")
			command = f"python3 {EffPath} {config_flow} {output_dir}/proj/proj_{suffix}_{iCutSets}.root -c {cent} -o {output_dir} -s {suffix}_{iCutSets} --batch"
			print(f"\033[32m{command}\033[0m")
			os.system(command)

		with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
			results_eff = list(executor.map(run_efficiency, range(nCutSets)))
	else:
		print("\033[33mWARNING: Efficiency will not be performed\033[0m")

	SimFitPath = "./../get_vn_vs_mass.py"
	def run_simfit(i):
		"""Run simulation fit for a given cutset index."""
		iCutSets = f"{i:02d}"
		command = f"python3 {SimFitPath} {config_flow} {cent} {output_dir}/proj/proj_{suffix}_{iCutSets}.root -o {output_dir}/ry -s _{suffix}_{iCutSets} -vn {vn_method} --batch"
		print(f"\033[32mProcessing cutset {iCutSets}...\033[0m")
		os.system(command)

	if vn:
		# Ensure the output directory exists
		check_dir(f"{output_dir}/ry")
		with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
			executor.map(run_simfit, range(nCutSets))
	else:
		print("\033[33mWARNING: vn extraction will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# Compute the fraction by cut variation method
	if frac_cut_var:
		check_dir(f"{output_dir}/CutVarFrac")
		CurVarFracPath = "./compute_frac_cut_var.py"

		print(f"\033[32mpython3 {CurVarFracPath} {config_flow} -dir {output_dir} -s {suffix}\033[0m")
		os.system(f"python3 {CurVarFracPath} {config_flow} -dir {output_dir} -s {suffix} --batch")
	else:
		print("\033[33mWARNING: Fraction by cut variation will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# Compute fraction by Data-driven method
	if data_driven_frac:
		check_dir(f"{output_dir}/DataDrivenFrac")
		DataDrivenFracPath = "./ComputeDataDriFrac_flow.py"

		print(f"\033[32mpython3 {DataDrivenFracPath} -dir {output_dir} -s {suffix}\033[0m")
		os.system(f"python3 {DataDrivenFracPath} -dir {output_dir} -s {suffix} --batch")
	else:
		print("\033[33mWARNING: Fraction by Data-driven method will not be performed\033[0m")

#___________________________________________________________________________________________________________________________
	# Compute v2 vs fraction
	if v2_vs_frac:
		check_dir(f"{output_dir}/V2VsFrac")
		v2vsFDFracPath = "./ComputeV2vsFDFrac.py"

		print(f"\033[32mpython3 {v2vsFDFracPath} {config_flow} -dir {output_dir} -s {suffix}\033[0m")
		os.system(f"python3 {v2vsFDFracPath} {config_flow} -dir {output_dir} -s {suffix}")
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
	parser.add_argument("--do_calc_weights", "-cw", action="store_true", help="skip calculation of weights")
	parser.add_argument("--do_make_yaml", "-my", action="store_true", help="skip make yaml")
	parser.add_argument("--do_proj_data", "-pd", action="store_true", help="skip projection for data")
	parser.add_argument("--do_proj_mc", "-pm", action="store_true", help="skip projection for MC")
	parser.add_argument("--do_efficiency", "-e", action="store_true", help="skip efficiency")
	parser.add_argument("--do_vn", "-vn", action="store_true", help="skip vn extraction")
	parser.add_argument("--do_frac_cut_var", "-f", action="store_true", help="skip fraction by cut variation")
	parser.add_argument("--do_data_driven_frac", "-ddf", action="store_true", help="skip fraction by data-driven method")
	parser.add_argument("--do_v2_vs_frac", "-v2fd", action="store_true", help="skip v2 vs FD fraction")
	parser.add_argument("--do_merge_images", "-mergeim", action="store_true", help="skip v2 vs FD fraction")
	args = parser.parse_args()

	start_time = time.time()
	run_full_cut_variation(args.flow_config, 
                           args.preprocessed,
						   args.do_calc_weights,
						   args.do_make_yaml, 
						   args.do_proj_data, 
						   args.do_proj_mc, 
						   args.do_efficiency, 
						   args.do_vn,
						   args.do_frac_cut_var, 
						   args.do_data_driven_frac, 
						   args.do_v2_vs_frac,
						   args.do_merge_images)

	end_time = time.time()
	execution_time = end_time - start_time
	print(f"\033[34mTotal execution time: {execution_time:.2f} seconds\033[0m")
 