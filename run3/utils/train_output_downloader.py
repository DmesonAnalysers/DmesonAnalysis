''' 
This script is used to download the output files from the grid,
merge them and convert them to parquet format.
Usage:
    python train_output_downloader.py <config_file> [--aod] [--analysis] [--parquet]
'''

import os
import json
import argparse
import sys
import pandas as pd
import ROOT
import uproot
from sklearn.model_selection import train_test_split
import xml.etree.ElementTree as ET
# pylint: disable=no-member

def split_dataset(df, config):
    """
    Split the dataset based on given selections and save the resulting dataframes as parquet files.

    Args:
        df (pandas.DataFrame): The input dataframe.
        config (dict): Configuration parameters.

    Returns:
        None
    """

    folder = "MC" if config["isMC"] else "Data"
    output_filename_no_suffix = config["output_directory"] + "/" + \
        folder + f"/Train{config['train_run']}/{config['suffix']}"

    for sel_name, selection in config["mc_selections"].items():
        df_sel = df.query(selection)
        df_sel.to_parquet(output_filename_no_suffix + f"_{sel_name}.parquet")
        if config["train_fraction"] < 1:
            df_train, df_eff = train_test_split(df_sel,
                train_size=config["train_fraction"],random_state=42)
            df_train.to_parquet(output_filename_no_suffix + f"_{sel_name}_Train.parquet")
            df_eff.to_parquet(output_filename_no_suffix + f"_{sel_name}_Eff.parquet")
            del df_train, df_eff
        del df_sel



def get_files_to_download(grid, config):
    """
    Retrieves a list of files to download based on the given grid and configuration.
    Args:
        grid (ROOT.TJAlien): The grid object used for file operations.
        config (dict): The configuration object containing input information.
    Returns:
        A list of files to download.
    """
    with open(config["input"], "r", encoding="utf8") as f:
        input_files = f.readlines()
        input_files = [x.strip() for x in input_files] # Remove leading/trailing whitespaces

    file_list = []
    if config["is_slim"]:
        file_list = input_files
    else:
        for file in input_files:
            grid_result = grid.Ls(file)
            if not grid_result:
                print(f"\033[93mWARINING\033[0m: File {file} not found on the grid.")
                continue
            file_info_list = grid_result.GetFileInfoList()

            for i in range(file_info_list.GetSize()):
                directory_name = grid_result.GetFileName(i)

                if "aod_collection" not in directory_name and "analysis" \
                    not in directory_name and "full_config" not in directory_name:
                    file_list.append(file + '/' + directory_name)

    return file_list

def download_analysis_results(input_list_filename, config):
    """
    Downloads analysis results from input file names and merges them into a single output file.

    Args:
        input_list_filename (list): List of input file names.
        config (dict): Configuration dictionary containing parameters.

    Returns:
        None
    """

    merge_analysis_results = ROOT.TFileMerger(False)
    count = 0

    for file_name in input_list_filename:
        count += 1
        if count < config["start_file"]:
            continue
        if count > config["max_files_to_download"] + config["start_file"]:
            break

        file_path = file_name.strip() # Remove leading/trailing whitespaces
        merge_analysis_results.AddFile(f"alien://{file_path}/AnalysisResults.root")

    # Create folder for output
    folder = "MC" if config["isMC"] else "Data"
    output_directory = config["output_directory"] + "/" + folder + f"/Train{config['train_run']}"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    analysis_output_filename = output_directory + f"/AnalysisResults_{config['suffix']}.root"

    merge_analysis_results.OutputFile(analysis_output_filename, "RECREATE")
    merge_analysis_results.Merge()

def download_aod(input_list_filename, config):
    """
    Downloads AOD files from a list of input file names and merges them into a single output file.
    Args:
        input_list_filename (str): The path to the file containing the list of input file names.
        config (dict): A dictionary containing configuration parameters.
    Returns:
        None
    """

    merge_aod = ROOT.TFileMerger(False)
    count = 0
    total_size = 0

    for file_name in input_list_filename:
        count += 1
        if count < config["start_file"]:
            continue
        if count > config["max_files_to_download"] + config["start_file"]:
            break

        file_path = file_name.strip() # Remove leading/trailing whitespaces
        file = ROOT.TFile.Open(f"alien://{file_path}/AO2D.root")
        total_size += file.GetSize()
        merge_aod.AddAdoptFile(file)

    total_size = total_size / 1024 / 1024
    print(f"Total size of files to download: {total_size:.2f} MB")
    if (total_size / 1024) > 3:
        print("\033[93mWARINING\033[0m: Total size of files to download "
            "is greater than 3 GB. Are you sure you want to continue? (y/n)")
        answer = input()
        if answer.lower() != 'y':
            print("Exiting.")
            sys.exit()

    # Create folder for output
    folder = "MC" if config["isMC"] else "Data"
    output_directory = config["output_directory"] + "/" + folder + f"/Train{config['train_run']}"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    analysis_output_filename = output_directory + f"/AO2D_{config['suffix']}.root"

    merge_aod.OutputFile(analysis_output_filename, "RECREATE")
    merge_aod.Merge()

def download_full_config(input_list_filename, config):
    """
    Downloads full_config.json file associated to train run
    Args:
        input_list_filename (str): The path to the file containing the list of input file names.
        config (dict): A dictionary containing configuration parameters.
    Returns:
        None
    """
    
    # Create folder for output
    folder = "MC" if config["isMC"] else "Data"
    output_directory = config["output_directory"] + "/" + folder + f"/Train{config['train_run']}"
    if not os.path.exists(output_directory):
        os.makedirs(output_directory)
    
    os.system(f"alien.py cp {input_list_filename[0]}/Stage_1.xml file:{output_directory}")

    # Load and parse the XML file
    xml_file = f"{output_directory}/Stage_1.xml"
    tree = ET.parse(xml_file)
    root = tree.getroot()

    # Find the specific event with name="1"
    event = root.find(".//event[@name='1']")  # XPath to find <event> with attribute name="1"
    if event is not None:
        file_element = event.find("file")
        if file_element is not None:
            lfn = os.path.dirname(file_element.get("lfn"))
            os.system(f"alien.py cp {lfn}/full_config.json file:{output_directory}")
            with open(f"{output_directory}/full_config.json", "r") as infile:
                data = json.load(infile)
            
            # Extract all dictionaries from the "configuration" fields
            task_cfg = {}
            for workflow in data["workflows"]:
                configuration = workflow.get("configuration", {})
                task_cfg.update(configuration)

            with open(f"{output_directory}/full_config.json", "w") as outfile:
                json.dump(task_cfg, outfile, indent=4)
            os.remove(xml_file)
            print('full_config.json loaded!')
            
def convert_aod_to_parquet(config): # pylint: disable=too-many-locals
    """
    Converts AOD (Analysis Object Data) file to Parquet format.
    Args:
        config (dict): Configuration parameters for the conversion process.
    Returns:
        None
    """

    folder = "MC" if config["isMC"] else "Data"
    input_filename = config["output_directory"] + "/" + folder + \
        f"/Train{config['train_run']}/AO2D_{config['suffix']}.root"
    output_filename = config["output_directory"] + "/" + folder + \
        f"/Train{config['train_run']}/{config['suffix']}.parquet"

    if os.path.exists(output_filename):
        print("File already exists. Removing it...")
        os.remove(output_filename)

    print("Converting to parquet...", end="\r")
    dfs = []
    with uproot.open(input_filename) as f:
        for folder_name, run_folder in f.items(recursive=False): # Loop over the run folders
            if "DF" in folder_name:
                dfs_folder = []
                if "*" in config["tree_name"]: # Loop over all trees in the folder
                    for obj_name, class_name in run_folder.classnames().items():
                        if "TTree" in class_name:
                            dfs_folder.append(run_folder[obj_name].arrays(library="pd"))
                else:
                    if not isinstance(config["tree_name"], list):
                        config["tree_name"] = [config["tree_name"]]
                    for tree_name in config["tree_name"]:
                        dfs_folder.append(run_folder[tree_name].arrays(library="pd"))
                dfs.append(pd.concat(dfs_folder, axis=1))
                del dfs_folder
                if config["selections"]:
                    dfs[-1] = dfs[-1].query(config["selections"])

    df = pd.concat(dfs)
    df.to_parquet(output_filename)
    print("Converting to parquet... Done!")

    if config["isMC"]:
        if config["mc_selections"]:
            split_dataset(df, config)

        if config["train_fraction"] < 1:
            df_train, df_eff = train_test_split(df,
                train_size=config["train_fraction"], random_state=42)
            df_train.to_parquet(output_filename.replace(".parquet", "_Train.parquet"))
            df_eff.to_parquet(output_filename.replace(".parquet", "_Eff.parquet"))
            del df

def download_files_from_grid(config, aod=False, analysis=False, parquet=False, full_config=False):
    """
    Downloads files from the grid based on the provided configuration.

    Args:
        config (object): The configuration object.
        aod (bool, optional): Flag indicating whether to download AOD files.
        analysis (bool, optional): Flag indicating whether to download analysis results.
        parquet (bool, optional): Flag indicating whether to convert AOD files to Parquet format.
        full_config (bool, optional): Flag indicating whether to download the full_config.json file.
    Note:
        If no flags are provided, all operations are performed.
    """

    if not aod and not analysis and not parquet:
        aod = True
        analysis = True
        parquet = True
    if not full_config:
        full_config = True

    if aod or analysis:
        # Get files from grid
        grid = ROOT.TGrid.Connect("alien://")
        if not grid:
            print("No grid connection available. Exiting.")
            sys.exit()
        files_to_download = get_files_to_download(grid, config)

        if aod:
            download_aod(files_to_download, config)
        if analysis:
            download_analysis_results(files_to_download, config)
    if full_config:
        download_full_config(files_to_download, config)
    if parquet:
        convert_aod_to_parquet(config)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Create list of files from grid and merge files.')
    parser.add_argument('config_file', type=str,
        help='Path to the JSON configuration file')
    parser.add_argument('--aod', action='store_true', default=False,
        help='Run only the AOD download and merge')
    parser.add_argument('--analysis', action='store_true', default=False,
        help='Run only the analysis results download and merge')
    parser.add_argument('--parquet', action='store_true', default=False,
        help='Run only the conversion to Parquet')
    parser.add_argument('--config', action='store_true', default=False,
        help='Download only the full_config.json')
    args = parser.parse_args()


    with open(args.config_file, encoding="utf8") as cfg_file:
        cfg = json.load(cfg_file)

    download_files_from_grid(cfg, aod=args.aod, analysis=args.analysis, parquet=args.parquet, full_config=args.config)
