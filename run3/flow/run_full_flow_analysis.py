"""
Script to run full analysis for D resonances
"""

import argparse
import os

def run_full_analysis(config,
                      an_res_file,
                      centrality,
                      suffix,
                      doEP,
                      skip_resolution,
                      skip_projection,
                      skip_rawyield
                      ):
    """
    function for full analysis

    Parameters
    ----------

    - config (str): path of directory with config files
    - an_res_file (str): path of directory with analysis results
    - centrality (str): centrality class
    - suffix (str): suffix for output files
    - doEP (bool): do EP resolution
    - skip_resolution (bool): skip resolution extraction
    - skip_projection (bool): skip projection extraction
    - skip_rawyield (bool): skip raw yield extraction
    """

    # get all parameters needed
    if suffix != "":
        suffix_withopt = f" -s {suffix}"
    else:
        suffix = an_res_file.split("AnalysisResults")[-1].replace(".root", "")
        suffix_withopt = f" -s {suffix}"
    if doEP:
        outputdir = "~/flowD/ep"
    else:
        outputdir = "~/flowD/sp"


    if not skip_resolution:
        # resolution extraction
        outputdir_reso = f"-o {outputdir}/resolution/"
        command_reso = f"python3 compute_reso.py {an_res_file} {suffix_withopt} {outputdir_reso}"
        if doEP:
            command_reso += " --doEP"
        print("\n\033[92m Starting resolution extraction\033[0m")
        print(f"\033[92m {command_reso}\033[0m")
        os.system(command_reso)

    if not skip_projection:
        # projection
        outputdir_proj = f"-o {outputdir}/proj"
        command_proj = f"python3 project_thnsparse.py {config} {an_res_file} {suffix_withopt} {outputdir_proj}"
        if doEP:
            command_proj += " --doEP"
        print("\n\033[92m Starting projection\033[0m")
        print(f"\033[92m {command_proj}\033[0m")
        os.system(command_proj)

    if not skip_rawyield:
        # raw yield
        outputdir_rawyield = f"-o {outputdir}/ry"
        proj_file = f"{outputdir}/proj/"
        if doEP:
            proj_file += f"proj{suffix}.root"
            command_sim_fit = f"python3 get_vn_vs_mass.py {config} {centrality} {proj_file} {outputdir_rawyield} {suffix_withopt} --doEP"
        else:
            proj_file += f"proj{suffix}.root"
            command_sim_fit = f"python3 get_vn_vs_mass.py {config} {centrality} {proj_file} {outputdir_rawyield} {suffix_withopt}"
        print("\n\033[92m Starting simultaneous fit\033[0m")
        print(f"\033[92m {command_sim_fit}\033[0m")
        os.system(command_sim_fit)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("an_res_file", metavar="text",
                        default="an_res.root", help="input ROOT file with anres")
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k3050", help="centrality class")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--doEP", action="store_true", default=False,
                        help="do EP resolution")
    parser.add_argument("--skip_resolution", action="store_true", default=False,
                        help="skip resolution extraction")
    parser.add_argument("--skip_projection", action="store_true", default=False,
                        help="skip projection extraction")
    parser.add_argument("--skip_rawyield", action="store_true", default=False,
                        help="skip raw yield extraction")
    args = parser.parse_args()

    run_full_analysis(
        args.config,
        args.an_res_file,
        args.centrality,
        args.suffix,
        args.doEP,
        args.skip_resolution,
        args.skip_projection,
        args.skip_rawyield
    )
