"""
Script to run full analysis for D resonances
"""

import argparse
import os

def run_full_analysis(config,
                      an_res_file,
                      centrality,
                      outputdir,
                      suffix,
                      vn_method,
                      wagon_id,
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
    - outputdir (str): output directory
    - suffix (str): suffix for output files
    - vn_method (str): vn technique (sp, ep, deltaphi)
    - wagon_id (str): wagon ID
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
    if vn_method not in ["sp", "ep", "deltaphi"]:
        print("\033[91m Invalid vn method. Only sp, ep, deltaphi implemented. Exit!\033[0m")
        return
    outputdir = f"{outputdir}/{vn_method}"
    vn_method_withopt = f" -vn {vn_method}"
    cent_withopt = f" -c {centrality}"
    wagon_id_withopt = f" -w {wagon_id}"

    if skip_resolution and skip_projection and skip_rawyield:
        print("\033[91m Nothing to do, all steps are skipped\033[0m")
        return

    if not os.path.exists(outputdir):
        print(f"\033[92m Creating output directory {outputdir}\033[0m")
        os.makedirs(outputdir)

    if not skip_resolution:
        # resolution extraction
        if not os.path.exists(f"{outputdir}/resolution"):
            os.makedirs(f"{outputdir}/resolution")
        outputdir_reso = f"-o {outputdir}/resolution/"
        command_reso = f"python3 compute_reso.py {an_res_file} {suffix_withopt} {outputdir_reso} {cent_withopt} {vn_method_withopt}"
        if wagon_id != "":
            command_reso += f" {wagon_id_withopt}"
        print("\n\033[92m Starting resolution extraction\033[0m")
        print(f"\033[92m {command_reso}\033[0m")
        os.system(command_reso)

    if not skip_projection:
        # projection
        if not os.path.exists(f"{outputdir}/proj"):
            os.makedirs(f"{outputdir}/proj")
        outputdir_proj = f"-o {outputdir}/proj"
        command_proj = f"python3 project_thnsparse.py {config} {an_res_file} {cent_withopt} {suffix_withopt} {outputdir_proj} {vn_method_withopt}"
        if wagon_id != "":
            command_proj += f" {wagon_id_withopt}"
        print("\n\033[92m Starting projection\033[0m")
        print(f"\033[92m {command_proj}\033[0m")
        os.system(command_proj)

    if not skip_rawyield:
        # raw yield
        if not os.path.exists(f"{outputdir}/ry"):
            os.makedirs(f"{outputdir}/ry")
        outputdir_rawyield = f"-o {outputdir}/ry"
        proj_file = f"{outputdir}/proj/"
        proj_file += f"proj{suffix}.root"
        command_vn = f"python3 get_vn_vs_mass.py {config} {centrality} {proj_file} {outputdir_rawyield} {suffix_withopt} {vn_method_withopt}"
        print("\n\033[92m Starting vn extraction\033[0m")
        print(f"\033[92m {command_vn}\033[0m")
        os.system(command_vn)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("config", metavar="text",
                        default="config.yaml", help="configuration file")
    parser.add_argument("an_res_file", metavar="text",
                        default="an_res.root", help="input ROOT file with anres")
    parser.add_argument("--centrality", "-c", metavar="text",
                        default="k3050", help="centrality class")
    parser.add_argument("--outputdir", "-o", metavar="text",
                        default=".", help="output directory")
    parser.add_argument("--suffix", "-s", metavar="text",
                        default="", help="suffix for output files")
    parser.add_argument("--vn_method", "-vn", metavar="text",
                        default="sp", help="vn technique (sp, ep, deltaphi)")
    parser.add_argument("--wagon_id", "-w", metavar="text",
                        default="", help="wagon ID", required=False)
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
        args.outputdir,
        args.suffix,
        args.vn_method,
        args.wagon_id,
        args.skip_resolution,
        args.skip_projection,
        args.skip_rawyield
    )
