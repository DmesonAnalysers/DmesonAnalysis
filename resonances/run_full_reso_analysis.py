"""
Script to run full analysis for D resonances
"""

import argparse
import os

def run_full_analysis(dir_config_proj, config_fit, trigger, pdg_d, pdg_v0):
    """
    function for full analysis

    Parameters
    ----------

    - dir_config_proj (str): path of directory with config files for projection step
    - config_fit (str): path config file for fit step
    - trigger (str): trigger class (MB or HM)
    - pdg_d (int): PDG code of the D meson
    - pdg_v0 (int): PDG code of the V0
    """

    for cfg_file in os.listdir(dir_config_proj):
        if "config_proj" not in cfg_file:
            continue

        # get all parameters needed
        cfg_file = os.path.join(dir_config_proj, cfg_file)
        suffix = cfg_file.split("config_proj_reso")[-1].replace(".yml", "")
        if suffix != "":
            suffix = f" -s {suffix}"

        if pdg_d == 411:
            name_d = "Dplus"
        elif pdg_d == 413:
            name_d = "Dstar"

        if pdg_v0 == 310:
            name_v0 = "K0S"
        elif pdg_v0 == 3122:
            name_v0 = "Lambda"

        # projection of data
        print("\n\033[92m Starting data projection\033[0m")
        command_proj = "python3 project_DV0reso.py"
        args_proj = f"{cfg_file} -t {trigger} -d {pdg_d} -v0 {pdg_v0} -o {dir_config_proj} {suffix}"
        os.system(f"{command_proj} {args_proj}")

        # raw yield extraction
        print("\n\033[92m Starting raw yield extraction\033[0m")
        input_4fit = os.path.join(dir_config_proj, f"{name_d}_{name_v0}_{trigger}{suffix}.parquet.gzip")
        command_fit = "python3 extract_rawyield_DV0reso.py"
        args_fit = f"{input_4fit} {config_fit} -o {dir_config_proj} {suffix}"
        os.system(f"{command_fit} {args_fit}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Arguments")
    parser.add_argument("dir_config_proj", metavar="text",
                        default=".", help="directory with config files for projection step")
    parser.add_argument("config_fit", metavar="text",
                        default="config_fit.yml", help="config file for fit")
    parser.add_argument("--trigger", "-t", metavar="text",
                        default="HM", required=True, help="trigger class (options: [HM, MB]")
    parser.add_argument("--pdg_D", "-d", type=int, required=True, default=413,
                        help="pdg code of the D meson (options: [411, 413])")
    parser.add_argument("--pdg_V0", "-v0", type=int, required=True, default=310,
                        help="pdg code of the V0 (options: [310, 3122])")
    args = parser.parse_args()

    run_full_analysis(args.dir_config_proj, args.config_fit, args.trigger, args.pdg_D, args.pdg_V0)
