import sys, os
import argparse, shutil
import MDAnalysis as mda

from config import Config
from analysis import Analysis

from helper import get_output_dir_name_from_config_file_args


parser = argparse.ArgumentParser(description="MD Analysis")

parser.add_argument("-f", help="Select a config file", type=str)

args = parser.parse_args()


def prepare_analysis():

    env = {}
    env["output_params"] = {}

    analysis_config = Config()
    env["input_params"] = analysis_config.get_analysis_config(args.f)

    config_file_name = get_output_dir_name_from_config_file_args(args.f)

    output_dir = env["input_params"].get("output_dir", "output/" + config_file_name)
    env["input_params"]["output_dir"] = output_dir

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    shutil.copy(args.f, output_dir)  # backup config file in output directory

    return env


def finalise_analysis():
    print("Completed analysis.")


def main():

    env = prepare_analysis()

    print("Starting analysis for: " + args.f + " ...\n")

    structure_analysis = Analysis(env["input_params"])

    structure_analysis.torsion_analysis()
    structure_analysis.namd_energy_analysis()
    structure_analysis.ring_pucker_analysis()

    finalise_analysis()


if __name__ == "__main__":
    main()
