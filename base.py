import sys, os
import argparse
import MDAnalysis as mda

from pmf_multi import scatter_without_pmf_contour

from config import Config
from analysis import Analysis
import scatter


parser = argparse.ArgumentParser(description="MD Analysis")

parser.add_argument("-f", help="Select a config file", type=str)

args = parser.parse_args()


def main():
    env = {}
    env["output_params"] = {}

    analysis_config = Config()
    env["input_params"] = analysis_config.get_analysis_config(args.f)

    output_dir = env["input_params"].get("output_dir", "output")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    structure_analysis = Analysis(env["input_params"])

    structure_analysis.torsion_analysis()


if __name__ == "__main__":
    main()
