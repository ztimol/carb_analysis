import sys, os
import argparse
import MDAnalysis as mda


from pmf_multi import scatter_without_pmf_contour

from config import get_analysis_config
from analysis import Analysis
import scatter


parser = argparse.ArgumentParser(description="MD Analysis")

parser.add_argument("-f", help="Select a config file", type=str)

args = parser.parse_args()


def main():
    env = {}
    env["output_params"] = {}

    env["input_params"] = get_analysis_config(args.f)

    dcd_file = env["input_params"].get("dcd_file", None)
    psf_file = env["input_params"].get("psf_file", None)

    output_dir = env["input_params"].get("output_dir", "output")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    structure_analysis = Analysis(dcd_file, psf_file)

    structure_analysis.torsion_analysis(env)


if __name__ == "__main__":
    main()
