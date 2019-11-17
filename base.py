import sys, os
import argparse
import MDAnalysis as mda

from torsion_angle import TorsionAngle
from pmf_multi import scatter_without_pmf_contour
from helper import confirm_critical_file_exists
from config import get_analysis_config
import scatter

parser = argparse.ArgumentParser(description="MD Analysis")

parser.add_argument("-f", help="Select a config file", type=str)

args = parser.parse_args()


class Base:
    def __init__(self, dcd_file_path, psf_file_path):
        self.dcd_file_path = dcd_file_path
        self.psf_file_path = psf_file_path

        confirm_critical_file_exists(dcd_file_path)
        confirm_critical_file_exists(psf_file_path)

        self.mda_universe = mda.Universe(psf_file_path, dcd_file_path)

    def torsion_analysis(self, env):

        start_frame = int(env["input_params"].get("start_frame", 0))
        output_dir = env["input_params"].get("output_dir", "output")

        try:
            input_torsion_params = env["input_params"]["torsions"]

            torsion_angles_dir = os.path.join(output_dir, "torsion_angles")
            if not os.path.exists(torsion_angles_dir):
                os.mkdir(torsion_angles_dir)
        except KeyError:
            return  # no torsions specified in config file. Don't run torsion analysis.

        torsion = TorsionAngle(torsion_angles_dir)
        torsion.torsion_trajectory_analysis(input_torsion_params, start_frame)

        import pdb

        pdb.set_trace()


def main():
    env = {}
    env["output_params"] = {}

    env["input_params"] = get_analysis_config(args.f)

    dcd_file = env["input_params"].get("dcd_file", None)
    psf_file = env["input_params"].get("psf_file", None)

    output_dir = env["input_params"].get("output_dir", "output")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    # we check for dcd and psf file existence here as opposed to a try,
    # except later in the script so that we catch the error earlier

    base_analysis = Base(dcd_file, psf_file)
    base_analysis.torsion_analysis(env)


if __name__ == "__main__":
    main()
