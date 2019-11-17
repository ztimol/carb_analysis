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

    def do_torsion_analysis(self, env):

        start_frame = int(env["input_params"].get("start_frame", 0))

        output_dir = env["input_params"].get("output_dir", "output")
        torsion_angles_dir = os.path.join(output_dir, "torsion_angles")

        if not os.path.exists(torsion_angles_dir):
            os.mkdir(torsion_angles_dir)

        input_torsion_params = env["input_params"]["torsions"]

        torsion = TorsionAngle(self.mda_universe, torsion_angles_dir)
        torsion.torsion_trajectory_analysis(input_torsion_params, start_frame)

        # self.create_torsion_plots(env, torsion_params)

    def create_torsion_plots(self, env, torsion_params):

        output_dir = os.path.join(
            os.getcwd(), env["input_params"].get("output_dir", None)
        )

        try:
            torsions_to_plot = env["input_params"]["torsions"]["plots"]
        except KeyError:
            return

        for torsion_name, torsion_values in torsions_to_plot.items():

            x_key = torsion_values.get("x_key", None)
            y_key = torsion_values.get("y_key", None)

            plot_file_path = os.path.join(
                output_dir, "torsion_angles", torsion_name, x_key + "_" + y_key
            )

            x_values = (
                torsion_params.get(torsion_name, {})
                .get(x_key, None)
                .get("trajectory_torsion_angles", [])
            )
            y_values = (
                torsion_params.get(torsion_name, {})
                .get(y_key, None)
                .get("trajectory_torsion_angles", [])
            )

            phi_list = [i[0] for i in x_values]
            psi_list = [i[0] for i in y_values]
            scatter_without_pmf_contour(
                phi_list, psi_list, plot_file_path, x_key, y_key
            )

            scatter.make_scatter(
                x_key,
                phi_list,
                [t for t in range(len(phi_list))],
                x_key,
                os.path.join(output_dir, "torsion_angles", torsion_name),
            )

            scatter.make_scatter(
                y_key,
                psi_list,
                [t for t in range(len(psi_list))],
                x_key,
                os.path.join(output_dir, "torsion_angles", torsion_name),
            )


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
    base_analysis.do_torsion_analysis(env)


if __name__ == "__main__":
    main()
