import sys, os
import argparse
import MDAnalysis as mda

from torsion_angle import TorsionAngle
from pmf_multi import scatter_without_pmf_contour
import scatter

from config import get_analysis_config

parser = argparse.ArgumentParser(description="Calculate NOE from an md simulation")

parser.add_argument("-f", help="Select a config file", type=str)

args = parser.parse_args()


class Base:
    def do_torsion_analysis(self, env):

        base_dir = os.getcwd()

        dcd_file_path = env["input_params"].get("dcd_file", None)
        psf_file_path = env["input_params"].get("psf_file", None)
        start_frame = int(env["input_params"].get("start_frame", 0))

        output_dir = os.path.join(base_dir, env["input_params"].get("output_dir", None))

        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

        torsion_angles_dir = os.path.join(output_dir, "torsion_angles")

        if not os.path.exists(torsion_angles_dir):
            os.mkdir(torsion_angles_dir)

        u = mda.Universe(psf_file_path, dcd_file_path)
        torsion_stats = {}
        torsion_params = {}
        for torsion_name, torsion_values in env["input_params"]["torsions"].items():
            if torsion_name == "plots":
                continue
            torsion_params[torsion_name] = {}
            torsion_stats[torsion_name] = {}
            for torsion_type, torsion_selection in torsion_values.items():
                torsion_params[torsion_name][torsion_type] = {}
                torsion_stats[torsion_name][torsion_type] = {}
                mda_atom_selection = u.select_atoms(torsion_selection)
                # coordinates_for_frame = s.positions
                torsion = TorsionAngle([mda_atom_selection], start_frame)

                torsion_name_dir = os.path.join(torsion_angles_dir, torsion_name)

                if not os.path.exists(torsion_name_dir):
                    os.mkdir(torsion_name_dir)

                torsion_type_file_path = os.path.join(
                    torsion_name_dir, torsion_type + ".dat"
                )

                torsion_angles = (
                    torsion.calculate_atom_selection_torsions_for_trajectory()
                )

                torsion_params[torsion_name][torsion_type][
                    "trajectory_torsion_angles"
                ] = torsion_angles

                torsion.write_trajectory_atom_selection_torsions_to_file(
                    torsion_angles, torsion_type_file_path
                )

                torsion_params[torsion_name][torsion_type][
                    "trajectory_torsion_angle_stats"
                ] = torsion.stats_analysis(torsion_angles)

                torsion_stats[torsion_name][torsion_type] = torsion_params[
                    torsion_name
                ][torsion_type]["trajectory_torsion_angle_stats"]

            torsion_stats_file_path = os.path.join(
                torsion_name_dir, torsion_name + "_torsion_angle_stats.txt",
            )

            torsion.write_stats_data(torsion_stats, torsion_stats_file_path)

        self.create_torsion_plots(env, torsion_params)

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
            import pdb

            pdb.set_trace()
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
                phi_list, psi_list, plot_file_path, x_key, y_key,
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

    base_analysis = Base()
    base_analysis.do_torsion_analysis(env)


if __name__ == "__main__":
    main()
