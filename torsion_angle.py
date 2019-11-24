import os
import numpy as np
from MDAnalysis.analysis.dihedrals import Dihedral
from trajectory import Trajectory
from torsion_angle_plot import TorsionAnglePlot


class TorsionAngle(TorsionAnglePlot, Trajectory):
    def __init__(self, env, mda_universe, torsion_angles_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.torsion_angles_dir = torsion_angles_dir

    def _set_torsion_angles_dir(self, torsion_name):
        torsion_name_dir = os.path.join(self.torsion_angles_dir, torsion_name)

        if not os.path.exists(torsion_name_dir):
            os.mkdir(torsion_name_dir)

        return torsion_name_dir

    def torsion_trajectory_analysis(self):
        torsion_stats = {}

        for torsion_name, torsion_values in self.env["torsions"]["vars"].items():
            torsion_stats[torsion_name] = {}

            for torsion_type, torsion_selection in torsion_values.items():
                torsion_stats[torsion_name][torsion_type] = {}

                try:
                    atom_indexes = torsion_selection.split()
                    mda_atom_selection = self.mda_universe.select_atoms(
                        "index " + atom_indexes[0],
                        "index " + atom_indexes[1],
                        "index " + atom_indexes[2],
                        "index " + atom_indexes[3],
                    )
                except:
                    mda_atom_selection = self.mda_universe.select_atoms(
                        torsion_selection
                    )

                mda_dihedral = Dihedral([mda_atom_selection])

                torsion_angles = self.atom_selection_torsions_for_trajectory(
                    mda_dihedral
                )

                torsion_name_dir = self._set_torsion_angles_dir(torsion_name)

                torsion_type_file_path = os.path.join(
                    torsion_name_dir, torsion_type + ".dat"
                )

                self._write_trajectory_atom_selection_torsions_to_file(
                    torsion_angles, torsion_type_file_path
                )

                torsion_stats[torsion_name][torsion_type] = self.stats_analysis(
                    torsion_angles
                )

                torsion_stats[torsion_name][torsion_type][
                    "torsion_angles"
                ] = torsion_angles

                time_series = np.arange(
                    0, self.get_trajectory_time_in_ns(), self.ns_per_frame()
                )

                self.time_series_scatter(
                    time_series, torsion_angles, torsion_name_dir, torsion_type
                )

                self.probability_histogram(
                    torsion_angles, torsion_name_dir, torsion_type
                )

            torsion_stats_file_path = os.path.join(
                torsion_name_dir, torsion_name + "_torsion_angle_stats.txt"
            )

            self._write_torsion_stats(torsion_stats, torsion_stats_file_path)

            self._create_additional_torsional_plots(
                torsion_name, torsion_name_dir, torsion_stats
            )

    def atom_selection_torsions_for_trajectory(self, mda_dihedral):
        return mda_dihedral.run().angles

    def _write_trajectory_atom_selection_torsions_to_file(
        self, torsion_angles, output_file_name
    ):
        with open(output_file_name, "w") as outf:
            for frame_num, torsion_angle in enumerate(torsion_angles):
                outf.write(str(frame_num))
                outf.write(" ")
                outf.write(str(torsion_angle[0]))
                outf.write("\n")

    def stats_analysis(self, torsion_angles):

        torsion_angles_from_start_frame = torsion_angles[self.get_start_frame() :]

        mean_average = np.mean(torsion_angles_from_start_frame)
        median_average = np.median(torsion_angles_from_start_frame)
        std_dev = np.std(torsion_angles_from_start_frame)
        max_torsion = max(torsion_angles_from_start_frame)
        min_torsion = min(torsion_angles_from_start_frame)

        torsion_stats = {
            "mean_average": mean_average,
            "median_average": median_average,
            "std_dev": std_dev,
            "max_torsion": max_torsion,
            "min_torsion": min_torsion,
        }

        return torsion_stats

    def _write_torsion_stats(self, torsion_stats, output_file_name):

        with open(output_file_name, "w") as outf:
            outf.write("TORSION STATS\n")
            outf.write("---------------\n")
            for torsion_name, torsion_values in torsion_stats.items():
                for torsion_type, torsion_stats_values in torsion_values.items():
                    outf.write(torsion_name)
                    outf.write(": ")
                    outf.write(torsion_type)
                    outf.write("\n")
                    for stat_name, stat_value in torsion_stats_values.items():
                        outf.write(stat_name)
                        outf.write(" ")
                        outf.write(str(stat_value))
                        outf.write("\n")
                outf.write("\n\n")

    def _create_additional_torsional_plots(
        self, torsion_name, torsion_name_dir, torsion_stats
    ):
        for torsion_plot_name, plot_values in self.env["torsions"]["plots"].items():
            if torsion_plot_name == torsion_name:
                x_key = plot_values["x_key"]
                y_key = plot_values["y_key"]

                x_series = torsion_stats[torsion_plot_name][x_key]["torsion_angles"]
                y_series = torsion_stats[torsion_plot_name][y_key]["torsion_angles"]
                self.torsion_angles_scatter(
                    x_series, y_series, x_key, y_key, torsion_name_dir
                )
