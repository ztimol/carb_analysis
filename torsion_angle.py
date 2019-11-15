import numpy as np
from MDAnalysis.analysis.dihedrals import Dihedral


class TorsionAngle:
    def __init__(self, atom_selection):
        self.torsion = Dihedral(atom_selection)

    def calculate_atom_selection_torsions_for_trajectory(self):
        atom_selection_torsions_for_trajectory = self.torsion.run()
        return atom_selection_torsions_for_trajectory.angles

    def write_trajectory_atom_selection_torsions_to_file(
        self, torsion_angles, output_file_name
    ):
        with open(output_file_name, "w") as outf:
            for frame_num, torsion_angle in enumerate(torsion_angles):
                outf.write(str(frame_num))
                outf.write(" ")
                outf.write(str(torsion_angle[0]))
                outf.write("\n")

    def stats_analysis(self, torsion_name, torsion_type, torsion_angles):
        mean_average = np.mean(torsion_angles)
        median_average = np.median(torsion_angles)
        std_dev = np.std(torsion_angles)
        max_torsion = max(torsion_angles)
        min_torsion = min(torsion_angles)

        stats_data = {
            torsion_name: {
                torsion_type: {
                    "mean_average": mean_average,
                    "median_average": median_average,
                    "std_dev": std_dev,
                    "max_torsion": max_torsion,
                    "min_torsion": min_torsion,
                }
            }
        }

        return stats_data

    def write_stats_data(self, torsion_stats, output_file_name):
        print(torsion_stats)
        with open(output_file_name, "a") as outf:
            outf.write("TORSION STATS\n")
            outf.write("---------------\n")
            for torsion_name, torsion_values in torsion_stats.items():
                outf.write(torsion_name)
                outf.write(":")
                for torsion_type, torsion_stats_values in torsion_values.items():
                    outf.write(torsion_type)
                    outf.write("\n")
                    for stat_name, stat_value in torsion_stats_values.items():
                        outf.write(stat_name)
                        outf.write(" ")
                        outf.write(str(stat_value))
                        outf.write("\n")
                outf.write("\n")
