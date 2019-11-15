import numpy as np
from MDAnalysis.analysis.dihedrals import Dihedral


class TorsionAngle:
    def __init__(self, atom_selection, start_frame):
        self.torsion = Dihedral(atom_selection)
        self.start_frame = start_frame

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

    def stats_analysis(self, torsion_angles):

        torsion_angles_from_start_frame = torsion_angles[self.start_frame :]

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

    def write_stats_data(self, torsion_stats, output_file_name):
        print(torsion_stats)
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
