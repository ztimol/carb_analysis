import sys, os
import argparse
import MDAnalysis as mda

from torsion_angle import TorsionAngle

from config import get_analysis_config

parser = argparse.ArgumentParser(description="Calculate NOE from an md simulation")

parser.add_argument("-f", help="Select a config file", type=str)

args = parser.parse_args()


def main():
    env = {}
    env = {"results": {}}

    env["input_params"] = get_analysis_config(args.f)
    do_analysis(env)


def do_analysis(env):
    dcd_file = "S_flexneri_7a_6RU_0-32ns_every100frms.dcd"
    u = mda.Universe("S_flexneri_7a_6RU.psf", dcd_file)
    base_dir = os.getcwd()
    output_dir = os.path.join(base_dir, "output")

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)

    torsion_angles_dir = os.path.join(output_dir, "torsion_angles")

    if not os.path.exists(torsion_angles_dir):
        os.mkdir(torsion_angles_dir)

    for torsion_name, torsion_values in env["input_params"]["torsions"].items():
        for torsion_type, torsion_selection in torsion_values.items():
            mda_atom_selection = u.select_atoms(torsion_selection)
            # coordinates_for_frame = s.positions
            torsion = TorsionAngle([mda_atom_selection])

            torsion_name_dir = os.path.join(torsion_angles_dir, torsion_name)

            if not os.path.exists(torsion_name_dir):
                os.mkdir(torsion_name_dir)

            torsion_type_file_path = os.path.join(
                torsion_name_dir, torsion_type + ".dat"
            )

            torsion_angles = torsion.calculate_atom_selection_torsions_for_trajectory()
            torsion.write_trajectory_atom_selection_torsions_to_file(
                torsion_angles, torsion_type_file_path
            )

            torsion_stats = torsion.stats_analysis(
                torsion_name, torsion_type, torsion_angles
            )

            torsion_stats_file_path = os.path.join(
                torsion_name_dir, torsion_name + "_torsion_angle_stats.txt"
            )

            torsion.write_stats_data(torsion_stats, torsion_stats_file_path)


if __name__ == "__main__":
    main()
