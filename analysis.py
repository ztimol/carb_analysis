import os
import MDAnalysis as mda
from trajectory import Trajectory
from torsion.torsion_angle import TorsionAngle
from atom_distance.atom_distance import AtomDistance


class Analysis(Trajectory):
    def __init__(self, env):
        self.env = env
        self.mda_universe = mda.Universe(
            self.get_psf_file_path(), self.get_dcd_file_path()
        )
        super().__init__(env)

    def torsion_analysis(self):
        output_dir = self.env.get("output_dir", "output")

        try:
            if self.env["torsions"]:
                torsion_angles_dir = os.path.join(output_dir, "torsion_angles")
                if not os.path.exists(torsion_angles_dir):
                    os.mkdir(torsion_angles_dir)
        except KeyError:
            print("no torsions specified in config file. Don't run torsion analysis.")
            return

        torsion = TorsionAngle(self.env, self.mda_universe, torsion_angles_dir)
        torsion.torsion_trajectory_analysis()

    def rmsd_analysis(self, env):
        output_dir = env["input_params"].get("output_dir", "output")

        try:
            if self.env["atom_rmsd"]:
                atom_rmsds_dir = os.path.join(output_dir, "atom_rmsds")
                if not os.path.exists(atom_rmsds_dir):
                    os.mkdir(atom_rmsds_dir)
        except KeyError:
            print(
                "no atom rmsds specified in config file. Don't run atom rmsd analysis."
            )
            return

        atom_rmsd = AtomRmsd(self.env, self.mda_universe, atom_rmsds_dir)
        atom_rmsd.atom_rmsd_trajectory_analysis()

    def distance_analysis(self, env):
        output_dir = env["input_params"].get("output_dir", "output")

        try:
            if self.env["atom_distance"]:
                atom_distances_dir = os.path.join(output_dir, "atom_distances")
                if not os.path.exists(atom_distances_dir):
                    os.mkdir(atom_distances_dir)
        except KeyError:
            print(
                "no atom distances specified in config file. Don't run atom distance analysis."
            )
            return

        atom_distance = AtomDistance(self.env, self.mda_universe, atom_distances_dir)
        atom_distance.atom_distance_trajectory_analysis()
