import os
import MDAnalysis as mda
from trajectory import Trajectory
from torsion.torsion_angle import TorsionAngle
from atom_distance.atom_distance import AtomDistance
from namd_energy.namd_energy import NAMDEnergy
from ring_pucker.cp_ring_pucker import CPRingPucker


class Analysis(Trajectory):
    def __init__(self, env):
        self.env = env
        self.mda_universe = mda.Universe(
            self.get_psf_file_path(), self.get_dcd_file_path()
        )
        super().__init__(env)

        self.output_dir = self.env.get("output_dir")

    def torsion_analysis(self):

        try:
            if self.env["torsions"]:
                torsion_angles_dir = os.path.join(self.output_dir, "torsion_angles")
                if not os.path.exists(torsion_angles_dir):
                    os.mkdir(torsion_angles_dir)
        except KeyError:
            print(
                "No torsions specified in config file. Will not perform torsion analysis.\n"
            )
            return

        print("Commencing torsion angle calculations...")

        torsion = TorsionAngle(self.env, self.mda_universe, torsion_angles_dir)
        torsion.torsion_trajectory_analysis()

        print("Completed torsion angle calculations.")

    def namd_energy_analysis(self):

        try:
            if self.env["namd_energies"]:
                namd_energy_dir = os.path.join(self.output_dir, "namd_energies")
                if not os.path.exists(namd_energy_dir):
                    os.mkdir(namd_energy_dir)
        except KeyError:
            print(
                "No NAMD energy params specified in config file. Will not perform NAMD energy analysis.\n"
            )
            return

        print("Commencing namd energy calculations...")

        # TO DO: check if namd2 executable path is included as namd_path field
        namd_energy = NAMDEnergy(self.env, self.mda_universe, namd_energy_dir)
        namd_energy.namd_single_point_energy_analysis()

        print("Completed namd energy calculations.")

    def ring_pucker_analysis(self):

        try:
            if self.env["ring_puckers"]:
                ring_pucker_dir = os.path.join(self.output_dir, "ring_pucker")
                if not os.path.exists(ring_pucker_dir):
                    os.mkdir(ring_pucker_dir)
        except KeyError:
            print(
                "No ring pucker params specified in config file. Will not perform ring pucker analysis.\n"
            )
            return

        print("Commencing puckering parameter calculations...")

        ring_pucker = CPRingPucker(self.env, self.mda_universe, ring_pucker_dir)
        ring_pucker.cp_ring_pucker_analysis()

        print("Completed puckering parameter calculations.")

    def distance_analysis(self):

        try:
            if self.env["atom_distances"]:
                atom_distances_dir = os.path.join(self.output_dir, "atom_distances")
                if not os.path.exists(atom_distances_dir):
                    os.mkdir(atom_distances_dir)
        except KeyError:
            print(
                "No ring Atom pair distances specified in config file. Will not perform ring distance analysis.\n"
            )
            return

        print("Commencing atom pair distance calculations...")

        atom_distance = AtomDistance(self.env, self.mda_universe, atom_distances_dir)
        atom_distance.atom_distance_trajectory_analysis()

        print("Completed atom pair distance calculations...")
