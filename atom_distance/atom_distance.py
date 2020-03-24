# import numpy as np
from trajectory import Trajectory

# from torsion_angle_plot import TorsionAnglePlot


class AtomDistance(Trajectory):
    def __init__(self, env, mda_universe, atom_distances_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.atom_distances_dir = atom_distances_dir

    def atom_distance_trajectory_analysis(self):
        for atom_distance_name, atom_distance_values in self.env["torsions"][
            "vars"
        ].items():
            torsion_stats[torsion_name] = {}
