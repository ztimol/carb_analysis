import os
import numpy as np
from trajectory import Trajectory
from atom_distance.atom_distance_plot import AtomDistancePlot


class AtomDistance(AtomDistancePlot, Trajectory):
    def __init__(self, env, mda_universe, atom_distances_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.atom_distances_dir = atom_distances_dir

    def atom_distance_trajectory_analysis(self):

        for atom_selection in self.env["atom_distances"]:
            atom_selection_out_dir = os.path.join(
                self.atom_distances_dir, atom_selection
            )

            if not os.path.exists(atom_selection_out_dir):
                os.mkdir(atom_selection_out_dir)

            data_out_file_path = os.path.join(
                atom_selection_out_dir, atom_selection + ".dat"
            )

            trajectory_atom_pair_distances = self._get_atom_pair_distances_for_trajectory(
                atom_selection
            )

            self._write_trajectory_atom_pair_distances(
                trajectory_atom_pair_distances, data_out_file_path
            )

            time_series = np.arange(
                0, self.get_trajectory_time_in_ns(), self.ns_per_frame()
            )

            self.time_series_scatter(
                time_series,
                trajectory_atom_pair_distances,
                atom_selection_out_dir,
                atom_selection,
            )

            self.probability_histogram(
                trajectory_atom_pair_distances, atom_selection_out_dir, atom_selection
            )

    def _get_atom_pair_distances_for_trajectory(self, atom_selection):
        first_atom = self.mda_universe.select_atoms(atom_selection)[0]
        second_atom = self.mda_universe.select_atoms(atom_selection)[-1]

        trajectory_atom_pair_distances = []
        for frame_number in self.mda_universe.trajectory:  # iterate through all frames
            r = (
                first_atom.position - second_atom.position
            )  # end-to-end vector from atom positions
            d = np.linalg.norm(r)  # end-to-end distance
            trajectory_atom_pair_distances.append(d)

        return trajectory_atom_pair_distances

    def _write_trajectory_atom_pair_distances(
        self, trajectory_atom_pair_distances, data_out_file_path
    ):

        with open(data_out_file_path, "w") as out_file:
            frame_number = 0
            for atom_pair_distance in trajectory_atom_pair_distances:
                out_file.write(str(frame_number))
                out_file.write(" ")
                out_file.write(str(atom_pair_distance))
                out_file.write("\n")
                frame_number += 1
