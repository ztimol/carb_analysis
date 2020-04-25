import os
import numpy as np
from plot.plot import Plot


class EnergyByPuckerAmplitude:

    # get cremer-pople puckering parameters from file
    def read_trajectory_cp_puckering_parameters(self, cp_trajectory_parameter_file):

        cp_trajectory_pucker_params = {}

        with open(cp_trajectory_parameter_file, "r") as infile:
            for line in infile:
                line = line.split()
                frame_number = eval(line[0])
                cp_phi = eval(line[1])
                cp_theta = eval(line[2])
                cp_q = eval(line[3])
                cp_trajectory_pucker_params[frame_number] = {
                    "cp_phi": cp_phi,
                    "cp_theta": cp_theta,
                    "cp_q": cp_q,
                }

        return cp_trajectory_pucker_params

    def read_trajectory_potential_energy(self, trajectory_pe_file):

        trajectory_energy_per_frame = []

        with open(trajectory_pe_file, "r") as infile:
            for line in infile:
                line = line.split()
                potential_energy = eval(line[1])
                trajectory_energy_per_frame.append(potential_energy)

        return trajectory_energy_per_frame

    def seperate_energy_by_cp_theta(
        self, cp_trajectory_pucker_params, trajectory_energy_per_frame
    ):

        pe_for_4c1_chair = []
        pe_for_1c4_chair = []
        chair_boat_twist_conformers = []
        other_pucker_conformers = []

        for frame_number, cp_params in cp_trajectory_pucker_params.items():
            cp_theta = cp_params["cp_theta"]
            potential_energy = trajectory_energy_per_frame[frame_number]
            if cp_theta < 30:
                pe_for_4c1_chair.append(potential_energy)
            elif cp_theta > 150:
                pe_for_1c4_chair.append(potential_energy)
            elif 75 < cp_theta < 105:
                chair_boat_twist_conformers.append(potential_energy)
            else:
                other_pucker_conformers.append(potential_energy)
            # else:
            #     raise Exception("Invalid theta value for frame: " + str(frame_number))

        pe_by_conformers = {
            "pe_for_4c1_chair": pe_for_4c1_chair,
            "pe_for_1c4_chair": pe_for_1c4_chair,
            "chair_boat_twist_conformers": chair_boat_twist_conformers,
            "other_pucker_conformers": other_pucker_conformers,
        }

        return pe_by_conformers

    def bin_energies_by_cp_theta(
        self, cp_trajectory_pucker_params, trajectory_energy_per_frame
    ):

        binned_energies = {}

        for frame_number, cp_params in cp_trajectory_pucker_params.items():

            cp_theta = cp_params["cp_theta"]
            potential_energy = trajectory_energy_per_frame[frame_number]

            for bin_value in range(0, 185, 5):
                if bin_value <= cp_theta < bin_value + 5:
                    try:
                        binned_energies[bin_value].append(potential_energy)
                    except:
                        binned_energies[bin_value] = []
                        binned_energies[bin_value].append(potential_energy)
                    break

        mean_energy_per_cp_theta_bin = {}

        for bin_value, energies in binned_energies.items():
            mean_energy_per_cp_theta_bin[bin_value] = np.mean(energies)

        return mean_energy_per_cp_theta_bin

    def plot_binned_energies_against_cp_theta(self, mean_energy_per_cp_theta_bin):

        scatter_params = {
            "y_label": "PE (KCal)",
            "x_label": "CP theta (deg)",
            "y_start": 290,
            "y_end": 320,
            "x_major_tick": 30,
            "x_end": 180,
            "color": "black",
        }

        plot_file_path = "generic_scripts/binned_pe_against_cp_theta.png"

        Plot().two_dimensional_scatter(
            mean_energy_per_cp_theta_bin.keys(),
            mean_energy_per_cp_theta_bin.values(),
            plot_file_path,
            scatter_params=scatter_params,
        )

    def plot_energy_against_cp_theta(
        self, cp_trajectory_pucker_params, trajectory_energy_per_frame
    ):

        trajectory_cp_theta = [
            cp_params["cp_theta"] for cp_params in cp_trajectory_pucker_params.values()
        ]

        scatter_params = {
            "y_label": "PE (KCal)",
            "x_label": "CP theta (deg)",
            "y_start": 250,
            "y_end": 350,
            "x_major_tick": 30,
            "x_end": 180,
            "color": "black",
        }

        plot_file_path = "generic_scripts/pe_against_cp_theta.png"

        Plot().two_dimensional_scatter(
            trajectory_cp_theta,
            trajectory_energy_per_frame,
            plot_file_path,
            scatter_params=scatter_params,
        )

    def calc_energy_stats_by_conformers(self, pe_by_conformers):

        for conformer_type, trajectory_energies in pe_by_conformers.items():
            print(conformer_type)
            print("number_of_frames_for_conformer", len(trajectory_energies))
            mean_average = np.mean(trajectory_energies)
            median_average = np.median(trajectory_energies)
            std_dev = np.std(trajectory_energies)
            max_energy = max(trajectory_energies)
            min_energy = min(trajectory_energies)

            print("mean", mean_average)
            print("median", median_average)
            print("std dev", std_dev)
            print("max_energy", max_energy)
            print("min_energy", min_energy)
            print()
