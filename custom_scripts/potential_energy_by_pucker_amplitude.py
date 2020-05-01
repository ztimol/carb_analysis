import os
import numpy as np
from collections import OrderedDict
from plot.plot import Plot


class PotentialEnergyByPuckerAmplitude:

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

    def read_cp_phi_and_cp_theta_binned_energies(self):
        data_file = "custom_scripts/output/cp_phi_and_cp_theta_binned_energies.dat"
        cp_phi_values = []
        cp_theta_values = []
        potential_energy_values = []
        with open(data_file, "r") as infile:
            for line in infile:
                line = line.split()
                cp_phi = eval(line[0])
                cp_theta = eval(line[1])
                pe = eval(line[2])
                cp_phi_values.append(cp_phi)
                cp_theta_values.append(cp_theta)
                potential_energy_values.append(pe)

        return cp_phi_values, cp_theta_values, potential_energy_values

    def plot_energy_against_cp_theta(
        self, cp_trajectory_pucker_params, trajectory_energy_per_frame
    ):

        trajectory_cp_theta = [
            cp_params["cp_theta"] for cp_params in cp_trajectory_pucker_params.values()
        ]

        scatter_params = {
            "y_label": "PE (kcal/mol)",
            "x_label": "CP theta (deg)",
            "y_start": 260,
            "y_end": 350,
            "x_major_tick": 30,
            "x_end": 180,
            "color": "black",
        }

        plot_file_path = "custom_scripts/output/pe_against_cp_theta.png"

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
