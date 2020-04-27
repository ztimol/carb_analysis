import os
import math
import numpy as np
from collections import OrderedDict
from plot.plot import Plot
from trajectory import Trajectory


class CpPuckerThetaFreeEnergy(Trajectory):
    def __init__(self):
        self.boltzmann_constant = 1.38064852e-23
        self.temp_kelvin = 300
        self.avogradros_constant = 6.02214086e23
        self.joues_per_kilocal = 4184
        self.cp_trajectory_pucker_params = {}
        # self.cp_phi_and_cp_theta_bin_values = OrderedDict()
        self.cp_theta_bin_values = OrderedDict()
        self.theta_with_max_count = {}
        self.theta_with_min_count = {}

    # get cremer-pople puckering parameters from file
    def read_trjectory_cp_puckering_parameters(self, cp_trajectory_parameter_file):

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

        self.cp_trajectory_pucker_params = cp_trajectory_pucker_params

    def count_cp_theta_by_bins(self):

        cp_theta_bin_values = OrderedDict()

        for theta_bin_value in range(0, 182, 2):
            cp_theta_bin_values[theta_bin_value] = {}
            cp_theta_bin_values[theta_bin_value]["count"] = 0

        for frame_number, cp_params in self.cp_trajectory_pucker_params.items():
            cp_theta = cp_params["cp_theta"]
            for theta_bin_value in range(0, 182, 2):
                if theta_bin_value <= cp_theta < theta_bin_value + 2:
                    cp_theta_bin_values[theta_bin_value]["count"] += 1

        self.cp_theta_bin_values = cp_theta_bin_values

    def _find_theta_bins_with_max_and_min_count(self):
        theta_with_max_count = {"count": 0, "theta_bin": None}

        for theta_bin, theta_bin_values in self.cp_theta_bin_values.items():
            count = theta_bin_values["count"]
            if count > theta_with_max_count["count"]:
                theta_with_max_count["count"] = count
                theta_with_max_count["theta_bin"] = theta_bin

        theta_with_min_count = {
            "count": theta_with_max_count["count"],
            "theta_bin": None,
        }

        for theta_bin, theta_bin_values in self.cp_theta_bin_values.items():
            count = theta_bin_values["count"]
            if count < theta_with_min_count["count"] and count != 0:
                theta_with_min_count["count"] = count
                theta_with_min_count["theta_bin"] = theta_bin

        self.theta_with_max_count = theta_with_max_count
        self.theta_with_min_count = theta_with_min_count

    def calc_free_energy_for_bins(self):

        cp_theta_bin_values = self.cp_theta_bin_values
        self._find_theta_bins_with_max_and_min_count()

        for theta_bin, bin_values in cp_theta_bin_values.items():
            count = bin_values["count"]
            max_count = self.theta_with_max_count["count"]
            try:
                free_energy_diff = (
                    -self.boltzmann_constant
                    * self.temp_kelvin
                    * math.log(count / max_count)
                    * self.avogradros_constant
                ) / self.joues_per_kilocal
                cp_theta_bin_values[theta_bin]["free_energy_diff"] = free_energy_diff
            except ValueError:
                cp_theta_bin_values[theta_bin]["free_energy_diff"] = -1
                cp_theta_bin_values[theta_bin]["count"] = -1

        cp_theta_bin_values = self._fix_free_energy(cp_theta_bin_values)
        self.cp_theta_bin_values = cp_theta_bin_values

    def write_theta_count_free_energy(self):

        data_file = "custom_scripts/output/theta_count_free_energy.dat"

        with open(data_file, "w") as out_file:
            for cp_theta, bin_values in self.cp_theta_bin_values.items():
                out_file.write(
                    "{cp_theta} {count} {free_energy}\n".format(
                        cp_theta=cp_theta,
                        count=bin_values["count"],
                        free_energy=bin_values["free_energy_diff"],
                    )
                )

    def _fix_free_energy(self, cp_theta_bin_values):

        min_theta_bin = self.theta_with_min_count["theta_bin"]
        max_energy = self.cp_theta_bin_values[min_theta_bin]["free_energy_diff"]

        for theta_bin in cp_theta_bin_values.keys():
            free_energy_diff = cp_theta_bin_values[theta_bin]["free_energy_diff"]
            if free_energy_diff == -1:
                cp_theta_bin_values[theta_bin]["free_energy_diff"] = None  # max_energy

        return cp_theta_bin_values

    def plot_free_energy_against_cp_theta(self):

        free_energies = [
            bin_value["free_energy_diff"]
            for bin_value in self.cp_theta_bin_values.values()
        ]

        scatter_params = {
            "y_label": "delta G (kcal/mol)",
            "x_label": "CP theta (deg)",
            "y_start": 0,
            "y_end": 8,
            "x_major_tick": 30,
            "x_end": 180,
            "color": "black",
        }

        plot_file_path = "custom_scripts/output/free_energy_against_cp_theta.png"

        Plot().two_dimensional_scatter(
            self.cp_theta_bin_values.keys(),
            free_energies,
            plot_file_path,
            scatter_params=scatter_params,
        )
