import os
import math
import numpy as np
from collections import OrderedDict
from plot.plot import Plot
from trajectory import Trajectory


class CpPuckerPhiThetaFreeEnergy(Trajectory):
    def __init__(self):
        self.boltzmann_constant = 1.38064852e-23
        self.temp_kelvin = 300
        self.avogradros_constant = 6.02214086e23
        self.joues_per_kilocal = 4184
        self.cp_trajectory_pucker_params = OrderedDict()
        self.cp_phi_theta_bin_values = OrderedDict()
        self.phi_theta_with_max_count = {}
        self.phi_theta_with_min_count = {}

    # get cremer-pople puckering parameters from file
    def read_trjectory_cp_puckering_parameters(self, cp_trajectory_parameter_file):

        cp_trajectory_pucker_params = OrderedDict()
        import math

        with open(cp_trajectory_parameter_file, "r") as infile:
            for line in infile:
                line = line.split()
                frame_number = eval(line[0])
                cp_phi = eval(line[1])
                cp_theta = eval(line[2])
                cp_q = eval(line[3])

                # plumed
                time = eval(line[0])
                cp_phi = (eval(line[4]) * 180) / math.pi
                cp_theta = (eval(line[5]) * 180) / math.pi
                cp_q = eval(line[6])
                cp_trajectory_pucker_params[time] = {
                    "cp_phi": cp_phi,
                    "cp_theta": cp_theta,
                    "cp_q": cp_q,
                }

        self.cp_trajectory_pucker_params = cp_trajectory_pucker_params

    def count_cp_phi_theta_by_bins(self):

        cp_phi_theta_bin_values = OrderedDict()

        for phi_bin_value in range(0, 362, 2):
            cp_phi_theta_bin_values[phi_bin_value] = OrderedDict()
            for theta_bin_value in range(0, 182, 2):
                cp_phi_theta_bin_values[phi_bin_value][theta_bin_value] = {}
                cp_phi_theta_bin_values[phi_bin_value][theta_bin_value]["count"] = 0

        for frame_number, cp_params in self.cp_trajectory_pucker_params.items():
            cp_phi = cp_params["cp_phi"]
            cp_theta = cp_params["cp_theta"]
            for phi_bin_value in range(0, 362, 2):
                if phi_bin_value <= cp_phi < phi_bin_value + 2:
                    for theta_bin_value in range(0, 182, 2):
                        if theta_bin_value <= cp_theta < theta_bin_value + 2:
                            cp_phi_theta_bin_values[phi_bin_value][theta_bin_value][
                                "count"
                            ] += 1

        self.cp_phi_theta_bin_values = cp_phi_theta_bin_values

    def calc_free_energy_for_bins(self):

        cp_phi_theta_bin_values = self.cp_phi_theta_bin_values
        self._find_phi_theta_bins_with_max_and_min_count()

        for phi_bin, phi_bin_values in cp_phi_theta_bin_values.items():
            for theta_bin, phi_theta_bin_values in phi_bin_values.items():
                count = phi_theta_bin_values["count"]
                max_count = self.phi_theta_with_max_count["count"]
                try:
                    free_energy_diff = (
                        -self.boltzmann_constant
                        * self.temp_kelvin
                        * math.log(count / max_count)
                        * self.avogradros_constant
                    ) / self.joues_per_kilocal
                    cp_phi_theta_bin_values[phi_bin][theta_bin][
                        "free_energy_diff"
                    ] = free_energy_diff
                except ValueError:
                    cp_phi_theta_bin_values[phi_bin][theta_bin]["free_energy_diff"] = -1
                    cp_phi_theta_bin_values[phi_bin][theta_bin]["count"] = -1

        cp_phi_theta_bin_values = self._fix_free_energy(cp_phi_theta_bin_values)
        self.cp_phi_theta_bin_values = cp_phi_theta_bin_values

    def _find_phi_theta_bins_with_max_and_min_count(self):
        phi_theta_with_max_count = {"count": 0, "theta_bin": None, "phi_bin": None}

        for phi_bin, phi_bin_values in self.cp_phi_theta_bin_values.items():
            for theta_bin, phi_theta_bin_values in phi_bin_values.items():
                count = phi_theta_bin_values["count"]
                if count > phi_theta_with_max_count["count"]:
                    phi_theta_with_max_count["count"] = count
                    phi_theta_with_max_count["phi_bin"] = phi_bin
                    phi_theta_with_max_count["theta_bin"] = theta_bin

        phi_theta_with_min_count = {
            "count": phi_theta_with_max_count["count"],
            "phi_bin": None,
            "theta_bin": None,
        }

        for phi_bin, phi_bin_values in self.cp_phi_theta_bin_values.items():
            for theta_bin, phi_theta_bin_values in phi_bin_values.items():
                count = phi_theta_bin_values["count"]
                if count < phi_theta_with_min_count["count"] and count != 0:
                    phi_theta_with_min_count["count"] = count
                    phi_theta_with_min_count["phi_bin"] = phi_bin
                    phi_theta_with_min_count["theta_bin"] = theta_bin

        self.phi_theta_with_max_count = phi_theta_with_max_count
        self.phi_theta_with_min_count = phi_theta_with_min_count
        print(self.phi_theta_with_max_count)

    def write_phi_theta_count_free_energy(self):

        data_file = "custom_scripts/output/cp_phi_theta_free_energy.dat"

        with open(data_file, "w") as out_file:
            for cp_phi, phi_bin_values in self.cp_phi_theta_bin_values.items():
                for cp_theta, phi_theta_bin_values in phi_bin_values.items():
                    out_file.write(
                        "{cp_phi} {cp_theta} {free_energy}\n".format(
                            cp_phi=cp_phi,
                            cp_theta=cp_theta,
                            #              count=phi_theta_bin_values["count"],
                            free_energy=phi_theta_bin_values["free_energy_diff"],
                        )
                    )

    def _fix_free_energy(self, cp_phi_theta_bin_values):

        min_phi_bin = self.phi_theta_with_min_count["phi_bin"]
        min_theta_bin = self.phi_theta_with_min_count["theta_bin"]
        max_energy = self.cp_phi_theta_bin_values[min_phi_bin][min_theta_bin][
            "free_energy_diff"
        ]

        for phi_bin, phi_bin_values in self.cp_phi_theta_bin_values.items():
            for theta_bin, phi_theta_bin_values in phi_bin_values.items():
                free_energy_diff = phi_theta_bin_values["free_energy_diff"]
                if free_energy_diff == -1:
                    cp_phi_theta_bin_values[phi_bin][theta_bin][
                        "free_energy_diff"
                    ] = max_energy

        return cp_phi_theta_bin_values

    def sort_cp_phi_cp_theta_binned_energies_for_contour(self):

        cp_phi_values = []
        cp_theta_values = []
        free_energy_diff_values = []
        count_values = []

        for phi_bin, phi_bin_values in self.cp_phi_theta_bin_values.items():
            cp_phi_values.append([])
            cp_theta_values.append([])
            free_energy_diff_values.append([])
            count_values.append([])
            for theta_bin, phi_theta_bin_values in phi_bin_values.items():
                cp_phi_values[-1].append(phi_bin)
                cp_theta_values[-1].append(theta_bin)
                free_energy_diff_values[-1].append(
                    phi_theta_bin_values["free_energy_diff"]
                )
                count_values[-1].append((phi_theta_bin_values["count"]))
                if phi_theta_bin_values["count"] == -1:
                    count_values[-1][-1] = 0

        return cp_phi_values, cp_theta_values, free_energy_diff_values, count_values

    def plot_binned_energies_against_cp_phi_and_cp_theta(self):

        (
            cp_phi_values,
            cp_theta_values,
            free_energy_values,
            count_values,
        ) = self.sort_cp_phi_cp_theta_binned_energies_for_contour()

        scatter_params = {
            "y_label": "Theta (deg)",
            "x_label": "Phi (deg)",
            "countour_levels": range(0, 6, 1),
        }

        plot_file_path = "custom_scripts/output/cp_phi_theta_free_energy_contour.png"

        Plot().contour_plot(
            cp_phi_values,
            cp_theta_values,
            free_energy_values,
            plot_file_path,
            scatter_params=scatter_params,
        )

        # cp_phi_1d_values = [
        #     cp_params["cp_phi"]
        #     for cp_params in self.cp_trajectory_pucker_params.values()
        # ]
        # cp_theta_1d_values = [
        #     cp_params["cp_theta"]
        #     for cp_params in self.cp_trajectory_pucker_params.values()
        # ]

        # polar_scatter_plot_file_path = (
        #     "custom_scripts/output/cp_phi_theta_free_energy_polar_scatter.png"
        # )

        # Plot().polar_scatter(
        #     cp_phi_1d_values,
        #     cp_theta_1d_values,
        #     polar_scatter_plot_file_path,
        #     scatter_params=scatter_params,
        # )

        polar_scatter_heatmap_file_path = (
            "custom_scripts/output/cp_phi_theta_count_heatmap.png"
        )

        Plot().polar_scatter_heatmap(
            cp_phi_values,
            cp_theta_values,
            count_values,
            polar_scatter_heatmap_file_path,
            scatter_params=scatter_params,
        )
