import os
import math
import numpy as np
from collections import OrderedDict
from plot.plot import Plot
from trajectory import Trajectory

# from ring_pucker.cp_ring_pucker import CPRingPucker


class CpPuckerEnergyPmf(Trajectory):
    def __init__(self):
        self.boltzmann_constant = 1.38064852e-23
        self.temp_kelvin = 300
        self.avogradros_constant = 6.02214086e23
        self.joues_per_kilocal = 4184
        self.cp_trajectory_pucker_params = {}
        self.cp_phi_and_cp_theta_bin_values = OrderedDict()
        self.cp_theta_bin_values = OrderedDict()
        self.theta_with_max_count = {}

    # get cremer-pople puckering parameters from file
    def read_trjectory_cp_puckering_parameters(self):

        cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlcNAc/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

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

    def count_cp_phi_and_cp_theta_by_bins(self):

        cp_phi_and_cp_theta_bin_values = OrderedDict()

        for phi_bin_value in range(0, 362, 2):
            cp_phi_and_cp_theta_bin_values[phi_bin_value] = OrderedDict()
            for theta_bin_value in range(0, 182, 2):
                cp_phi_and_cp_theta_bin_values[phi_bin_value][
                    theta_bin_value
                ] = OrderedDict()
                cp_phi_and_cp_theta_bin_values[phi_bin_value][theta_bin_value][
                    "count"
                ] = 0

        for cp_params in self.cp_trajectory_pucker_params.values():
            cp_phi = cp_params["cp_phi"]
            cp_theta = cp_params["cp_theta"]
            for phi_bin_value in range(0, 362, 2):
                if phi_bin_value <= cp_phi < phi_bin_value + 2:
                    for theta_bin_value in range(0, 182, 2):
                        if theta_bin_value <= cp_theta < theta_bin_value + 2:
                            cp_phi_and_cp_theta_bin_values[phi_bin_value][
                                theta_bin_value
                            ]["count"] += 1

        self.cp_phi_and_cp_theta_bin_values = cp_phi_and_cp_theta_bin_values

    def count_cp_theta_by_bins(self):

        cp_theta_bin_values = OrderedDict()

        for theta_bin_value in range(0, 182, 2):
            cp_theta_bin_values[theta_bin_value] = {}
            cp_theta_bin_values[theta_bin_value]["count"] = 0

        for cp_params in self.cp_trajectory_pucker_params.values():
            cp_theta = cp_params["cp_theta"]
            for theta_bin_value in range(0, 182, 2):
                if theta_bin_value <= cp_theta < theta_bin_value + 2:
                    cp_theta_bin_values[theta_bin_value]["count"] += 1

        self.cp_theta_bin_values = cp_theta_bin_values

    def _find_theta_bin_with_max_count(self):
        theta_with_max_count = {"count": 0, "theta_bin": None}

        for theta_bin, theta_bin_values in self.cp_theta_bin_values.items():
            count = theta_bin_values["count"]
            if count > theta_with_max_count["count"]:
                theta_with_max_count["count"] = count
                theta_with_max_count["theta_bin"] = theta_bin

        self.theta_with_max_count = theta_with_max_count

    def calc_free_energy_for_bins(self):

        cp_theta_bin_values = self.cp_theta_bin_values
        self._find_theta_bin_with_max_count()

        for theta_bin, theta_bin_values in cp_theta_bin_values.items():
            count = theta_bin_values["count"]
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
        import pdb

        pdb.set_trace()
        self.cp_theta_bin_values = cp_theta_bin_values

    # def fix_free_energy(self):

    #     max_theta_bin = self.theta_with_max_count["theta_bin"]
    #     max_energy = self.cp_theta_bin_values[max_theta_bin]["free_energy_diff"]

    #     for theta_bin, cp_theta_bin_values in self.cp_theta_bin_values.items():
    #         free_energy_diff = cp_theta_bin_values["free_energy_diff"]
    #         if free_energy_diff == -1:
    #             cp_theta_bin_values[theta_bin]["free_energy_diff"] = max_energy

    #     import pdb

    #     pdb.set_trace()
