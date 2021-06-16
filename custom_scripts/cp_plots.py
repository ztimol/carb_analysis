import os
import math
import numpy as np
from ring_pucker.cp_ring_pucker_plot import CPRingPuckerPlot

# from trajectory import Trajectory


class CpPlots:
    def __init__(self):
        self.cp_trajectory_pucker_params = {}

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

    # def plot_cp_theta_from_file_values(self):

    #     trajectory_cp_theta = [
    #         pucker_params["cp_theta"]
    #         for pucker_params in self.cp_trajectory_pucker_params.values()
    #     ]

    #     time_series = range(0, len(trajectory_cp_theta))  # [
    #     # frame_num / 25000 for frame_num in range(0, len(trajectory_cp_theta))
    #     # ]

    #     CPRingPuckerPlot().cp_theta_time_series_scatter(
    #         time_series,
    #         trajectory_cp_theta,
    #         "output/aLRha13bDGlcNAc_3_10-aDGlc14bDGlcNAc_-31_-45_attempt_2/ring_pucker/GlcNAc/",
    #         "trajectory_cp_theta",
    #     )
