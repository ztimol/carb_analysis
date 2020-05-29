import os
import math
import numpy as np
from plot.plot import Plot


class SinglePointEnergy:
    def __init__(self):
        self.trajectory_pe_values = []
        self.pe_periodic_mean_values = []

    # get cremer-pople puckering parameters from file
    def read_trjectory_pe_values(self, namd_energy_file):

        trajectory_pe_values = []

        with open(namd_energy_file, "r") as infile:
            for line in infile:
                line = line.split()
                pe = eval(line[1])
                trajectory_pe_values.append(pe)

        self.trajectory_pe_values = trajectory_pe_values

    def calc_periodic_mean_of_energy(self):

        pe_periodic_mean_values = []
        temp_values = []

        for index, pe in enumerate(self.trajectory_pe_values):

            temp_values.append(pe)
            if (index + 1) % 1000 == 0:
                pe_periodic_mean_values.append(np.mean(temp_values))
                temp_values = []

        self.pe_periodic_mean_values = pe_periodic_mean_values

    def plot_energy_with_average(self):

        scatter_params = {
            "y_label": "PE (kcal/mol)",
            "x_label": "time (ns)",
            "y_start": 230,
            "y_end": 310,
            "x_major_tick": 200,
            "x_end": 2000,
            "color": "black",
        }

        plot_file_path = "custom_scripts/output/trajectory_single_point_potential_energy_with_average.png"

        # time_series = np.arange(0, 1500, 1500 / 60000)
        time_series = np.arange(0, 2000, 2000 / 80000)

        Plot().two_dimensional_scatter_plot_with_average(
            time_series,
            self.trajectory_pe_values,
            np.arange(12.5, 2000, 25),
            self.pe_periodic_mean_values,
            plot_file_path,
            scatter_params=scatter_params,
        )
