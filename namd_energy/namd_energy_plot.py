import os
import numpy as np
from plot.plot import Plot
from trajectory import Trajectory


class NAMDEnergyPlot(Plot):
    def time_series_scatter(
        self, time_series, energies, namd_energy_name_dir, namd_energy_name
    ):

        scatter_params = {
            "y_label": "Potential Energy (kcal)",
            "x_label": "time (ns)",
            "x_major_tick": 50,
            # "y_end": 180,
            "y_start": 230,
            "x_end": self.get_trajectory_time_in_ns(),
        }

        plot_file_path = os.path.join(
            namd_energy_name_dir, namd_energy_name + "_time_series.png"
        )
        self.two_dimensional_scatter(
            time_series, energies, plot_file_path, scatter_params=scatter_params
        )

    def two_dimensional_scatter(self, x_values, y_values, plot_file_path, **kwargs):
        scatter_params = kwargs.get("scatter_params", {})
        super().two_dimensional_scatter(
            x_values, y_values, plot_file_path, scatter_params=scatter_params
        )

    def probability_histogram(self, data_list, namd_energy_name_dir, namd_energy_type):
        histogram_params = {
            "x_label": "Potential energy (kcal)",
            "x_major_tick": 5,
            "bins": list(np.arange(240, 290, 5)),
        }
        plot_file_path = os.path.join(
            namd_energy_name_dir, namd_energy_type + "_histogram.png"
        )
        self.histogram(data_list, plot_file_path, histogram_params=histogram_params)

    def histogram(self, data_list, plot_file_path, **kwargs):
        histogram_params = kwargs.get("histogram_params", {})
        super().histogram(data_list, plot_file_path, histogram_params=histogram_params)
