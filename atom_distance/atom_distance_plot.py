import os
import numpy as np
from plot.plot import Plot
from trajectory import Trajectory


class AtomDistancePlot(Plot):
    def time_series_scatter(
        self, time_series, atom_distances, atom_selection_out_dir, atom_selection
    ):

        scatter_params = {
            "y_label": "Atom Pair Distance (A)",
            "x_label": "time (ns)",
            "x_major_tick": 200,
            # "y_end": 180,
            "y_start": 0,
            "x_end": self.get_trajectory_time_in_ns(),
        }

        plot_file_path = os.path.join(
            atom_selection_out_dir, atom_selection + "_time_series.png"
        )
        self.two_dimensional_scatter(
            time_series, atom_distances, plot_file_path, scatter_params=scatter_params
        )

    def two_dimensional_scatter(self, x_values, y_values, plot_file_path, **kwargs):
        scatter_params = kwargs.get("scatter_params", {})
        super().two_dimensional_scatter(
            x_values, y_values, plot_file_path, scatter_params=scatter_params
        )

    def probability_histogram(self, data_list, atom_selection_out_dir, atom_selection):
        histogram_params = {
            "x_label": "Atom Pair Distance (A)",
            "x_major_tick": 1,
            "bins": list(np.arange(0, 5, 0.5)),
        }
        plot_file_path = os.path.join(
            atom_selection_out_dir, atom_selection + "_histogram.png"
        )
        self.histogram(data_list, plot_file_path, histogram_params=histogram_params)

    def histogram(self, data_list, plot_file_path, **kwargs):
        histogram_params = kwargs.get("histogram_params", {})
        super().histogram(data_list, plot_file_path, histogram_params=histogram_params)
