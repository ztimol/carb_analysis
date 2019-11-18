import os
from plot import Plot


class TorsionAnglePlot(Plot):
    def time_series_scatter(
        self, time_series, torsion_angles, torsion_name_dir, torsion_type,
    ):
        scatter_params = {
            "y_label": "$\\" + torsion_type + "$",
            "x_label": "time (ns)",
            "x_major_tick": 5,
            "y_end": 180,
            "y_start": -180,
            "torsion_type": torsion_type,
            "x_end": self.get_trajectory_time_in_ns(),
        }

        plot_file_path = os.path.join(
            torsion_name_dir, torsion_type + "_time_series.png"
        )
        self.two_dimensional_scatter(
            time_series,
            torsion_angles,  # [i[0] for i in torsion_angles],
            plot_file_path,
            scatter_params=scatter_params,
        )

    def torsion_angles_scatter(
        self, x_series, y_series, x_key, y_key, torsion_name_dir
    ):
        scatter_params = {
            "y_label": "$\\" + y_key + "$",
            "x_label": "$\\" + x_key + "$",
            "x_major_tick": 60,
            "x_end": 180,
            "x_start": -180,
            "y_end": 180,
            "y_start": -180,
        }

        plot_file_path = os.path.join(torsion_name_dir, x_key + "_" + y_key + ".png")
        self.two_dimensional_scatter(
            x_series, y_series, plot_file_path, scatter_params=scatter_params
        )

    def two_dimensional_scatter(self, x_values, y_values, plot_file_path, **kwargs):
        scatter_params = kwargs.get("scatter_params", {})
        super().two_dimensional_scatter(
            x_values, y_values, plot_file_path, scatter_params=scatter_params,
        )
