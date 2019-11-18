import os
from plot import Plot


class TorsionAnglePlot(Plot):
    def two_dimensional_scatter(
        self, time_series, torsion_angles, torsion_name_dir, torsion_type,
    ):

        scatter_params = {
            "y_label": "$\\" + torsion_type + "$",
            "x_major_tick": 5,
            "torsion_type": torsion_type,
            "x_end": self.get_trajectory_time_in_ns(),
        }

        plot_file_path = os.path.join(
            torsion_name_dir, torsion_type + "_time_series.png"
        )
        super().two_dimensional_scatter(
            time_series,
            [i[0] for i in torsion_angles],
            plot_file_path,
            scatter_params=scatter_params,
        )
