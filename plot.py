import matplotlib.pyplot as plt
import numpy as np
from trajectory import Trajectory


class Plot:

    # def _getOutfileName(fName, extra_name=None):

    # try:
    #     x = fName.split(".")[:-1]
    #     outfile_name = "".join(x) + extra_name + ".png"
    # except:
    #     outfile_name = fName + extra_name + ".png"

    # return outfile_name

    def two_dimensional_scatter(self, x_values, y_values, plot_file_path, **kwargs):

        scatter_params = kwargs.get("scatter_params", {})

        x_label = scatter_params.get("x_label")
        x_start = scatter_params.get("x_start", 0)
        x_end = scatter_params.get("x_end")

        x_major_tick = scatter_params.get("x_major_tick", 10)

        y_label = scatter_params.get("y_label", "y")
        y_start = scatter_params.get("y_start", 0)
        y_end = scatter_params.get("y_end")

        fig = plt.figure()
        ax = fig.gca()
        # font = {"size": 40}
        # plt.rc("font", **font)
        plt.scatter(x_values, y_values, s=5, color="g")

        ax.xaxis.set_ticks(np.arange(x_start, x_end + x_major_tick, x_major_tick))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim([x_start, x_end])
        plt.ylim([y_start, y_end])
        ax.tick_params(axis="x")
        ax.tick_params(axis="y")
        plt.grid()

        # fig.savefig(outfileName, dpi=400, format="png")

        plt.savefig(plot_file_path, dpi=400, format="png")
        plt.close()
