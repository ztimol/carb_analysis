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
        y_label = scatter_params.get("y_label", "y")

        xend = scatter_params.get("x_end")
        xstart = scatter_params.get("x_start", 0)
        x_major_tick = scatter_params.get("x_major_tick", 10)

        fig = plt.figure(figsize=(24, 12.76))
        ax = fig.gca()
        font = {"size": 40}
        plt.rc("font", **font)
        plt.scatter(x_values, y_values, s=5, color="g")

        ax.xaxis.set_ticks(np.arange(xstart, xend + x_major_tick, x_major_tick))
        plt.xlabel("time (ns)", fontsize=40)
        plt.ylabel(y_label, fontsize=40)
        plt.xlim([xstart, xend])
        plt.ylim([-200, 200])
        ax.tick_params(axis="x", labelsize=35)
        ax.tick_params(axis="y", labelsize=35)
        plt.grid()

        # fig.savefig(outfileName, dpi=400, format="png")

        plt.savefig(plot_file_path, dpi=400, format="png")
        plt.close()
