import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
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

    def histogram(data_list, data_file_name, **kwargs):

        histogram_params = kwargs.get("histogram_params", {})

        x_label = histogram_params.get("x_label")
        x_start = histogram_params.get("x_start", 0)
        x_end = histogram_params.get("x_end")

        x_major_tick = histogram_params.get("x_major_tick", 10)

        y_label = histogram_params.get("y_label", "y")
        y_start = histogram_params.get("y_start", 0)
        y_end = histogram_params.get("y_end")

        # fig = plt.figure(figsize=(8, 5))
        fig = plt.figure()
        ax = fig.gca()

        # font = {"size": 40}
        plt.hist(
            data_list,
            #             density=True,
            bins=list(np.arange(0, 20, 0.2)),
            color="#000000",
            edgecolor="white",
            weights=np.ones(len(data_list)) / len(data_list),
        )
        plt.ylabel(ylabel)
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        plt.xlabel(xlabel)
        # plt.xlim([0, 10])
        plt.ylim([0, 0.05])
        ax.xaxis.set_ticks(np.arange(0, 21, 1))
        ax.yaxis.set_ticks(np.arange(0, 0.06, 0.01))
        plt.rc("font", **font)

        fig.savefig(outfile_name, dpi=400, format="png")
        plt.close()
