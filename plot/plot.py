import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
from scipy.stats import gaussian_kde
from trajectory import Trajectory


class Plot:
    def two_dimensional_scatter(self, x_values, y_values, plot_file_path, **kwargs):

        scatter_params = kwargs.get("scatter_params", {})

        x_label = scatter_params.get("x_label")
        x_start = scatter_params.get("x_start", 0)
        x_end = scatter_params.get("x_end")

        x_major_tick = scatter_params.get("x_major_tick", 10)

        y_label = scatter_params.get("y_label", "y")
        y_start = scatter_params.get("y_start", 0)
        y_end = scatter_params.get("y_end")

        color = scatter_params.get("color", "black")

        fig = plt.figure()
        ax = fig.gca()
        # font = {"size": 40}
        # plt.rc("font", **font)

        plt.scatter(x_values, y_values, s=5, color=color)

        ax.xaxis.set_ticks(np.arange(x_start, x_end + x_major_tick, x_major_tick))
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.xlim([x_start, x_end])
        plt.ylim([y_start, y_end])
        ax.tick_params(axis="x")
        ax.tick_params(axis="y")
        plt.grid()

        plt.savefig(plot_file_path, dpi=400, format="png")
        plt.close()

    def histogram(self, data_list, plot_file_path, **kwargs):

        histogram_params = kwargs.get("histogram_params", {})

        bins = histogram_params.get("bins")

        x_label = histogram_params.get("x_label")
        x_start = histogram_params.get("x_start", min(bins))
        x_end = histogram_params.get("x_end", max(bins))

        x_major_tick = histogram_params.get("x_major_tick", 10)

        y_label = histogram_params.get("y_label", "Probability (%)")
        y_start = histogram_params.get("y_start", 0)
        y_end = histogram_params.get("y_end")

        fig = plt.figure()
        ax = fig.gca()

        plt.hist(
            data_list,
            bins=bins,
            color="#000000",
            edgecolor="white",
            weights=np.ones(len(data_list)) / len(data_list),
        )

        plt.ylabel(y_label)
        plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
        plt.xlabel(x_label)
        # plt.xlim([0, 10])
        # plt.ylim([0, 0.05])
        ax.xaxis.set_ticks(np.arange(x_start, x_end + x_major_tick, x_major_tick))
        # ax.yaxis.set_ticks(np.arange(y_start, 0.06, 0.01))

        fig.savefig(plot_file_path, dpi=400, format="png")
        plt.close()

    def probability_density(x_list, y_list, plot_file_path, **kwargs):
        # fig = plt.figure(figsize=(48, 25))
        fig = plt.figure()
        ax = fig.gca()

        gaussian_density = gaussian_kde(y_list)
        xs = np.linspace(-180, 180, 8000)
        gaussian_density.covariance_factor = lambda: 0.05
        gaussian_density._compute_covariance()

        i = list(densityA(xs)).index(max(densityA(xs)))

        plt.plot(xs, gaussian_density(xs), color="blue")
        ax.fill_between(xs, gaussian_density(xs), interpolate=True, color="blue")

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        # font = {"size": 30}
        # plt.rc("font", **font)

        fig.savefig(plot_file_path, dpi=400, format="png")
