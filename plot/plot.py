import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.ticker import PercentFormatter, MaxNLocator
from matplotlib import cm
import numpy as np
import scipy
import helper
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
        plt.xlabel(x_label, fontsize=20)
        plt.ylabel(y_label, fontsize=20)
        plt.xlim([x_start, x_end])
        plt.ylim([y_start, y_end])
        ax.tick_params(axis="x")
        ax.tick_params(axis="y")
        plt.grid()

        plt.savefig(plot_file_path, dpi=400, format="png")
        plt.close()

    def two_dimensional_scatter_plot_with_average(
        self,
        x_values,
        y_values,
        average_x_values,
        average_y_values,
        plot_file_path,
        **kwargs
    ):

        scatter_params = kwargs.get("scatter_params", {})

        x_label = scatter_params.get("x_label")
        x_start = scatter_params.get("x_start", 0)
        x_end = scatter_params.get("x_end")

        x_major_tick = scatter_params.get("x_major_tick", 10)

        y_label = scatter_params.get("y_label", "y")
        y_start = scatter_params.get("y_start", 0)
        y_end = scatter_params.get("y_end")

        color = scatter_params.get("color", "black")
        average_line_color = scatter_params.get("average_line_color", "yellow")

        fig = plt.figure()
        ax = fig.gca()

        plt.scatter(x_values, y_values, s=5, color=color)

        plt.plot(average_x_values, average_y_values, color=average_line_color)

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

    def two_dimensional_line_plot(self, x_values, y_values, plot_file_path, **kwargs):
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

        plt.plot(x_values, y_values, color=color)

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
        plt.xlim([x_start, x_end])
        plt.ylim([y_start, y_end])
        ax.xaxis.set_ticks(np.arange(x_start, x_end + x_major_tick, x_major_tick))
        # ax.yaxis.set_ticks(np.arange(y_start, 0.06, 0.01))

        fig.savefig(plot_file_path, dpi=400, format="png")
        plt.close()

    # def probability_density(self, x_list, y_list, plot_file_path, **kwargs):
    #     # fig = plt.figure(figsize=(48, 25))
    #     fig = plt.figure()
    #     ax = fig.gca()

    #     gaussian_density = gaussian_kde(y_list)
    #     xs = np.linspace(-180, 180, 8000)
    #     gaussian_density.covariance_factor = lambda: 0.05
    #     gaussian_density._compute_covariance()

    #     i = list(densityA(xs)).index(max(densityA(xs)))

    #     plt.plot(xs, gaussian_density(xs), color="blue")
    #     ax.fill_between(xs, gaussian_density(xs), interpolate=True, color="blue")

    #     plt.xlabel(x_label)
    #     plt.ylabel(y_label)
    #     # font = {"size": 30}
    #     # plt.rc("font", **font)

    #     fig.savefig(plot_file_path, dpi=400, format="png")

    def contour_plot(self, x_variables, y_variables, z_variables, out_file, **kwargs):

        scatter_params = kwargs.get("scatter_params", {})

        x_label = scatter_params.get("x_label")
        y_label = scatter_params.get("y_label")
        countour_levels = scatter_params.get("countour_levels")

        if not countour_levels:
            countour_levels = range(1, 11)  # levels to draw contours at

        fig = plt.figure()
        a = fig.add_subplot(1, 1, 1)
        smoothed_z_variables = scipy.ndimage.gaussian_filter(
            z_variables, sigma=1
        )  # bigger sigma = more smoothing; can go <1

        CS = a.contour(
            x_variables,
            y_variables,
            smoothed_z_variables,
            countour_levels,
            cmap=cm.copper_r,
            linewidths=2,
        )

        plt.clabel(CS, CS.levels, inline=True, fmt="%r ", fontsize=8)

        a.set_xlabel(x_label, fontsize=20)
        a.set_ylabel(y_label, fontsize=20)
        a.set_xticks((0, 60, 120, 180, 240, 300, 360))
        a.set_yticks((0, 30, 90, 60, 120, 150, 180))
        a.tick_params(axis="both", labelsize=10)

        fig.savefig(out_file, dpi=400, format="png")

        plt.close()

    def polar_contour(self, phi_values, theta_values, z_values, out_file, **kwargs):

        scatter_params = kwargs.get("scatter_params", {})

        countour_levels = scatter_params.get("countour_levels")

        if not countour_levels:
            countour_levels = range(0, 16)  # levels to draw contours at

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_rlim(top=180, bottom=0)
        contour_plot = ax.contour(
            phi_values,
            theta_values,
            z_values,
            countour_levels,
            cmap=cm.binary,
            linewidths=2,
        )

        plt.clabel(
            contour_plot, contour_plot.levels, inline=True, fmt="%r ", fontsize=8
        )

        fig.savefig(out_file, dpi=400, format="png")

    def polar_scatter(self, phi, theta, out_file, **kwargs):

        scatter_params = kwargs.get("scatter_params", {})

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_rlim(top=180, bottom=0)
        ax.scatter(phi, theta, c=theta, s=1, cmap="hsv", alpha=0.75)

        fig.savefig(out_file, dpi=400, format="png")

    def polar_scatter_heatmap(
        self, phi_values, theta_values, count_values, out_file, **kwargs
    ):

        scatter_params = kwargs.get("scatter_params", {})

        smoothed_count_values = scipy.ndimage.gaussian_filter(
            count_values, sigma=1
        )  # bigger sigma = more smoothing; can go <1

        # count_values = np.array(count_values)  # [:-1, :-1]
        # levels = MaxNLocator(nbins=30).tick_values(
        #     count_values.min(), count_values.max()
        # )
        levels = range(0, 100)

        cmap = plt.get_cmap("bone_r")

        fig = plt.figure()
        ax = fig.add_subplot(111, projection="polar")
        ax.set_theta_zero_location("N")
        ax.set_theta_direction(-1)
        ax.set_rlim(top=180, bottom=0)

        cf = ax.contourf(
            np.array(phi_values),  # [:-1, :-1] + 2 / 2.0,
            np.array(theta_values),  # [:-1, :-1] + 2 / 2.0,
            np.array(count_values),
            levels=levels,
            cmap=cmap,
            extend="max"
            # norm=cm.colors.Normalize(),
        )

        # heatmap = np.histogram2d(phi_values, theta_values, bins=(360, 180))
        # extent = (0, 360, 0, 180)

        # plt.clf()

        # plt.imshow(heatmap[0].T, extent=extent, origin="lower")
        fig.colorbar(cf)

        fig.savefig(out_file, dpi=400, format="png")
