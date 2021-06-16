import matplotlib.pyplot as plt
from matplotlib import cm
import math
import numpy as np


INFILE_PATH_1 = "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlc/torsion_angles/aDGlc13bDGlc/phi.dat"
# INFILE_PATH_2 = "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlc/torsion_angles/aDGlc13bDGlc/psi.dat"

OUTPUT_PATH = "./"

MERCATOR_PLOT_NAME = "out"
POLAR_PLOT_NAME = "cp_polar_heat_map.png"

PLOT_COLOUR = "red"
y_axis_label = r"$\psi$"
x_axis_label = r"$\phi$"

LABEL_FONT_SIZE = 30
TICK_FONT_SIZE = 20

# Y_TICKS = (0, 30, 60, 90, 120, 150, 180)
Y_TICKS = (-180, -120, -60, 0, 60, 120, 180)
X_TICKS = (-180, -120, -60, 0, 60, 120, 180)

AXIS_RANGE = (-180, 180, -180, 180)

# Y_VALUE_INDEX = 1

FACE_COLOUR = "#D3D3D3"
FACE_COLOUR = "#FFFFFF"

CMAP = "inferno_r"


def calc_bins():

    phi_values = []
    psi_values = []
    with open(INFILE_PATH_1) as fp:
        for line in fp:
            line = line.split()
            phi_values.append(float(line[1]))

    # with open(INFILE_PATH_2) as fp:
    #     for line in fp:
    #         line = line.split()
    #         psi_values.append(float(line[1]))

    bin_counts = {}
    bin_count_list = []

    for i in range(min(X_TICKS), max(X_TICKS), 1):
        bin_counts[i] = {}
        for j in range(min(Y_TICKS), max(Y_TICKS), 1):
            bin_counts[i][j] = 0

    for frame_number, phi in enumerate(phi_values):
        floored_phi = math.floor(phi / 1) * 1
        floored_theta = math.floor(psi_values[frame_number] / 1) * 1
        bin_counts[floored_phi][floored_theta] += 1

    for frame_number, phi in enumerate(phi_values):
        floored_phi = math.floor(phi / 1) * 1
        floored_theta = math.floor(psi_values[frame_number] / 1) * 1
        count = bin_counts[floored_phi][floored_theta]
        bin_count_list.append(math.log(count))

    # with open("cp_phi_theta_bins.txt", "w") as outfile:
    #     outfile.write("#phi theta count\n")
    #     for i in range(-180, 180, 1):
    #         for j in range(-180, 180, 1):
    #             count = bin_counts[i][j]
    #             outfile.write(f"{i} {j} {count}\n")
    #         outfile.write("\n")

    return bin_counts, phi_values, psi_values, bin_count_list


def plot(phi_values, theta_values, bin_count_list):

    x = []
    y = []
    z = []
    max_z = 0
    with open("cp_phi_theta_bins.txt", "r") as f:
        lines = f.readlines()
        xline = []
        yline = []
        zline = []
        for line in lines:
            line = line.strip()
            tmp = line.split()
            if len(tmp) and line[0] != "#" and len(line) >= 3:
                # print("Line is:'", tmp)
                nextX = float(line.split()[0])
                if xline and (xline[-1] != nextX):  # end of row
                    x.append(xline)
                    y.append(yline)
                    z.append(zline)
                    xline = []
                    yline = []
                    zline = []
                nextY = float(line.split()[1])
                nextZ = float(line.split()[2])

                xline.append(nextX)
                yline.append(nextY)
                zline.append(nextZ)
                if nextZ > max_z:
                    max_z = nextZ

        x.append(xline)  # do last append
        y.append(yline)
        z.append(zline)

    fig = plt.figure()
    a = fig.add_subplot(1, 1, 1)

    colors1 = plt.cm.Blues(np.linspace(0.0, 1, 128))
    # colors2 = plt.cm.jet(np.linspace(0, 1, 128))
    colors2 = plt.cm.copper_r(np.linspace(0, 1, 128))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    cmap = cm.colors.LinearSegmentedColormap.from_list("my_colormap", colors)

    cmap = plt.get_cmap(CMAP)

    # cf = a.contourf(
    #     np.array(x),  # [:-1, :-1] + 2 / 2.0,
    #     np.array(y),  # [:-1, :-1] + 2 / 2.0,
    #     np.array(z),
    #     levels=range(0, 100),
    #     cmap=cmap,
    #     extend="max"
    #     # norm=cm.colors.Normalize(),
    # )

    cf = a.scatter(phi_values, theta_values, c=bin_count_list, cmap=cmap, s=1)

    # heatmap = np.histogram2d(phi_values, theta_values, bins=(360, 180))
    # extent = (0, 360, 0, 180)

    # plt.clf()
    # plt.imshow(
    #     heatmap[0].T, extent=extent, origin="lower", cmap=cmap, interpolation="nearest"
    # )

    # plt.title(PLOT_NAME)
    # plt.gca().set_aspect("equal", "box")
    # plt.axis((0, 4, 0, 180))
    plt.axis(AXIS_RANGE)

    a.set_xlabel(x_axis_label, fontsize=LABEL_FONT_SIZE)
    h = a.set_ylabel(y_axis_label, fontsize=LABEL_FONT_SIZE)
    # h.set_rotation(0)
    a.set_xticks(X_TICKS)
    a.set_yticks(Y_TICKS)
    a.set_facecolor(FACE_COLOUR)
    # cbar = fig.colorbar(cf)
    # cbar.ax.get_yaxis().set_ticks([])

    for tick in a.xaxis.get_major_ticks():
        tick.label.set_fontsize(TICK_FONT_SIZE)
    for tick in a.yaxis.get_major_ticks():
        tick.label.set_fontsize(TICK_FONT_SIZE)

    plt.savefig(OUTPUT_PATH + MERCATOR_PLOT_NAME + ".png", dpi=400, bbox_inches="tight")
    return x, y, z


def plot_polor_heatmap(x, y, z, phi_values, theta_values, bin_count_list):

    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(top=180, bottom=0)

    ax.set_rgrids((0, 30, 60, 90, 120, 150, 180))
    ax.set_facecolor(FACE_COLOUR)

    colors1 = plt.cm.Blues(np.linspace(0.0, 1, 128))
    colors2 = plt.cm.copper_r(np.linspace(0, 1, 128))

    # combine them and build a new colormap
    colors = np.vstack((colors1, colors2))
    cmap = cm.colors.LinearSegmentedColormap.from_list("my_colormap", colors)

    cmap = plt.get_cmap(CMAP)
    phi_list = []
    theta_list = []
    count_list = []

    for frame_number, phi in enumerate(phi_values):
        phi_list.append(math.radians(phi))

    cf = ax.scatter(phi_list, theta_values, c=bin_count_list, cmap=cmap, s=1)

    # heatmap = np.histogram2d(phi_values, theta_values, bins=(360, 180))
    # extent = (0, 360, 0, 180)

    # ax.clf()
    # ax.imshow(heatmap[0].T, cmap=cmap)

    # cbar = fig.colorbar(cf)
    # cbar.ax.get_yaxis().set_ticks([])

    plt.savefig(OUTPUT_PATH + POLAR_PLOT_NAME, dpi=400, bbox_inches="tight")
    # plt.show()


def main():
    bin_counts, phi_values, theta_values, logged_bin_count_list = calc_bins()
    x, y, z = plot(phi_values, theta_values, logged_bin_count_list)
    plot_polor_heatmap(x, y, z, phi_values, theta_values, logged_bin_count_list)


main()
