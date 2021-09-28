import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
import math
import numpy as np

STRUCTURE_NAME = "bDGlcNAc12aLRha"

INFILE_PATH_1 = "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/torsion_angles/{structure_name}_torsion_ru2-ru5_4C1_and_boat_only.dat".format(
    structure_name=STRUCTURE_NAME
)

OUTPUT_PATH = "./"

MERCATOR_PLOT_NAME = "{structure_name}_torsion_ru2-ru5_4C1_and_boat_only".format(
    structure_name=STRUCTURE_NAME
)

PLOT_COLOUR = "red"
y_axis_label = r"$\psi$"
x_axis_label = r"$\phi$"

LABEL_FONT_SIZE = 40
TICK_FONT_SIZE = 40

Y_TICKS = (-180, -120, -60, 0, 60, 120, 180)
X_TICKS = (-180, -120, -60, 0, 60, 120, 180)

AXIS_RANGE = (-180, 180, -180, 180)

# Y_VALUE_INDEX = 1

FACE_COLOUR = "#D3D3D3"
FACE_COLOUR = "#FFFFFF"

CMAP = "inferno_r"


def calc_bins():

    x_values = []
    y_values = []
    with open(INFILE_PATH_1) as fp:
        for line in fp:
            line = line.split()
            x_values.append(float(line[1]))
            y_values.append(float(line[2]))

    bin_counts = {}
    bin_count_list = []

    for i in range(min(X_TICKS), max(X_TICKS), 1):
        bin_counts[i] = {}
        for j in range(min(Y_TICKS), max(Y_TICKS), 1):
            bin_counts[i][j] = 0

    for frame_number, x in enumerate(x_values):
        floored_x = math.floor(x / 1) * 1
        floored_y = math.floor(y_values[frame_number] / 1) * 1
        bin_counts[floored_x][floored_y] += 1

    for frame_number, x in enumerate(x_values):
        floored_x = math.floor(x / 1) * 1
        floored_y = math.floor(y_values[frame_number] / 1) * 1
        count = bin_counts[floored_x][floored_y]
        bin_count_list.append(math.log(count))

    return bin_counts, x_values, y_values, bin_count_list


def plot(x_values, y_values, bin_count_list):

    x = []
    y = []
    z = []

    fig = plt.figure(figsize=(12, 6.5))
    a = fig.add_subplot(1, 1, 1)

    colors = [(1, 1, 1), (0, 0, 1), (1, 0.75, 0), (1, 0, 0), (0.75, 0, 0)]

    newcmp = cm.colors.LinearSegmentedColormap.from_list("my_colormap", colors)

    cf = a.scatter(x_values, y_values, c=bin_count_list, cmap=newcmp, s=1)

    plt.axis(AXIS_RANGE)

    # a.set_xlabel(x_axis_label, fontsize=LABEL_FONT_SIZE)
    # h = a.set_ylabel(y_axis_label, fontsize=LABEL_FONT_SIZE)
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

    a.grid(b=True, which="major", axis="both")
    a.set_xticklabels([])
    a.set_yticklabels([])
    # plt.gca().set_position([0, 0, 1, 1])
    plt.savefig(
        OUTPUT_PATH + MERCATOR_PLOT_NAME + ".svg", bbox_inches="tight", format="svg",
    )
    return x, y, z


def main():
    bin_counts, x_values, y_values, logged_bin_count_list = calc_bins()
    x, y, z = plot(x_values, y_values, logged_bin_count_list)


main()
