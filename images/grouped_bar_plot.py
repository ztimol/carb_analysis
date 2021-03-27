import matplotlib.pyplot as plt
from matplotlib import cm
import math
import numpy as np
from matplotlib.ticker import MultipleLocator, AutoMinorLocator


# INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_bDGlc14_bDGlcNAc/torsion_angles/aDGlc13bDGlcNAc/phi.dat"

# INFILE_PATHS = [
#     "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlc/ring_pucker/Glc/trajectory_cp_phi_theta_Q.dat",
#     "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlcNAc/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat",
#     "/home/timol/C6W/Studies/structure_analysis/output/bDGlc13_aDGlc14_bDGlcNAc/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat",
#     "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_bDGlc14_bDGlcNAc/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat",
#     "/home/timol/C6W/Studies/structure_analysis/output/bDGlc13_bDGlc14_bDGlcNAc/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat",
# ]

INFILE_PATHS = [
    "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlc_glycam/ring_pucker/Glc/trajectory_cp_phi_theta_Q.dat",
    "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlcNAc_glycam/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat",
    "/home/timol/C6W/Studies/structure_analysis/output/bDGlc13_aDGlc14_bDGlcNAc_glycam/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat",
    "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_bDGlc14_bDGlcNAc_glycam/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat",
    "/home/timol/C6W/Studies/structure_analysis/output/bDGlc13_bDGlc14_bDGlcNAc_glycam/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat",
]

# OUTPUT_PATH = "/home/timol/C6W/Studies/Dynamics/NAMD/MD/solution/general_structures/aDGlc13_bDGlc14_bDGlcNAc/images/ring_pucker/"
OUTPUT_PATH = "./"
PLOT_NAME = "out_bar"

PLOT_COLOUR = "orange"
y_axis_label = r"$\theta$"
x_axis_label = r"$\phi$"

LABEL_FONT_SIZE = 20
TICK_FONT_SIZE = 15

# Y_TICKS = (0, 30, 60, 90, 120, 150, 180)
Y_TICKS = (0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5)
X_TICKS = (0, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)
# AXIS_RANGE = (0, 360, 0, 5)

# Y_TICKS = (-180, -120, -60, 0, 60, 120, 180)

# Y_VALUE_INDEX = 1

FACE_COLOUR = "#D3D3D3"
# FACE_COLOUR = "#FFE4E1"

CMAP = "inferno_r"


def get_conformers_list():

    return [
        # {
        #     "label": r"$^{4}C_{1}$",
        #     "phi_lb": 0,
        #     "phi_ub": 360,
        #     "theta_lb": 0,
        #     "theta_ub": 25,
        #     "count": 0,
        # },
        # {
        #     "label": r"$^{2}H_{3}$",
        #     "phi_lb": 110,
        #     "phi_ub": 130,
        #     "theta_lb": 45,
        #     "theta_ub": 65,
        #     "count": 0,
        # },
        # {
        #     "label": r"$E_{3}$",
        #     "phi_lb": 170,
        #     "phi_ub": 190,
        #     "theta_lb": 45,
        #     "theta_ub": 65,
        #     "count": 0,
        # },
        # {
        #     "label": r"$^{4}H_{3}$",
        #     "phi_lb": 200,
        #     "phi_ub": 220,
        #     "theta_lb": 45,
        #     "theta_ub": 65,
        #     "count": 0,
        # },
        # {
        #     "label": r"$^{4}E$",
        #     "phi_lb": 230,
        #     "phi_ub": 250,
        #     "theta_lb": 45,
        #     "theta_ub": 65,
        #     "count": 0,
        # },
        # {
        #     "label": r"$^{4}H_{5}$",
        #     "phi_lb": 260,
        #     "phi_ub": 280,
        #     "theta_lb": 45,
        #     "theta_ub": 65,
        #     "count": 0,
        # },
        # {
        #     "label": r"$^{O}H_{5}$",
        #     "phi_lb": 320,
        #     "phi_ub": 340,
        #     "theta_lb": 45,
        #     "theta_ub": 65,
        #     "count": 0,
        # },
        {
            "label": r"$^{3,O}B$",
            "phi_lb": 0,
            "phi_ub": 10,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$^{3}S_{1}$",
            "phi_lb": 20,
            "phi_ub": 40,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$B_{1,4}$",
            "phi_lb": 50,
            "phi_ub": 70,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$^{5}S_{1}$",
            "phi_lb": 80,
            "phi_ub": 100,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$^{2,5}B$",
            "phi_lb": 110,
            "phi_ub": 130,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$^{2}S_{O}$",
            "phi_lb": 140,
            "phi_ub": 160,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$B_{3,O}$",
            "phi_lb": 170,
            "phi_ub": 190,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$^{1}S_{3}$",
            "phi_lb": 200,
            "phi_ub": 220,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$^{1,4}B$",
            "phi_lb": 230,
            "phi_ub": 250,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$^{1}S_{5}$",
            "phi_lb": 260,
            "phi_ub": 280,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$B_{2,5}$",
            "phi_lb": 290,
            "phi_ub": 310,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        {
            "label": r"$^{O}S_{2}$",
            "phi_lb": 320,
            "phi_ub": 340,
            "theta_lb": 80,
            "theta_ub": 100,
            "count": 0,
        },
        # {
        #     "label": r"$^{1}C_{4}$",
        #     "phi_lb": 0,
        #     "phi_ub": 360,
        #     "theta_lb": 155,
        #     "theta_ub": 180,
        #     "count": 0,
        # },
    ]


def conformer_count_per_structure(infile):

    phi_values = []
    theta_values = []
    with open(infile) as fp:
        for line in fp:
            line = line.split()
            phi_values.append(float(line[1]))
            theta_values.append(float(line[2]))

    conformer_type_and_counts = get_conformers_list()

    for phi, theta in zip(phi_values, theta_values):
        for conformer in conformer_type_and_counts:
            phi_lb = conformer["phi_lb"]
            phi_ub = conformer["phi_ub"]
            if phi_lb <= phi <= phi_ub:
                theta_lb = conformer["theta_lb"]
                theta_ub = conformer["theta_ub"]
                if theta_lb <= theta <= theta_ub:
                    conformer["count"] += 1
                    break

            if phi_lb == 0 and phi_ub == 10:
                if 350 <= phi <= 360:
                    theta_lb = conformer["theta_lb"]
                    theta_ub = conformer["theta_ub"]
                    if theta_lb <= theta <= theta_ub:
                        conformer["count"] += 1
                        break

    number_of_frames = len(phi_values)
    for conformer in conformer_type_and_counts:
        conformer["percentage_count"] = (conformer["count"] / number_of_frames) * 100

    return conformer_type_and_counts


def prepare_values_for_bar_chart(conformers_with_counts_for_all_structures):

    label_list = [x["label"] for x in get_conformers_list()]

    percentage_population_list = []

    for conformers_of_structure in conformers_with_counts_for_all_structures:
        percentage_population_list.append([])
        for conformer in conformers_of_structure:
            percentage_population_list[-1].append(conformer["percentage_count"])

    return label_list, percentage_population_list


def plot(labels, percentage_population_list):

    x = np.arange(len(labels))  # the label locations
    width = 0.15  # the width of the bars

    fig, ax = plt.subplots(figsize=(12, 6))

    rects1 = ax.barh(x - 2 * width, percentage_population_list[0], width, label="aaG")
    rects2 = ax.barh(x - width, percentage_population_list[1], width, label="aaGN")
    rects3 = ax.barh(x, percentage_population_list[2], width, label="baGN")
    rects4 = ax.barh(x + width, percentage_population_list[3], width, label="abGN")
    rects5 = ax.barh(x + 2 * width, percentage_population_list[4], width, label="bbGN")

    # Add some text for labels, title and custom x-axis tick labels, etc.
    # ax.set_xlabel("% Occupancy", fontsize=LABEL_FONT_SIZE)

    ax.set_xticks(Y_TICKS)

    # ax.set_xticklabels(labels, fontsize=TICK_FONT_SIZE)
    # ax.legend()
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax.set_yticks(x)
    ax.set_yticklabels(labels, fontsize=TICK_FONT_SIZE)

    for tick in ax.xaxis.get_major_ticks():
        tick.label.set_fontsize(TICK_FONT_SIZE)
    for tick in ax.yaxis.get_major_ticks():
        tick.label.set_fontsize(TICK_FONT_SIZE)

    # autolabel(rects1, ax)
    # autolabel(rects2, ax)
    # autolabel(rects3, ax)
    # autolabel(rects4, ax)
    # autolabel(rects5, ax)
    ax.set_facecolor(FACE_COLOUR)

    plt.savefig(OUTPUT_PATH + PLOT_NAME + ".png", dpi=400, bbox_inches="tight")


def autolabel(rects, ax):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for rect in rects:
        height = round(rect.get_height(), 1)
        if height > 2:
            ax.annotate(
                "{}".format(height),
                xy=(rect.get_x() + rect.get_width() / 2, height),
                xytext=(0, 3),  # 3 points vertical offset
                textcoords="offset points",
                ha="center",
                va="bottom",
            )


def main():
    conformers_with_counts_for_all_structures = []
    for infile in INFILE_PATHS:
        counts = conformer_count_per_structure(infile)
        conformers_with_counts_for_all_structures.append(counts)

    labels, percentage_population_list = prepare_values_for_bar_chart(
        conformers_with_counts_for_all_structures
    )

    plot(labels, percentage_population_list)


#        count_by_conformer(bin_counts)
# x, y, z = plot(phi_values, theta_values, logged_bin_count_list)


main()
