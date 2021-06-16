import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math
from scipy.stats import binned_statistic
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import statistics


OUTPUT_PATH = "./"
PLOT_NAME = "out"

# INFILE_PATH = #"/home/timol/C6W/Studies/structure_analysis/output/bDGlc13_bDGlc14_bDGlcNAc_glycam/torsion_angles/bDGlc14bDGlcNAc/psi.dat"
# INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/y_s_flexneri_6ru/atom_distances/index 6 508/index 6 508.dat"
# INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/ring_pucker/ru5/trajectory_cp_phi_theta_Q.dat"

# INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/y_s_flexneri_6ru/atom_distances/index 33 488/index 33 488.dat"
INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/atom_distances/index 75 740/index 75 740.dat"


PLOT_COLOUR = "orange"
y_axis_label = r"$\it{r}$ ($\AA$)"
# y_axis_label = r"$\theta$"
x_axis_label = r"t (ns)"

# FACE_COLOUR = "#D3D3D3"
FACE_COLOUR = "#FFFFFF"

LABEL_FONT_SIZE = 30
TICK_FONT_SIZE = 20

# Y_TICKS = (0, 30, 60, 90, 120, 150, 180)
# Y_TICKS = (0, 60, 120, 180, 240, 300, 360)
# Y_TICKS = (-180, -120, -60, 0, 60, 120, 180)
Y_TICKS = (0, 10, 20, 30, 40, 50, 60, 70, 80)
X_TICKS = (0, 400, 800, 1200, 1600, 2000)
# X_TICKS = (0, 0.5, 1, 1.5, 2)

# Y_TICKS = (1.7, 1.9, 2.1, 2.3, 2.5, 2.7, 2.9, 3.1, 3.3, 3.5)


AXIS_RANGE = (0, 1000, 1.7, 3.5)
Y_VALUE_INDEX = 1

x_values = []
y_values = []


def binned_stats(data, bins, data_range, bin_size):
    bin_stats = binned_statistic(
        y_values, y_values, bins=bins, range=data_range, statistic="count"
    )

    perc_bins = [(x / len(data)) * 100 for x in bin_stats[0]]

    full_bins = [x for x in range(data_range[0], data_range[1], bin_size)]

    i = 0
    for bin_start in full_bins:
        print(bin_start, perc_bins[i])
        i += 1


with open(INFILE_PATH) as fp:
    for line in fp:
        line = line.split()
        x_values.append((float(line[0]) * 25000) / 1000000)
        y_values.append(float(line[Y_VALUE_INDEX]))

import pdb

pdb.set_trace()

# --------------
# binned_stats(y_values, 36, (-180, 180), 10)
# binned_stats(y_values, 1, (-20, 20), 40)
# binned_stats(y_values, 1, (-180, -140), 40)
# binned_stats(y_values, 1, (140, 180), 40)

# print(np.average(y_values[: math.ceil(len(y_values))]))
fig = plt.figure(figsize=(12, 6.5))
a = fig.add_subplot(1, 1, 1)

a.scatter(x_values, y_values, s=0.3, color=PLOT_COLOUR)  # , marker="+")

# plt.title(PLOT_NAME)
# plt.gca().set_aspect("equal", "box")
plt.axis(AXIS_RANGE)

a.set_xlabel(x_axis_label, fontsize=LABEL_FONT_SIZE)
h = a.set_ylabel(y_axis_label, fontsize=LABEL_FONT_SIZE)
# h.set_rotation(0)
a.set_xticks(X_TICKS)
a.set_yticks(Y_TICKS)
a.set_facecolor(FACE_COLOUR)

for tick in a.xaxis.get_major_ticks():
    tick.label.set_fontsize(TICK_FONT_SIZE)
for tick in a.yaxis.get_major_ticks():
    tick.label.set_fontsize(TICK_FONT_SIZE)

a.xaxis.set_minor_locator(MultipleLocator(50))

plt.savefig(OUTPUT_PATH + PLOT_NAME + ".png", dpi=300, bbox_inches="tight")
