import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math
from scipy.stats import binned_statistic
from matplotlib.ticker import PercentFormatter


OUTPUT_PATH = "./"
PLOT_NAME = "out"

# INFILE_PATH = #"/home/timol/C6W/Studies/structure_analysis/output/bDGlc13_bDGlc14_bDGlcNAc_glycam/torsion_angles/bDGlc14bDGlcNAc/psi.dat"

# INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/y_s_flexneri_6ru/atom_distances/index 6 508/index 6 508.dat"

INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/ring_pucker/ru4/trajectory_cp_phi_theta_Q.dat"

PLOT_COLOUR = "gray"
# y_axis_label = r"$\it{r}$ ($\AA$)"
# x_axis_label = r"t (ns)"

y_axis_label = r"$\it{P}$"
x_axis_label = r"$\theta$"

X_TICKS = (0, 30, 60, 90, 120, 150, 180)
Y_TICKS = (0, 0.05, 0.1, 0.15, 0.2, 0.25)

# FACE_COLOUR = "#D3D3D3"
FACE_COLOUR = "#FFFFFF"

LABEL_FONT_SIZE = 30
TICK_FONT_SIZE = 20

Y_VALUE_INDEX = 2

x_values = []
y_values = []

with open(INFILE_PATH) as fp:
    for line in fp:
        line = line.split()
        x_values.append((float(line[0]) * 25000) / 1000000)
        y_values.append(float(line[Y_VALUE_INDEX]))

fig = plt.figure(figsize=(12, 6.5))
a = fig.add_subplot(1, 1, 1)

plt.gca().yaxis.set_major_formatter(PercentFormatter(1))

a.hist(
    y_values,
    color=PLOT_COLOUR,
    edgecolor="black",
    bins=36,
    weights=np.ones(len(y_values)) / len(y_values),
)

a.set_xlabel(x_axis_label, fontsize=LABEL_FONT_SIZE)
h = a.set_ylabel(y_axis_label, fontsize=LABEL_FONT_SIZE)

a.set_facecolor(FACE_COLOUR)

a.set_xticks(X_TICKS)
a.set_yticks(Y_TICKS)

for tick in a.xaxis.get_major_ticks():
    tick.label.set_fontsize(TICK_FONT_SIZE)
for tick in a.yaxis.get_major_ticks():
    tick.label.set_fontsize(TICK_FONT_SIZE)

plt.savefig(OUTPUT_PATH + PLOT_NAME + ".png", dpi=300, bbox_inches="tight")
