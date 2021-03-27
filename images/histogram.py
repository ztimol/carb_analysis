import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math
from scipy.stats import binned_statistic


OUTPUT_PATH = "./"
PLOT_NAME = "out"

# INFILE_PATH = #"/home/timol/C6W/Studies/structure_analysis/output/bDGlc13_bDGlc14_bDGlcNAc_glycam/torsion_angles/bDGlc14bDGlcNAc/psi.dat"

# INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/atom_distances/index 6 760/index 6 760.dat"
INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/7b_s_flexneri_6ru/atom_distances/index 6 790/index 6 790.dat"
# INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/y_s_flexneri_6ru/atom_distances/index 6 508/index 6 508.dat"

PLOT_COLOUR = "green"
# y_axis_label = r"$\it{r}$ ($\AA$)"
# x_axis_label = r"t (ns)"

y_axis_label = r"$\it{P}$"
x_axis_label = r"$\it{r}$ ($\AA$)"

# FACE_COLOUR = "#D3D3D3"
FACE_COLOUR = "#FFFFFF"

LABEL_FONT_SIZE = 30
TICK_FONT_SIZE = 20

Y_VALUE_INDEX = 1

x_values = []
y_values = []

with open(INFILE_PATH) as fp:
    for line in fp:
        line = line.split()
        x_values.append((float(line[0]) * 25000) / 1000000)
        y_values.append(float(line[Y_VALUE_INDEX]))

fig = plt.figure(figsize=(12, 6.5))
a = fig.add_subplot(1, 1, 1)

a.hist(y_values, color=PLOT_COLOUR, edgecolor="black", bins=45, density=True)

a.set_xlabel(x_axis_label, fontsize=LABEL_FONT_SIZE)
h = a.set_ylabel(y_axis_label, fontsize=LABEL_FONT_SIZE)

a.set_facecolor(FACE_COLOUR)

for tick in a.xaxis.get_major_ticks():
    tick.label.set_fontsize(TICK_FONT_SIZE)
for tick in a.yaxis.get_major_ticks():
    tick.label.set_fontsize(TICK_FONT_SIZE)

plt.savefig(OUTPUT_PATH + PLOT_NAME + ".png", dpi=300, bbox_inches="tight")
