import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
import math
import pandas as pd
from scipy.stats import binned_statistic
from matplotlib.ticker import MultipleLocator, AutoMinorLocator
import statistics


OUTPUT_PATH = "./"
PLOT_NAME = "out_test"

INFILE_PATH = "test.dat"


PLOT_COLOUR = "blue"
x_axis_label = r"$\phi$"
y_axis_label = r"$\psi$"

# FACE_COLOUR = "#D3D3D3"
FACE_COLOUR = "#FFFFFF"

LABEL_FONT_SIZE = 50
TICK_FONT_SIZE = 35

Y_TICKS = (-180, -120, -60, 0, 60, 120, 180)
X_TICKS = (-180, -120, -60, 0, 60, 120, 180)

AXIS_RANGE = (-180, 180, -180, 180)

Y_VALUE_INDEX = 1
X_VALUE_INDEX = 0

x_values = []
y_values = []

with open(INFILE_PATH) as fp:
    for line in fp:
        line = line.split()
        x_values.append(float(line[X_VALUE_INDEX]))
        y_values.append(float(line[Y_VALUE_INDEX]))

# print(np.average(y_values[: math.ceil(len(y_values))]))
# fig = plt.figure(figsize=(12, 6.5))
fig = plt.figure(figsize=(12, 9))
# fig = plt.figure()
a = fig.add_subplot(1, 1, 1)

a.plot(x_values, y_values, color=PLOT_COLOUR)  # , marker="+")

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
