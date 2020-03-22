import os
import numpy as np
import matplotlib.pyplot as plt


def getOutfileName(fName, extra_name=None):

    try:
        x = fName.split(".")[:-1]
        outfile_name = "".join(x) + extra_name + ".png"
    except:
        outfile_name = fName + extra_name + ".png"

    return outfile_name


def make_scatter(fName, yList, xList, ylabel, out_dir):
    xmajortick = 10000  # scatter_params["xmajor"]
    # ylabel =  # scatter_params["ylabel"]

    outfileName = getOutfileName(fName, extra_name="time_series")
    outfileName = os.path.join(out_dir, fName + "_" + outfileName)

    fig = plt.figure(figsize=(24, 12.76))
    ax = fig.gca()
    font = {"size": 40}
    plt.rc("font", **font)
    plt.scatter(xList, yList, s=5, color="r")

    xend = len(xList)  # scatter_params["xend"]  # or start
    xstart = 0  # scatter_params["xstart"]  # or end

    ax.xaxis.set_ticks(np.arange(xstart, xend + xmajortick, xmajortick))
    plt.xlabel("time (ns)", fontsize=40)
    plt.ylabel(ylabel, fontsize=40)
    plt.xlim([xstart, xend])
    plt.ylim([-200, 200])
    ax.tick_params(axis="x", labelsize=35)
    ax.tick_params(axis="y", labelsize=35)
    plt.grid()

    # fig.savefig(outfileName, dpi=400, format="png")

    plt.savefig(outfileName, dpi=400, format="png")
    plt.close()
