import sys, os
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import subprocess as sp
import argparse

parser = argparse.ArgumentParser(description="Create a histogram")

parser.add_argument("-f", help="Select data file.")

args = parser.parse_args()


def getOutfileName(fName):

    try:
        x = fName.split(".")[:-1]
        outfileName = "".join(x) + "_histogram.png"
    except:
        outfileName = fName + "_histogram.png"

    return outfileName


def makeHistogram(data_file_name, data_list, histogram_params):

    outfile_name = getOutfileName(data_file_name)
    ylabel = histogram_params["ylabel"]
    xlabel = histogram_params["xlabel"]

    fig = plt.figure(figsize=(8, 5))
    ax = fig.gca()

    font = {"size": 40}
    plt.hist(
        data_list,
        #             density=True,
        bins=list(np.arange(0, 20, 0.2)),
        color="#000000",
        edgecolor="white",
        weights=np.ones(len(data_list)) / len(data_list),
    )
    plt.ylabel(ylabel)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel(xlabel)
    # plt.xlim([0, 10])
    plt.ylim([0, 0.08])
    ax.xaxis.set_ticks(np.arange(0, 21, 1))
    ax.yaxis.set_ticks(np.arange(0, 0.06, 0.01))
    plt.rc("font", **font)

    fig.savefig(outfile_name, dpi=400, format="png")
    plt.close()


def main():

    histogram_params = {"dcd_freq": 250, "ylabel": "Probability", "xlabel": "RMSD"}
    data_list, time_list = getData(args.f, histogram_params)

    makeHistogram(args.f, data_list, histogram_params)


main()
