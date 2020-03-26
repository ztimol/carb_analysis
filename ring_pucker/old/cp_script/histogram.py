import sys, os
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import numpy as np
import subprocess as sp
import argparse
import statistics

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


def getData(infile, histogram_params):

    yList = []
    xList = []

    inf = open(infile, "r")
    for line in inf:
        time_in_ns = eval(line.split(",")[0])
        if time_in_ns >= 200:
            xList.append(time_in_ns)
            y_value = eval(line.split(",")[2])
            #            print(frame_num, time_in_ns, y_value)
            yList.append(y_value)

    inf.close()

    return yList, xList


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
        bins=list(range(0, 180, 1)),
        color="g",
        edgecolor="white",
        weights=np.ones(len(data_list)) / len(data_list),
    )
    plt.ylabel(ylabel, fontsize=16)
    plt.gca().yaxis.set_major_formatter(PercentFormatter(1))
    plt.xlabel(xlabel, fontsize=16)
    plt.xlim([0, 180])
    ax.xaxis.set_ticks(np.arange(0, 185, 20))
    # plt.rc("font", **font)
    ax.tick_params(labelsize=14)

    fig.savefig(outfile_name, dpi=400, format="png")
    plt.close()


def main():

    histogram_params = {
        "dcd_freq": 250,
        "ylabel": "Probability",
        "xlabel": r"${\Theta}$[deg]",
    }
    data_list, time_list = getData(args.f, histogram_params)

    # print("mean:", statistics.mean(data_list))
    # print("median:", statistics.median(data_list))
    # print("standard deviation:", statistics.stdev(data_list))
    makeHistogram(args.f, data_list, histogram_params)


main()
