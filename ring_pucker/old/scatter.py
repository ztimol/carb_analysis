import numpy as np
import matplotlib.pyplot as plt
import subprocess
import argparse
import statistics

parser = argparse.ArgumentParser(description="Create a scatter plot")

parser.add_argument("-f", help="Select data file.")
parser.add_argument("-xstart", help="x-axis start value", type=int)
parser.add_argument("-xend", help="x-axis end value", type=int)
parser.add_argument("-xmajor", help="x-axis major tick interval", type=int)
parser.add_argument("-xlabel", help="x-axis label", type=str)
parser.add_argument("-xsym", help="x-axis symbol", type=str)
parser.add_argument("-ylabel", help="y-axis label", type=str)
parser.add_argument("-ysym", help="y-axis symbol", type=str)

args = parser.parse_args()


def validateParser():

    if args.ylabel and args.ysym:
        input("Do not specify both ylabel and ysym flags.")
        raise SystemExit(0)

    if args.xlabel and args.xsym:
        input("Do not specify both ylabel and ysym flags.")
        raise SystemExit(0)


def getOutfileName(fName):

    try:
        x = fName.split(".")[:-1]
        outfileName = "".join(x) + "_time_series.png"
    except:
        outfileName = fName + "_time_series.png"

    return outfileName


def getData(infile, scatter_params):

    yList = []
    xList = []

    inf = open(infile, "r")
    for line in inf:
        frame_num = eval(line.split()[0])
        if frame_num >= 0:
            time_in_ns = (frame_num * scatter_params["dcd_freq"] * 100) / 1000000.0
            xList.append(time_in_ns)
            y_value = eval(line.split()[1])
            yList.append(y_value)

    inf.close()

    return yList, xList


def makeScatter(fName, yList, xList, scatter_params):

    xmajortick = scatter_params["xmajor"]
    ylabel = scatter_params["ylabel"]
    xlabel = scatter_params["xlabel"]

    outfileName = getOutfileName(fName)

    fig = plt.figure(figsize=(24, 12.75))
    ax = fig.gca()
    font = {"size": 40}
    # plt.rc("font", **font)
    plt.scatter(xList, yList, s=7, color="black")

    #    start, end = ax.get_xlim()
    xend = scatter_params["xend"]  # or start
    xstart = scatter_params["xstart"]  # or end

    ax.xaxis.set_ticks(np.arange(xstart, xend + xmajortick, xmajortick))
    plt.xlabel(xlabel, fontsize=45)
    plt.ylabel(ylabel, fontsize=45)
    plt.xlim([xstart, xend])
    plt.ylim([0, 1])
    ax.tick_params(axis="x", labelsize=40)
    ax.tick_params(axis="y", labelsize=40)
    plt.grid()

    fig.savefig(outfileName, dpi=400, format="png")
    plt.close()
    print("mean:", statistics.mean(yList))
    print("median:", statistics.median(yList))
    print("standard deviation:", statistics.stdev(yList))
    # plt.show()


def main():

    scatter_params = {
        "xmajor": 10,
        "dcd_freq": 250,
        "xstart": 0,
        "xend": 150,
        "ylabel": "Q",
        "xlabel": "t (ns)",
    }

    # validateParser()
    yList, xList = getData(args.f, scatter_params)
    makeScatter(args.f, yList, xList, scatter_params)


main()


# def makeScatterDbl(yListA, xListA, yListB, xListB, yListC, xListC,fName):

#     fig = plt.figure(figsize=(24, 12.5))
#     ax = fig.gca()
#     plt.scatter(xListA, yListA, s=9, color='green')
#     plt.scatter(xListB, yListB, s=9, color='black')
#     plt.scatter(xListC, yListC, s=9, color='#dd1c77')
#     start, end = ax.get_xlim()
#     ax.xaxis.set_ticks(np.arange(start, 310, 10))
#     ax.set_ylim([-200,200])
#     # add some text for labels, title and axes ticks
#     plt.xlabel('t (ns)')
#     plt.ylabel(r'$\phi$')
#     #plt.tight_layout()
#     plt.xlim([0,300])
#     plt.grid()
#     font = {'size':50}
#     plt.rc('font', **font)

#     fig.savefig(fName + ".eps", dpi=1000,format="eps")
#     #plt.show()
