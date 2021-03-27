import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.ndimage
import helper
import math


def scatter_with_pmf_plumed_contour(pmf_file, xList=[], yList=[]):
    """Plot single PMF and output pdf
    """
    x = []
    y = []
    z = []

    # process PMF file into 2D arrays
    with open(pmf_file) as f:
        lines = f.readlines()
        xline = []
        yline = []
        zline = []
        for line in lines:
            if line[0] != "#" and len(line) >= 3:
                nextY = float(line.split()[1]) * (180 / math.pi)
                if yline and (yline[-1] != nextY):  # end of row
                    x.append(xline)
                    y.append(yline)
                    z.append(zline)
                    xline = []
                    yline = []
                    zline = []
                yline.append(nextY)
                xline.append(float(line.split()[0]) * (180 / math.pi))
                zline.append(float(line.split()[2]))
        x.append(xline)  # do last append
        y.append(yline)
        z.append(zline)

    # smoothing data
    z = scipy.ndimage.gaussian_filter(
        z, sigma=1
    )  # bigger sigma = more smoothing; can go <1

    fig = plt.figure()
    a = fig.add_subplot(1, 1, 1)

    levs = range(1, 8)  # levels to draw contours at

    CS = a.contour(x, y, z, levs, cmap=cm.copper_r, linewidths=2)
    # plt.clabel(CS, inline=1, fontsize=10)

    # Recast levels to new class
    CS.levels = [helper.nf(val) for val in CS.levels]

    plt.clabel(CS, CS.levels, inline=True, fmt="%r ", fontsize=8)

    if xList and yList:
        plt.scatter(xList, yList, s=5, color="g")

    a.set_xlabel(r"$\phi$", fontsize=20)
    a.set_ylabel(r"$\theta$", fontsize=20)
    a.set_xticks((0, 60, 120, 180, 240, 300, 360))
    a.set_yticks(range(0, 190, 10))
    a.tick_params(axis="both", labelsize=10)

    fig.savefig("out_pmf.png", dpi=400, format="png")

    plt.close()


# def scatter_without_pmf_contour(xList, yList, outfile_name, x_key, y_key):

#     # plt.scatter(xList, yList, s=5, color="g")

#     # a.set_xlabel("$\phi$", fontsize=20)
#     # a.set_ylabel("$\psi$", fontsize=20)
#     # a.set_xticks((-120, -60, 0, 60, 120))
#     # a.set_yticks((-120, -60, 0, 60, 120))
#     # a.tick_params(axis="both", labelsize=10)

#     plt.scatter(xList, yList, s=5, color="g")
#     plt.xlabel("$\\" + x_key + "$", fontsize=20)
#     plt.ylabel("$\\" + y_key + "$", fontsize=20)
#     # plt.ylabel("$\psi$", fontsize=20)
#     # plt.show()

#     plt.savefig(outfile_name + ".png", dpi=800, format="png")

#     plt.close()
