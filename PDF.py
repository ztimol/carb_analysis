import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
import seaborn as sns


def getData(fName):

    infile = fName + ".dat"

    yList = []
    xList = []

    inf = open(infile, "r")

    for line in inf:
        xList.append(eval(line.split()[0]))
        yList.append(eval(line.split()[1]))

    inf.close()

    omgList = []

    for omg in yList:
        if omg < 0:
            omg += 360
        omgList.append(omg)

    for omg in yList:
        if omg > 0:
            omg -= 360
        omgList.append(omg)

    return omgList, xList

def getDataDbl(fNameA, fNameB):

    infileA = fNameA + ".dat"
    infileB = fNameB + ".dat"

    yListA = []
    xListA = []

    yListB = []
    xListB = []

    infA = open(infileA, "r")

    for line in infA:
        # xListA.append(eval(line.split()[0]))
        yListA.append(eval(line.split()[1]))

    infA.close()

    infB = open(infileB, "r")

    for line in infB:
        xListB.append(eval(line.split()[0]))
        yListB.append(eval(line.split()[1]))

    infB.close()

    return yListA, xListA, yListB, xListB


def makeplot(yList, xList, fName):
    fig = plt.figure(figsize=(48, 25))
    ax = fig.gca()

    densityA = gaussian_kde(yList)
    xs = np.linspace(-180,180,8000)
    densityA.covariance_factor = lambda : .05
    densityA._compute_covariance()

    i = list(densityA(xs)).index(max(densityA(xs)))
    print(max(densityA(xs)))
    print(list(xs)[i])
    

    plt.plot(xs, densityA(xs), color='blue')
    ax.fill_between(xs,densityA(xs), interpolate=True, color='blue')

    #ax.xaxis.set_ticks(np.arange(-200, 200, 10))
    # add some text for labels, title and axes ticks
    plt.xlabel(r'$\psi$')
    #plt.ylabel("Probability")
    #plt.tight_layout()
    #plt.xlim([0,200])
    #plt.grid()
    font = {'size':30}
    plt.rc('font', **font)



    fig.savefig(fName + "_PDF.eps", dpi=1000,format="eps")

def main():

    fNameList = ["MenY3RU_PsiL4_0-200ns_every100frms"]

    for fName in fNameList:
        yList, xList = getData(fName)
        makeplot(yList, xList, fName)

main()