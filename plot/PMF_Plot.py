"""
 Simple python3 script to create a dihedral angle plot overlaid with a PMF contour map
 Adapted from contour map script written by Michelle Kuttel (2019)
 Written by Ryan Lazar (LZRRYA001) December 2019
"""
#
#
#
import matplotlib.pyplot as plt
from matplotlib import cm
import scipy.ndimage

#### !!! EDIT THIS PART !!! ####

# Input and output file paths
PMF_PATH = "/Users/RyanLazar/Dropbox/UCT/Masters/MetaDynamics/MetaD_aLRha12aLRha_GLY.pmf"
# The paths to your phi and psi angle files, extracted using VMD, corresponding to the structure you ran the pmf on here
PHI_ANGLES_PATH = "/Users/RyanLazar/Dropbox/UCT/Masters/Molecular_Dynamics/Dimers/aLRha12aLRha_phi.txt"
PSI_ANGLES_PATH = "/Users/RyanLazar/Dropbox/UCT/Masters/Molecular_Dynamics/Dimers/aLRha12aLRha_phi.txt"
# Saves the output figure to the below file path NB! dont forget the end '/'
OUTPUT_PATH = "/Users/RyanLazar/Dropbox/UCT/Masters/CHARMM_vs_GLYCAM/2_Simulations/Molecular_Dynamics/Dimers/Outputs/"
PLOT_NAME = "aLRha12aLRha_GLY_PMF" # Create a name for the plot and the output file

# Plot colouring
PLOT_COLOUR = 'purple'  # see other valid colours on matplotlib's guide (https://matplotlib.org/2.0.2/examples/color/named_colors.html)
COLOUR_MAP_GRADIENT = cm.coolwarm  # cm.gray looks better for overlay see more options at https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html

# Font sizes and tick frequency
LABEL_FONT_SIZE = 30
TICK_FONT_SIZE = 15
Y_TICKS = (-120, -60, 0, 60, 120)
X_TICKS = (-120, -60, 0, 60, 120)

#### SCRIPT ####

linkage = [PHI_ANGLES_PATH, PSI_ANGLES_PATH ]

x = []
y = []
z = []

#
phi = []
psi = []

phiInput = linkage[0]
psiInput = linkage[1]
    # Extract all phi, line by line
with open(phiInput) as fp:
    line = fp.readline()
    while line:
        sline = line.split()
        phi.append(float(sline[1]))
        line = fp.readline()

    # Extract all psi, line by line
with open(psiInput) as fp:
    line = fp.readline()
    while line:
        sline = line.split()
        psi.append(float(sline[1]))
        line = fp.readline()
    # Plot all points phi and psi

# process PMF file into 2D arrays
with open(PMF_PATH) as f:
    lines = f.readlines()
    xline = []
    yline = []
    zline = []
    for line in lines:
        line = line.strip()
        tmp = line.split()
        if len(tmp) and line[0] != '#' and len(line) >= 3:
            print("Line is:'", tmp)
            nextX = float(line.split()[0])
            if xline and (xline[-1] != nextX):  # end of row
                x.append(xline)
                y.append(yline)
                z.append(zline)
                xline = []
                yline = []
                zline = []
            nextY = float(line.split()[1])
            nextZ = float(line.split()[2])

            xline.append(nextX)
            yline.append(nextY)
            zline.append(nextZ)

    x.append(xline)  # do last append
    y.append(yline)
    z.append(zline)

# print(x,y,z)
# smoothing data
z = scipy.ndimage.gaussian_filter(z, sigma=1)  # bigger sigma = more smoothing; can go <1

fig = plt.figure()
a = fig.add_subplot(1, 1, 1)

levs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]  # levels to draw contours at
CS = a.contour(x, y, z, levs, cmap=COLOUR_MAP_GRADIENT)
a.scatter(phi, psi, s=0.3, color=PLOT_COLOUR)

# plt.clabel(CS, inline=1, fontsize=10)


###labels of contour lines better
# Define a class that forces representation of float to look a certain way
# This remove trailing zero so '1.0' becomes '1'
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()


# Recast levels to new class
CS.levels = [nf(val) for val in CS.levels]
# only every second level
# CS.levels =CS.levels[1::2]

# Label levels with specially formatted floats
# if plt.rcParams["text.usetex"]:
#    fmt = r'%r \%%'
# else:
#    fmt = '%r %%'
plt.clabel(CS, CS.levels, inline=True, fmt='%r ', fontsize=8)

# make a colorbar for the contour lines
CB = plt.colorbar(CS, shrink=0.9, extend='both')

# smoothing of contour lines


# Label for graph and axes
plt.title(PLOT_NAME)
plt.gca().set_aspect('equal', 'box')
plt.axis((-180,180,-180,180))
a.set_xlabel('$\phi$', fontsize=LABEL_FONT_SIZE)
h = a.set_ylabel('$\psi$', fontsize=LABEL_FONT_SIZE)
h.set_rotation(0)
a.set_xticks(X_TICKS)
a.set_yticks(Y_TICKS)

for tick in a.xaxis.get_major_ticks():
            tick.label.set_fontsize(TICK_FONT_SIZE)
for tick in a.yaxis.get_major_ticks():
            tick.label.set_fontsize(TICK_FONT_SIZE)


# a.set_zlabel('E ($kcal.mol^{-1}$)')

plt.savefig(OUTPUT_PATH + PLOT_NAME+'.png', dpi=300, bbox_inches='tight')
# pp = PdfPages(sys.argv[2]+'.png')
# pp.savefig()
# pp.close()

# mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
plt.show()
# save to file


