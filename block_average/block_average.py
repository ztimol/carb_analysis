import os
import numpy as np
from MDAnalysis.analysis import rms
from trajectory import Trajectory
import math

# from torsion.torsion_angle_plot import TorsionAnglePlot

import matplotlib

matplotlib.use("agg")  # no interactive plotting, only save figures
from pylab import errorbar, subplot, xlabel, ylabel, savefig, plot


class BlockAverage(Trajectory):
    def __init__(self, env, mda_universe, block_ave_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.block_ave_dir = block_ave_dir

    def _set_block_ave_dir(self, torsion_name):
        block_name_dir = os.path.join(self.block_ave_dir, torsion_name)

        if not os.path.exists(block_name_dir):
            os.mkdir(block_name_dir)

        return block_name_dir

    def block_average_analysis(self):
        results = []
        for nblocks in range(1, 10):
            print(nblocks)
            results.append(self.blocked(nblocks))

        result = np.array(results)
        self.write_block_averages(result)
        self.plot(result)

    def blocked(self, nblocks):
        size = int(self.mda_universe.trajectory.n_frames / nblocks)
        blocks = []

        # self.mda_universe.trajectory[0]  # first frame
        # first_frame_positions = self.mda_universe.select_atoms("all").positions

        for block in range(nblocks):
            a = []
            for frame in self.mda_universe.trajectory[
                block * size : (block + 1) * size
            ]:
                # positions = self.mda_universe.select_atoms("all").positions
                # a.append(rms.rmsd(first_frame_positions, positions))
                a.append(self.mda_universe.select_atoms("all").radius_of_gyration())
            blocks.append(np.average(a))

        blockaverage = np.average(blocks)
        blockstd = np.std(blocks) / math.sqrt(len(blocks))
        return nblocks, size, blockaverage, blockstd

    def write_block_averages(self, result):
        pass

    def plot(self, result):
        # subplot(211)

        # errorbar(result[:, 0], result[:, 2], yerr=result[:, 3])
        plot(result[:, 0], result[:, 3])
        xlabel("block length")
        ylabel(r"$\langle R_{\rm{gyr}} \rangle$ ($\AA$)")
        savefig("./block_average/figures/blocks.png")

        print("Wrote ./block_average/figures/blocks.{{pdf,png}}".format(*vars()))


# def main():
#     u = MDAnalysis.Universe(PSF, DCD)
#     results = []
#     for num_of_blocks in xrange(2, 10):
#         results.append(blocked(u, num_of_blocks, rgyr))
#     r = np.array(results)
