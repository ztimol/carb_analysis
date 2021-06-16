import os
import numpy as np
from MDAnalysis.analysis import rms
from trajectory import Trajectory
import math
import time

# from torsion.torsion_angle_plot import TorsionAnglePlot


class BlockAverage(Trajectory):
    def __init__(self, env, mda_universe, block_ave_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.block_ave_dir = block_ave_dir

    def block_average_analysis(self):
        block_average_data = []  # initialise with origin
        print("start")
        for nblocks in range(1, 20):
            print(nblocks)
            start = time.time()
            block_average_data.append(self.blocked(nblocks))
            end = time.time()
            print(end - start)

        block_calc_data = np.array(block_average_data)

        data_out_file_path = os.path.join(self.block_ave_dir, "block_average.dat")

        self._write_block_averages(block_calc_data, data_out_file_path)
        self.plot(block_calc_data)

    def blocked(self, nblocks):
        size = int(self.mda_universe.trajectory.n_frames / nblocks)
        blocks = []

        for block in range(nblocks):
            a = []
            for frame in self.mda_universe.trajectory[
                block * size : (block + 1) * size
            ]:
                a.append(self.mda_universe.select_atoms("all").radius_of_gyration())

            blocks.append(np.average(a))

        blockaverage = np.average(blocks)
        blockstd = np.std(blocks)  # / math.sqrt(len(blocks))
        block_length = size / self.get_frames_per_ns()
        return block_length, size, blockaverage, blockstd

    # def blocked(universe, nblocks, analyze):
    #     size = universe.trajectory.numframes/nblocks
    #     blocks = []
    #     for block in xrange(nblocks):
    #         a = []

    #         for ts in u.trajectory[block*size:(block+1)*size]:
    #             a.append(analyze(universe))
    #         blocks.append(numpy.average(a))

    #         blockaverage = numpy.average(blocks)
    #         blockstd = numpy.std(blocks)
    #     return nblocks, size, blockaverage, blockstd

    def _write_block_averages(self, block_calc_data, data_out_file_path):
        block_lengths, block_sizes, block_averages, block_std_errors = (
            np.flip(block_calc_data[:, 0]),
            np.flip(block_calc_data[:, 1]),
            block_calc_data[:, 2],
            block_calc_data[:, 3],
        )

        with open(data_out_file_path, "w") as out_file:
            for block_length, block_size, block_average, block_std_error in zip(
                block_lengths, block_sizes, block_averages, block_std_errors
            ):
                line = "{block_length} {block_size} {block_average} {block_std_error}\n".format(
                    block_length=block_length,
                    block_size=block_size,
                    block_average=block_average,
                    block_std_error=block_std_error,
                )
                out_file.write(line)

    # def plot(self, result):
    #     plot(np.flip(result[:, 0]), result[:, 3])
    #     xlabel("block length")
    #     ylabel(r"Block Standard Error $\angle R_{\rm{gyr}} \rangle$ ($\AA$)")
    #     savefig("./block_average/figures/blocks.png")

    #     print("Wrote ./block_average/figures/blocks.{{pdf,png}}".format(*vars()))
