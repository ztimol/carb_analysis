import os
import numpy as np
from MDAnalysis.analysis.dihedrals import Dihedral
from trajectory import Trajectory
from torsion.torsion_angle_plot import TorsionAnglePlot


class BlockAverage(Trajectory):
    def __init__(self, env, mda_universe, torsion_angles_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.torsion_angles_dir = torsion_angles_dir

    def _set_torsion_angles_dir(self, torsion_name):
        torsion_name_dir = os.path.join(self.torsion_angles_dir, torsion_name)

        if not os.path.exists(torsion_name_dir):
            os.mkdir(torsion_name_dir)

        return torsion_name_dir

    def torsion_trajectory_analysis(self):



        
def blocked(universe, nblocks, analyze):
    size = universe.trajectory.numframes / nblocks
    blocks = []
    for block in xrange(nblocks):
        a = []
        for ts in u.trajectory[block * size:(block + 1) * size]:
            a.append(analyze(universe))
        blocks.append(np.average(a))
    blockaverage = np.average(blocks)
    blockstd = np.std(blocks)

    return nblocks, size, blockaverage, blockstd


def rgyr(universe):
    return universe.select_atoms('protein').radius_of_gyration()


def main():    
    u = MDAnalysis.Universe(PSF, DCD)
    results = []
    for num_of_blocks in xrange(2, 10):
        results.append(blocked(u, num_of_blocks, rgyr))
    r = np.array(results)
