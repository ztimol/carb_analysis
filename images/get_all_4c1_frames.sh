import MDAnalysis as mda
import pandas as pd


PSF_FILE = "/home/timol/C6W/Studies/Dynamics/NAMD/Shigella/MD/Solution/7a_s_flexneri/7a_s_flexneri_6ru/ring_restrained/trajectories/7a_s_flexneri_6ru.psf"

DCD_FILE = "/home/timol/C6W/Studies/Dynamics/NAMD/Shigella/MD/Solution/7a_s_flexneri/7a_s_flexneri_6ru/ring_restrained/trajectories/7a_s_flexneri_6ru_ring_restrained_0-1000ns.dcd"

OUTFILE = "/home/timol/C6W/Studies/Dynamics/NAMD/Shigella/MD/Solution/7a_s_flexneri/7a_s_flexneri_6ru/ring_restrained/trajectories/7a_s_flexneri_6ru_ring_restrained_4C1_only_0-1000ns.dcd"

INFILE_PATH = "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru_ring_restrained/ring_pucker/all_theta.dat"

df = pd.read_csv(INFILE_PATH,sep="\t")

del df["Frame"]

df4c1 = df[df < 45].dropna()

frames_4c1 = df4c1[4000:].index

universe = mda.Universe(PSF_FILE, DCD_FILE)

carb_all = universe.select_atoms("all")
with mda.Writer(OUTFILE, carb_all.n_atoms) as W:
    for ts in universe.trajectory:
        if ts.frame in frames_4c1:
            W.write(carb_all)


