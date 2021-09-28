import pandas as pd
import os

BASE_DIR = (
    "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/ring_pucker/"
)

STRUCTURE_NAME = "bDGlcNAc"

FILE_NAMES = [
    "{structure_name}_ru2/trajectory_cp_phi_theta_Q.dat".format(
        structure_name=STRUCTURE_NAME
    ),
    "{structure_name}_ru3/trajectory_cp_phi_theta_Q.dat".format(
        structure_name=STRUCTURE_NAME
    ),
    "{structure_name}_ru4/trajectory_cp_phi_theta_Q.dat".format(
        structure_name=STRUCTURE_NAME
    ),
    "{structure_name}_ru5/trajectory_cp_phi_theta_Q.dat".format(
        structure_name=STRUCTURE_NAME
    ),
]

OUTFILE = os.path.join(BASE_DIR, "bDGlcNAc_ring_pucker_ru2-ru5.dat")

initial_file = FILE_NAMES[0]

df = pd.read_csv(
    os.path.join(BASE_DIR, initial_file), sep=" ", names=["frame", "phi", "theta", "Q"]
)
df = df[df.index > 3999]

for file_name in FILE_NAMES[1:]:
    file_path = os.path.join(BASE_DIR, file_name)
    df1 = pd.read_csv(
        os.path.join(BASE_DIR, file_path), sep=" ", names=["frame", "phi", "theta", "Q"]
    )
    df1 = df1[df1.index > 3999]
    df = df.append(df1, ignore_index=True)

del df["frame"]
df.to_csv(OUTFILE, sep=" ", header=False, index_label=False)
