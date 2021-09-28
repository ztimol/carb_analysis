import pandas as pd
import os

BASE_DIR = "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/torsion_angles/"

STRUCTURE_NAME = "aDGlc12aDGlc"

FILE_NAMES = [
    "{structure_name}_ru2/phi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru3/phi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru4/phi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru5/phi.dat".format(structure_name=STRUCTURE_NAME),
]

OUTFILE = os.path.join(BASE_DIR, "aDGlc12aDGlc_torsion_ru2-ru5.dat")

initial_file = FILE_NAMES[0]

phi_df = pd.read_csv(
    os.path.join(BASE_DIR, initial_file), sep=" ", names=["frame", "phi"]
)["phi"]
phi_df = phi_df[phi_df.index > 3999]

for file_name in FILE_NAMES[1:]:
    file_path = os.path.join(BASE_DIR, file_name)
    df1 = pd.read_csv(file_path, sep=" ", names=["frame", "phi"])["phi"]
    df1 = df1[df1.index > 3999]
    phi_df = phi_df.append(df1, ignore_index=True)

FILE_NAMES = [
    "{structure_name}_ru2/psi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru3/psi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru4/psi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru5/psi.dat".format(structure_name=STRUCTURE_NAME),
]

initial_file = FILE_NAMES[0]

psi_df = pd.read_csv(
    os.path.join(BASE_DIR, initial_file), sep=" ", names=["frame", "psi"]
)["psi"]

psi_df = psi_df[psi_df.index > 3999]

for file_name in FILE_NAMES[1:]:
    file_path = os.path.join(BASE_DIR, file_name)
    df1 = pd.read_csv(file_path, sep=" ", names=["frame", "psi"])["psi"]
    df1 = df1[df1.index > 3999]
    psi_df = psi_df.append(df1, ignore_index=True)

df = pd.DataFrame()
df["phi"] = phi_df
df["psi"] = psi_df

df.to_csv(OUTFILE, sep=" ", header=False, index_label=False)
