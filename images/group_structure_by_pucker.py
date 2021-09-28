import pandas as pd
import os

BASE_DIR = "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/torsion_angles/"

PUCKER_DIR = (
    "/home/timol/C6W/Studies/structure_analysis/output/7a_s_flexneri_6ru/ring_pucker/"
)

STRUCTURE_NAME = "bDGlcNAc12aLRha"

PUCKER_FILE_NAMES = [
    "GlcNAc_ru2/trajectory_cp_phi_theta_Q.dat".format(structure_name=STRUCTURE_NAME),
    "GlcNAc_ru3/trajectory_cp_phi_theta_Q.dat".format(structure_name=STRUCTURE_NAME),
    "GlcNAc_ru4/trajectory_cp_phi_theta_Q.dat".format(structure_name=STRUCTURE_NAME),
    "GlcNAc_ru5/trajectory_cp_phi_theta_Q.dat".format(structure_name=STRUCTURE_NAME),
]


FILE_NAMES = [
    "{structure_name}_ru2/phi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru3/phi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru4/phi.dat".format(structure_name=STRUCTURE_NAME),
    "{structure_name}_ru5/phi.dat".format(structure_name=STRUCTURE_NAME),
]

OUTFILE = os.path.join(
    BASE_DIR,
    "{structure_name}_torsion_ru2-ru5_4C1_and_boat_only.dat".format(
        structure_name=STRUCTURE_NAME
    ),
)

initial_file = FILE_NAMES[0]


# read torsions
phi_df = pd.read_csv(
    os.path.join(BASE_DIR, initial_file), sep=" ", names=["frame", "phi"]
)["phi"]
phi_df = phi_df[phi_df.index > 3999]

# read pucker
pucker_df = pd.read_csv(
    os.path.join(PUCKER_DIR, PUCKER_FILE_NAMES[0]),
    sep=" ",
    names=["frame", "phi", "theta", "Q"],
)["theta"]

df = pd.DataFrame()
df["phi"] = phi_df
df["theta"] = pucker_df

phi_df = df.query("theta < 120")  # .drop("theta", axis=1)


for file_name, pucker_file_name in zip(FILE_NAMES[1:], PUCKER_FILE_NAMES[1:]):
    file_path = os.path.join(BASE_DIR, file_name)
    df1 = pd.read_csv(file_path, sep=" ", names=["frame", "phi"])["phi"]
    df1 = df1[df1.index > 3999]

    pucker_df = pd.read_csv(
        os.path.join(PUCKER_DIR, pucker_file_name),
        sep=" ",
        names=["frame", "phi", "theta", "Q"],
    )["theta"]
    df = pd.DataFrame()
    df["phi"] = df1
    df["theta"] = pucker_df
    df1 = df.query("theta < 60")  # .drop("theta", axis=1)

    phi_df = phi_df.append(df1, ignore_index=True)


# ---------

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

pucker_df = pd.read_csv(
    os.path.join(PUCKER_DIR, PUCKER_FILE_NAMES[0]),
    sep=" ",
    names=["frame", "phi", "theta", "Q"],
)["theta"]

df = pd.DataFrame()
df["psi"] = psi_df
df["theta"] = pucker_df

psi_df = df.query("theta < 120")  # .drop("theta", axis=1)

for file_name, pucker_file_name in zip(FILE_NAMES[1:], PUCKER_FILE_NAMES[1:]):
    file_path = os.path.join(BASE_DIR, file_name)
    df1 = pd.read_csv(file_path, sep=" ", names=["frame", "psi"])["psi"]
    df1 = df1[df1.index > 3999]

    pucker_df = pd.read_csv(
        os.path.join(PUCKER_DIR, pucker_file_name),
        sep=" ",
        names=["frame", "phi", "theta", "Q"],
    )["theta"]
    df = pd.DataFrame()
    df["psi"] = df1
    df["theta"] = pucker_df
    df1 = df.query("theta < 60")  # .drop("theta", axis=1)

    psi_df = psi_df.append(df1, ignore_index=True)

# ------

df = pd.DataFrame()
df["phi"] = phi_df["phi"]
df["psi"] = psi_df["psi"]
# df["ring_theta"] = phi_df["theta"]

df.to_csv(OUTFILE, sep=" ", header=False, index_label=False)
