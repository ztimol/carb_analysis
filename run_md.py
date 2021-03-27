import sys, os
import argparse
import subprocess as sp
import shutil

parser = argparse.ArgumentParser(description="Run MD")

parser.add_argument("-d", help="Select run directory.", type=str)
parser.add_argument(
    "-s",
    help="Create 'setup_files' directory with included parameter file.",
    default=False,
    type=str,
)
parser.add_argument("-em", help="Energy Minimisation", default=False, type=bool)
parser.add_argument("-eq", help="MD Equilibration", default=False, type=bool)
parser.add_argument("-solvate", help="Solvation", default=False, type=bool)
parser.add_argument(
    "-sol_box_size", help="Solvation Box Size", default=False, type=float
)

args = parser.parse_args()


def create_em_conf_file_from_template(em_dir, structure_name):

    current_dir = os.getcwd()
    namd_em_template_conf_file_path = os.path.join(
        current_dir, "templates", "namd", "em_template.conf"
    )

    psf_file_relative_path = "../structures/" + structure_name + ".psf"
    pdb_file_relative_path = "../structures/" + structure_name + ".pdb"

    namd_em_conf_file_path = os.path.join(em_dir, "em.conf")

    with open(namd_em_template_conf_file_path, "r") as namd_em_template_file:
        x = namd_em_template_file.read()

        x = x.replace("__PSF_FILE_PATH_TOKEN__", psf_file_relative_path)

        x = x.replace("__PDB_FILE_PATH_TOKEN__", pdb_file_relative_path)

    with open(namd_em_conf_file_path, "w") as namd_em_file:
        namd_em_file.write(x)


def _do_em(em_dir):

    log_file_path = os.path.join(em_dir, "em.log")
    conf_file_path = os.path.join(em_dir, "em.conf")

    namd_energy_calc_command = [
        "/home/timol/.NAMD_2.13_Linux-x86_64-multicore/namd2",
        "+p8",
        conf_file_path,
    ]

    f = open(log_file_path, "w")
    sp.call(namd_energy_calc_command, stdout=f)
    f.close()


def energy_minimisation(base_dir, structure_name):

    setup_files_dir = os.path.join(base_dir, "setup_files")
    em_dir = os.path.join(setup_files_dir, "em")

    if not os.path.exists(em_dir):
        os.mkdir(em_dir)

    # structure_dir = os.path.join(setup_files_dir, "structure")

    create_em_conf_file_from_template(em_dir, structure_name)
    _do_em(em_dir)

    shutil.copy(
        os.path.join(em_dir, "run_output", "em.coor"),
        os.path.join(em_dir, "run_output", "em.pdb"),
    )


def create_solvation_file_from_template(solvation_dir, structure_name):

    current_dir = os.getcwd()
    solvate_template_file_path = os.path.join(
        current_dir, "templates", "namd", "solvate_template.tcl"
    )

    solvation_script_file_path = os.path.join(solvation_dir, "solvate.tcl")

    psf_file_relative_path = "../structures/" + structure_name + ".psf"
    pdb_file_relative_path = "../em/run_output/em.pdb"

    MIN = "-30"
    MAX = "30"

    with open(solvate_template_file_path, "r") as solvation_template:
        x = solvation_template.read()

        x = x.replace("__SOLVATION_DIR_TOKEN__", solvation_dir)
        x = x.replace("__PSF_FILE_PATH_TOKEN__", psf_file_relative_path)
        x = x.replace("__PDB_FILE_PATH_TOKEN__", pdb_file_relative_path)
        x = x.replace("__LOWER_BOUND_TOKEN__", MIN)
        x = x.replace("__UPPER_BOUND_TOKEN__", MAX)

    with open(solvation_script_file_path, "w") as solvation_script:
        solvation_script.write(x)


def _do_solvation(solvation_dir):

    solvation_data_file_path = os.path.join(solvation_dir, "solvate.dat")
    solvation_script = os.path.join(solvation_dir, "solvate.tcl")

    vmd_solvate_command = [
        "vmd",
        "-dispdev",
        "text",
        "-e",
        solvation_script,
    ]

    f = open(solvation_data_file_path, "w")
    sp.call(vmd_solvate_command, stdout=f)
    f.close()


def solvate_minimised_structure(base_dir, structure_name):

    setup_files_dir = os.path.join(base_dir, "setup_files")
    solvation_dir = os.path.join(setup_files_dir, "solvated")

    if not os.path.exists(solvation_dir):
        os.mkdir(solvation_dir)

    create_solvation_file_from_template(solvation_dir, structure_name)
    _do_solvation(solvation_dir)


def main():

    base_dir = "/home/timol/C6W/Studies/Dynamics/NAMD/Shigella/MD/Solution/7b_s_flexneri/7b_s_flexneri_3ru/simulations/"

    structure_name = "7b_s_flexneri_3ru"

    energy_minimisation(base_dir, structure_name)
    solvate_minimised_structure(base_dir, structure_name)


if __name__ == "__main__":
    main()
