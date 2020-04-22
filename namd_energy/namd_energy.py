import os
import numpy as np
import subprocess as sp
import shutil
from trajectory import Trajectory
from namd_energy.namd_energy_plot import NAMDEnergyPlot


class NAMDEnergy(NAMDEnergyPlot, Trajectory):
    def __init__(self, env, mda_universe, namd_energies_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.namd_energies_dir = namd_energies_dir

    def namd_single_point_energy_analysis(self):
        for namd_energy_name, energy_type in self.env["namd_energies"].items():

            print(
                "Performing single point energy calculation for "
                + namd_energy_name
                + " ..."
            )

            namd_energy_name_dir = os.path.join(
                self.namd_energies_dir, namd_energy_name
            )

            if not os.path.exists(namd_energy_name_dir):
                os.mkdir(namd_energy_name_dir)

            self._create_namd_energy_config_file()

            potential_energies = self._calculate_namd_single_point_energy(
                namd_energy_name_dir, namd_energy_name
            )

            time_series = np.arange(
                0, self.get_trajectory_time_in_ns(), self.ns_per_frame()
            )

            self.time_series_scatter(
                time_series, potential_energies, namd_energy_name_dir, namd_energy_name
            )

            self.probability_histogram(
                potential_energies, namd_energy_name_dir, namd_energy_name
            )

            print(
                "Completed single point energy calculation for "
                + namd_energy_name
                + "."
            )

    def _calculate_namd_single_point_energy(
        self, namd_energy_name_dir, namd_energy_name
    ):

        log_file_path = os.path.join(namd_energy_name_dir, namd_energy_name + ".log")
        data_out_file_path = os.path.join(
            namd_energy_name_dir, namd_energy_name + ".dat"
        )

        namd_energy_calc_command = [
            self.env["namd_path"],
            "+p2",
            "./namd_energy/md_energy.conf",
        ]

        f = open(log_file_path, "w")
        sp.call(namd_energy_calc_command, stdout=f)
        f.close()

        os.unlink("./namd_energy/md_energy.conf")

        potential_energies_per_frame = self._get_potential_energies_from_log_file(
            log_file_path
        )
        self._write_potential_energies(potential_energies_per_frame, data_out_file_path)
        return potential_energies_per_frame

    def _get_potential_energies_from_log_file(self, log_file_path):

        potential_energies_per_frame = []

        with open(log_file_path, "r") as lf:
            for line in lf:
                line = line.split()
                if line and line[0] == "ENERGY:":
                    potential_energy = eval(line[-3])
                    potential_energies_per_frame.append(potential_energy)

        return potential_energies_per_frame[:-1]  # last energy output is a duplicate

    def _write_potential_energies(
        self, potential_energies_per_frame, data_out_file_path
    ):

        with open(data_out_file_path, "w") as out_file:
            frame_number = 0
            for potential_energy in potential_energies_per_frame:
                out_file.write(str(frame_number))
                out_file.write(" ")
                out_file.write(str(potential_energy))
                out_file.write("\n")
                frame_number += 1

    def _create_namd_energy_config_file(self):

        if self.is_amber_mode_enabled():
            self._create_namd_energy_amber_config_file()
        else:
            self._create_namd_energy_charmm_config_file()

    def _create_namd_energy_charmm_config_file(self):

        with open(
            "./namd_energy/md_energy_template.conf", "r"
        ) as namd_energy_template_file:

            namd_energy_config = namd_energy_template_file.read()

            namd_energy_config = namd_energy_config.replace(
                "__PSF_FILE_TOKEN__", self.get_psf_file_path()
            )
            namd_energy_config = namd_energy_config.replace(
                "__PDB_FILE_TOKEN__", self.get_pdb_file_path()
            )
            namd_energy_config = namd_energy_config.replace(
                "__DCD_TRAJCTORY_FILE_TOKEN__", self.get_dcd_file_path()
            )

        with open("./namd_energy/md_energy.conf", "w") as namd_energy_config_file:
            namd_energy_config_file.write(namd_energy_config)

            namd_energy_config_file_path = os.path.join(
                self.namd_energies_dir, "md_energy.conf"
            )

        shutil.copy("./namd_energy/md_energy.conf", namd_energy_config_file_path)

    def _create_namd_energy_amber_config_file(self):
        with open(
            "./namd_energy/md_energy_amber_template.conf"
        ) as namd_energy_amber_template_file:

            namd_energy_config = namd_energy_amber_template_file.read()

            namd_energy_config = namd_energy_config.replace(
                "__PARM7_FILE_TOKEN__", self.get_parm7_file_path()
            )
            namd_energy_config = namd_energy_config.replace(
                "__RST7_FILE_TOKEN__", self.get_rst7_file_path()
            )

            namd_energy_config = namd_energy_config.replace(
                "__DCD_TRAJCTORY_FILE_TOKEN__", self.get_dcd_file_path()
            )

        with open("./namd_energy/md_energy.conf", "w") as namd_energy_config_file:
            namd_energy_config_file.write(namd_energy_config)

            namd_energy_config_file_path = os.path.join(
                self.namd_energies_dir, "md_energy.conf"
            )

        shutil.copy("./namd_energy/md_energy.conf", namd_energy_config_file_path)
