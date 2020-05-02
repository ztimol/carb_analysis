import os
from helper import is_whole_string_comment, is_empty_string, clean_string


class Config:
    def __init__(self, env=None):
        self.env = env

    def get_analysis_config(self, config_file):

        is_config_file_real = self._check_if_config_file_exists(config_file)

        if is_config_file_real:
            return self._read_config_params_from_file(config_file)

        input(config_file + " is not a valid file. Press enter to exit the program.")
        return None

    def _read_config_params_from_file(self, config_file):

        config_params = {}
        with (open(config_file, "r")) as infile:
            for line in infile:
                if not is_whole_string_comment(line) and not is_empty_string(line):
                    cleaned_line = clean_string(line)
                    field_name = cleaned_line.split()[0]
                    if field_name == "torsion":
                        config_params = self._handle_torsion_field(
                            cleaned_line, config_params
                        )
                    elif field_name == "namd_energy":
                        config_params = self._handle_namd_energy_field(
                            cleaned_line, config_params
                        )
                    elif field_name == "ring_pucker":
                        config_params = self._handle_ring_pucker_field(
                            cleaned_line, config_params
                        )
                    elif field_name == "atom_distance":
                        config_params = self._handle_atom_distance_field(
                            cleaned_line, config_params
                        )
                    elif field_name == "frames_per_ns":
                        config_params[field_name] = eval(line.split()[1])
                    else:
                        config_params[field_name] = cleaned_line.split()[1]
        return config_params

    def _check_if_config_file_exists(self, config_file):
        if os.path.isfile(config_file):
            return True
        else:
            return False

    def _handle_torsion_field(self, cleaned_line, config_params):
        try:
            config_params["torsions"] = self._get_torsion_params(
                cleaned_line, config_params["torsions"]
            )
        except KeyError:
            config_params["torsions"] = {}
            config_params["torsions"] = self._get_torsion_params(
                cleaned_line, config_params["torsions"]
            )
        return config_params

    def _get_torsion_params(self, line, torsion_params):

        torsion_input_values = line.split()
        torsion_name = torsion_input_values[1].strip()
        torsion_type = torsion_input_values[2].strip()

        if "vars" not in torsion_params:
            torsion_params["vars"] = {}
        if "plots" not in torsion_params:
            torsion_params["plots"] = {}

        if torsion_type == "scatter":
            try:
                torsion_params["plots"][torsion_name] = {
                    "x_key": torsion_input_values[3].strip(),
                    "y_key": torsion_input_values[4].strip(),
                }

            except KeyError:
                torsion_params["plots"][torsion_name] = {
                    "x_key": torsion_input_values[3].strip(),
                    "y_key": torsion_input_values[4].strip(),
                }
        else:
            index_of_first_flag = -1  # line.find("-")
            index_of_torsion_type = line.find(torsion_type) + len(torsion_type)
            torsion_selection = (
                line[index_of_torsion_type:]
                .strip()[1:index_of_first_flag]  # - index_of_torsion_type]
                .replace("-", "")
                .replace("'", "")
                .replace('"', "")
                .strip()
            )  # extract atom selection from the line

            try:
                torsion_params["vars"][torsion_name][torsion_type] = torsion_selection
            except KeyError:
                torsion_params["vars"][torsion_name] = {}
                torsion_params["vars"][torsion_name][torsion_type] = torsion_selection
        return torsion_params

    def get_dcd_file_path(self):
        return self.env.get("dcd_file", None)

    def get_psf_file_path(self):
        return self.env.get("psf_file", None)

    def get_pdb_file_path(self):
        return self.env.get("pdb_file", None)

    def get_parm7_file_path(self):
        return self.env.get("parm7_file", None)

    def get_rst7_file_path(self):
        return self.env.get("rst7_file", None)

    def is_amber_mode_enabled(self):
        try:
            if self.env["amber"] == "yes":
                return True
            return False
        except KeyError:
            return False

    def get_start_frame(self):
        try:
            return int(self.env["start_frame"])
        except KeyError:
            return 0

    # trajectory frames per nanoseconds
    def get_frames_per_ns(self):
        return int(self.env["frames_per_ns"])

    def get_torsion_x_axis_key(torsion_name):
        return

    def _handle_namd_energy_field(self, cleaned_line, config_params):
        try:
            config_params["namd_energies"] = self._get_namd_energy_params(
                cleaned_line, config_params["namd_energies"]
            )
        except KeyError:
            config_params["namd_energies"] = {}
            config_params["namd_energies"] = self._get_namd_energy_params(
                cleaned_line, config_params["namd_energies"]
            )
        return config_params

    def _get_namd_energy_params(self, line, namd_energy_params):
        namd_energy_params[line.split()[1].strip()] = line.split()[2].strip()
        return namd_energy_params

    def _handle_ring_pucker_field(self, cleaned_line, config_params):
        try:
            config_params["ring_puckers"] = self._get_ring_pucker_params(
                cleaned_line, config_params["ring_puckers"]
            )
        except KeyError:
            config_params["ring_puckers"] = {}
            config_params["ring_puckers"] = self._get_ring_pucker_params(
                cleaned_line, config_params["ring_puckers"]
            )
        return config_params

    def _get_ring_pucker_params(self, line, ring_pucker_params):
        ring_pucker_name = line.split()[1].strip()
        ring_pucker_selection = line[
            line.find(ring_pucker_name) + len(ring_pucker_name) :
        ].strip()[1:-1]
        ring_pucker_params[ring_pucker_name] = ring_pucker_selection
        return ring_pucker_params

    def _handle_atom_distance_field(self, cleaned_line, config_params):
        try:
            config_params["atom_distances"].append(
                self._get_atom_distance_params(cleaned_line)
            )
        except KeyError:
            config_params["atom_distances"] = []
            config_params["atom_distances"].append(
                self._get_atom_distance_params(cleaned_line)
            )
        return config_params

    def _get_atom_distance_params(self, line):
        field_name = line.split()[0].strip()
        atom_distance_selection = line[
            line.find(field_name) + len(field_name) :
        ].strip()[1:-1]
        return atom_distance_selection
