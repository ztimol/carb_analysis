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
                        try:
                            config_params["torsions"] = self._get_torsion_params(
                                cleaned_line, config_params["torsions"]
                            )
                        except KeyError:
                            config_params["torsions"] = {}
                            config_params["torsions"] = self._get_torsion_params(
                                cleaned_line, config_params["torsions"]
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

    def _get_torsion_params(self, line, torsion_params):

        torsion_input_values = line.split()
        torsion_name = torsion_input_values[1].strip()
        torsion_type = torsion_input_values[2].strip()

        if torsion_type == "scatter":
            try:
                torsion_params["plots"][torsion_name] = {
                    "x_key": line[3].strip(),
                    "y_key": line[4].strip(),
                }
            except KeyError:
                torsion_params["plots"] = {}
                torsion_params["plots"][torsion_name] = {
                    "x_key": torsion_input_values[3].strip(),
                    "y_key": torsion_input_values[4].strip(),
                }
        else:
            torsion_selection = line[
                line.find(torsion_type) + len(torsion_type) :
            ].strip()[1:-1]
            try:
                torsion_params["vars"] = {}
                torsion_params[torsion_name][torsion_type] = torsion_selection
            except KeyError:
                torsion_params["vars"][torsion_name] = {}
                torsion_params["vars"][torsion_name][torsion_type] = torsion_selection

        return torsion_params

    def get_dcd_file_path(self):
        return self.env.get("dcd_file", None)

    def get_psf_file_path(self):
        return self.env.get("psf_file", None)

    def get_start_frame(self):
        return int(self.env["start_frame"])

    # trajectory frames per nanoseconds
    def get_frames_per_ns(self):
        return int(self.env["frames_per_ns"])

    def get_torsion_x_axis_key(torsion_name):
        return
