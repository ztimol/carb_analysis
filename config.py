import os
import re
from helper import is_whole_string_comment, is_empty_string, clean_string


def get_analysis_config(config_file):

    is_config_file_real = _check_if_config_file_exists(config_file)

    if is_config_file_real:
        return _read_config_params_from_file(config_file)
    else:
        input(config_file + " is not a valid file. Press enter to exit the program.")


def _read_config_params_from_file(config_file):

    config_params = {}
    config_params["torsions"] = {}
    with (open(config_file, "r")) as infile:
        for line in infile:
            if not is_whole_string_comment(line) and not is_empty_string(line):
                cleaned_line = clean_string(line)
                field_name = cleaned_line.split()[0]
                if field_name == "torsion":
                    config_params["torsions"] = _get_torsion_params(
                        cleaned_line, config_params["torsions"]
                    )
                elif field_name == "frames_per_ns":
                    config_params[field_name] = eval(line.split()[1])
                else:
                    config_params[field_name] = cleaned_line.split()[1]
    return config_params


def _check_if_config_file_exists(config_file):
    if os.path.isfile(config_file):
        return True
    else:
        return False


def _get_torsion_params(line, torsion_params):

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
        torsion_selection = line[line.find(torsion_type) + len(torsion_type) :].strip()[
            1:-1
        ]
        try:
            torsion_params["vars"] = {}
            torsion_params[torsion_name][torsion_type] = torsion_selection
        except KeyError:
            torsion_params["vars"][torsion_name] = {}
            torsion_params["vars"][torsion_name][torsion_type] = torsion_selection

    return torsion_params
