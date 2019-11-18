import MDAnalysis as mda
from helper import confirm_critical_file_exists


class Trajectory:
    def __init__(self, dcd_file_path, psf_file_path):
        self.dcd_file_path = dcd_file_path
        self.psf_file_path = psf_file_path

        # we check for dcd and psf file existence here as opposed to a try,
        # except later in the script so that we catch the error earlier
        confirm_critical_file_exists(dcd_file_path)
        confirm_critical_file_exists(psf_file_path)
