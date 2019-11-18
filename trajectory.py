import MDAnalysis as mda
from helper import confirm_critical_file_exists
from config import Config


class Trajectory(Config):
    def __init__(self, env):
        super().__init__(env)
        self.dcd_file_path = self.get_dcd_file_path()
        self.psf_file_path = self.get_psf_file_path()

        # we check for dcd and psf file existence here as opposed to a try,
        # except later in the script so that we catch the error earlier
        # confirm_critical_file_exists(dcd_file_path)
        # confirm_critical_file_exists(psf_file_path)

    def get_trajectory_time_in_ns(self):
        return self.mda_universe.trajectory.n_frames / self.get_frames_per_ns()

    # nanoseconds per frame
    def ns_per_frame(self):
        return self.get_trajectory_time_in_ns() / self.mda_universe.trajectory.n_frames
