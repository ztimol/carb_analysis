import os
import MDAnalysis as mda
from trajectory import Trajectory
from torsion_angle import TorsionAngle


class Analysis(Trajectory):
    def __init__(self, env):
        self.env = env
        self.mda_universe = mda.Universe(
            self.get_psf_file_path(), self.get_dcd_file_path()
        )
        super().__init__(env)

    def torsion_analysis(self, env):
        start_frame = int(env["input_params"].get("start_frame", 0))
        output_dir = env["input_params"].get("output_dir", "output")

        try:
            input_torsion_params = env["input_params"]["torsions"]

            torsion_angles_dir = os.path.join(output_dir, "torsion_angles")
            if not os.path.exists(torsion_angles_dir):
                os.mkdir(torsion_angles_dir)
        except KeyError:
            return  # no torsions specified in config file. Don't run torsion analysis.

        torsion = TorsionAngle(self.env, self.mda_universe, torsion_angles_dir)
        torsion.torsion_trajectory_analysis(input_torsion_params, start_frame)
