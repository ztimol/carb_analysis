import os
import math
from trajectory import Trajectory

# from ring_pucker.cp_ring_pucker import CPRingPucker


class CpPuckerEnergyPmf(Trajectory):
    # def __init__(self, env, mda_universe, ring_pucker_dir):
    def __init__(self):
        # self.mda_universe = mda_universe
        # self.env = env
        # self.ring_pucker_dir = ring_pucker_dir
        self.gas_constant_calories = 1.987
        # self.gas_constant_calories = 1.380649e9 - 23
        self.temp_kelvin = 300

    # def calc_pmf(self):
    #     ring_pucker = CPRingPucker(self.env, self.mda_universe, ring_pucker_dir)
    #     ring_pucker.cp_ring_pucker_analysis()

    # get cremer-pople puckering parameters from file
    def _read_trjectory_cp_puckering_parameters(self):

        cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13_aDGlc14_bDGlcNAc_glycam/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

        cp_trajectory_pucker_params = {}

        with open(cp_trajectory_parameter_file, "r") as infile:
            for line in infile:
                line = line.split()
                frame_number = eval(line[0])
                cp_phi = eval(line[1])
                cp_theta = eval(line[2])
                cp_q = eval(line[3])
                cp_trajectory_pucker_params[frame_number] = {
                    "cp_phi": cp_phi,
                    "cp_theta": cp_theta,
                    "cp_q": cp_q,
                }

        return cp_trajectory_pucker_params

    def _calc_theta_probability(self, cp_trajectory_pucker_params):

        rounded_trajectory_cp_theta_values = [
            round(cp_params["cp_theta"])
            for cp_params in cp_trajectory_pucker_params.values()
        ]

        number_of_frames = len(
            rounded_trajectory_cp_theta_values
        )  # self.mda_universe.trajectory.n_frames

        theta_probability = {}
        free_energy_per_theta_value = {}

        for theta_value in range(0, 181):
            theta_count = rounded_trajectory_cp_theta_values.count(theta_value)
            probability_of_theta = theta_count / number_of_frames
            theta_probability[theta_value] = probability_of_theta
            try:
                free_energy_per_theta_value[theta_value] = (
                    -self.gas_constant_calories
                    * self.temp_kelvin
                    * math.log(probability_of_theta)
                )
            except ValueError:
                free_energy_per_theta_value[theta_value] = -1
            import pdb

            pdb.set_trace()
