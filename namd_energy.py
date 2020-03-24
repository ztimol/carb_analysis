class NAMDEnergy:
    def __init__(self, env, mda_universe, torsion_angles_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.torsion_angles_dir = torsion_angles_dir

    def calculate_namd_single_point_energy():
        pass

    def get_potential_energies_from_log_file():

        potential_energies = []

        with open(args.log, "r") as log_file:
            for line in log_file:
                line = line.split()
                if line and line[0] == "ENERGY:":
                    potential_energy = eval(line[-3])
                    potential_energies.append(potential_energy)

        return potential_energies

    def write_potential_energies(potential_energies):

        with open(args.out, "w") as out_file:
            frame_number = 0
            for potential_energy in potential_energies:
                out_file.write(str(frame_number))
                out_file.write(" ")
                out_file.write(str(potential_energy))
                out_file.write("\n")
                frame_number += 1
