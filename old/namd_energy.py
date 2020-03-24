import argparse

parser = argparse.ArgumentParser(description="Extract potential energies from a namd log file")

parser.add_argument("-out", help="Select log file.")  # output data file
parser.add_argument("-log", help="Select log file.")  # input namd log file

args = parser.parse_args()


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


def main():

    potential_energies = get_potential_energies_from_log_file()
    write_potential_energies(potential_energies)


if __name__ == "__main__":
    main()
