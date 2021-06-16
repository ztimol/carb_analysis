import numpy as np
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Create a scatter plot")

parser.add_argument("-f", help="Select data file.")

args = parser.parse_args()


def get_outfile_name(fName):

    try:
        x = fName.split(".")[:-1]
        outfileName = "".join(x) + ".png"
    except:
        outfileName = fName + ".png"

    return outfileName


def get_cremer_pople_parameters(env_params):

    trajectory_times_in_ns = []
    theta_per_ns = []
    phi_per_ns = []
    q_per_ns = []

    inf = open(args.f, "r")
    for line in inf:
        line = line.split(",")
        time_in_ns = eval(line[0])
        if time_in_ns >= 0:
            q = eval(line[1])
            theta = eval(line[2])
            phi = eval(line[3])
            trajectory_times_in_ns.append(time_in_ns)
            q_per_ns.append(q)
            theta_per_ns.append(theta)
            phi_per_ns.append(phi)
            print(theta, time_in_ns)

    inf.close()

    env_params["trajectory_times_in_ns"] = trajectory_times_in_ns
    env_params["theta_per_ns"] = theta_per_ns
    env_params["phi_per_ns"] = phi_per_ns
    env_params["q_per_ns"] = q_per_ns

    return


def make_polar_scatter(phi, r, file_name):

    fig = plt.figure()
    ax = fig.add_subplot(111, projection="polar")
    # ax.scatter(phi, r, c=colors, s=area, cmap="hsv", alpha=0.75)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_rlim(top=180, bottom=0)
    ax.scatter(phi, r, c=r, s=1, cmap="hsv", alpha=0.75)

    fig.savefig(get_outfile_name(file_name), dpi=400, format="png")


def main():
    env_params = {}
    get_cremer_pople_parameters(env_params)

    make_polar_scatter(env_params["theta_per_ns"], env_params["phi_per_ns"], args.f)


if __name__ == "__main__":
    main()
