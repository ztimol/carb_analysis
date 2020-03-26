import numpy as np
import argparse
import statistics
import math
import cp_polar_scatter

parser = argparse.ArgumentParser(description="Create a scatter plot")

parser.add_argument("-f", help="Select data file.", type=str)
parser.add_argument("-cp_calc", help="create plot", type=bool, default=True)
parser.add_argument("-polar_plot", help="create plot.", type=bool, default=False)

args = parser.parse_args()


def validateParser():

    if args.ylabel and args.ysym:
        input("Do not specify both ylabel and ysym flags.")
        raise SystemExit(0)

    if args.xlabel and args.xsym:
        input("Do not specify both ylabel and ysym flags.")
        raise SystemExit(0)


def get_outfile_name(fName):

    try:
        x = fName.split(".")[:-1]
        outfileName = "".join(x) + "_cp_params_time_series.dat"
    except:
        outfileName = fName + "_cp_params_time_series.png"

    return outfileName


def get_ring_atom_coordinates_per_frame(infile, env_params):

    trajectory_times_in_ns = []
    ring_atom_coordinates_per_frame = []
    ring_atom_x_coordinates_per_frame = []
    ring_atom_y_coordinates_per_frame = []
    ring_atom_z_coordinates_per_frame = []

    inf = open(infile, "r")
    i = 0
    for line in inf:
        line = line.split(":")
        frame_num = eval(line[0])
        ring_atom_coordinates_per_frame.append([])
        ring_atom_x_coordinates_per_frame.append([])
        ring_atom_y_coordinates_per_frame.append([])
        ring_atom_z_coordinates_per_frame.append([])

        if frame_num >= 0:
            time_in_ns = (frame_num * env_params["dcd_freq"] * 100) / 1000000.0
            trajectory_times_in_ns.append(time_in_ns)
            coordinates = line[1].split(",")

            for coordinate in coordinates:  # ignore new line character
                ring_atom_coordinates_per_frame[i].append(coordinate.split())
                x_coordinates = eval(coordinate.split()[0])
                y_coordinates = eval(coordinate.split()[1])
                z_coordinates = eval(coordinate.split()[2])

                ring_atom_x_coordinates_per_frame[i].append(x_coordinates)
                ring_atom_y_coordinates_per_frame[i].append(y_coordinates)
                ring_atom_z_coordinates_per_frame[i].append(z_coordinates)
            # assert (
            #     len(ring_atom_x_coordinates_per_frame[i])
            #     == len(ring_atom_y_coordinates_per_frame[i])
            #     == len(ring_atom_z_coordinates_per_frame[i]),
            #     "Coordinate error for frame: " + str(frame_num),
            # )

            i += 1

    inf.close()

    env_params["trajectory_times_in_ns"] = trajectory_times_in_ns
    env_params["ring_atom_coordinates_per_frame"] = ring_atom_coordinates_per_frame
    env_params["ring_atom_x_coordinates_per_frame"] = ring_atom_x_coordinates_per_frame
    env_params["ring_atom_y_coordinates_per_frame"] = ring_atom_y_coordinates_per_frame
    env_params["ring_atom_z_coordinates_per_frame"] = ring_atom_z_coordinates_per_frame


def write_cremer_pople_values(env_params):

    with open(get_outfile_name(args.f), "w") as outfile:
        for t_ns, Q, theta2, phi2 in zip(
            env_params["trajectory_times_in_ns"],
            env_params["Q_values"],
            env_params["theta2_deg"],
            env_params["phi2_deg"],
        ):
            outfile.write(
                "{t_ns}, {Q}, {theta2}, {phi2}\n".format(
                    t_ns=t_ns, Q=Q, theta2=theta2, phi2=phi2
                )
            )
    return


def calc_ring_geometric_center(env_params, frame_number):

    x_coors = env_params["ring_atom_x_coordinates_per_frame"][frame_number]
    y_coors = env_params["ring_atom_y_coordinates_per_frame"][frame_number]
    z_coors = env_params["ring_atom_z_coordinates_per_frame"][frame_number]
    x_coor_ave = statistics.mean(x_coors)
    y_coor_ave = statistics.mean(y_coors)
    z_coor_ave = statistics.mean(z_coors)
    ring_geometric_center_coord = [x_coor_ave, y_coor_ave, z_coor_ave]

    return ring_geometric_center_coord


def calc_atom_position_vectors_from_geometric_center(
    env_params, frame_number, ring_geometric_center_coord
):

    x_coors = env_params["ring_atom_x_coordinates_per_frame"][frame_number]
    y_coors = env_params["ring_atom_y_coordinates_per_frame"][frame_number]
    z_coors = env_params["ring_atom_z_coordinates_per_frame"][frame_number]

    number_of_atoms_in_ring = len(x_coors)

    ring_atom_position_vectors = []

    for atom_number in range(number_of_atoms_in_ring):
        x_coor = x_coors[atom_number]
        y_coor = y_coors[atom_number]
        z_coor = z_coors[atom_number]
        x_vector = x_coor - ring_geometric_center_coord[0]
        y_vector = y_coor - ring_geometric_center_coord[1]
        z_vector = z_coor - ring_geometric_center_coord[2]
        ring_atom_position_vectors.append([x_vector, y_vector, z_vector])

    return ring_atom_position_vectors


def calc_projection_value(atom_number, number_of_atoms_in_ring, m=1, method="cos"):

    a = (2 * math.pi * m * atom_number) / number_of_atoms_in_ring

    if method == "cos":
        return math.cos(a)
    elif method == "sin":
        return math.sin(a)

    return None


def calc_position_vector_projections(ring_atom_position_vectors, projection_type="x"):

    number_of_atoms_in_ring = len(ring_atom_position_vectors)
    projected_vectors = []

    for atom_number in range(number_of_atoms_in_ring):
        if projection_type == "x":
            projection_value = calc_projection_value(
                atom_number, number_of_atoms_in_ring
            )
        elif projection_type == "y":
            projection_value = calc_projection_value(
                atom_number, number_of_atoms_in_ring, method="sin"
            )

        x_coor = ring_atom_position_vectors[atom_number][0]
        y_coor = ring_atom_position_vectors[atom_number][1]
        z_coor = ring_atom_position_vectors[atom_number][2]

        projected_vectors.append(
            [
                x_coor * projection_value,
                y_coor * projection_value,
                z_coor * projection_value,
            ]
        )

    return projected_vectors


def sum_vectors(vectors):

    summed_x_coordinates_vectors = sum([vector[0] for vector in vectors])
    summed_y_coordinates_vectors = sum([vector[1] for vector in vectors])
    summed_z_coordinates_vectors = sum([vector[2] for vector in vectors])

    return [
        summed_x_coordinates_vectors,
        summed_y_coordinates_vectors,
        summed_z_coordinates_vectors,
    ]


def calc_cross_product_unit_vector(x_projected_vectors, y_projected_vectors):
    summed_x_projected_vector = sum_vectors(x_projected_vectors)
    summed_y_projected_vector = sum_vectors(y_projected_vectors)

    cross_product_of_x_against_y_projected_vectors = np.cross(
        summed_x_projected_vector, summed_y_projected_vector
    )

    maginitude_of_crossed_vectors = np.linalg.norm(
        cross_product_of_x_against_y_projected_vectors
    )

    unit_vector_of_cross_product = (
        cross_product_of_x_against_y_projected_vectors / maginitude_of_crossed_vectors
    )

    return unit_vector_of_cross_product


def calc_ring_atom_z_values(
    ring_atom_position_vectors, x_projected_vectors, y_projected_vectors
):

    unit_vector_of_cross_product = calc_cross_product_unit_vector(
        x_projected_vectors, y_projected_vectors
    )

    number_of_atoms_in_ring = len(ring_atom_position_vectors)

    ring_atom_z_values = []

    for atom_number in range(number_of_atoms_in_ring):
        position_vector = ring_atom_position_vectors[atom_number]
        atom_z_value = np.dot(position_vector, unit_vector_of_cross_product)
        ring_atom_z_values.append(atom_z_value)

    return ring_atom_z_values


def calc_cremer_pople_Q_value(ring_atom_z_values):
    return math.sqrt(sum([z_value ** 2 for z_value in ring_atom_z_values]))


def calc_cremer_pople_q2_projection(
    atom_number, number_of_atoms_in_ring, projection_type="x"
):
    if projection_type == "x":
        return calc_projection_value(atom_number, number_of_atoms_in_ring, m=2)
    elif projection_type == "y":
        return calc_projection_value(
            atom_number, number_of_atoms_in_ring, m=2, method="sin"
        )


def calc_cremer_pople_q2_and_phi(ring_atom_z_values):

    number_of_atoms_in_ring = len(ring_atom_z_values)
    ring_atoms_zj_cos_phi = []
    ring_atoms_zj_sin_phi = []

    for atom_number in range(number_of_atoms_in_ring):
        x_projection_value = calc_cremer_pople_q2_projection(
            atom_number, number_of_atoms_in_ring, projection_type="x"
        )
        y_projection_value = calc_cremer_pople_q2_projection(
            atom_number, number_of_atoms_in_ring, projection_type="y"
        )
        ring_atoms_zj_cos_phi.append(
            ring_atom_z_values[atom_number] * x_projection_value
        )
        ring_atoms_zj_sin_phi.append(
            ring_atom_z_values[atom_number] * y_projection_value
        )

    summed_ring_atoms_zj_cos_phi = sum(ring_atoms_zj_cos_phi)
    summed_ring_atoms_zj_sin_phi = sum(ring_atoms_zj_sin_phi)

    q2_cos_phi2 = math.sqrt(2 / number_of_atoms_in_ring) * summed_ring_atoms_zj_cos_phi
    q2_sin_phi2 = (
        -1 * math.sqrt(2 / number_of_atoms_in_ring) * summed_ring_atoms_zj_sin_phi
    )

    q2 = calc_cremer_pople_q2(q2_cos_phi2, q2_sin_phi2)
    phi2 = calc_cremer_pople_phi(q2_cos_phi2, q2_sin_phi2)
    return q2, phi2


def calc_cremer_pople_q2(q2_cos_phi2, q2_sin_phi2):
    return math.sqrt(q2_cos_phi2 ** 2 + q2_sin_phi2 ** 2)


def calc_cremer_pople_phi(q2_cos_phi2, q2_sin_phi2):
    tan_phi2 = q2_sin_phi2 / q2_cos_phi2

    # not full sure why the '+ 180' is needed but it is.
    return math.atan(tan_phi2) + math.pi


def calc_cremer_pople_q3(ring_atom_z_values):

    number_of_atoms_in_ring = len(ring_atom_z_values)
    c = sum(
        [
            ring_atom_z_values[atom_number] * math.cos(math.pi * atom_number)
            for atom_number in range(number_of_atoms_in_ring)
        ]
    )

    q3 = (number_of_atoms_in_ring ** -0.5) * c

    return q3


def calc_cremer_pople_theta(q2, q3):
    tan_theta = q2 / q3

    if tan_theta < 0:
        return abs(math.atan(round(tan_theta, 5)))
    elif tan_theta > 0:
        return math.pi - math.atan(round(tan_theta, 5))


def calc_cremer_pople_ring_values_per_frame(env_params, frame_number):

    ring_geometric_center_coord = calc_ring_geometric_center(env_params, frame_number)

    ring_atom_position_vectors = calc_atom_position_vectors_from_geometric_center(
        env_params, frame_number, ring_geometric_center_coord
    )

    x_projected_vectors = calc_position_vector_projections(
        ring_atom_position_vectors, projection_type="x"
    )

    y_projected_vectors = calc_position_vector_projections(
        ring_atom_position_vectors, projection_type="y"
    )

    ring_atom_z_values = calc_ring_atom_z_values(
        ring_atom_position_vectors, x_projected_vectors, y_projected_vectors
    )

    cremer_pople_Q_value = calc_cremer_pople_Q_value(ring_atom_z_values)

    cremer_pople_q2, cremer_pople_phi2 = calc_cremer_pople_q2_and_phi(
        ring_atom_z_values
    )

    q3 = calc_cremer_pople_q3(ring_atom_z_values)

    cremer_pople_theta = calc_cremer_pople_theta(cremer_pople_q2, q3)

    cremer_pople_phi2_deg = (cremer_pople_phi2 * 180) / math.pi
    cremer_pople_theta_deg = (cremer_pople_theta * 180) / math.pi

    env_params["Q_values"].append(cremer_pople_Q_value)
    env_params["q2"].append(cremer_pople_q2)
    env_params["phi2_deg"].append(cremer_pople_phi2_deg)
    env_params["theta2_deg"].append(cremer_pople_theta_deg)


def calc_all_cremer_pople_values(env_params):

    for frame_number in range(env_params["number_of_frames"]):
        calc_cremer_pople_ring_values_per_frame(env_params, frame_number)


def main():

    env_params = {"dcd_freq": 250}

    env_params["Q_values"] = []
    env_params["q2"] = []
    env_params["q3"] = []
    env_params["theta2_deg"] = []
    env_params["phi2_deg"] = []

    # validateParser()
    get_ring_atom_coordinates_per_frame(args.f, env_params)

    env_params["number_of_frames"] = len(env_params["trajectory_times_in_ns"])

    if args.cp_calc:
        calc_all_cremer_pople_values(env_params)
        write_cremer_pople_values(env_params)

    if args.polar_plot:
        cp_polar_scatter.make_polar_scatter(
            env_params["phi2_deg"], env_params["theta2_deg"], args.f
        )


if __name__ == "__main__":
    main()
