import os
import math, statistics
import numpy as np
from trajectory import Trajectory


class CPRingPucker(Trajectory):
    def __init__(self, env, mda_universe, ring_pucker_dir):
        self.mda_universe = mda_universe
        self.env = env
        self.ring_pucker_dir = ring_pucker_dir

    def cp_ring_pucker_analysis(self):

        for ring_pucker_name, ring_pucker_atom_selection in self.env[
            "ring_puckers"
        ].items():
            self._calc_cremer_pople_ring_values_per_frame(
                ring_pucker_name, ring_pucker_atom_selection
            )

    def _calc_cremer_pople_ring_values_per_frame(
        self, ring_pucker_name, ring_pucker_atom_selection
    ):

        cremer_pople_pucker_values_per_frame = {}

        for frame in self.mda_universe.trajectory:
            ring_center_of_geo = self.mda_universe.select_atoms(
                ring_pucker_atom_selection
            ).center_of_geometry()

            ring_atom_position_vectors = self._calc_ring_atom_position_vectors_from_geometric_center(
                ring_center_of_geo, ring_pucker_atom_selection
            )

            x_projected_vectors = self._calc_ring_atom_position_vector_projections(
                ring_atom_position_vectors
            )

            y_projected_vectors = self._calc_ring_atom_position_vector_projections(
                ring_atom_position_vectors, projection_type="y"
            )

            ring_atom_z_values = self._calc_ring_atom_z_values(
                ring_atom_position_vectors, x_projected_vectors, y_projected_vectors
            )

            cremer_pople_Q_value = self._calc_cremer_pople_Q_value(ring_atom_z_values)

            cremer_pople_q2, cremer_pople_phi2 = self._calc_cremer_pople_q2_and_phi(
                ring_atom_z_values
            )

            q3 = self._calc_cremer_pople_q3(ring_atom_z_values)

            cremer_pople_theta = self._calc_cremer_pople_theta(cremer_pople_q2, q3)

            cremer_pople_phi2_deg = (cremer_pople_phi2 * 180) / math.pi
            cremer_pople_theta_deg = (cremer_pople_theta * 180) / math.pi

            cremer_pople_pucker_values_per_frame[frame.frame] = {
                "cremer_pople_Q_value": cremer_pople_Q_value,
                "cremer_pople_phi2_deg": cremer_pople_phi2_deg,
                "cremer_pople_theta_deg": cremer_pople_theta_deg,
            }

        import pdb

        pdb.set_trace()

    def _calc_ring_atom_position_vectors_from_geometric_center(
        self, ring_center_of_geo, ring_pucker_atom_selection
    ):
        ring_atom_coordinates = self.mda_universe.select_atoms(
            ring_pucker_atom_selection
        ).positions

        ring_x_coors = [i[0] for i in ring_atom_coordinates]
        ring_y_coors = [i[1] for i in ring_atom_coordinates]
        ring_z_coors = [i[2] for i in ring_atom_coordinates]

        ring_atom_position_vectors = []

        for atom_number in range(len(ring_atom_coordinates)):
            x_vector = ring_x_coors[atom_number] - ring_center_of_geo[0]
            y_vector = ring_y_coors[atom_number] - ring_center_of_geo[1]
            z_vector = ring_z_coors[atom_number] - ring_center_of_geo[2]
            ring_atom_position_vectors.append([x_vector, y_vector, z_vector])

        return ring_atom_position_vectors

    def _calc_ring_atom_position_vector_projections(
        self, ring_atom_position_vectors, projection_type="x"
    ):

        number_of_atoms_in_ring = len(ring_atom_position_vectors)
        projected_vectors = []

        for atom_number in range(number_of_atoms_in_ring):
            if projection_type == "x":
                projection_value = self._calc_projection_value(
                    atom_number, number_of_atoms_in_ring
                )
            elif projection_type == "y":
                projection_value = self._calc_projection_value(
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

    def _calc_projection_value(
        self, atom_number, number_of_atoms_in_ring, m=1, method="cos"
    ):

        a = (2 * math.pi * m * atom_number) / number_of_atoms_in_ring

        if method == "cos":
            return math.cos(a)
        elif method == "sin":
            return math.sin(a)

        return None

    def _calc_ring_atom_z_values(
        self, ring_atom_position_vectors, x_projected_vectors, y_projected_vectors
    ):

        unit_vector_of_cross_product = self._calc_cross_product_unit_vector(
            x_projected_vectors, y_projected_vectors
        )

        number_of_atoms_in_ring = len(ring_atom_position_vectors)

        ring_atom_z_values = []

        for atom_number in range(number_of_atoms_in_ring):
            position_vector = ring_atom_position_vectors[atom_number]
            atom_z_value = np.dot(position_vector, unit_vector_of_cross_product)
            ring_atom_z_values.append(atom_z_value)

        return ring_atom_z_values

    def _calc_cross_product_unit_vector(self, x_projected_vectors, y_projected_vectors):
        summed_x_projected_vector = self._sum_vectors(x_projected_vectors)
        summed_y_projected_vector = self._sum_vectors(y_projected_vectors)

        cross_product_of_x_against_y_projected_vectors = np.cross(
            summed_x_projected_vector, summed_y_projected_vector
        )

        maginitude_of_crossed_vectors = np.linalg.norm(
            cross_product_of_x_against_y_projected_vectors
        )

        unit_vector_of_cross_product = (
            cross_product_of_x_against_y_projected_vectors
            / maginitude_of_crossed_vectors
        )

        return unit_vector_of_cross_product

    def _sum_vectors(self, vectors):

        summed_x_coordinates_vectors = sum([vector[0] for vector in vectors])
        summed_y_coordinates_vectors = sum([vector[1] for vector in vectors])
        summed_z_coordinates_vectors = sum([vector[2] for vector in vectors])

        return [
            summed_x_coordinates_vectors,
            summed_y_coordinates_vectors,
            summed_z_coordinates_vectors,
        ]

    def _calc_cremer_pople_Q_value(self, ring_atom_z_values):
        return math.sqrt(sum([z_value ** 2 for z_value in ring_atom_z_values]))

    def _calc_cremer_pople_q2_and_phi(self, ring_atom_z_values):

        number_of_atoms_in_ring = len(ring_atom_z_values)
        ring_atoms_zj_cos_phi = []
        ring_atoms_zj_sin_phi = []

        for atom_number in range(number_of_atoms_in_ring):
            x_projection_value = self._calc_cremer_pople_q2_projection(
                atom_number, number_of_atoms_in_ring, projection_type="x"
            )
            y_projection_value = self._calc_cremer_pople_q2_projection(
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

        q2_cos_phi2 = (
            math.sqrt(2 / number_of_atoms_in_ring) * summed_ring_atoms_zj_cos_phi
        )
        q2_sin_phi2 = (
            -1 * math.sqrt(2 / number_of_atoms_in_ring) * summed_ring_atoms_zj_sin_phi
        )

        q2 = self._calc_cremer_pople_q2(q2_cos_phi2, q2_sin_phi2)
        phi2 = self._calc_cremer_pople_phi(q2_cos_phi2, q2_sin_phi2)
        return q2, phi2

    def _calc_cremer_pople_q2_projection(
        self, atom_number, number_of_atoms_in_ring, projection_type="x"
    ):
        if projection_type == "x":
            return self._calc_projection_value(
                atom_number, number_of_atoms_in_ring, m=2
            )
        elif projection_type == "y":
            return self._calc_projection_value(
                atom_number, number_of_atoms_in_ring, m=2, method="sin"
            )

    def _calc_cremer_pople_q2(self, q2_cos_phi2, q2_sin_phi2):
        return math.sqrt(q2_cos_phi2 ** 2 + q2_sin_phi2 ** 2)

    def _calc_cremer_pople_phi(self, q2_cos_phi2, q2_sin_phi2):
        tan_phi2 = q2_sin_phi2 / q2_cos_phi2

        # not full sure why the '+ 180' is needed but it is.
        return math.atan(tan_phi2) + math.pi

    def _calc_cremer_pople_q3(self, ring_atom_z_values):

        number_of_atoms_in_ring = len(ring_atom_z_values)
        c = sum(
            [
                ring_atom_z_values[atom_number] * math.cos(math.pi * atom_number)
                for atom_number in range(number_of_atoms_in_ring)
            ]
        )

        q3 = (number_of_atoms_in_ring ** -0.5) * c

        return q3

    def _calc_cremer_pople_theta(self, q2, q3):
        tan_theta = q2 / q3

        if tan_theta < 0:
            return abs(math.atan(round(tan_theta, 5)))
        elif tan_theta > 0:
            return math.pi - math.atan(round(tan_theta, 5))
