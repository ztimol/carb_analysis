from generic_scripts.energy_by_pucker_amplitude import EnergyByPuckerAmplitude


def calc_cp_theta_energy_stats():

    cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13_aDGlc14_bDGlcNAc_glycam/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

    trajectory_pe_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13_aDGlc14_bDGlcNAc_glycam/namd_energies/trisaccharide_PE/trisaccharide_PE.dat"

    # cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13bDGlcNAc_3_10-aDGlc14bDGlcNAc_-31_-45_attempt_2/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

    # trajectory_pe_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13bDGlcNAc_3_10-aDGlc14bDGlcNAc_-31_-45_attempt_2/namd_energies/trisaccharide_PE/trisaccharide_PE.dat"

    x = EnergyByPuckerAmplitude()
    cp_trajectory_pucker_params = x.read_trajectory_cp_puckering_parameters(
        cp_trajectory_parameter_file
    )

    trajectory_energy_per_frame = x.read_trajectory_potential_energy(trajectory_pe_file)
    pe_by_conformers = x.seperate_energy_by_cp_theta(
        cp_trajectory_pucker_params, trajectory_energy_per_frame
    )

    x.plot_energy_against_cp_theta(
        cp_trajectory_pucker_params, trajectory_energy_per_frame
    )

    mean_energy_per_cp_theta_bin = x.bin_energies_by_cp_theta(
        cp_trajectory_pucker_params, trajectory_energy_per_frame
    )

    x.plot_binned_energies_against_cp_theta(mean_energy_per_cp_theta_bin)

    x.calc_energy_stats_by_conformers(pe_by_conformers)


def main():
    calc_cp_theta_energy_stats()


main()