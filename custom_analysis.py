from custom_scripts.potential_energy_by_pucker_amplitude import (
    PotentialEnergyByPuckerAmplitude,
)
from custom_scripts.cp_pucker_theta_free_energy import CpPuckerThetaFreeEnergy
from custom_scripts.cp_pucker_phi_theta_free_energy import CpPuckerPhiThetaFreeEnergy

# cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13_aDGlc14_bDGlcNAc_glycam/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

# trajectory_pe_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13_aDGlc14_bDGlcNAc_glycam/namd_energies/trisaccharide_PE/trisaccharide_PE.dat"

# cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13bDGlcNAc_3_10-aDGlc14bDGlcNAc_-31_-45_attempt_2/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

# trajectory_pe_file = "/home/timol/C6W/Studies/structure_analysis/output/aLRha13bDGlcNAc_3_10-aDGlc14bDGlcNAc_-31_-45_attempt_2/namd_energies/trisaccharide_PE/trisaccharide_PE.dat"


def calc_potential_energy_by_pucker_amplitude():

    cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlcNAc/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

    trajectory_pe_file = "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlcNAc/namd_energies/trisaccharide_PE/trisaccharide_PE.dat"

    x = PotentialEnergyByPuckerAmplitude()
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

    x.calc_energy_stats_by_conformers(pe_by_conformers)


def calc_cp_pucker_amplitude_free_energy():

    cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlcNAc/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

    x = CpPuckerThetaFreeEnergy()
    x.read_trjectory_cp_puckering_parameters(cp_trajectory_parameter_file)
    x.count_cp_theta_by_bins()
    x.calc_free_energy_for_bins()
    x.write_theta_count_free_energy()
    x.plot_free_energy_against_cp_theta()


def calc_cp_pucker_phi_theta_free_energy():
    cp_trajectory_parameter_file = "/home/timol/C6W/Studies/structure_analysis/output/aDGlc13_aDGlc14_bDGlcNAc/ring_pucker/GlcNAc/trajectory_cp_phi_theta_Q.dat"

    x = CpPuckerPhiThetaFreeEnergy()
    x.read_trjectory_cp_puckering_parameters(cp_trajectory_parameter_file)
    x.count_cp_phi_theta_by_bins()
    x.calc_free_energy_for_bins()
    x.write_phi_theta_count_free_energy()


def main():
    # calc_potential_energy_by_pucker_amplitude()

    # calc_cp_pucker_amplitude_free_energy()

    calc_cp_pucker_phi_theta_free_energy()


main()
