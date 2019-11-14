#from MDAnalysis import Universe
import MDAnalysis as mda

class PdbReader():
    
    def get_all_atom_coordinates_frames_from_trajectory(pdb_file_path):

        dcd_file = 'S_flexneri_7a_6RU_0-32ns_every100frms.dcd'

        u = mda.Universe('S_flexneri_7a_6RU.psf', dcd_file)
        s = u.select_atoms('name H1 H2 H3 H4 H5 H6')
        print(u.atoms)
        ring_proton_coordinates_for_frame = s.positions

        def get_single_frames_from_trajectory(pdb_file_path, frame_number):
            pass


def main():
    dcd_file = 'S_flexneri_7a_6RU_0-32ns_every100frms.dcd'

    u = mda.Universe('S_flexneri_7a_6RU.psf', dcd_file)
    s = u.select_atoms('name H1 H2 H3 H4 H5 H6')
    for atom in s.atoms:
        print(atom.index)
    ring_proton_coordinates_for_frame = s.positions


if __name__ == '__main__':
    main()
