/**
 * Parm7File Class: stores and manipulates sections of a parm7 file
 *
 * @author Ryan Lazar (LZRRYA001)
 * July 2019
 */
import java.io.*;
import java.util.HashMap;

public class Parm7File {

    private HashMap<String,String> parmFile;
    private String fileName;
    /**
     * Constructs a Parm7File object by reading and storing each section and storing in a HashMap.
     *
     * @param inputParmFile the name of the parm7 file, in the same director as the invoked application, to be converted. 
     */
    Parm7File(String inputParmFile) throws IOException {
        fileName = inputParmFile;
        //Create hashmap representing the parm7 file
        parmFile = new HashMap<String, String>();
        InputStream is = new FileInputStream(fileName);
        BufferedReader buf = new BufferedReader(new InputStreamReader(is));
        String line = buf.readLine();
        //Store version of the parm7 file (first line) separately
        parmFile.put("VERSION",line);
        //IStores the current section being read in
        String current = "";
        line = buf.readLine();
        //Read input parm7 file in line by line
        while(line != null)
        {
            line.concat("\n");
            //Indicates a new section has been reached
            if (line.startsWith("%FLAG")){
                //Store the section as a key in the parmFile hashmap
                String key = line.split("\\s+")[1];
                parmFile.put(key,"\n"+line);
                //Indicate the all subsequent lines belong in this key's value store
                current = key;
            //For each non-flag, non-blank line (i.e. the lines of a section)    
            }else if (!current.equalsIgnoreCase("")) {
                //Concatenate to the value store of the current section being read
                parmFile.put(current,parmFile.get(current).concat("\n").concat(line));
            }
            line = buf.readLine();
        }
    }

    /**
     * Reorder and write Parm7 output in the correct order as expected by NAMD
     *
     * Output is stored in a parm7 file called "reordered.parm7"
     */
    public void reorder() throws IOException {
        //Open "reordered.parm7" for writing
        FileWriter fw=new FileWriter("reordered.parm7");
        //Write each section of the parm7 into an order that NAMD expects
        fw.write(parmFile.get("VERSION"));
        fw.write(parmFile.get("TITLE"));
        fw.write(parmFile.get("POINTERS"));
        fw.write(parmFile.get("ATOM_NAME"));
        fw.write(parmFile.get("CHARGE"));
        fw.write(parmFile.get("MASS"));
        fw.write(parmFile.get("ATOM_TYPE_INDEX"));
        fw.write(parmFile.get("NUMBER_EXCLUDED_ATOMS"));
        fw.write(parmFile.get("NONBONDED_PARM_INDEX"));
        fw.write(parmFile.get("RESIDUE_LABEL"));
        fw.write(parmFile.get("RESIDUE_POINTER"));
        fw.write(parmFile.get("BOND_FORCE_CONSTANT"));
        fw.write(parmFile.get("BOND_EQUIL_VALUE"));
        fw.write(parmFile.get("ANGLE_FORCE_CONSTANT"));
        fw.write(parmFile.get("ANGLE_EQUIL_VALUE"));
        fw.write(parmFile.get("DIHEDRAL_FORCE_CONSTANT"));
        fw.write(parmFile.get("DIHEDRAL_PERIODICITY"));
        fw.write(parmFile.get("DIHEDRAL_PHASE"));
        //Empty sections must be trimmed
        fw.write(trimEnd(parmFile.get("SOLTY")));
        fw.write(parmFile.get("SCEE_SCALE_FACTOR"));
        fw.write(parmFile.get("LENNARD_JONES_ACOEF"));
        fw.write(parmFile.get("LENNARD_JONES_BCOEF"));
        fw.write(parmFile.get("BONDS_INC_HYDROGEN"));
        fw.write(parmFile.get("BONDS_WITHOUT_HYDROGEN"));
        fw.write(parmFile.get("ANGLES_INC_HYDROGEN"));
        fw.write(parmFile.get("ANGLES_WITHOUT_HYDROGEN"));
        fw.write(parmFile.get("DIHEDRALS_INC_HYDROGEN"));
        fw.write(parmFile.get("DIHEDRALS_WITHOUT_HYDROGEN"));
        fw.write(parmFile.get("EXCLUDED_ATOMS_LIST"));
        fw.write(parmFile.get("HBOND_ACOEF"));
        fw.write(parmFile.get("HBOND_BCOEF"));
        fw.write(parmFile.get("HBCUT"));
        fw.write(parmFile.get("AMBER_ATOM_TYPE"));
        fw.write(parmFile.get("TREE_CHAIN_CLASSIFICATION"));
        fw.write(parmFile.get("JOIN_ARRAY"));
        fw.write(parmFile.get("IROTAT"));
        fw.write(parmFile.get("SOLVENT_POINTERS"));
        fw.write(parmFile.get("ATOMS_PER_MOLECULE"));
        fw.write(parmFile.get("BOX_DIMENSIONS"));
        fw.write(parmFile.get("RADIUS_SET"));
        fw.write(parmFile.get("RADII"));
        fw.write(parmFile.get("SCREEN"));
        fw.write(parmFile.get("SCNB_SCALE_FACTOR"));
        fw.write(parmFile.get("ATOMIC_NUMBER"));
        fw.close();
    }

    /**
     * Remove white space from the end of input string
     * @param value the string to be trimmed of whitespace
     */
    private String trimEnd(String value) {
        int len = value.length();
        int st = 0;
        while ((st < len) && value.charAt(len - 1) == '\n') {
            len--;
        }
        return value.substring(0, len);
    }
}
