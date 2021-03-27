/**
 * Parm7Prepper Application: a utility to reorder Parm7 files prepared by Glycam web's Carbohydrate builder unitlity into an ordering that is 
 * expected by NAMD. 
 *
 * @author Ryan Lazar (LZRRYA001)
 * July 2019
 */
import java.io.IOException;

public class Parm7Prepper {
	/**
	 * Application main method. 
	 * Expects command line argument 0, the name of the parm7 file to be converted.
	 *
	 * Example Usage:
	 * java -cp /path/to/Parm7Prepper/bin Parm7Prepper input.parm7
	 *
	 * where "input.parm7" is the .parm7 file for a solvated structure generated by glycam.com's Carbohydrate builder utility. 
	 */
    public static void main(String[] args) {
    	//Extract input file name
        String fileName = args[0];
        try {
        	//Represent internally
            Parm7File parm = new Parm7File(fileName);
            //Reorder and write output in order expected by NAMD
            parm.reorder();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }
}