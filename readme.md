THIS README IS INCOMPLETE. INFORMATION WILL BE ADDED AND CHANGED CONTINUOUSLY.

# Molecular Dynamics Structure Analysis module


## Multiprocessing

The underlying MDAnalysis libraries do not employ multiprocessing - see the MDAnalysis docs. The way trajectories are handled by MDAnalysis also do not make implementing multiprocessing on top of it trivial. There is a good reason for this - again google and the MDAnalysis docs provide more information.


## Config file 

### Config Parameters

* **frames_per_ns**
  * Optional: False
  * Acecptable value(s): numeric (integer)
  * Default value: *n/a*
  * Description: Used to calculate the trajectory length in time from the number of frames. Value must be specific in config_file.
  
* **start_frame**:
  * Optional: True
  * Acceptable value(s): numeric (integer)
  * Default value: *0*
  * Description: Used as the trajectory start point in the analysis.
  
* dcd_file:
  * Optional: False
  * Acceptable value(s): string
  * Default value: *n/a*
  * Description: Absolute or relative file path to trajectory DCD file.
    
* psf_file:
  * Optional: False
  * Acceptable value(s): string
  * Default value: *n/a*
  * Description: Absolute or relative file path to structure file.
  
  
* amber:
  * Optional: True
  * Acceptable value(s): *yes* or *no*
  * Default value: *no*
  * Description: Absolute or relative file path to structure file.
  
 
### Config file example

Below is an example config file that does the following:

* Measured the torsion angles for two sets of glycisidic linkages - phi and psi for each;
* Calculates the single point potential energies for the trajectory with the NAMD Amber force field parameters enabled;
* Calculates the Cremer-Pople ring pucker polar coordiates for two six membered rings.

```
start_frame 0
frames_per_ns 40

cores 4


dcd_file /home/[user]/C6W/Studies/Dynamics/NAMD_glycam/MD/solution/general_structures/aLRha13_aDGlc14_bDGlcNAc/trajectories/aLRha13_aDGlc14_bDGlcNAc_glycam_0-1000ns.dcd

psf_file /home/[user]/C6W/Studies/Dynamics/NAMD_glycam/MD/solution/general_structures/aLRha13_aDGlc14_bDGlcNAc/trajectories/1_noWAT.psf

amber yes

parm7_file /home/[user]/C6W/Studies/Dynamics/NAMD_glycam/MD/solution/general_structures/aLRha13_aDGlc14_bDGlcNAc/trajectories/1_noWAT.parm7

rst7_file /home/[user]/C6W/Studies/Dynamics/NAMD_glycam/MD/solution/general_structures/aLRha13_aDGlc14_bDGlcNAc/trajectories/1_noWAT.rst7

torsion aLRha13bDGlcNAc phi "51 50 27 15"
torsion aLRha13bDGlcNAc psi "50 27 15 16"
torsion aLRha13bDGlcNAc scatter phi psi

torsion aDGlc14bDGlcNAc phi "29 28 14 12"
torsion aDGlc14bDGlcNAc psi "28 14 12 13"
torsion aDGlc14bDGlcNAc scatter phi psi

namd_path  /home/timol/.NAMD_2.13_Linux-x86_64-multicore/namd2
namd_energy trisaccharide_PE potential_energy

ring_pucker GlcNAc "resid 2 and name O5 C1 C2 C3 C4 C5"
ring_pucker Glc "resid 3 and name O5 C1 C2 C3 C4 C5"
```
