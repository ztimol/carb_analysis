cd __SOLVATION_DIR_TOKEN__
package require solvate
solvate __PSF_FILE_PATH_TOKEN__ __PDB_FILE_PATH_TOKEN__ -o solvate -s WT -minmax {{__LOWER_BOUND_TOKEN__ __LOWER_BOUND_TOKEN__ __LOWER_BOUND_TOKEN__} {__UPPER_BOUND_TOKEN__ __UPPER_BOUND_TOKEN__ __UPPER_BOUND_TOKEN__}} -rotate -rotsel { all } -rotinc 10 -x 0 -y 0 -z 0 +x 0 +y 0 +z 0 -b 0.5
set everyone [atomselect top all]
set minmax_value [measure minmax $everyone]
#puts $minmax_value
set center_value [measure center $everyone]
#puts $center_value 
exit
