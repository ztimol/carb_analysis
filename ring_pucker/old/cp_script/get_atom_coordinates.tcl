

set num_frames [molinfo top get numframes]

#set sel [atomselect top "resid 1 and name O5 C1 C2 C3 C4 C5"]

set resid "1"
set linkage_number "L1"

set filename "ring_coordinates_GlcNAc_${linkage_number}_resid$resid.dat"


set o5 [atomselect top "resid $resid and name O5"]
set c1 [atomselect top "resid $resid and name C1"]
set c2 [atomselect top "resid $resid and name C2"]
set c3 [atomselect top "resid $resid and name C3"]
set c4 [atomselect top "resid $resid and name C4"]
set c5 [atomselect top "resid $resid and name C5"]

set fileId [open $filename "w"]

for {set i 0} {$i < $num_frames} {incr i} {
animate goto $i
#set atom_coordinates [$sel get {x y z}]
puts -nonewline $fileId "$i: "

set o5_coor [lindex [$o5 get {x y z}] 0]
set c1_coor [lindex [$c1 get {x y z}] 0]
set c2_coor [lindex [$c2 get {x y z}] 0]
set c3_coor [lindex [$c3 get {x y z}] 0]
set c4_coor [lindex [$c4 get {x y z}] 0]
set c5_coor [lindex [$c5 get {x y z}] 0]
puts -nonewline $fileId "$o5_coor,"
puts -nonewline $fileId "$c1_coor,"
puts -nonewline $fileId "$c2_coor,"
puts -nonewline $fileId "$c3_coor,"
puts -nonewline $fileId "$c4_coor,"
puts -nonewline $fileId "$c5_coor\n"
}
