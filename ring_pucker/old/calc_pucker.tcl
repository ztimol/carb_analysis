

set num_frames [molinfo top get numframes]

set filename "pucker_GlcNAc_L1.dat"

set fileId [open $filename "w"]

for {set i 0} {$i < $num_frames} {incr i} {
    animate goto $i
    set sel [atomselect top "index 6"]
    set pucker_value [$sel get pucker]
    puts -nonewline $fileId "$i $pucker_value \n"
}
