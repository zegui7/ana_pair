set coord [lindex $argv 0]
set traj [lindex $argv 1]
set filename [lindex $argv 2]
set inter [lindex $argv 3]
mol new $coord
mol addfile $traj waitfor all molid $coord
file mkdir [format %s/all_frames $inter]
set allsel [atomselect top protein]
puts $tcl_version
#for {set i 0} {$i < [molinfo top get numframes]} {incr i} {
for {set i 0} {$i < 3} {incr i} {
	catch {
		set n [expr $i + 1]
		set framename [format %s/all_frames/%s_%s.pdb $inter $filename $n]
		molinfo top set frame $i
		$allsel writepdb $framename
		puts Frame:$n
	}
	}
exit
