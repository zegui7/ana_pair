package require saltbr
set pdb [lindex $argv 0]
set pre_1 [split $pdb "."]
set pre_2 [lindex $pre_1 0]
set pre_3 [split $pre_2 "/"]
set temp_index [expr [llength $pre_3] - 1]
set fold_name [lindex $pre_3 $temp_index]
set pre_4 [lreplace $pre_3 $temp_index $temp_index]
lappend pre_4 "saltbr_output" $fold_name
set temp_dir [join $pre_4 "/"]
file mkdir $temp_dir
set pre_log [lappend pre_4 $fold_name.log]
set log_name [join $pre_log "/"]
mol new $pdb
set allatoms [atomselect top "all"]
saltbr -sel $allatoms -ondist 5.0 -writefiles yes -outdir $temp_dir -log $log_name
exit
