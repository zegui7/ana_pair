package require saltbr
set curr_dir [pwd]
set pdb [lindex $argv 0]
set pre_1 [split $pdb "."]
set pre_2 [lindex $pre_1 0]
set pre_3 [split $pre_2 "/"]
set temp_index [expr [llength $pre_3] - 1]
set fold_name [lindex $pre_3 $temp_index][]
set pre_4 [lreplace $pre_3 $temp_index $temp_index]
set pre_5 [lappend pre_4 "saltbr_output" $fold_name]
set temp_dir [join $pre_5 "/"]
file mkdir $temp_dir
set pre_log [lappend pre_5 $fold_name.log]
set log_name [join $pre_log "/"]
append curr_dir "/" $log_name
mol new $pdb
set allatoms [atomselect top "all"]
saltbr -sel $allatoms -ondist 5.0 -writefiles yes -outdir $temp_dir -log $curr_dir
exit