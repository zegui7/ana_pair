package require hbonds
set pdb [lindex $argv 0]
set pre_1 [split $pdb "."]
set pre_2 [lindex $pre_1 0]
set pre_3 [split $pre_2 "/"]
set temp_index [expr [llength $pre_3] - 1]
set fold_name [lindex $pre_3 $temp_index]
puts $fold_name
set pre_4 [lreplace $pre_3 $temp_index $temp_index]
lappend pre_4 "sasa_hb"
set pre_5 $pre_4
set temp_dir [join $pre_4 "/"]
file mkdir $temp_dir
lappend pre_4 sasa_$fold_name.log
lappend pre_5 hb_$fold_name.log
set log_name [join $pre_4 "/"]
set log_name_hb [join $pre_5 "/"]
if {[file exists $log_name] == 0} {
      mol new $pdb
      set output [open $log_name w]
      set allatoms [atomselect top "all"]
      set residlist [lsort -unique -dictionary [$allatoms get [list chain resid resname]]]
      set residuelist [lsort -unique -dictionary [$allatoms get residue]]
      hbonds -sel1 $allatoms -writefile yes -detailout $log_name_hb -type all
      foreach r $residuelist {
            set current [lindex $residlist $r]
            set curr_chain [lindex $current 0]
            set curr_res [lindex $current 1] 
            set curr_name [lindex $current 2]
            set sasa_temp [measure sasa 1.0 $allatoms -restrict [atomselect top "resid $curr_res ch $curr_chain"]]
            puts $output "$curr_res $curr_chain $curr_name $sasa_temp"
      }
      close $output
}
exit