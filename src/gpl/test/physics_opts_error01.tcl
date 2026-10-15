# -physics_force_ratio out of range must fail with GPL-0434
source helpers.tcl
read_lef ./nangate45.lef
read_def ./error01.def

catch { global_placement -physics_driven -physics_force_ratio -1 } msg
puts $msg
if { [string match "*GPL-0434*" $msg] } {
  puts pass
} else {
  puts "fail: expected GPL-0434"
}
