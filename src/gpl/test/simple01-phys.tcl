# physics-driven global placement on simple01 (copy of simple01-td)
source helpers.tcl
set test_name simple01-phys
read_liberty ./library/nangate45/NangateOpenCellLibrary_typical.lib
read_lef ./nangate45.lef
read_def ./simple01-td.def
create_clock -name core_clock -period 2 clk
set_wire_rc -signal -layer metal3
set_wire_rc -clock -layer metal5

global_placement -physics_driven \
  -physics_driven_overflows {79 64 40 20}

report_physics_placement

set def_file [make_result_file $test_name.def]
write_def $def_file
diff_files $def_file $test_name.defok
