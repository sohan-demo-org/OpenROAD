# Physics-driven global placement demo on gcd / Nangate45.
#
# Headless:  openroad -exit physics_demo_gcd.tcl
# GUI:       openroad -gui physics_demo_gcd.tcl
#
# Images are written to $PHYSICS_DEMO_OUT (default: results/physics_demo next
# to this script) when the GUI is up.  Not registered as a regression test.

set script_dir [file dirname [file normalize [info script]]]
set or_test_dir [file normalize [file join $script_dir .. .. .. test]]
if { [info exists ::env(PHYSICS_DEMO_OUT)] } {
  set out_dir [file normalize $::env(PHYSICS_DEMO_OUT)]
} else {
  set out_dir [file join $script_dir results physics_demo]
}
file mkdir $out_dir

proc section { title } {
  puts ""
  puts "=================================================================="
  puts "== $title"
  puts "=================================================================="
}

set gui_up [expr { [info commands gui::enabled] ne "" && [gui::enabled] }]

section "1. Read Nangate45 libraries and the gcd netlist"
cd $or_test_dir
source flow_helpers.tcl
source Nangate45/Nangate45.vars
read_libraries
read_verilog gcd_nangate45.v
link_design gcd
source flow_helpers.tcl ;# the SDC calls set_all_input_output_delays
read_sdc gcd_nangate45.sdc

section "2. Floorplan"
initialize_floorplan -site $site \
  -die_area {0 0 100.13 100.8} \
  -core_area {10.07 11.2 90.25 91}
source $tracks_file
place_pins -hor_layers metal3 -ver_layers metal2
set_wire_rc -signal -layer metal3
set_wire_rc -clock -layer metal5

if { $gui_up } {
  gui::fit
  save_image -area {0 0 100.13 100.8} -resolution 0.08 [file join $out_dir 01_floorplan.png]
}

section "3. global_placement -physics_driven"
global_placement -density 0.6 \
  -physics_driven \
  -physics_driven_overflows {79 64 40 20}

section "4. report_physics_placement"
report_physics_placement

section "5. Heat maps"
if { $gui_up } {
  gui::fit
  save_image -area {0 0 100.13 100.8} -resolution 0.08 [file join $out_dir 02_placement.png]
  foreach {label short_name name} {
    temperature PhysicsTemperature "Physics: Temperature"
    ir_drop PhysicsIRDrop "Physics: IR Drop"
    stress PhysicsStress "Physics: Reliability Stress"
    power_density PhysicsPowerDensity "Physics: Power Density"
  } {
    puts "-- $name"
    set control "Heat Maps/$name"
    gui::set_display_controls $control visible true
    gui::set_heatmap $short_name ShowLegend 1
    gui::set_heatmap $short_name rebuild
    save_image -area {0 0 100.13 100.8} -resolution 0.08 [file join $out_dir 03_heatmap_$label.png]
    gui::set_display_controls $control visible false
  }
  gui::set_display_controls "Heat Maps/Physics: Temperature" visible true
  puts "Images written to $out_dir"
} else {
  puts "GUI not running; skipping heat map images."
  puts "Re-run with 'openroad -gui physics_demo_gcd.tcl' to render them."
}

section "Done"
