# Copyright (C) 2010,2011,2012 The ESPResSo project
#  
# This file is part of ESPResSo.
#  
# ESPResSo is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#  
# ESPResSo is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#  
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>. 
#
source "tests_common.tcl"

require_feature ELECTROSTATICS 
require_feature MASS
require_feature EXTERNAL_FORCES

proc veccompare { a b } {
#  puts "$a $b"
  if { [ llength $a ]  != [ llength $b ] } {
    return 0
  }
  for { set i 0 } { $i < [ llength $a ] } { incr i } {
    set A [ lindex $a $i ] 
    set B [ lindex $b $i ] 
    if { !( $A == $B || [ expr abs($A - $B) < 1e-6 ]) } {
      return 0
    }
  }
  return 1
}


puts "---------------------------------------------------------------"
puts "- Testcase observable.tcl running on [format %02d [setmd n_nodes]] nodes"
puts "---------------------------------------------------------------"

setmd box_l 8 8 8
setmd skin 0.5
setmd time_step 0.01
thermostat off

part 0 pos 1 0 0 v 2 0 0 type 0 q 1 mass 1 ext_force 1 0 0 
part 1 pos 2 0 0 v 4 0 0 type 0 q -1 ext_force -2 0 0
part 2 pos 1 1 0 v 2 2 0 type 1
part 3 pos 2 1 0 v 4 2 0 type 1
part 4 pos 3 0 0 v 6 0 0 type 0 mass 2
integrate 0


############### Observable particle_position ##########################
## Here we are checking if particle specifications are working
# position particle 0
set p0 [ observable new particle_positions id 0 ]
if { ![ veccompare [ observable $p0 print ] { 1 0 0 } ] }  {
  error "P0 is not working"
}
# position particle 0, 3
set p1 [ observable new particle_positions id { 2 3 } ]
if { ![ veccompare [ observable $p1 print ] { 1 1 0 2 1 0} ] }  {
  error "P1 is not working"
}
# positions particle 0 1
set p2 [ observable new particle_positions type { 0 } ]
if { ![ veccompare [ observable $p2 print ] { 1 0 0 2 0 0 3 0 0 } ] }  {
  error "P2 is not working"
}
# position particle 0-3
set p3 [ observable new particle_positions type { 0 1 } ]
if { ![ veccompare [ observable $p3 print ] { 1 0 0 2 0 0 1 1 0 2 1 0  3 0 0 } ] }  {
  error "P3 is not working"
}


############## Observable particle velocities #######################

set v0 [ observable new particle_velocities id 0 ]
if { ![ veccompare [ observable $v0 print ] { 2 0 0 } ] }  {
  error "V0 is not working"
}
# position particle 0, 3
set v1 [ observable new particle_velocities id { 2 3 } ]
if { ![ veccompare [ observable $v1 print ] { 2 2 0 4 2 0} ] }  {
  error "V1 is not working"
}
# positions particle 0 1
set v2 [ observable new particle_velocities type { 0 } ]
if { ![ veccompare [ observable $v2 print ] { 2 0 0 4 0 0 6 0 0 } ] }  {
  error "V2 is not working"
}
# position particle 0-3
set v3 [ observable new particle_velocities type { 0 1 } ]
if { ![ veccompare [ observable $v3 print ] { 2 0 0 4 0 0 2 2 0 4 2 0  6 0 0 } ] }  {
  error "V3 is not working"
}
# position all particles
set v4 [ observable new particle_velocities all ]
if { ![ veccompare [ observable $v4 print ] { 2 0 0 4 0 0 2 2 0 4 2 0  6 0 0 } ] }  {
  error "V3 is not working"
}

############# Observable com_position ##########

set com_pos1 [ observable new com_position id { 0 1 } ]
if { ![ veccompare [ observable $com_pos1 print ] { 1.5 0 0 } ] }  {
  error "com_pos1 is not working"
}
set com_pos2 [ observable new com_position type 0 ]
if { ![ veccompare [ observable $com_pos2 print ] { 2.25 0 0 } ] }  {
  error "com_pos2 is not working"
}
set temp [ observable new particle_positions id { 0 1 2 3 } ]
set block_average [ observable new block_average $temp blocksize 2 stride 3 ]
puts [ observable $block_average print ]

############# Observable com_position ##########
set com_vel1 [ observable new com_velocity id { 0 1 } ]
if { ![ veccompare [ observable $com_vel1 print ] { 3 0 0 } ] }  {
  error "com_vel1 is not working"
}
set com_vel2 [ observable new com_velocity type 0 ]
if { ![ veccompare [ observable $com_vel2 print ] { 4.5 0 0 } ] }  {
  error "com_vel2 is not working"
}


############# Observable dipole_moment #####################
set dipm [ observable new dipole_moment id { 0 1 } ]
if { ![ veccompare [ observable $dipm print ] { -1 0 0 } ] }  {
  error "dipm is not working"
}

############# Observable current       #####################
set current [ observable new current id { 0 1 } ]
if { ![ veccompare [ observable $current print ] { -2 0 0 } ] }  {
  error "current is not working"
}

############# Observable particle_currents #####################
set particle_currents [ observable new particle_currents id { 0 1 } ]
if { ![ veccompare [ observable $particle_currents print ] { 2 0 0 -4 0 0 } ] }  {
  error "particle_currents is not working"
}

############# Observable stress tensor ##########
set stress_tensor [ observable new stress_tensor ] 
if { ![ veccompare [ observable $stress_tensor print ] [ lreplace [ lindex [ analyze stress ] 0 ] 0 0 ] ] } {
  error "stress_tensor is not working"
}

############# Observable com force and particle force ###### 
set particle_forces [ observable new particle_forces id { 0 1 } ] 
if { ![ veccompare [ observable $particle_forces print ] { 1 0 0 -2 0 0 } ] } {
  error "particle_forces is not working"
}

set com_force [ observable new com_force id { 0 1 } ] 
if { ![ veccompare [ observable $com_force print ] { -1 0 0 } ] } {
  error "particle_forces is not working"
}

############# Observable interacts_with and interaction_lifetimes ###### 
set pid [setmd n_part];
part $pid           pos 0.0 1.0 1.0 type 0
part [expr $pid+1]  pos -0.5 1.0 1.0 type 1
set cut 1.0;
set pos_mask 1;
set neg_mask -1;
set inter_id [observable new interacts_with id $pid id [expr $pid+1] $cut];
set lft_pos_id [observable new interaction_lifetimes id $pid id [expr $pid+1] $cut  $pos_mask];
set lft_neg_id [observable new interaction_lifetimes id $pid id [expr $pid+1] $cut  $neg_mask];
setmd time 0;
set inter_expect 1.0;
if { [observable $inter_id print] != $inter_expect } {
    error "inter failed, should be \"$inter_expect\", is \"[observable $inter_id print]\" ";
}
for {set i 0} {$i<10} {incr i} {
    setmd time [expr [setmd time]+$i]
    part [expr $pid+1] pos [expr (0.5+$i%2)*$cut] 1.0 1.0;
    observable $lft_pos_id update;
    observable $lft_neg_id update;
}
set inter_expect 0.0;
if { [observable $inter_id print] != $inter_expect } {
    error "inter failed, should be \"$inter_expect\", is \"[observable $inter_id print]\" ";
}
set lft_pos_expect "1.0 3.0 5.0 7.0 9.0 "
if { [observable $lft_pos_id print] != $lft_pos_expect } {
    error "lft_pos failed, should be \"$lft_pos_expect\", is \"[observable $lft_pos_id print]\" ";
}
set lft_neg_expect "2.0 4.0 6.0 8.0 "
if { [observable $lft_neg_id print] != $lft_neg_expect } {
    error "lft_neg failed, should be \"$lft_neg_expect\", is \"[observable $lft_neg_id print]\" ";
}

############# Observable tcl command #######################
proc p1 {} {
  return { 1 2 3 }
}
set tclcommand [ observable new tclcommand 3 "p1" ] 
puts [ p1 ]
puts "This is the output [ observable $tclcommand print ]"
if { ![ veccompare [ observable $tclcommand print ] { 1 2 3 } ] } {
  error "tclcommand is not working"
}

part 0 pos 0 0 0
set p0 [ observable new particle_positions id 0 ]
set av [ observable new average $p0 ]
observable $av update
part 0 pos 2 2 2
observable $av update
if { ![ veccompare [ observable $av print ] { 1 1 1 } ] } {
  error "average is not working"
}

set stddev_o [ observable new stddev $p0 ]
set variance_o [ observable new variance $p0 ]
set valuelist [ list ]
for { set i 0  } { $i < 10 } { incr i } {
  part 0 pos $i $i $i
  observable $stddev_o update
  observable $variance_o update
  lappend valuelist $i
}

proc mean { l } {
  set sum 0.
  foreach x $l {
    set sum [ expr $sum + $x ]
  }
  return [ expr $sum / [ llength $l ] ]
}

proc stdev { l } {
  set sum 0.
  set themean [ mean $l ]
  foreach x $l {
    set sum [ expr $sum  + ($x-$themean)*($x-$themean) ]
  }
  return [ expr sqrt($sum / [ llength $l] ) ]
}

set stddev1 [ lindex [ observable $stddev_o print ] 0 ]
set stddev2 [ stdev $valuelist ]
set variance1 [ lindex [ observable $variance_o print ] 0 ]
set variance2 [ expr $stddev2*$stddev2 ]

if { [ expr abs($stddev1 - $stddev2) > 1e-5 ] } {
  error "stddev does not work"
}


if { [ expr abs($variance1 - $variance2) > 1e-5 ] } {
  error "stddev does not work"
}


