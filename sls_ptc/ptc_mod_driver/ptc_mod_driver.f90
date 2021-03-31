program ptc_mod_driver

use ptc_mod

implicit none

character(100) ptc_flat_file
type(layout), pointer :: ptc_layout
type(probe_8) one_turn_ray8
type(probe) co_pr

call alloc(one_turn_ray8)

call getarg(1,ptc_flat_file)

!parse lattice
call my_ptc_parse(ptc_flat_file,ptc_layout,3)

!get 1-turn map and closed orbit at element 1
call ptc_one_turn_taylor(ptc_layout,1,one_turn_ray8,co_pr)

!calculat terms of 1-turn map
call ptc_populate_terms(ptc_layout, one_turn_ray8, co_pr)

!display the terms
call ptc_print_terms()

call kill(one_turn_ray8)

end program
