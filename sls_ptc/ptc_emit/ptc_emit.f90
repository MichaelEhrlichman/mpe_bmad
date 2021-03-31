program ptc_emit

use bmad
use ptc_layout_mod

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: bmad_co(:)
type(coord_struct) ptc_co
type (normal_modes_struct) bmad_mode
type (normal_modes_struct) ptc_mode
real(rp) sigma_mat(6,6)
integer i, status

call getarg(1, lat_file)

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)
bmad_com%radiation_damping_on = .true.
call lat_to_ptc_layout(lat)

call ptc_emit_calc(lat%ele(0), ptc_mode, sigma_mat, ptc_co)
ptc_mode%sig_z = SQRT(sigma_mat(5,5))
ptc_mode%sigE_E = SQRT(sigma_mat(6,6))

write(*,*) "ptc emit a:   ", ptc_mode%a%emittance
write(*,*) "ptc emit b:   ", ptc_mode%b%emittance
write(*,*) "ptc emit z:   ", ptc_mode%z%emittance
write(*,*) "ptc sig z:    ", ptc_mode%sig_z
write(*,*) "ptc sigE_E:   ", ptc_mode%sigE_E
write(*,*) "ptc nua:      ", ptc_mode%a%tune*twopi
write(*,*) "ptc nub:      ", ptc_mode%b%tune*twopi
write(*,*) "ptc nus:      ", ptc_mode%z%tune*twopi
write(*,*) "ptc taua:     ", lat%param%total_length / c_light / ptc_mode%a%alpha_damp
write(*,*) "ptc taub:     ", lat%param%total_length / c_light / ptc_mode%b%alpha_damp
write(*,*) "ptc tauc:     ", lat%param%total_length / c_light / ptc_mode%z%alpha_damp

!call twiss_and_track(lat, bmad_co)
call closed_orbit_calc(lat,bmad_co,6)
call twiss_and_track(lat,bmad_co,status)
call radiation_integrals(lat, bmad_co, bmad_mode)
call calc_z_tune(lat)
bmad_mode%a%tune = mod(lat%ele(lat%n_ele_track)%a%phi,twopi)
bmad_mode%b%tune = mod(lat%ele(lat%n_ele_track)%b%phi,twopi)
bmad_mode%z%tune = lat%z%tune

write(*,*) 
write(*,*) "bmad emit a:  ", bmad_mode%a%emittance
write(*,*) "bmad emit b:  ", bmad_mode%b%emittance
write(*,*) "bmad emit z:  ", bmad_mode%z%emittance
write(*,*) "bmad sig z:   ", bmad_mode%sig_z
write(*,*) "bmad sigE_E:  ", bmad_mode%sigE_E
write(*,*) "bmad nua:     ", bmad_mode%a%tune
write(*,*) "bmad nub:     ", bmad_mode%b%tune
write(*,*) "bmad nus:     ", bmad_mode%z%tune
write(*,*) "bmad taua:    ", lat%param%total_length / c_light / bmad_mode%a%alpha_damp
write(*,*) "bmad taub:    ", lat%param%total_length / c_light / bmad_mode%b%alpha_damp
write(*,*) "bmad tauc:    ", lat%param%total_length / c_light / bmad_mode%z%alpha_damp

end program ptc_emit












