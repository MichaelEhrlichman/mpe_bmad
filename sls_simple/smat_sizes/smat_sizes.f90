program smat_sizes

use bmad
use mode3_mod

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(normal_modes_struct) mode

integer i
integer status

logical err

real(rp) t6(6,6), sigma_mat(6,6)

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)
call calc_z_tune(lat)

mode%a%emittance = 1.0e-9
mode%b%emittance = 0.1e-9
mode%z%emittance = 1.0e-6
mode%a%tune = lat%a%tune
mode%b%tune = lat%b%tune
mode%z%tune = lat%z%tune

open(45,file='smat_sizes.out')

write(45,*) "!   a-mode emittance:          ", mode%a%emittance
write(45,*) "!   b-mode emittance:          ", mode%b%emittance
write(45,*) "!   z-mode emittance:          ", mode%z%emittance

write(45,'(2a11,3a25)') '!        ix', 's (m)', 'sqrt(sigma_mat(1,1))', 'sqrt(sigma_mat(3,3))', 'sqrt(sigma_mat(5,5))'
do i=1,lat%n_ele_track
  call transfer_matrix_calc (lat, t6, ix1=i, one_turn=.true.)
  call make_smat_from_abc(t6, mode, sigma_mat, err)
  write(45,'(i11,f11.5,3es25.5)') i, lat%ele(i)%s, sqrt(sigma_mat(1,1)), sqrt(sigma_mat(3,3)), sqrt(sigma_mat(5,5))
enddo
close(45)


end program












