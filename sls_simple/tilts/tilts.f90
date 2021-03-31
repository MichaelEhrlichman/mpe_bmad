program tilts

use bmad
use mode3_mod
use sls_lib

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(normal_modes_struct) mode

integer i
integer status

real(rp) t6(6,6), sigma_mat(6,6)
real(rp) txz

logical err

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)
call radiation_integrals(lat,orb,mode)

lat%a%tune = mod(lat%ele(lat%n_ele_track)%a%phi,twopi)
lat%b%tune = mod(lat%ele(lat%n_ele_track)%b%phi,twopi)
call calc_z_tune(lat)
mode%a%tune = lat%a%tune
mode%b%tune = lat%b%tune
mode%z%tune = lat%z%tune
mode%a%emittance = 38.0e-9
mode%b%emittance = 0.190e-9
mode%z%emittance = 4.4e-4 * 7.0594e-3
mode%sigE_E = 4.4e-4

open(45,file='tilts.out')
write(45,'(a11,a11,3a14)') '!        ix', 's (m)', 'theta_xz', 'sig11', 'linear x size'
do i=1,lat%n_ele_track
  call transfer_matrix_calc (lat, t6, ix1=i, one_turn=.TRUE.)
  call make_smat_from_abc(t6, mode, sigma_mat, err)
  txz = 0.5*atan2(2.0*sigma_mat(1,5), sigma_mat(5,5)-sigma_mat(1,1))
  write(45,'(i11,f11.4,3es14.5)') i, lat%ele(i)%s, txz, sqrt(sigma_mat(1,1)), sqrt(lat%ele(i)%a%beta*mode%a%emittance+(lat%ele(i)%x%eta*mode%sigE_E)**2)
enddo
close(45)

end program












