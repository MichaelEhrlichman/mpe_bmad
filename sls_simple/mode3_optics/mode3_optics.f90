program mode3_optics

use bmad
use mode3_mod
!use sls_lib

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s
real(rp), allocatable :: sbend_k2_state(:)
real(rp), allocatable :: sextupole_state(:)
real(rp), allocatable :: multipole_state(:)

integer i,j
integer status

logical error

real(rp) Ha, Hb
real(rp) nat_chrom_x, nat_chrom_y
real(rp) cor_chrom_x, cor_chrom_y
real(rp) s, delta_s
real(rp) oor
real(rp) total_bend_angle

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)
call twiss3_at_start(lat,error)
call twiss3_propagate_all(lat)

!
! Twiss by ele
!

open(45,file='twiss3.out')

write(45,'(2a)') '# ', lat_file
write(45,'(a11,2a11,22a14,"   ",a17)') '!        ix', 's (m)', 'l (m)', 'beta_a', 'alpha_a', 'gamma_a', 'eta_a', 'etap_a', 'phi_a', &
                                                                        'beta_b', 'alpha_b', 'gamma_b', 'eta_b', 'etap_b', 'phi_b', &
                                                                        'beta_c', 'alpha_c', 'gamma_c', 'eta_c', 'etap_c', 'phi_c', 'name'
do i=1,lat%n_ele_track
  write(45,'(i11,f11.5,f11.4,18es14.4,"   ",a14)') i, lat%ele(i)%s, lat%ele(i)%value(l$), &
        lat%ele(i)%mode3%a%beta, lat%ele(i)%mode3%a%alpha, lat%ele(i)%mode3%a%gamma, lat%ele(i)%mode3%a%eta, lat%ele(i)%mode3%a%etap, lat%ele(i)%mode3%a%phi, &
        lat%ele(i)%mode3%b%beta, lat%ele(i)%mode3%b%alpha, lat%ele(i)%mode3%b%gamma, lat%ele(i)%mode3%b%eta, lat%ele(i)%mode3%b%etap, lat%ele(i)%mode3%b%phi, &
        lat%ele(i)%mode3%c%beta, lat%ele(i)%mode3%c%alpha, lat%ele(i)%mode3%c%gamma, lat%ele(i)%mode3%c%eta, lat%ele(i)%mode3%c%etap, lat%ele(i)%mode3%c%phi, &
        lat%ele(i)%name
enddo
close(45)

end program












