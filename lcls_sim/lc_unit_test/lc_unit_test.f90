program lc_unit_test

use bmad
use bmad_parser_mod

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s

integer i,j
integer ix
integer status

real(rp) rmat(6,6), vec0(6)
real(rp) bi, bf

bp_com%always_parse = .true.
bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

allocate(orb(0:lat%n_ele_track))
orb(0)%vec(6) = 0.0d0
call twiss_and_track(lat,orb,status)

! Twiss by ele end

open(45,file='twiss.out')

write(45,*) "!   a-mode phase advance:    ", lat%ele(lat%n_ele_track)%a%phi /2./pi
write(45,*) "!   b-mode phase advance:    ", lat%ele(lat%n_ele_track)%b%phi /2./pi
write(45,*) "!   total length:            ", lat%param%total_length

write(45,'(a11,a15,a11,21a16,"   ",a14)') '!        ix', 's (m)', 'l (m)', 'beta_a', 'alpha_a', 'gamma_a', 'eta_a', 'etap_a', 'eta_x', 'etap_x', 'phi_x', &
                      'beta_b', 'alpha_b', 'gamma_b', 'eta_b', 'etap_b', 'eta_y', 'etap_y', 'phi_y', 'energy', 'r55', 'r56', 'r65', 'r66', 'name'

open(46,file='rmat_mad8.out')
write(46,'(a11,a15,a11,21a16,"   ",a14)') '!        ix', 's (m)', 'l (m)', 'r55', 'r56', 'r65', 'r66', 'name'

do i=1,lat%n_ele_track
  call transfer_matrix_calc(lat,rmat,vec0,i-1,i)
  bi = orb(i-1)%beta
  bf = orb(i)%beta
  write(45,'(i11,f15.8,f11.4,21es16.7,"   ",a14)') i, &
        lat%ele(i)%s, lat%ele(i)%value(l$), &
        lat%ele(i)%a%beta, lat%ele(i)%a%alpha, lat%ele(i)%a%gamma, lat%ele(i)%a%eta, lat%ele(i)%a%etap, &
        lat%ele(i)%x%eta, lat%ele(i)%x%etap, lat%ele(i)%a%phi, &
        lat%ele(i)%b%beta, lat%ele(i)%b%alpha, lat%ele(i)%b%gamma, lat%ele(i)%b%eta, lat%ele(i)%b%etap, &
        lat%ele(i)%y%eta, lat%ele(i)%y%etap, lat%ele(i)%b%phi, value_of_attribute(lat%ele(i),'e_tot'), &
        rmat(5,5), rmat(5,6), rmat(6,5), rmat(6,6), &
        lat%ele(i)%name
  write(46,'(i11,f15.8,f11.4,4es16.7,"   ",a14)') i, lat%ele(i)%s, lat%ele(i)%value(l$), &
        bi/bf*rmat(5,5), rmat(5,6)/bi/bf, bi*bf*rmat(6,5), bi/bf*rmat(6,6), &
        lat%ele(i)%name
enddo
close(45)
close(46)

end program












