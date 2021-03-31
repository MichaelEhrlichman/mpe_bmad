program off_optics

use bmad
!use sls_lib

implicit none

character lat_file*200
character*7 de_str

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s
type(coord_struct) orb_at_s

integer i
integer status

real(rp) de
real(rp) t6(6,6), v6(6)
real(rp) symp_check
real(rp) s, delta_s

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call getarg(2, de_str)
read(de_str,*) de

call bmad_parser(lat_file, lat)
call set_on_off(rfcavity$, lat, off$)

allocate(orb(0:lat%n_ele_track))
orb(0)%vec(:) = 0.0d0
orb(0)%vec(6) = de
call twiss_and_track(lat,orb,status)

if(status /= ok$) then
  write(*,*) "Lattice is unstable."
  stop
endif

call transfer_matrix_calc(lat, t6, v6, ix1=0, one_turn=.true.)
symp_check = mat_symp_error(t6)

open(45,file='off_twiss_'//adjustl(trim(de_str))//'.out')

write(45,*) "!   a-mode tune:          ", lat%ele(lat%n_ele_track)%a%phi /2./pi
write(45,*) "!   b-mode tune:          ", lat%ele(lat%n_ele_track)%b%phi /2./pi
write(45,*) "!   de:                   ", de
write(45,*) "!   symplectic:           ", symp_check

write(45,'(a11,2a11,17a14,"   ",a14)') '!        ix', 's (m)', 'l (m)', 'beta_a', 'alpha_a', 'gamma_a', 'eta_a', 'etap_a', 'eta_x', 'etap_x', 'phi_x', &
                      'beta_b', 'alpha_b', 'gamma_b', 'eta_b', 'etap_b', 'eta_y', 'etap_y', 'phi_y', 'x_co', 'name'
s=0.0
delta_s = 0.01
i=1
do while (s .lt. lat%param%total_length)
  call twiss_and_track_at_s(lat, s, ele_at_s, orb, orb_at_s)
  write(45,'(i11,2f11.3,17es14.4,"   ",a14)') i, &
        ele_at_s%s, ele_at_s%value(l$), &
        ele_at_s%a%beta, ele_at_s%a%alpha, ele_at_s%a%gamma, ele_at_s%a%eta, ele_at_s%a%etap, &
        ele_at_s%x%eta, ele_at_s%x%etap, ele_at_s%a%phi, &
        ele_at_s%b%beta, ele_at_s%b%alpha, ele_at_s%b%gamma, ele_at_s%b%eta, ele_at_s%b%etap, &
        ele_at_s%y%eta, ele_at_s%y%etap, ele_at_s%b%phi, orb_at_s%vec(1), &
        ele_at_s%name
  s = s + delta_s
  i = i + 1
enddo
close(45)

end program


