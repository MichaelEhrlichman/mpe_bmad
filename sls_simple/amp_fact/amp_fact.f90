program amp_fact

use bmad
!use sls_lib

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s
type(ele_struct), pointer :: ele_at_kick

integer i,j
integer ikick
integer status

real(rp) s, delta_s
real(rp) k1l_at_kick
real(rp) nu_x, dphi_x
real(rp) x0, x_offset

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)

!
! Twiss by ele
!

ikick = 6
x_offset = 10.0d-6

ele_at_kick => lat%ele(ikick)
k1l_at_kick = ele_at_kick%value(k1$) * ele_at_kick%value(l$)
nu_x = lat%ele(lat%n_ele_track)%a%phi / twopi

open(45,file='amp_fact.out')

write(45,'(a14,a14)') "# loc (s)", "Dx"

s=0.0
delta_s = 0.10
do while (s .lt. lat%param%total_length)
  call twiss_and_track_at_s(lat, s, ele_at_s)
  !thin quad amplification factors
  dphi_x = ele_at_s%a%phi - ele_at_kick%a%phi
  x0 = sqrt(ele_at_s%a%beta*ele_at_kick%a%beta)*cos(pi*nu_x-abs(dphi_x))/2.0d0/sin(pi*nu_x) * k1l_at_kick * x_offset
  write(45,'(f14.5,es14.5)') s, x0
  s = s + delta_s
enddo
close(45)

end program












