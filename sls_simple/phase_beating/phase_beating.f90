program phase_beating

use bmad

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)

integer i, ix
integer status

real(rp) dKl, B_ix, P_ix, dPhi
real(rp) C0, nu_x

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)

ix = 449
dKl = -0.00005
B_ix = lat%ele(ix)%a%beta
P_ix = lat%ele(ix)%a%phi
C0 = B_ix * dKl / 2.0d0
nu_x = lat%ele(lat%n_ele_track)%a%phi

open(45,file='phase_beating.out')
do i=1,lat%n_ele_track
  if(lat%ele(i)%name(1:4) == 'DET_') then
    if(i > ix) then
      dPhi = C0*(1 + cos(2*P_ix-lat%ele(i)%a%phi)*sin(lat%ele(i)%a%phi-nu_x)/sin(nu_x))
    else
      dPhi = C0*(sin(lat%ele(i)%a%phi)*cos(lat%ele(i)%a%phi-2*P_ix+nu_x)/sin(nu_x))
    endif
    write(45,'(i6,f11.5,es14.4,"   ",a)') i, lat%ele(i)%s, dPhi, lat%ele(i)%name
  endif
enddo
close(45)

end program





