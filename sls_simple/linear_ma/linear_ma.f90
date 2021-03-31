program linear_ma_program

use bmad
use sls_lib
use touschek_mod
use linear_aperture_mod

implicit none

character lat_file*200
character use_line*20

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(momentum_aperture_struct), allocatable :: ma(:)

integer i, n_ma
integer status

real(rp) delta_s

call getarg(1, lat_file)
call getarg(2, use_line)
call bmad_parser(lat_file, lat, use_line=use_line)

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.
call set_on_off(rfcavity$,lat,off$)

call twiss_and_track(lat,orb,status)

delta_s = 0.01
n_ma = floor(lat%param%total_length/delta_s) + 1
allocate(ma(n_ma))
do i=1,n_ma
  ma(i)%s = delta_s*(i-1)
enddo

call linear_ma(lat,ma)

open(55,file='linear_ma.dat')

do i=1,n_ma
  write(55,'(i6,f14.5, 2es14.4)') i, ma(i)%s, ma(i)%neg, ma(i)%pos
enddo

close(55)

end program












