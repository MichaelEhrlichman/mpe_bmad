program one_turn_mat

use bmad
use sls_lib

implicit none

character lat_file*200
character ix_str*6

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)

real(rp) t6(6,6), v6(6)

integer ix
integer i, status

bmad_com%radiation_damping_on = .true.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call getarg(2, ix_str)

call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)

read(ix_str,*) ix
call transfer_matrix_calc(lat, t6, v6, ix1=ix, one_turn=.true.)

write(*,*) "One-turn transfer matrix:"
do i=1,6
  write(*,'(6es20.10)') t6(i,:)
enddo

write(*,*) 
write(*,*) "One-turn transfer vector:"
write(*,'(6es20.10)') v6(:)

end program












