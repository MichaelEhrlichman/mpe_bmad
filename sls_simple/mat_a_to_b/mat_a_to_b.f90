program mat_a_to_b

use bmad

implicit none

character lat_file*200
character ix_a_str*6
character ix_b_str*6

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)

real(rp) t6(6,6), v6(6)

integer ix_a, ix_b
integer i, status

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call getarg(2, ix_a_str)
call getarg(3, ix_b_str)

call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)

read(ix_a_str,*) ix_a
read(ix_b_str,*) ix_b
call transfer_matrix_calc(lat, t6, v6, ix_a, ix_b)

write(*,*) "transfer matrix:"
do i=1,6
  write(*,'(6es20.10)') t6(i,:)
enddo

write(*,*) 
write(*,*) "transfer vector:"
write(*,'(6es20.10)') v6(:)

end program












