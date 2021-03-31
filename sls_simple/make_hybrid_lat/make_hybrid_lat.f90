program hybridize

use bmad

implicit none

character lat_file*200
type(lat_struct) lat, hybrid_lat
type(coord_struct), allocatable :: co(:), hybrid_co(:)
integer i

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

do i=1,lat%n_ele_track
  if( any(lat%ele(i)%key == [sextupole$,multipole$,wiggler$] )) then
    lat%ele(i)%select = .true.
  else
    lat%ele(i)%select = .false.
  endif
enddo

call twiss_and_track(lat,co)
do i=1,lat%n_ele_track
  write(400,*) lat%ele(i)%s, lat%ele(i)%a%beta, lat%ele(i)%b%beta
enddo
call make_hybrid_lat(lat, hybrid_lat, .true.)

call twiss_and_track(hybrid_lat,hybrid_co)
do i=1,hybrid_lat%n_ele_track
  write(500,*) hybrid_lat%ele(i)%s, hybrid_lat%ele(i)%a%beta, hybrid_lat%ele(i)%b%beta
enddo


call write_digested_bmad_file ('hybrid.lat', hybrid_lat)

end program












