PROGRAM co

use bmad

implicit none

character*100 lat_file
character*6 pz_str
integer i
real(rp) pz

type(lat_struct) lat
type(coord_struct), allocatable :: clo(:)

call getarg(1,lat_file)
call getarg(2,pz_str)
read(pz_str,*) pz

call bmad_parser(lat_file, lat)

allocate(clo(0:lat%n_ele_track))

clo(0)%vec = 0.0d0
clo(0)%vec(6) = pz

call set_on_off(rfcavity$,lat,off$)
bmad_com%radiation_damping_on=.false.
call closed_orbit_calc(lat,clo,4)
!call set_on_off(rfcavity$,lat,on$)
!bmad_com%radiation_damping_on=.true.
!call closed_orbit_calc(lat,clo,6)

open(20,file='co_'//trim(pz_str)//'.dat')
write(20,'(a11,6a14)') "# loc (m)", "x", "px", "y", "py", "z", "pz"
do i=1,lat%n_ele_track
  write(20,'(f11.4,6es14.5)') lat%ele(i)%s, clo(i)%vec(1:6)
enddo
close(20)

end program












