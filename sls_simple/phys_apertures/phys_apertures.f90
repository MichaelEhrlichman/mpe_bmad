program phys_apertures

use bmad

implicit none

character lat_file*200
type(lat_struct) lat
integer i

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

open(45,file='phys_apertures.dat')
write(45,'(a7,a12,4a14)') '# ix', 's', 'x1', 'x2', 'y1', 'y2'
do i=0,lat%n_ele_track
  write(45,'(i7,f12.4,4es14.5,2a)') i, lat%ele(i)%s, lat%ele(i)%value(x1_limit$), lat%ele(i)%value(x2_limit$), &
                                                    lat%ele(i)%value(y1_limit$), lat%ele(i)%value(y2_limit$), &
                                                    "   ", lat%ele(i)%name
enddo
close(45)

end program












