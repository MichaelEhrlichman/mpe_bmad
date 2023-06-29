program element_list

use bmad

implicit none

character lat_file*200
type(lat_struct) lat
integer i

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

open(45,file='element_list.dat')
do i=0,lat%n_ele_track
  if((lat%ele(i)%name(1:1) == 'Q' .or. lat%ele(i)%name(1:1) == 'B') .or. lat%ele(i)%name(1:1) == 'S') then
    write(45,'(i7,f14.6,f15.6,a,a,a,i3)') i, lat%ele(i)%s, lat%ele(i)%value(l$), '    ', lat%ele(i)%name, lat%ele(i)%alias, lat%ele(i)%slave_status
  endif
enddo

!super elements
write(45,*) "#"
write(45,*) "# Super Elements"
write(45,*) "#"
do i=lat%n_ele_track+1, lat%n_ele_max
  write(45,'(i7,f12.4,f13.4,a,a,i3)') i, lat%ele(i)%s, lat%ele(i)%value(l$), '    ', lat%ele(i)%name, lat%ele(i)%slave_status
enddo
close(45)

open(45,file='element_list.all')
do i=0,lat%n_ele_track
  write(45,'(i7,f14.6,f15.6,a,a,a,i3)') i, lat%ele(i)%s, lat%ele(i)%value(l$), '    ', lat%ele(i)%name, lat%ele(i)%alias, lat%ele(i)%slave_status
enddo
close(45)

end program












