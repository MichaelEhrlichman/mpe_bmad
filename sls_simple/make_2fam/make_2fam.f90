program make_2fam

use bmad

implicit none

character lat_file*200
character(10) leftNumber
type(lat_struct) lat
type(coord_struct), allocatable :: co(:)
integer i, status

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,co,status)

open(45,file='2fam.dat')
write(45,'(a)') "! horizontal"
do i=0,lat%n_ele_track
  if(lat%ele(i)%key == sextupole$) then
    if(lat%ele(i)%a%beta .gt. lat%ele(i)%b%beta) then
      write(leftNumber,'(f10.5)') lat%ele(i)%a%beta
      leftNumber = adjustl(leftNumber)
      !write(45,'(a,a,a)') trim(lat%ele(i)%name), '[k2]:', trim(leftNumber)//"*a"
      write(45,'(a,a,a)') trim(lat%ele(i)%name), '[k2]:', trim(leftNumber)
    endif
  endif
enddo

write(45,'(a)') "! vertical"
do i=0,lat%n_ele_track
  if(lat%ele(i)%key == sextupole$) then
    if(lat%ele(i)%a%beta .lt. lat%ele(i)%b%beta) then
      write(leftNumber,'(f10.5)') lat%ele(i)%b%beta
      leftNumber = adjustl(leftNumber)
      !write(45,'(a,a,a)') trim(lat%ele(i)%name), '[k2]:', trim(leftNumber)//"*a"
      write(45,'(a,a,a)') trim(lat%ele(i)%name), '[k2]:', trim(leftNumber)
    endif
  endif
enddo

close(45)

end program












