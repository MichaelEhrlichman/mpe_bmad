program write_quads
  use bmad

  implicit none

  type(lat_struct) ring
  character*100 lat_file
  integer i, counter
  character(40) typed(50)

  typed = ''

  call getarg(1,lat_file)
  call bmad_parser(lat_file,ring)

  open(10,file='quads.bmad')
  counter = 1
  do i=1,ring%n_ele_track
    if(ring%ele(i)%key == quadrupole$) then
      if ( .not. any(typed == ring%ele(i)%name) ) then
        typed(counter) = trim(ring%ele(i)%name)
        counter = counter + 1
        write(10,'(a,es14.6)') trim(ring%ele(i)%name)//'[k1]=', ring%ele(i)%value(k1$)
      endif
    endif
  enddo
  close(10)

end program
