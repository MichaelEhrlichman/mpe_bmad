program toss
  use bmad
  use bmad_parser_mod, only: bp_com
  use em_field_mod

  implicit none

  type(lat_struct) lat
  type(coord_struct) orb
  type(em_field_struct) field
  type (ele_pointer_struct), allocatable :: eles(:)

  integer i,j,k
  integer, parameter :: nsteps=2000
  integer n_loc

  logical err

  real(rp) s, delta_s

  character*100 lat_file

  bp_com%always_parse = .true.
  call getarg(1,lat_file)
  call bmad_parser(lat_file, lat)

  orb%vec(:) = 0.0d0
  orb%vec(1) = 0.00
  orb%vec(3) = 0.00

  open(20,file='field.dat')
  do i=1,nsteps
    s = 0.8/(nsteps-1)*(i-1)
    call lat_ele_locator('lgb_map1', lat, eles, n_loc, err)
    call em_field_calc(eles(1)%ele, lat%param, s, 0.0d0, orb, .true., field)
    write(20,'(i6,f11.5,3es14.4)') i, s, field%B(1:3)
  enddo
  close(20)

end program





