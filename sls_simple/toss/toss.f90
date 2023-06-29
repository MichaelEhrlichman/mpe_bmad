program toss
  use bmad
  use bmad_parser_mod, only: bp_com

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:)
  type(coord_struct) orb1

  integer i,j,k
  integer, parameter :: nsteps=2000

  real(rp) vec_offset(6)
  real(rp) s, delta_s

  character*100 lat_file

  bp_com%always_parse = .true.

  call getarg(1,lat_file)
  call bmad_parser(lat_file, lat)

  allocate(orb(0:lat%n_ele_track))

  vec_offset = (/ 0.0000, 0.0, 0.0000, 0.0, 0.0, 0.015 /)

  open(20,file='toss.dat')
  orb(0)%vec = vec_offset
  call init_coord(orb1,vec_offset,lat%ele(0),element_end=upstream_end$)
  call track_all(lat,orb)
  do i=0,lat%n_ele_track
    write(20,'(i8,f14.6,6es15.6,a,a)') i, lat%ele(i)%s, orb(i)%vec, "   ", lat%ele(i)%name
  enddo

  call init_coord(orb1,vec_offset,lat%ele(0),element_end=upstream_end$)
  orb1%vec = vec_offset
  s = 0.0d0
  delta_s = lat%param%total_length/nsteps
  open(21,file='toss_s.dat')
  do i=1,nsteps
    write(21,'(f14.5,6es15.6)') s, orb1%vec
    call track_from_s_to_s(lat,s,s+delta_s,orb1,orb1)
    s=s+delta_s
  enddo
  write(21,'(f14.5,6es15.6)') s, orb1%vec
  close(21)
end program





