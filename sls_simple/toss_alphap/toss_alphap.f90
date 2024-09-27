program toss_alphap
  use bmad
  use bmad_parser_mod, only: bp_com

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:), co(:)

  integer i
  integer, parameter :: npart=100

  real(rp) vec_offset(6)
  real(rp), parameter :: pzm=0.100
  real(rp) pz

  character*100 lat_file

  bp_com%always_parse = .true.

  call getarg(1,lat_file)
  call bmad_parser(lat_file, lat)
  call twiss_and_track(lat,co)

  allocate(orb(0:lat%n_ele_track))

  open(20,file='toss_alphap.dat')
  do i=1, npart
    pz = -pzm + 2*(i-1)*pzm/(npart-1)
    vec_offset = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, pz /)
    call init_coord(orb(0),co(0)%vec+vec_offset,lat%ele(0),element_end=upstream_end$)
    call track_all(lat,orb)
    write(20,'(i8,es15.6,6es15.6)') i, pz, orb(lat%n_ele_track)%vec
  enddo
  close(20)

end program





