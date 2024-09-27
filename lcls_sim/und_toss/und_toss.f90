program und_toss
  use bmad

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:)

  integer i,j,k
  integer track_state

  logical err_flag

  real(rp) vec_offset(6)
  real(rp) xmin,ymin,xmax,ymax
  real(rp) xi,yi,xf,yf
  real(rp) dx, dy

  integer nx, ny

  character*100 lat_file

  call getarg(1,lat_file)
  call bmad_parser(lat_file, lat)

  allocate(orb(0:lat%n_ele_track))

  nx = 10
  ny = 10
  xmin = -0.001
  xmax =  0.001
  ymin = -0.001
  ymax =  0.001
  dx = (xmax-xmin)/(nx-1.0)
  dy = (ymax-ymin)/(ny-1.0)

  open(20,file='und.dat')
  write(20,'(4a14)') "# xi", "yi", "xf", "yf"
  do i=1,nx
    xi = xmin + dx*(i-1.0)
    do j=1,ny
      yi = ymin + dy*(j-1.0)
      vec_offset = (/ xi, 0.0d0, yi, 0.0d0, 0.0d0, 0.0d0 /)
      call init_coord(orb(0),vec_offset,lat%ele(0),element_end=upstream_end$)
      call track_all(lat,orb, track_state=track_state)
      if(track_state .ne. moving_forward$) then
        write(*,*) "particle lost at xi, yi = ", xi, yi
      else
        xf = orb(lat%n_ele_track)%vec(1)
        yf = orb(lat%n_ele_track)%vec(3)
        write(20,'(6es15.6)') xi, yi, xf, yf
      endif
    enddo
  enddo
  close(20)

end program





