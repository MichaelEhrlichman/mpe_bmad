program bend_maker
  use bmad
 
  implicit none
 
  integer, parameter :: n_intra_steps = 10

  type(lat_struct) lat
  type(coord_struct), allocatable :: orbit(:)
  type(coord_struct) trajectory(0:n_intra_steps)

  integer i, j

  real(rp) gap0, tanf, bl0
  real(rp) g0tanf, xa, xb
  real(rp) bl_kick_avg, xold
  real(rp) vec(6)

  character(100) lat_file
  character(25) set_str

  logical err_flag

  call getarg(1,lat_file)

  call bmad_parser(lat_file, lat)
  
  bmad_com%auto_bookkeeper = .false.
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  allocate(orbit(0:lat%n_ele_track))

  vec(:) = 0.0d0
  call init_coord(orbit(0), vec, ele=lat%ele(1), element_end=upstream_end$,  particle=electron$)

  target_bend_angle = 0.1745329
  gap0 = 0.02
  tanf = tan(0.1745329)  ! 10 degrees
  g0tanf = gap0*tanf
  bl0 = value_of_attribute(lat%ele(1), 'bl_kick')
  do while(.true.)
    do i=1,lat%n_ele_track-1
      do while(.true.)
        xold = orbit(i)%vec(1)
        call intra_ele_trajectory(lat%ele(i), lat%param, orbit(i-1), trajectory)
        orbit(i) = trajectory(n_intra_steps)
        if(abs(xold-orbit(i)%vec(1)) .lt. 1e-12) exit
        bl_kick_avg = 0.0d0
        do j=1,n_intra_steps
          xa = trajectory(i-1)%vec(1)
          xb = trajectory(i)%vec(1)
          bl_kick_avg = bl_kick_avg + AvgB(xa,xb)
        enddo
        bl_kick_avg = bl_kick_avg / n_intra_steps
        write(set_str,'(a,f10.8)') "bl_kick=", bl_kick_avg
        write(*,*) "FOO set_str: ", set_str
        call set_ele_attribute(lat%ele(i), set_str, err_flag)
        call lattice_bookkeeper(lat)
      enddo
    enddo
    write(*,*) "orbit out: ", orbit(lat%n_ele_track-1)%vec(1:6)
  endif
  

  contains
    function AvgB(xa,xb) result(a)
      real(rp) xa, xb, a
      a = bl0*g0tanf/(xb-xa)*Log((g0tanf+xb)/(g0tanf+xa))
    end function

    subroutine intra_ele_trajectory(ele, param, coord0, trajectory)
      type(ele_struct) ele
      type(lat_param_struct) param
      type(coord_struct) coord0
      type(coord_struct) trajectory(0:n_intra_steps)
      real(rp) la, lb, dl
      integer i

      dl = ele%value(l$)/n_intra_steps

      trajectory(0) = coord0
      la = 0.0d0
      do i=1,n_intra_steps
        lb = la + dl 
        call twiss_and_track_intra_ele(ele,param,la,lb,.false.,.false.,trajectory(i-1),trajectory(i))
        la = lb
      enddo
    end subroutine

end program






