program portrait
  use bmad
  use sls_lib
 
  implicit none

  type(lat_struct) lat
  type(coord_struct), allocatable :: co(:)
  type(coord_struct), allocatable :: orbit(:)

  integer i, k

  real(rp) x_max, y_max

  character(100) lat_file
  character(100) in_file

  !parameters from .in
  integer ix_inj
  real(rp) init_vec(6)
  integer n_turns
  integer tracking_method
  integer dims
  integer n_steps

  namelist / portrait_params / lat_file, ix_inj, x_max, y_max, n_turns, n_steps, tracking_method, dims

  call getarg(1,in_file)

  tracking_method = -1

  open (unit = 10, file = in_file)
  read (10, nml = portrait_params)
  close (10)

  write(*,*) "Preparing lattice..."

  call bmad_parser(lat_file, lat)
  
  if(dims==4) then
    CALL set_on_off(rfcavity$, lat, off$)
  elseif (dims==6) then
    CALL set_on_off(rfcavity$, lat, on$)
  endif
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .true.
  
  if(tracking_method > -1) then
    do i = 1, lat%n_ele_max
      lat%ele(i)%tracking_method = tracking_method
    enddo
  endif

  call closed_orbit_calc(lat,co,dims)
  call lat_make_mat6(lat, -1, co)
  call twiss_at_start(lat)
  call twiss_propagate_all(lat)

  allocate(orbit(0:lat%n_ele_track))

  write(*,*) "Tracking..."

  open(100,file='portrait_x.dat')
  write(100,'(A6,4A14)') "# turn", "x", "px", "y", "py"
  open(200,file='portrait_x_Anorm.dat')
  write(200,'(A6,4A14)') "# turn", "Jx", "Jx'", "Jy", "Jy'"
  do k = 1,n_steps
    init_vec = 0.0
    init_vec(1) = x_max/(1.0d0*n_steps)*k
    init_vec(3) = 0.00001
    write(*,'(a,i6,a,i6)') "x step ", k, " of ", n_steps
    call track_turns()
    write(100,*)
    write(100,*)
    write(200,*)
    write(200,*)
  enddo
  close(100)
  close(200)

  open(100,file='portrait_y.dat')
  write(100,'(A6,4A14)') "# turn", "x", "px", "y", "py"
  open(200,file='portrait_y_Anorm.dat')
  write(200,'(A6,4A14)') "# turn", "Jx", "Jx'", "Jy", "Jy'"
  do k = 1,n_steps
    init_vec = 0.0
    init_vec(3) = y_max/(1.0d0*n_steps)*k
    init_vec(1) = 0.00001
    write(*,'(a,i6,a,i6)') "y step ", k, " of ", n_steps
    call track_turns()
    write(100,*)
    write(100,*)
    write(200,*)
    write(200,*)
  enddo
  close(100)
  close(200)

  contains

  subroutine track_turns()
    implicit none

    integer i
    integer track_state
    real(rp) tbt_coords(6)
    real(rp) norm_coords(6)
    logical ok

    orbit(ix_inj)%vec = co(ix_inj)%vec + init_vec
    do i = 1, n_turns
      call track_many(lat, orbit, ix_inj, ix_inj, 1, track_state=track_state)
      if(track_state /= moving_forward$) then
        exit
      endif
      tbt_coords(1:6) = orbit(ix_inj)%vec(1:6) - co(ix_inj)%vec(1:6)
      write(100,'(I6,4ES14.5)') i, tbt_coords(1:4)

      call xy_to_action(lat, ix_inj, tbt_coords, norm_coords, ok)
      if (.not. ok) then
        write(*,*) "Error from xy_to_action"
        stop
      endif
      write(200,'(I6,4ES14.5)') i, norm_coords(1:4)
    enddo
  end subroutine
end program






