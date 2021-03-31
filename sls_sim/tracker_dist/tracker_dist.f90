program tracker_dist_program
  use bmad
  use sls_lib
 
  implicit none

  type(lat_struct) lat
  type(coord_struct), allocatable :: co(:)
  type(coord_struct), allocatable :: orbit(:)
  type(ele_struct) eles
  type(em_field_struct) field

  integer i, j, septum_ix
  integer track_state
  integer status, ios
  integer n_part, n_survive

  logical ok, err

  character(100) in_file, lat_file, dist_file

  !parameters from .in
  integer ix_inj
  real(rp) init_vec(6), junk(2), t0, test_vec(6)
  integer n_turns

  namelist / tracker_dist / lat_file, ix_inj, dist_file, n_turns, septum_ix, t0, test_vec

  call getarg(1,in_file)

  ix_inj = -1

  open (unit = 10, file = in_file)
  read (10, nml = tracker_dist)
  close (10)

  write(*,*) "Preparing lattice..."

  call bmad_parser(lat_file, lat)

  call set_on_off(rfcavity$, lat, off$)
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  call twiss_and_track(lat,co,status)

  call set_ele_real_attribute (lat%ele(59), 'hkick', 0.660d-3, err, .true.)
  call set_ele_real_attribute (lat%ele(475), 'hkick', 1.0750d-03, err, .true.)
  call lattice_bookkeeper(lat)

  allocate(orbit(0:lat%n_ele_track))

  !call init_coord(orbit(ix_inj), co(ix_inj), lat%ele(ix_inj), downstream_end$, t_offset=t0) 
  call init_coord(orbit(ix_inj), co(ix_inj))
  orbit(ix_inj)%t = t0
  open(300,file='timing.dat')
  do i=1,n_turns
    call track_many(lat, orbit, ix_inj, ix_inj, 1, track_state=track_state)
    call em_field_calc(lat%ele(septum_ix), lat%param, 0.0d0, orbit(septum_ix), .true., field)
    write(300,'(i6,es15.7,es14.5)') i, orbit(septum_ix)%t, 1000.0d0*field%B(2)/6.671d0*1.8
  enddo
  close(300)

  !call init_coord(orbit(ix_inj), co(ix_inj)%vec+test_vec, lat%ele(ix_inj), downstream_end$, t_offset=t0) 
  call init_coord(orbit(ix_inj), co(ix_inj)%vec+test_vec)
  orbit(ix_inj)%t = t0
  open(300,file='first_10_turns.dat')
  do i=1,10
    call track_many(lat, orbit, ix_inj, ix_inj, 1, track_state=track_state)
    do j=ix_inj+1,lat%n_ele_track
      write(300,'(i6,7es14.5)') j, orbit(j)%t, orbit(j)%vec(1:6)
    enddo
    do j=1,ix_inj
      write(300,'(i6,7es14.5)') j, orbit(j)%t, orbit(j)%vec(1:6)
    enddo
    if(track_state /= moving_forward$) then
      write(*,*) "(first 10) Test particle lost at turn ", i
      exit
    endif
  enddo
  close(300)

  !call init_coord(orbit(ix_inj), co(ix_inj)%vec+test_vec, lat%ele(ix_inj), downstream_end$, t_offset=t0) 
  call init_coord(orbit(ix_inj), co(ix_inj)%vec+test_vec)
  orbit(ix_inj)%t = t0
  open(300,file='test_trajectory.dat')
  do i=1,n_turns
    call track_many(lat, orbit, ix_inj, ix_inj, 1, track_state=track_state)
    write(300,'(i6,7es14.5)') i, orbit(septum_ix)%t, orbit(septum_ix)%vec(1:6)
      if(track_state /= moving_forward$) then
      write(*,*) "Test particle lost at turn ", i
      exit
    endif
  enddo
  close(300)

  write(*,*) "Tracking..."
  n_part = 0
  n_survive = 0
  open(200,file=dist_file)
  do while(.true.)
    read(200,*, iostat=ios) init_vec(1:6), junk(2)
    if(ios .ne. 0) exit
    n_part = n_part + 1
    write(*,'(i6,a,6es14.5)') n_part, ": initial vector: ", init_vec
    !orbit(ix_inj)%vec = co(ix_inj)%vec + init_vec
    call init_coord(orbit(ix_inj), co(ix_inj)%vec+init_vec, lat%ele(ix_inj), downstream_end$, t_offset=t0)
    do i = 1, n_turns
      call track_many(lat, orbit, ix_inj, ix_inj, 1, track_state=track_state)
      if(track_state /= moving_forward$) then
        write(*,'(a,i6)') "Particle lost at turn ", i
        exit
      endif
    enddo
    if(track_state == moving_forward$) n_survive = n_survive + 1
  enddo
  close(200)

  write(*,'(a,i6)')    "Total particles tracked: ", n_part
  write(*,'(a,i6,a,f5.2,a)') "Total survived:          ", n_survive, " (", (100.0*n_survive)/n_part, "%)"
end program






