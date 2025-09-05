program tracker_icer
  use bmad
 
  implicit none

  type(lat_struct) lat
  type(coord_struct), allocatable :: co(:), orbit(:)

  procedure(track1_custom_def) :: track1_custom

  integer i, n_turns, track_state, status

  character(100) lat_file

  track1_custom_ptr => track1_custom

  call getarg(1,lat_file)

  n_turns = 500000

  call bmad_parser(lat_file, lat)

  !call set_on_off(rfcavity$, lat, on$)
  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .true.

  call twiss_and_track(lat,co,status)

  allocate(orbit(0:lat%n_ele_track))

  open(100,file='tracker_simple.dat')
  write(100,'(a8,6a14,a14)') "# turn", "x", "px", "y", "py", "z", "pz", "track state"

  orbit(0)%vec = co(0)%vec + (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
  !orbit(0)%vec = (/-3.82749E-05, -2.10461E-05,  8.11201E-07, -1.08117E-06,  1.55813E-02, -7.74824E-04/)
  write(*,'(a,6es14.5)') "Closed Orbit:   ", co(0)%vec
  write(*,'(a,6es14.5)') "Initial Vector: ", orbit(0)%vec

  do i = 1, n_turns
    call track_all(lat, orbit, track_state=track_state)
    write(100,'(i8,6es14.5,i14)') i, orbit(0)%vec(1:6), track_state
    if(track_state /= moving_forward$) then
      write(*,*) "Particle lost at turn ", i
      exit
    endif
    orbit(0) = orbit(lat%n_ele_track)
    if (mod(i,10000) == 0) then
      write(*,'(a,i8,a,i8,a)') "Turn ", i, " of ", n_turns, " complete."
    endif
  enddo

  close(100)
end program






