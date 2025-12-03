program tracker_icer
  use bmad
  use bmad_parser_mod, only: bp_com
 
  implicit none

  type(lat_struct) lat
  type(coord_struct), allocatable :: co(:), orbit(:)

  procedure(track1_custom_def) :: track1_custom
  procedure(make_mat6_custom_def) :: make_mat6_custom

  integer i, n_turns, track_state, status
  integer part_ix

  character(100) lat_file, out_file
  character(4) part_ix_str

  real(rp) vec_offset(6), initial_vector(6)

  track1_custom_ptr => track1_custom
  make_mat6_custom_ptr => make_mat6_custom

  call getarg(1,lat_file)
  call getarg(2,part_ix_str)

  read(part_ix_str,*) part_ix

  n_turns = 20e6

  bp_com%always_parse = .true.

  call bmad_parser(lat_file, lat)

  !call set_on_off(rfcavity$, lat, on$)
  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .true.

  call twiss_and_track(lat,co,status)
  !call closed_orbit_calc(lat,co,5)
  !call twiss_at_start(lat)
  !call twiss_propagate_all(lat)

  allocate(orbit(0:lat%n_ele_track))

  write(out_file,'(a,i4.4,a)') 'tracker_icer_',part_ix,'.dat'
  write(*,*) out_file
  open(100,file=out_file)
  write(100,'(a8,6a14,a14)') "# turn", "x", "px", "y", "py", "z", "pz", "track state"

  !vec_offset = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
  !initial_vector = co(0)%vec + vec_offset
  initial_vector = 0.0d0
  call init_coord(orbit(0),initial_vector,lat%ele(0),element_end=upstream_end$)

  !write(*,'(a,6es14.5)') "Closed Orbit:   ", co(0)%vec
  write(*,'(a,6es14.5)') "Initial Vector: ", initial_vector

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






