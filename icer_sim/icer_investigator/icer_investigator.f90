program icer_investigator
  use bmad
  use bmad_parser_mod, only: bp_com

  implicit none

  procedure(track1_custom_def) :: track1_custom
  procedure(make_mat6_custom_def) :: make_mat6_custom

  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:), co(:)
  type(normal_modes_struct) mode
  type(all_pointer_struct) v_ptr
  type(ele_pointer_struct), allocatable :: eles(:)

  integer i, j, n_loc
  integer, parameter :: npart=101

  real(rp) vec_offset(6), initial_vector(6)
  real(rp), parameter :: pzm=0.05
  real(rp) pz

  character*100 lat_file

  logical err_flag

  track1_custom_ptr => track1_custom
  make_mat6_custom_ptr => make_mat6_custom

  bp_com%always_parse = .true.

  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .false.
  !bmad_com%rel_tol_tracking = 1d-9  !Relative error. Default = 1d-8
  !bmad_com%abs_tol_tracking = 1d-11 !Absolute error. Default = 1d-10

  call getarg(1,lat_file)
  call bmad_parser(lat_file, lat)
  !call twiss_and_track(lat,co)
  call closed_orbit_calc(lat,co,5)
  call twiss_at_start(lat)
  call twiss_propagate_all(lat)

  open(20,file='closed_orbit.out')
  do i=1,lat%n_ele_track
    write(20,'(i6,6es14.4)') i, co(i)%vec
  enddo
  close(20)

  !call lat_ele_locator('INDUCTION_CELL', lat, eles, n_loc)
  !write(*,*) "number of induction cells found: ", n_loc
  !call pointer_to_attribute(eles(1)%ele,'VOLTAGE', .true., v_ptr, err_flag)

  call radiation_integrals (lat, co, mode)
  write(*,*) "energy loss per turn: ", mode%e_loss

  !write(*,*) "a: ", v_ptr%r
  !v_ptr%r = mode%e_loss
  !call set_flags_for_changed_attribute(eles(1)%ele, v_ptr%r)
  !call lattice_bookkeeper(lat)
  !call lat_make_mat6(lat, eles(1)%ele%ix_ele)
  !call closed_orbit_calc(lat,co,5)
  !call twiss_at_start(lat)
  !call twiss_propagate_all(lat)
  !write(*,*) "b: ", v_ptr%r

  allocate(orb(0:lat%n_ele_track))

  vec_offset = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0 /)
  initial_vector = co(0)%vec + vec_offset
  initial_vector(6) = 1.0e-6
  !initial_vector(5) = 1.0d-6
  call init_coord(orb(0),initial_vector,lat%ele(0),element_end=upstream_end$)
  write(*,'(a,6es15.6)') "closed orbit: ", co(0)%vec
  write(*,'(a,6es15.6)') "initial:      ", orb(0)%vec
  do i=1,20000
    call track_all(lat,orb)
    if(i==1) then
      open(20,file='turn_1.out')
      do j=1,lat%n_ele_track
        write(20,'(i6,6es14.4)') j, orb(j)%vec
      enddo
      close(20)
    endif
    write(*,'(i6,6es15.6)') i, orb(lat%n_ele_track)%vec
    orb(0) = orb(lat%n_ele_track)
  enddo

end program





