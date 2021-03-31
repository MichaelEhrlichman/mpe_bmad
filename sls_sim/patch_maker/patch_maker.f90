program patch_maker
  use bmad
  use sls_lib
  use bmad_parser_mod, only: bp_com

  implicit none

  type (lat_struct) ring
  type (coord_struct), allocatable :: orb(:)
  type (ele_pointer_struct), allocatable :: eles(:)

  integer i,j,k
  integer n_loc

  logical over, err

  character*28 var_str
  character*48 set_str
  character*100 lat_file

  real(rp) field_scale, accuracy
  real(rp) tOffset

  !bmad_com%convert_to_kinetic_momentum = .true.
  bp_com%always_parse = .true.

  call getarg(1,lat_file)
  call bmad_parser(lat_file, ring)
  write(*,*) "lattice parsing complete ..."

  allocate(orb(0:ring%n_ele_track))

  !for finding patch angle and offset

  call binary_search(.false.,0._rp,accuracy,.true.) !calling binary_search with last arg .true. resets its state variables
  field_scale = 1.1 ! initial guess
  write(*,*) "field_scale found: ", field_scale, accuracy
  call lat_ele_locator('lgb_map1', ring, eles, n_loc, err)
  !call lat_ele_locator('lgb_grid1', ring, eles, n_loc, err)
  do while (accuracy .gt. 1.0d-15)
    write(var_str,'(f28.18)') field_scale
    set_str = 'cartesian_map(1)%field_scale='//trim(adjustl(var_str))
    !set_str = 'grid_field(1)%field_scale='//trim(adjustl(var_str))
    do i=1,n_loc
      call set_ele_attribute (eles(i)%ele, set_str, err)
    enddo
    call lattice_bookkeeper(ring)

    orb(0)%vec = 0.0d0
    call track_all(ring,orb)

    over = (orb(eles(1)%ele%ix_ele+1)%vec(2) .lt. 0.0d0)

    call binary_search(over, field_scale, accuracy, .false.) !updates delta_m according to one iteration of binary search
    write(*,'(a,f20.15,es14.4)') "field_scale found: ", field_scale, accuracy
  enddo
  deallocate(eles)

  stop !FOO

  call binary_search(.false.,0._rp,accuracy,.true.) !calling binary_search with last arg .true. resets its state variables
  tOffset = 1.0d-10 ! initial guess
  write(*,*) "tOffset found: ", tOffset, accuracy
  call lat_ele_locator('pa_out', ring, eles, n_loc, err)
  do while (accuracy .gt. 1.0d-15)
    write(var_str,'(es28.18)') -1*tOffset
    set_str = 't_offset='//trim(adjustl(var_str))
    do i=1,n_loc
      call set_ele_attribute (eles(i)%ele, set_str, err)
    enddo
    call lattice_bookkeeper(ring)

    orb(0)%vec = 0.0d0
    call track_all(ring,orb)

    over = (orb(eles(1)%ele%ix_ele)%vec(5) .lt. 0.0d0)

    call binary_search(over, tOffset, accuracy, .false.) !updates delta_m according to one iteration of binary search
    write(*,'(a,es25.18,es14.4)') "tOffset found: ", tOffset, accuracy
  enddo

  write(*,*) "Done!"

end program





