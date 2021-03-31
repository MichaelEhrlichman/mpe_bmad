program tracker_adts
  use bmad
  use dynap_mod
  use sls_lib
  use linear_aperture_mod
  use adts_mod
  use custom_dynamic_aperture_mod

  use namelist_general !general: lat_file, use_hybrid
  use namelist_adts !adts: x_min, x_max, y_min, y_max, use_linear_bounds, n_steps

  implicit none

  type (lat_struct) ring
  type (lat_struct) ring_use
  type (coord_struct), allocatable :: co(:)
  type (coord_struct), allocatable :: orb(:)
  type (custom_aperture_scan_struct) da_config

  real(rp) tune_x, tune_y
  real(rp) nu_x, nu_y

  integer i
  integer iargs
  integer track_state
  integer halfn
  integer status

  character*100 in_file
  character*100 lat_file_override

  iargs = iargc()
  lat_file_override = ''
  if( iargs == 1 ) then
    call getarg(1,in_file)
  elseif (iargs == 2 ) then
    call getarg(1,in_file)
    call getarg(2,lat_file_override)
  endif

  use_linear_bounds = .false. !default

  open (unit = 10, file = in_file)
  read (10, nml = general)
  read (10, nml = adts)
  close (10)

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  call bmad_parser(lat_file, ring)
  call set_on_off(rfcavity$, ring, off$)
  call closed_orbit_calc(ring,co,4)
  call lat_make_mat6(ring, -1, co)

  if(use_hybrid) then
    do i=0,ring%n_ele_track
      if( (ring%ele(i)%key == sbend$) .or. &
          (ring%ele(i)%key == sextupole$) .or. &
          (ring%ele(i)%key == rfcavity$) .or. &
          (ring%ele(i)%key == multipole$) ) then
        ring%ele(i)%select = .true.
      else
        ring%ele(i)%select = .false.
      endif
    enddo
    call make_hybrid_lat(ring, ring_use)
  else
    ring_use = ring
  endif

  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.
  bmad_com%aperture_limit_on = .false.

  allocate(orb(0:ring_use%n_ele_track))

  co(0)%vec = 0.0d0
  co(0)%vec(6) = 0.0d0

  call twiss_and_track(ring_use,co,status)

  if(use_linear_bounds) then
    da_config%min_angle = 0.0001
    da_config%max_angle = pi-0.0001
    da_config%n_angle = 3
    allocate(da_config%aperture(da_config%n_angle))
    da_config%param%closed_orbit%vec = co(0)%vec
    call linear_aperture(ring_use,da_config)
    x_max = da_config%aperture(1)%x
    x_min = da_config%aperture(3)%x
    y_min = 0.00001d0
    y_max = da_config%aperture(2)%y
  endif

  if( mod(n_steps,2)==0 ) then
    n_steps = n_steps + 1
  endif

  halfn = (n_steps-1)/2

  open(45,file='tracker_adts.out')
  do i=1,halfn
    orb(0)%vec = 0.0d0
    orb(0)%vec(1) = x_min - x_min/halfn*(i-1.0)
    orb(0)%vec(3) = 1.0e-6
    call tune_ele_by_ele(ring_use,orb(0),200,track_state,nu_x,nu_y)
    write(45,'(4es17.8)') orb(0)%vec(1), orb(0)%vec(3), nu_x, nu_y
  enddo
  orb(0)%vec = 0.0d0
  orb(0)%vec(1) = 1.0d-6
  orb(0)%vec(3) = 1.0e-6
  call tune_ele_by_ele(ring_use,orb(0),200,track_state,nu_x,nu_y)
  write(45,'(4es17.8)') orb(0)%vec(1), orb(0)%vec(3), nu_x, nu_y
  do i=1,halfn
    orb(0)%vec = 0.0d0
    orb(0)%vec(1) = x_max/halfn*i
    orb(0)%vec(3) = 1.0e-6
    call tune_ele_by_ele(ring_use,orb(0),200,track_state,nu_x,nu_y)
    write(45,'(4es17.8)') orb(0)%vec(1), orb(0)%vec(3), nu_x, nu_y
  enddo
  write(45,*)
  write(45,*)

  do i=1,n_steps
    orb(0)%vec = 0.0d0
    orb(0)%vec(1) = 1.0e-6
    orb(0)%vec(3) = y_min + (y_max-y_min)/(n_steps-1.0)*(i-1.0)
    call tune_ele_by_ele(ring_use,orb(0),200,track_state,nu_x,nu_y)
    write(45,'(4es17.8)') orb(0)%vec(1), orb(0)%vec(3), nu_x, nu_y
  enddo
  write(45,*)
  write(45,*)

  do i=1,n_steps
    orb(0)%vec = 0.0d0
    orb(0)%vec(1) = x_min + (x_max-x_min)/(n_steps-1.0)*(i-1.0)
    orb(0)%vec(3) = y_min + (y_max-y_min)/(n_steps-1.0)*(i-1.0)
    call tune_ele_by_ele(ring_use,orb(0),200,track_state,nu_x,nu_y)
    write(45,'(4es17.8)') orb(0)%vec(1), orb(0)%vec(3), nu_x, nu_y
  enddo
  close(45)

end program
                                                            










