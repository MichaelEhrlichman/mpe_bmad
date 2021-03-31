program dynap_slim
  use bmad
  use custom_dynamic_aperture_mod
  use sls_lib
  !use calc_ring_mod
  use linear_aperture_mod

  use namelist_general !general: lat_file, use_hybrid

  implicit none

  type (lat_struct) ring
  type (lat_struct) ring0
  type (lat_struct) ring_hybrid
  type (custom_aperture_scan_struct) da_config
  type (custom_aperture_scan_struct) da_config_linear
  type (coord_struct), allocatable :: co(:)

  integer i,j
  integer n_dE
  integer iargs
  integer status
  integer n_aperture_test

  real(rp) linear_vec, da_vec, delta, metric, area

  logical err

  integer, parameter :: max_dE=3
  integer tracking_method
  integer n_adts
  integer n_turn
  integer n_angle
  integer track_dims
  real(rp) dE(max_dE)
  real(rp) max_angle
  real(rp) init_len
  real(rp) adts_x_min, adts_x_max
  real(rp) adts_y_min, adts_y_max

  namelist / slim /   tracking_method, &
                      max_angle, &
                      n_adts, &
                      n_turn, &          !Number of turns particle must survive
                      n_angle, &
                      track_dims, &            !either 4 or 6
                      dE, &
                      adts_x_min, &
                      adts_x_max, &
                      adts_y_min, &
                      adts_y_max, &
                      init_len

  character*60 in_file
  character*100 lat_file_override

  dE(:) = -999.
  n_turn = 100
  use_hybrid = .false.
  max_angle = pi

  iargs = iargc()
  lat_file_override = ''
  if( iargs == 1 ) then
    call getarg(1,in_file)
  elseif (iargs == 2 ) then
    call getarg(1,in_file)
    call getarg(2,lat_file_override)
  endif
    
  open (unit = 10, file = in_file, action='read')
  read (10, nml = general)
  read (10, nml = slim)
  close (10)
  use_hybrid = .false. !broke for very misaligned lattices

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  do i=1,max_dE
    if( dE(i) .lt. -100 ) then
      exit
    endif
  enddo
  n_dE = i-1

  call bmad_parser(lat_file, ring)
  IF( track_dims == 4 ) then
    call set_on_off(rfcavity$, ring, off$)
    bmad_com%radiation_damping_on = .false.
  elseif( track_dims == 6 ) then
    call set_on_off(rfcavity$, ring, on$)
    bmad_com%radiation_damping_on = .true.
  else
    write(*,*) "error setting track_dims.  terminate."
    stop
  endif
  bmad_com%radiation_fluctuations_on = .false.
  bmad_com%aperture_limit_on = .true.

  allocate(co(0:ring%n_ele_track))

  do i=1,ring%n_ele_track
    if(ring%ele(i)%key == wiggler$) then
      ring%ele(i)%value(x1_limit$) = 1.0
      ring%ele(i)%value(x2_limit$) = 1.0
      ring%ele(i)%value(y1_limit$) = 1.0
      ring%ele(i)%value(y2_limit$) = 1.0
      ring%ele(i)%aperture_type = elliptical$
    endif
  enddo

  n_aperture_test = 0
  do i=1,ring%n_ele_track
    if( ring%ele(i)%value(x1_limit$) .gt. 1e-4 ) n_aperture_test = n_aperture_test + 1
  enddo
  if( (1.0*n_aperture_test)/ring%n_ele_track .lt. 0.2 ) then
    write(*,*) "Less than 20% the elements have a x1 physical aperture."
    write(*,*) "Probably something is wrong.  Check that lattice file defines aperture."
    write(*,*) "Aborting"
    stop
  endif

  if(tracking_method .gt. 0) then
    do i = 1, ring%n_ele_max
      if(ring%ele(i)%key .ne. wiggler$) then
        ring%ele(i)%tracking_method = tracking_method
      endif
    enddo
  endif

  da_config%param%n_turn = n_turn
  da_config%param%accuracy = 0.00001d0
  da_config%param%init_len = init_len
  da_config%param%step_len = 0.00005d0
  da_config%min_angle = 0.01d0
  da_config%max_angle = max_angle-0.01d0
  da_config%n_angle = n_angle
  da_config%param%adts_x_min = adts_x_min
  da_config%param%adts_x_max = adts_x_max
  da_config%param%adts_y_min = adts_y_min
  da_config%param%adts_y_max = adts_y_max
  allocate(da_config%aperture(1:n_angle))

  da_config_linear%min_angle = da_config%min_angle
  da_config_linear%max_angle = da_config%max_angle
  da_config_linear%n_angle = da_config%n_angle
  allocate(da_config_linear%aperture(1:n_angle))

  if(use_hybrid) then
    do i=1,ring%n_ele_track
      if( (ring%ele(i)%key == sbend$) .or. &
          (ring%ele(i)%key == sextupole$) .or. &
          (ring%ele(i)%key == rfcavity$) .or. &
          (ring%ele(i)%key == multipole$) ) then
        ring%ele(i)%select = .true.
      else
        ring%ele(i)%select = .false.
      endif
    enddo
  endif

  write(*,*) "Starting DA calculation..."
  open(20,file='da.out')
  open(21,file='da_linear.out')
  ring0 = ring
  do i=1,n_dE
    write(*,*) "Calculating ", dE(i)
    da_config%param%n_adts = -1
    if(i==1) da_config%param%n_adts = n_adts

    !
    ! calculate dynamic aperture
    !
    ring = ring0
    co(0)%vec = 0.0d0
    if(track_dims .eq. 4) co(0)%vec(6) = dE(i)
    !call calc_ring(ring,track_dims,co,err)
    call twiss_and_track(ring,co,status)
    if(status .ne. ok$) then
      da_config%aperture(:)%x = 0.0d0
      da_config%aperture(:)%y = 0.0d0
      da_config%aperture(:)%i_turn = 0
      da_config_linear%aperture(:)%x = 0.0d0
      da_config_linear%aperture(:)%y = 0.0d0
    else
      da_config%param%closed_orbit%vec = co(0)%vec
      if(track_dims .eq. 6) da_config%param%closed_orbit%vec(6) = da_config%param%closed_orbit%vec(6) + dE(i)
      if(use_hybrid) then
        call make_hybrid_lat(ring, ring_hybrid)
        call custom_dynamic_aperture_scan(ring_hybrid, da_config)
      else
        write(*,*) "Starting da scan."
        call custom_dynamic_aperture_scan(ring, da_config)
        write(*,*) "DA scan complete."
      endif

      !
      ! calculate linear aperture
      !
      ring = ring0
      do j = 1, ring%n_ele_max
        if(ring%ele(j)%key .ne. wiggler$) then
          ring%ele(j)%tracking_method = 4    ! tracking method 4 is linear tracking.
        endif
      enddo
      co(0)%vec = 0.0d0
      co(0)%vec(6) = dE(i)
      !call calc_ring(ring,track_dims,co,err)
      call twiss_and_track(ring,co,status)

      da_config_linear%param%closed_orbit%vec = co(0)%vec
      if(use_hybrid) then
        call make_hybrid_lat(ring, ring_hybrid)
        call linear_aperture(ring_hybrid,da_config_linear)
      else
        call linear_aperture(ring,da_config_linear)
      endif
    endif

    ! calculate binary search metric used by MOGA
    metric = 0.0d0
    area = 0.0d0
    do j=1,n_angle
      linear_vec = sqrt( (da_config_linear%aperture(j)%x-da_config_linear%param%closed_orbit%vec(1))**2 + &
                         (da_config_linear%aperture(j)%y-da_config_linear%param%closed_orbit%vec(3))**2 )
      da_vec = sqrt( (da_config%aperture(j)%x-da_config%param%closed_orbit%vec(1))**2 + &
                     (da_config%aperture(j)%y-da_config%param%closed_orbit%vec(3))**2 )
      delta = (linear_vec - da_vec)/linear_vec/sqrt(1.0d0*n_angle)
      if(j .lt. n_angle-1) then
        area = area + 0.5*abs(da_config%aperture(j)%x*da_config%aperture(j+1)%y-da_config%aperture(j+1)%x*da_config%aperture(j)%y)
      endif
      if(delta .lt. 0) then
        delta = 0.0d0  !no contribution from exceeding physical aperture
      endif
      metric = metric + delta**2
    enddo
    write(20,'(a,f10.4,a,f15.7,es14.4)') "Optimizer metric and area (", dE(i), "%): ", metric, area

    ! output DA
    do j=1, n_angle
      write(20,'(3ES17.8,I11)') dE(i), da_config%aperture(j)%x, da_config%aperture(j)%y, da_config%aperture(j)%i_turn
      write(21,'(3ES17.8,I11)') dE(i), da_config_linear%aperture(j)%x, da_config_linear%aperture(j)%y, da_config_linear%aperture(j)%i_turn
    enddo 
    write(20,*)
    write(20,*)
    write(21,*)
    write(21,*)
  enddo
  close(20)
  close(21)

end program
                                                            










