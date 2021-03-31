program dynap_linear

  use bmad
  use bmad_routine_interface
  use custom_dynamic_aperture_mod
  use linear_aperture_mod

  use namelist_general !general: lat_file, use_hybrid

  implicit none

  type (lat_struct) ring
  type (lat_struct) ring0
  !type (aperture_scan_struct) da_config
  type (custom_aperture_scan_struct) da_config
  type (coord_struct), allocatable :: co(:)

  integer i,j
  integer n_dE
  integer iargs
  integer status

  real(rp) area
  real(rp) pz

  logical err

  !namelist slim
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
                      n_turn, & 
                      n_angle, &
                      track_dims, & 
                      dE, &   !from namelist slim, dynap_linear only uses dE
                      adts_x_min, &
                      adts_x_max, &
                      adts_y_min, &
                      adts_y_max, &
                      init_len

  character*60 in_file
  character*100 lat_file_override

  dE(:) = -999.
  max_angle = pi/2.0d0

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

  n_angle = 200

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  do i=1,max_dE
    if( dE(i) .lt. -100 ) THEN
      EXIT
    endif
  enddo
  n_dE = i-1

  call bmad_parser(lat_file, ring)

  call set_on_off(rfcavity$, ring, off$)
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  allocate(co(0:ring%n_ele_track))

  do i=1,ring%n_ele_track
    if(ring%ele(i)%key == wiggler$) then
      ring%ele(i)%value(x1_limit$) = 1.0
      ring%ele(i)%value(x2_limit$) = 1.0
      ring%ele(i)%value(y1_limit$) = 1.0
      ring%ele(i)%value(y2_limit$) = 1.0
    endif
  enddo

  bmad_com%aperture_limit_on = .true.

  !da_config%param%accuracy = 0.00001
  !da_config%param%step_len = 0.001
  da_config%min_angle = 0.0001
  da_config%max_angle = max_angle-0.0001
  da_config%n_angle = n_angle

  allocate(da_config%aperture(1:da_config%n_angle))

  write(*,*) "Starting LA calculation..."
  open(20,file='da.out')
  open(33,file='linear_boundary.dat')
  ring0 = ring
  do i=1,n_dE
    ring = ring0
    co(0)%vec = 0.0d0
    co(0)%vec(6) = dE(i)
    call closed_orbit_calc(ring,co,4,err_flag=err)
    if(err) then
      write(*,*) "Unable to calculate closed orbit.  Moving on."
    else
      da_config%param%closed_orbit%vec = co(0)%vec
      call lat_make_mat6(ring, -1, co,err_flag=err)
      if(err) then
        write(*,*) "Unable to calculate transfer matrices.  Moving on."
      else
        call twiss_at_start(ring,status)
        if( status .ne. ok$ ) then
          write(*,*) "Unable to calculate twiss at start.  Moving on."
          err = .true.
        else
          call twiss_propagate_all(ring,err_flag=err)
          if(err) then
            write(*,*) "Unable to propagate twiss.  Moving on."
          else
            write(*,*) "Calling linear_aperture"
            call linear_aperture(ring,da_config)
            write(*,*) "call complete"
          endif
        endif
      endif
    endif

    !make linear aperture bounds file
    write(33,'(3es14.5)') da_config%aperture(1)%x, da_config%aperture(da_config%n_angle)%x, da_config%aperture((da_config%n_angle+1)/2)%y

    if(err) then
      da_config%aperture(:)%x = 0.0d0
      da_config%aperture(:)%y = 0.0d0
      da_config%aperture(:)%i_turn = 0
    endif

    area = 0.0d0
    do j=1,da_config%n_angle-1
      area = area + 0.5d0*abs( da_config%aperture(j)%x*da_config%aperture(j+1)%y - da_config%aperture(j+1)%x*da_config%aperture(j)%y)
    enddo
    write(20,'(a,f6.3,a,es14.4,a)') "# Area(", dE(i)*100,"%) = ", area, " m^2"
    do j=1, da_config%n_angle
      write(20,'(3ES17.8,I11)') dE(i), da_config%aperture(j)%x, da_config%aperture(j)%y, da_config%aperture(j)%i_turn
    enddo 
    write(20,*)
    write(20,*)

    if(i==1) then
      open(22,file='bounds_for_adts.out')
      write(22,'(3es14.5)') da_config%aperture(1)%x, da_config%aperture(da_config%n_angle)%x, da_config%aperture((da_config%n_angle+1)/2)%y
      close(22)
    endif

  enddo
  close(20)
  close(33)

  ! open(20,file='da_pz.out')
  ! da_config%min_angle = 0.0001
  ! da_config%max_angle = pi-0.0001
  ! da_config%n_angle = 3
  ! deallocate(da_config%aperture)
  ! allocate(da_config%aperture(1:da_config%n_angle))
  ! ring0 = ring
  ! do i=1,100
  !   ring = ring0
  !   co(0)%vec = 0.0d0
  !   pz = -0.05 + 0.10*(i-1.0)/99.0
  !   co(0)%vec(6) = pz
  !   call closed_orbit_calc(ring,co,4,err_flag=err)
  !   if(.not. err) then
  !     da_config%param%closed_orbit%vec = co(0)%vec
  !     call lat_make_mat6(ring, -1, co,err_flag=err)
  !     if(.not. err) then
  !       call twiss_at_start(ring,status)
  !       err = (status .ne. ok$)
  !       if(.not. err) then
  !         call twiss_propagate_all(ring,err_flag=err)
  !         if(.not. err) then
  !           call linear_aperture(ring,da_config)
  !         endif
  !       endif
  !     endif
  !   endif

  !   if(.not. err) then
  !     write(20,'(4es14.5)') pz, da_config%aperture(1)%x, da_config%aperture(3)%x, da_config%aperture(2)%y
  !   else
  !     write(20,'(4es14.5)') pz, 0.0, 0.0, 0.0
  !   endif
  ! enddo
  ! close(20)
end program
                                                            










