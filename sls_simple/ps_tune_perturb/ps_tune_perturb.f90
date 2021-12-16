program ps_tune_perturb
  use bmad

  implicit none

  type (lat_struct) ring
  type (coord_struct), allocatable :: co(:)
  type(ele_struct) ele
  type (ele_pointer_struct), allocatable :: ps_eles(:)

  integer n_eles, n_seeds, n_time
  integer i,j,k
  integer status

  character(35) set_str

  real(rp) omega_t, ideal_val, rnum
  real(rp) magnet_ps_phase, net_pert
  real(rp) ppm, Qx0, Qy0, g
  real(rp) Qx, Qy
  real(rp) max_Qx, max_Qy, max_x, max_pz

  logical err

  character(60) lat_name

  call getarg(1,lat_name)
    
  call bmad_parser(lat_name, ring)
  
  call twiss_and_track(ring,co,status)
  Qx0 = ring%ele(ring%n_ele_track)%a%phi/twopi
  Qy0 = ring%ele(ring%n_ele_track)%b%phi/twopi

  call set_custom_attribute_name('PS_PHASE', err)
  call set_custom_attribute_name('K1_0', err)
  if(err) write(*,*) "Error in set_custom_attribute_name."

  call lat_ele_locator("quadrupole::* sbend::*", ring, ps_eles, n_eles)

  do j=1,n_eles
    write(set_str,'(a,es18.8)') 'K1_0=', value_of_attribute(ps_eles(j)%ele, 'K1', err, err_print_flag = .true.)
    call set_ele_attribute(ps_eles(j)%ele, set_str, err, .true.)
    if(err) write(*,*) "Error in set_ele_attribute."
  enddo

  call ran_seed_put(0)
  !n_seeds = 10000
  !n_time = 25
  n_seeds = 5
  n_time = 50
  ppm = 1.0e-6
  open(45,file='slow_errors.out')
  do k=1,n_seeds
    max_Qx = 0.0d0
    max_Qy = 0.0d0
    max_x = 0.0d0
    max_pz = 0.0d0
    do j=1,n_eles
      ! Random Distribution
      call ran_uniform(rnum)
      !rnum = 0
      write(set_str,'(a,es14.5)') 'PS_PHASE=',rnum*twopi

      ! Super Bad K1 Distribution
      !ideal_val = value_of_attribute(ps_eles(j)%ele, 'K1', err, err_print_flag = .true.)
      !write(set_str,'(a,es14.5)') 'PS_PHASE=', pi*(sign(1.0d0,ideal_val)+1.0) / 2.0

      ! Super Bad g Distribution
      !if( has_attribute(ps_eles(j)%ele, 'G')) then
      !  ideal_val = value_of_attribute(ps_eles(j)%ele, 'G', err, err_print_flag = .true.)
      !  write(set_str,'(a,es14.5)') 'PS_PHASE=', pi*(sign(1.0d0,ideal_val)+1.0) / 2.0
      !endif

      call set_ele_attribute(ps_eles(j)%ele, set_str, err, .true.)
      if(err) write(*,*) "Error in set_ele_attribute."
    enddo
    do i=1,n_time
      omega_t = twopi/(n_time-1.)*(i-1.)
      do j=1,n_eles
        ideal_val = value_of_attribute(ps_eles(j)%ele, 'K1_0', err, err_print_flag = .true.)
        if( .not. err ) then
          magnet_ps_phase = value_of_attribute(ps_eles(j)%ele, 'PS_PHASE' , err_print_flag = .true.)
          net_pert = ppm * sin(omega_t + magnet_ps_phase)
          write(set_str,'(a,es16.9)') 'K1=', ideal_val * (1.0d0 + net_pert)
          call set_ele_attribute(ps_eles(j)%ele, set_str, err, .true.)
          if(err) write(*,*) "Error from set_ele_attribute for K1."
          if( has_attribute(ps_eles(j)%ele, 'G_ERR')) then
            g = value_of_attribute(ps_eles(j)%ele, 'G', err, err_print_flag = .true.)
            if(err) write(*,*) "Error from set_ele_attribute for g."
            write(set_str,'(a,es16.9)') 'G_ERR=', net_pert*g
            call set_ele_attribute(ps_eles(j)%ele, set_str, err, .true.)
          endif
        else
          write(*,*) "ele is missing attribute"
        endif
      enddo
      call lattice_bookkeeper(ring)
      call twiss_and_track(ring,co,status)
      !write(55,'(es18.9,4es18.9)') (i-1.)/(n_time-1.), ring%ele(ring%n_ele_track)%a%phi/twopi, ring%ele(ring%n_ele_track)%b%phi/twopi, co(0)%vec(1), co(0)%vec(6)
      max_x = max(max_x,abs(co(0)%vec(1)))
      max_pz = max(max_pz,abs(co(0)%vec(6)))
      Qx = ring%ele(ring%n_ele_track)%a%phi/twopi
      Qy = ring%ele(ring%n_ele_track)%b%phi/twopi
      max_Qx = max(max_Qx,abs(Qx-Qx0)/Qx0)
      max_Qy = max(max_Qy,abs(Qy-Qy0)/Qy0)
    enddo
    write(55,*)
    write(55,*)
    write(*,'(i6,4es14.7)')  k, max_Qx, max_Qy, max_x, max_pz
    write(45,'(i6,4es14.7)') k, max_Qx, max_Qy, max_x, max_pz
  enddo
  close(45)
end program

