program cor_perturb

  use bmad

  implicit none

  type (lat_struct) ring
  type (coord_struct), allocatable :: co(:)
  type(ele_struct) ele
  type (ele_pointer_struct), allocatable :: cor_eles(:)

  integer n_eles, n_seeds, n_time
  integer i,j,k
  integer status

  character(35) set_str

  real(rp) omega_t, ideal_val, rnum
  real(rp) x_cor_phase, y_cor_phase, x_net_pert, y_net_pert
  real(rp) ppm
  real(rp) max_x, max_y

  logical err

  character(60) lat_name

  call getarg(1,lat_name)
    
  call bmad_parser(lat_name, ring)
  
  call twiss_and_track(ring,co,status)

  call set_custom_attribute_name('X_COR_PHASE', err)
  if(err) write(*,*) "Error in set_custom_attribute_name."
  call set_custom_attribute_name('Y_COR_PHASE', err)
  if(err) write(*,*) "Error in set_custom_attribute_name."

  call lat_ele_locator("ALIAS::CORR", ring, cor_eles, n_eles)

  call ran_seed_put(0)
  n_seeds = 10000
  n_time = 50
  ppm = 1.0e-6
  open(45,file='slow_errors.out')
  do k=1,n_seeds
    max_x = 0.0d0
    max_y = 0.0d0
    do j=1,n_eles
      ! Random Distribution
      !rnum = 0

      ! Super Bad K1 Distribution
      !ideal_val = value_of_attribute(cor_eles(j)%ele, 'K1', err, err_print_flag = .true.)
      !write(set_str,'(a,es14.5)') 'COR_PHASE=', pi*(sign(1.0d0,ideal_val)+1.0) / 2.0

      ! Super Bad g Distribution
      !if( has_attribute(cor_eles(j)%ele, 'G')) then
      !  ideal_val = value_of_attribute(cor_eles(j)%ele, 'G', err, err_print_flag = .true.)
      !  write(set_str,'(a,es14.5)') 'COR_PHASE=', pi*(sign(1.0d0,ideal_val)+1.0) / 2.0
      !endif

      call ran_uniform(rnum)
      write(set_str,'(a,es14.5)') 'X_COR_PHASE=',rnum*twopi
      call set_ele_attribute(cor_eles(j)%ele, set_str, err, .true.)
      if(err) write(*,*) "Error in set_ele_attribute."
      call ran_uniform(rnum)
      write(set_str,'(a,es14.5)') 'Y_COR_PHASE=',rnum*twopi
      call set_ele_attribute(cor_eles(j)%ele, set_str, err, .true.)
      if(err) write(*,*) "Error in set_ele_attribute."
    enddo
    do i=1,n_time
      omega_t = twopi/(n_time-1.)*(i-1.)
      do j=1,n_eles
        x_cor_phase = value_of_attribute(cor_eles(j)%ele, 'X_COR_PHASE' , err_print_flag = .true.)
        y_cor_phase = value_of_attribute(cor_eles(j)%ele, 'Y_COR_PHASE' , err_print_flag = .true.)
        x_net_pert = ppm * sin(omega_t + x_cor_phase)
        y_net_pert = ppm * sin(omega_t + y_cor_phase)
        write(set_str,'(a,es16.9)') 'HKICK=', x_net_pert
        call set_ele_attribute(cor_eles(j)%ele, set_str, err, .true.)
        if(err) write(*,*) "Error from set_ele_attribute for HKICK."
        write(set_str,'(a,es16.9)') 'VKICK=', y_net_pert
        call set_ele_attribute(cor_eles(j)%ele, set_str, err, .true.)
        if(err) write(*,*) "Error from set_ele_attribute for VKICK."
      enddo
      call lattice_bookkeeper(ring)
      call twiss_and_track(ring,co,status)
      max_x = max(max_x,abs(co(0)%vec(1)))
      max_y = max(max_y,abs(co(0)%vec(3)))
    enddo
    write(*,'(i6,2es14.7)')  k, max_x, max_y
    write(45,'(i6,2es14.7)') k, max_x, max_y
  enddo
  close(45)
end program

