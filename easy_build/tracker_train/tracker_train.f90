program tracker_train
  use bmad
  use bmad_parser_mod, only: bp_com
  use beam_mod
  use random_mod
  use mode3_mod
 
  implicit none

  character(100) lat_file
  character(100) in_file
  character(180) line
  character(20) filename
  character(40) TFB_pickup_name

  type(lat_struct) lat, hybrid
  type(coord_struct), allocatable :: co(:), coh(:)
  type(coord_struct), allocatable :: orbit(:)
  type(beam_init_struct) beam_init
  type(beam_struct) beam
  type(ele_pointer_struct), allocatable :: rf_eles(:)
  type(ele_struct), pointer :: wake_ele => null()
  type(ele_pointer_struct), allocatable :: pickup_ele(:)

  integer i, j, k
  integer status
  integer taylor_order
  integer n_rf, n_loc
  integer TFB_mask_x(500)
  integer TFB_mask_y(500)

  logical ok, radiation, err

  real(rp) phase
  real(rp) harmon
  real(rp) ax(20), TFB_kick_x, max_kick_x
  real(rp) ay(20), TFB_kick_y, max_kick_y
  real(rp) ax_use(20), ay_use(20)
  real(rp), allocatable :: x(:,:), xnew(:)
  real(rp), allocatable :: y(:,:), ynew(:)
  real(rp) t6(6,6), N(6,6), abz_tunes(3)

  !parameters from .in
  logical hybridize, tfb_on, alpha_damp_on
  real(rp) gain_x, gain_y, alpha_damp(3), TFB_saturation
  integer ix_inj
  integer n_turns
  integer n_part, n_bunches
  integer hyb_TFB_pickup
  integer init_mode
  integer Ntaps_x, Ntaps_y
  real(rp) init_vec(6), r, init_J(3)
  real(rp) bunch_charge
  real(rp) Gx, Qx, phase_x, use_phase_x
  real(rp) Gy, Qy, phase_y, use_phase_y
  character(4) train_init ! 'vec', 'ran', 'mod'
  character(9) coeffs
  

  namelist / parameters / lat_file, hybridize, tfb_on, ix_inj, init_vec, n_turns, &
                          radiation, n_part, n_bunches, bunch_charge, &
                          gain_x, gain_y, &
                          train_init, init_mode, init_J, &
                          TFB_pickup_name, TFB_mask_x, TFB_mask_y, alpha_damp_on, alpha_damp, &
                          ax, ay, Gx, Qx, phase_x, Ntaps_x, Gy, Qy, phase_y, Ntaps_y, coeffs, TFB_saturation


  ! Gain=1.000; Phase=0; Freq=0.22; Taps=11
  !ax(:) = 0.0d0
  !ax = FIR_coeffs(1.0d0,0.22d0,0.0d0,11)*32767 
  !ax(1:11) = (/-2465.0,27358.0,8711.0,-28100.0,-23249.0,15380.0,25006.0,-10015.0,-32767.0,-6270.0,26410.0/)
  !ax = ax/32767.0

  ! Gain=1.000; Phase=0; Freq=0.33; Taps=11
  !ay(:) = 0.0d0
  !ax = FIR_coeffs(1.0d0,0.3285d0,0.0d0,11)*32767 
  !ay(1:11) = (/-2784.0,28877.0,-32767.0,-6054.0,30295.0,-30838.0,-9298.0,31437.0,-28677.0,-12487.0,32296.0/)
  !ay = ay/32767.0
      
  call getarg(1,in_file)

  coeffs = 'specified'
  Gx = 0
  Qx = 0
  phase_x = 0
  Ntaps_x = 0
  Gy = 0
  Qy = 0
  phase_y = 0
  Ntaps_y = 0
  ax = 0
  ay = 0
  gain_x = 0.0d0
  gain_y = 0.0d0
  alpha_damp = 0.0d0
  alpha_damp_on = .false.
  TFB_pickup_name = 'none'
  lat_file = ''
  hybridize = .true.
  tfb_on = .true.
  ix_inj = 0
  init_vec = 0.0d0
  init_J = 0.0d0
  init_mode = 0
  n_turns = 10
  radiation = .false.
  train_init = 'ran'
  n_part = 10
  n_bunches = 1
  bunch_charge = 1.0e-9
  TFB_mask_x = -1
  TFB_mask_y = -1
  TFB_saturation = 1.0

  open(unit=10, file=in_file)
  read(10, nml = parameters)
  close(10)

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, lat)
  bmad_com%radiation_damping_on = radiation
  bmad_com%radiation_fluctuations_on = radiation
  
  call twiss_and_track(lat,co)

  beam_init%distribution_type(1:3) = 'ran_gauss'
  beam_init%n_particle = n_part  !number of particles per bunch
  beam_init%n_bunch = n_bunches
  beam_init%bunch_charge = bunch_charge
  beam_init%a_emit = 1.806E-09
  beam_init%b_emit = 1.0e-12
  beam_init%sig_e = 8.527E-04
  beam_init%sig_z = 0.00571
  beam_init%species = 'electron'
  beam_init%full_6D_coupling_calc = .true.
  beam_init%center = co(ix_inj)%vec

  call lat_ele_locator('rfcav', lat, rf_eles, n_rf)
  beam_init%dt_bunch = 1.0d0 / value_of_attribute(rf_eles(1)%ele,'RF_FREQUENCY')
  harmon = value_of_attribute(rf_eles(1)%ele,'HARMON')
  write(*,'(a,es14.6,a)') "Bunch spacing: ", beam_init%dt_bunch, ' seconds'
  write(*,'(a,es14.6)') "Harmonic Number: ", harmon

  open(45,file='~bunch_vs_number.out')
  write(45,'(a12,a6,12a14)') "# turn", "bunch", "x0", "x0'", "y0", "y0'", "z0", "z0'", "r<xx>", "r<x'x'>", "r<yy>", "r<y'y'>", "r<zz>", "r<z'z'>"
  open(46,file='~bunch_vs_number.env')
  write(46,'(a12,a6,12a14)') "# turn", "bunch", "a0", "a0'", "b0", "b0'", "c0", "c0'"

  open(50,file='~max_kick.out')
  write(50,'(a)') "Max kick in radians for kicker in center of injection straight."
  write(50,'(a8,2a14)') "# turn", "max_kick_x", "max_kick_y"

  do i=1,n_bunches
    write(filename,'(a,i0.4,a)') '~bunch_',i,'.out'
    open(1000+i,file=filename)
    write(filename,'(a,i0.4,a)') '~bunch_',i,'.env'
    open(3000+i,file=filename)
  enddo

  lat%ele(:)%select = .false.
  !mark pickup locations before hybridizing
  if(tfb_on) then
    call lat_ele_locator(TFB_pickup_name,lat,pickup_ele,n_loc)
    if(n_loc .ne. 1) then
      write(*,*) "Error locating TFB_pickup_name. n_loc: ", n_loc
      stop
    endif
    pickup_ele(1)%ele%ix_pointer = 1
    pickup_ele(1)%ele%select = .true.
  else
    call lat_ele_locator('BEGINNING',lat,pickup_ele,n_loc)
    if(n_loc .ne. 1) then
      write(*,*) "Error locating BEGINNING. n_loc: ", n_loc
      stop
    endif
    pickup_ele(1)%ele%ix_pointer = 1
    pickup_ele(1)%ele%select = .true.
  endif

  do i=1,lat%n_ele_track
    if(associated(lat%ele(i)%wake)) then
      wake_ele => lat%ele(i)
      write(*,'(a,a,a)') "Element ", lat%ele(i)%name, " has wake."
      lat%ele(i)%select = .true.
    elseif (any(lat%ele(i)%key==(/sbend$,rbend$/)) .and. radiation) then
      lat%ele(i)%select = .true.
    endif
  enddo
  if(hybridize) then
    call make_hybrid_lat(lat, hybrid, .true.)
  else
    hybrid = lat
  endif
  call twiss_and_track(hybrid,coh)
  write(*,'(a,f14.5)') "βx(0) = ", hybrid%ele(0)%a%beta 
  call calc_z_tune(hybrid)
  call transfer_matrix_calc (hybrid, t6, ix1=0, one_turn=.true.)
  abz_tunes(1) = hybrid%a%tune
  abz_tunes(2) = hybrid%b%tune
  abz_tunes(3) = hybrid%z%tune
  call make_N(t6, N, err, abz_tunes)
  if(err) then
    write(*,*) "make_N error"
    stop
  endif

  if(hybridize) then
    taylor_order = -1
    do i=1,hybrid%n_ele_track
      if(associated(hybrid%ele(i)%taylor(1)%term)) then
        taylor_order = 0
        do k=1,6
          do j=1, size(hybrid%ele(i)%taylor(k)%term)
            taylor_order = max(taylor_order, sum(hybrid%ele(i)%taylor(k)%term(j)%expn))
          enddo
        enddo
        exit
      endif
    enddo
    write(*,*) "taylor_order: ", taylor_order
  endif

  do i=1,hybrid%n_ele_track
    if(hybrid%ele(i)%ix_pointer == 1) exit
  enddo
  hyb_TFB_pickup = i

  write(*,*) "Horizontal pickup-to-kicker phase Δ (units 1): ", (hybrid%ele(hybrid%n_ele_track)%a%phi - hybrid%ele(hyb_TFB_pickup)%a%phi)/twopi
  write(*,*) "Vertical pickup-to-kicker phase Δ (units 1):   ", (hybrid%ele(hybrid%n_ele_track)%b%phi - hybrid%ele(hyb_TFB_pickup)%b%phi)/twopi

  do i=1,hybrid%n_ele_track
    if(associated(hybrid%ele(i)%wake)) then
      wake_ele => hybrid%ele(i)
    endif
  enddo

  call init_beam_distribution(hybrid%ele(ix_inj), lat%param, beam_init, beam)

  if(train_init == 'vec') then
    do i=1,n_bunches
      do j=1,n_part
        beam%bunch(i)%particle(j)%vec(:) = beam%bunch(i)%particle(j)%vec(:) + init_vec
      enddo
    enddo
  elseif(train_init == 'ran') then
    call ran_seed_put(seed=0)
    do i=1,n_bunches
      do j=1,n_part
        call ran_gauss(r)
        beam%bunch(i)%particle(j)%vec(:) = beam%bunch(i)%particle(j)%vec(:) + r*init_vec
      enddo
    enddo
  elseif(train_init == 'mod') then
    do i=1,n_bunches
      phase = pi/2.0d0
      !init_vec(1) = init_vec(1) + sqrt(init_J(1)*hybrid%ele(ix_inj)%a%beta)*cos(phase+twopi*init_mode*i/(n_bunches*1.0d0))
      !init_vec(2) = init_vec(2) -1.0d0*sqrt(init_J(1)/hybrid%ele(ix_inj)%a%beta)* &
      !             (hybrid%ele(ix_inj)%a%alpha*cos(phase+twopi*init_mode*i/(n_bunches*1.0d0)) + &
      !                                         sin(phase+twopi*init_mode*i/(n_bunches*1.0d0))) 
      !init_vec(3) = init_vec(3) + sqrt(init_J(2)*hybrid%ele(ix_inj)%b%beta)*cos(phase+twopi*init_mode*i/(n_bunches*1.0d0))
      !init_vec(4) = init_vec(4) -1.0d0*sqrt(init_J(2)/hybrid%ele(ix_inj)%b%beta)* &
      !             (hybrid%ele(ix_inj)%b%alpha*cos(phase+twopi*init_mode*i/(n_bunches*1.0d0)) + &
      !                                         sin(phase+twopi*init_mode*i/(n_bunches*1.0d0)))
      init_vec(1) = sqrt(init_J(1)*hybrid%ele(ix_inj)%a%beta)*cos(phase+twopi*init_mode*i/(n_bunches*1.0d0))
      init_vec(2) = -1.0d0*sqrt(init_J(1)/hybrid%ele(ix_inj)%a%beta)* &
                   (hybrid%ele(ix_inj)%a%alpha*cos(phase+twopi*init_mode*i/(n_bunches*1.0d0)) + &
                                               sin(phase+twopi*init_mode*i/(n_bunches*1.0d0))) 
      init_vec(3) = sqrt(init_J(2)*hybrid%ele(ix_inj)%b%beta)*cos(phase+twopi*init_mode*i/(n_bunches*1.0d0))
      init_vec(4) = -1.0d0*sqrt(init_J(2)/hybrid%ele(ix_inj)%b%beta)* &
                   (hybrid%ele(ix_inj)%b%alpha*cos(phase+twopi*init_mode*i/(n_bunches*1.0d0)) + &
                                               sin(phase+twopi*init_mode*i/(n_bunches*1.0d0)))
      init_vec(5:6) = 0 !longitudinal not yet implemented.
      do j=1,n_part
        beam%bunch(i)%particle(j)%vec(:) = beam%bunch(i)%particle(j)%vec(:) + init_vec
      enddo
    enddo
  endif

  if(associated(wake_ele)) then
    do i=1,size(wake_ele%wake%lr%mode)
      write(filename,'(a,i0.4,a)') '~wake_',i,'.out'
      open(2000+i,file=filename)
    enddo
    open(2999,file='~wake_sum.out')
  endif

  allocate(xnew(n_bunches))
  allocate(ynew(n_bunches))
  allocate(x(n_bunches,20))
  allocate(y(n_bunches,20))

  call report(lat,beam,0)
  if(associated(wake_ele)) call report_wake(wake_ele,0)

  if(alpha_damp_on .and. radiation) then
    write(*,*) "Inconsistent .in file: alpha_damp_on and radiation should not both be true"
    stop
  endif

  if(coeffs == 'specified') then
    ax_use = ax/32767.0
    ay_use = ay/32767.0
  elseif(coeffs == 'generated') then
    use_phase_x = (hybrid%ele(hybrid%n_ele_track)%a%phi-hybrid%ele(hyb_TFB_pickup)%a%phi) / twopi*360.0d0 + 90.0d0 + phase_x
    use_phase_y = (hybrid%ele(hybrid%n_ele_track)%b%phi-hybrid%ele(hyb_TFB_pickup)%b%phi) / twopi*360.0d0 + 90.0d0 + phase_y
    ax_use = FIR_coeffs(Gx,Qx,use_phase_x,Ntaps_x)
    ay_use = FIR_coeffs(Gy,Qy,use_phase_y,Ntaps_y)
  else
    write(*,*) "coeffs setting invalid"
    stop
  endif

  open(100,file='coefficients.out')
  write(100,'(a)') "# ax"
  do i=1,20
    write(100,'(f16.4)') ax_use(i)
  enddo
  write(100,*)
  write(100,*)
  write(100,'(a)') "# ay"
  do i=1,20
    write(100,'(f16.4)') ay_use(i)
  enddo
  close(100)

  x(:,:) = 0.0d0
  y(:,:) = 0.0d0
  do j=1,n_turns
    if(mod(j,100)==0) write(*,'(a,i8)') "Tracking turn ", j
    if(tfb_on) then
      call track_beam(hybrid,beam,hybrid%ele(0),hybrid%ele(hyb_TFB_pickup),err=err)
      do i=1,n_bunches
        xnew(i) = beam%bunch(i)%particle(1)%vec(1)
        ynew(i) = beam%bunch(i)%particle(1)%vec(3)
      enddo
      call track_beam(hybrid,beam,hybrid%ele(hyb_TFB_pickup),hybrid%ele(hybrid%n_ele_track),err=err)
    else
      call track_beam(hybrid,beam,err=err)
    endif
    call report(lat,beam,j)
    if(associated(wake_ele)) call report_wake(wake_ele,j)
    if(any(beam%bunch(:)%n_live .ne. n_part)) then
      do i=1,n_bunches
        if(beam%bunch(i)%n_live .ne. n_part) then
          write(*,'(a,i6)') "bunch lost: ", i
          write(*,'(a,i6,a)') "  lost at element: ", beam%bunch(i)%particle(1)%state, hybrid%ele(beam%bunch(i)%particle(1)%state)%name
        endif
      enddo
      exit
    endif
    if(tfb_on) then
      max_kick_x = 0.0d0
      max_kick_y = 0.0d0
      do i=1,n_bunches
        TFB_kick_x = gain_x*FIR(xnew(i), x(i,:), ax_use)
        TFB_kick_y = gain_y*FIR(ynew(i), y(i,:), ay_use)
        TFB_kick_x = sign(1.0d0,TFB_kick_x)*min(abs(TFB_kick_x),TFB_saturation)
        TFB_kick_y = sign(1.0d0,TFB_kick_y)*min(abs(TFB_kick_y),TFB_saturation)
        if(.not. any(i==TFB_mask_x)) then
          max_kick_x = max(max_kick_x,abs(TFB_kick_x))
          beam%bunch(i)%particle(1)%vec(1) = beam%bunch(i)%particle(1)%vec(1) + TFB_kick_x
        endif
        if(.not. any(i==TFB_mask_y)) then
          max_kick_y = max(max_kick_y,abs(TFB_kick_y))
          beam%bunch(i)%particle(1)%vec(3) = beam%bunch(i)%particle(1)%vec(3) + TFB_kick_y
        endif
      enddo
      write(50,'(i8,2es14.4)') j, max_kick_x, max_kick_y
    endif
    if(alpha_damp_on) then
      do i=1,n_bunches
        call rad_damp(hybrid,beam%bunch(i),alpha_damp)
      enddo
    endif
  enddo

  close(45)
  close(46)
  close(50)
  do i=1,n_bunches
    close(1000+i)
    close(2000+i)
    close(2999)
    close(3000+i)
  enddo

  contains
    subroutine rad_damp(lat,bunch,alpha_damp)
      !use mode3_mod
      type(lat_struct) lat
      type(bunch_struct) bunch
      real(rp) alpha_damp(3), J(6)
      real(rp) T
      integer i
      logical error

      do i=1,size(bunch%particle)
        !call xyz_to_action(lat,0,bunch%particle(i)%vec,J,error)
        J = matmul(mat_symp_conj(N), bunch%particle(i)%vec)
        !write(*,*) "Jx: ", (lat%ele(0)%a%gamma*bunch%particle(i)%vec(1)**2+ &
        !                    2.0d0*lat%ele(0)%a%alpha*bunch%particle(i)%vec(1)*bunch%particle(i)%vec(2) + &
        !                    lat%ele(0)%a%beta*bunch%particle(i)%vec(2)**2)/2.0d0, (J(1)**2 + J(2)**2)/2.0d0
        if(error) then
          write(*,*) "xyz_to_action error"
          stop
        endif
        J(1) = exp(-2.0d0*alpha_damp(1))*J(1)
        J(2) = exp(-2.0d0*alpha_damp(1))*J(2)
        J(3) = exp(-2.0d0*alpha_damp(2))*J(3)
        J(4) = exp(-2.0d0*alpha_damp(2))*J(4)
        J(5) = exp(-2.0d0*alpha_damp(3))*J(5)
        J(6) = exp(-2.0d0*alpha_damp(3))*J(6)
        !call action_to_xyz(lat,0,J,bunch%particle(i)%vec,error)
        bunch%particle(i)%vec = matmul(N, J)
        if(error) then
          write(*,*) "action_to_xyz error"
          stop
        endif
      enddo
    end subroutine

    function FIR(xnew,x,a) result(y)
      real(rp) xnew, x(20)
      real(rp) a(20)
      real(rp) y

      x(2:20) = x(1:19)
      x(1) = xnew

      y = dot_product(x,a)
    end function

    function FIR_coeffs(G,Q,phi,N) result(a)
      real(rp) G, Q, phi
      real(rp) a(20)
      real(rp) f(20)
      real(rp) offset
      integer i, N

      a = 0.0d0
      f = 0.0d0

      do i=1,N
        f(i) = sin(twopi*Q*(i-1) + phi/360.d0*twopi)
      enddo
      offset = -1.0d0*sum(f(1:N))/N
      a(1:N) = offset + G*f(1:N)
      a(1:N) = G*a(1:N)/maxval(abs(a(1:N)))
    end function

    subroutine report_wake(ele,turn)
      type(ele_struct) ele
      integer i, turn
      real(rp) a_sin_sum, a_cos_sum, b_sin_sum, b_cos_sum

      a_sin_sum = 0.0d0
      a_cos_sum = 0.0d0
      b_sin_sum = 0.0d0
      b_cos_sum = 0.0d0
      do i=1, size(ele%wake%lr%mode)
        write(2000+i,'(i8,4es14.5)') turn, ele%wake%lr%mode(i)%a_sin, ele%wake%lr%mode(i)%a_cos, ele%wake%lr%mode(i)%b_sin, ele%wake%lr%mode(i)%b_cos
        a_sin_sum = a_sin_sum + ele%wake%lr%mode(i)%a_sin
        a_cos_sum = a_cos_sum + ele%wake%lr%mode(i)%a_cos
        b_sin_sum = b_sin_sum + ele%wake%lr%mode(i)%b_sin
        b_cos_sum = b_cos_sum + ele%wake%lr%mode(i)%b_cos
      enddo
      write(2999,'(i8,4es14.5)') turn, a_sin_sum, a_cos_sum, b_sin_sum, b_cos_sum
    end subroutine

    subroutine report(lat,beam,turn)
      !use mode3_mod
      type(lat_struct) lat
      type(beam_struct) beam
      integer i, turn
      real(rp) centroid(6), moments(6), J(6)
      logical ok, error
      type(coord_struct) norm_coords

      do i=1,size(beam%bunch(:))
        !call xy_to_action(lat, i, beam%bunch(i)%particle(1)%vec, norm_coords%vec, ok)
        !if(.not. ok) then
        !  write(*,*) "xy_to_action is not ok."
        !  error stop
        !endif

        call make_moments(beam%bunch(i)%particle(:), centroid, moments)
        J = matmul(mat_symp_conj(N), centroid)
        write(45,'(i12,i6,6es14.5)') turn, i, centroid(1:6)
        write(46,'(i12,i6,6es14.5)') turn, i, J(1:6)
        write(1000+i,'(i12,i6,6es14.5)') turn, i, centroid(1:6)
        write(3000+i,'(i12,i6,6es14.5)') turn, i, J(1:6)
      enddo
      write(45,*)
      write(45,*)
      write(46,*)
      write(46,*)
    end subroutine

    subroutine make_moments(particles,centroid,moments)
      implicit none

      type(coord_struct) particles(:)
      real(rp) particles_cent(size(particles(:)),6)
      real(rp) moments(6)
      real(rp) centroid(6)

      integer j
      integer nparts

      nparts = size(particles(:))

      do j=1,6
        centroid(j) = sum(particles(:)%vec(j))/nparts
        particles_cent(:,j) = particles(:)%vec(j) - centroid(j)
        moments(j) = sqrt(sum(particles_cent(:,j)*particles_cent(:,j)) / nparts)
      enddo
    end subroutine

    subroutine make_sig_mat(particles,sigma_mat)
      implicit none

      real(rp) particles(:,:)
      real(rp) particles_cent(size(particles(:,1)),size(particles(1,:)))
      real(rp) sigma_mat(6,6)
      real(rp) moments(6)

      real(rp) centroid(6)
      integer j, k
      integer nparts

      nparts = size(particles(:,1))

      do j=1,6
        centroid(j) = sum(particles(:,j))/nparts
        particles_cent(:,j) = particles(:,j) - centroid(j)
      enddo
      
      do j=1,6
        do k=j,6
          sigma_mat(j,k) = sum(particles_cent(:,j)*particles_cent(:,k)) / nparts
          if( j .ne. k ) then
            sigma_mat(k,j) = sigma_mat(j,k)
          endif
        enddo
      enddo
    end subroutine
end program






