program tracker_train
  use bmad
  use bmad_parser_mod, only: bp_com
  use beam_mod
 
  implicit none

  character(100) lat_file
  character(100) in_file
  character(180) line
  character(20) filename

  type(lat_struct) lat, hybrid
  type(coord_struct), allocatable :: co(:), coh(:)
  type(coord_struct), allocatable :: orbit(:)
  type(beam_init_struct) beam_init
  type(beam_struct) beam
  type(ele_pointer_struct), allocatable :: rf_eles(:)
  type(ele_struct), pointer :: wake_ele => null()

  integer i, j, k
  integer status
  integer n_rf

  logical ok, radiation, err

  real(rp) prev

  !parameters from .in
  integer ix_inj
  integer n_turns
  integer n_part, n_bunches
  real(rp) init_vec(6)
  real(rp) bunch_charge
  integer MBF_pickup_ix, MBF_kicker_ix
  

  namelist / parameters / lat_file, ix_inj, init_vec, n_turns, &
                             radiation, n_part, n_bunches, bunch_charge, &
                             MBF_pickup_ix, MBF_kicker_ix

  call getarg(1,in_file)

  lat_file = ''
  ix_inj = 0
  init_vec = 0.0d0
  n_turns = 10
  radiation = .false.
  n_part = 10
  n_bunches = 1
  bunch_charge = 1.0e-9

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
  beam_init%center = co(ix_inj)%vec + init_vec
  beam_init%full_6D_coupling_calc = .true.

  call lat_ele_locator('rfcav', lat, rf_eles, n_rf)
  beam_init%dt_bunch = 1.0d0 / value_of_attribute(rf_eles(1)%ele,'RF_FREQUENCY')
  write(*,'(a,es14.6,a)') "Bunch spacing: ", beam_init%dt_bunch, ' seconds'

  call init_beam_distribution(lat%ele(ix_inj), lat%param, beam_init, beam)

  open(45,file='~bunch_vs_number.out')
  write(45,'(a12,a6,12a14)') "# turn", "bunch", "x0", "x0'", "y0", "y0'", "z0", "z0'", "r<xx>", "r<x'x'>", "r<yy>", "r<y'y'>", "r<zz>", "r<z'z'>"

  do i=1,n_bunches
    write(filename,'(a,i0.4,a)') '~bunch_',i,'.out'
    open(1000+i,file=filename)
  enddo

  do i=1,lat%n_ele_track
    if(associated(lat%ele(i)%wake)) then
      wake_ele => lat%ele(i)
      write(*,'(a,a,a)') "Element ", lat%ele(i)%name, " has wake."
      lat%ele(i)%select = .true.
    elseif (any(lat%ele(i)%key==(/sbend$,rbend$/)) .and. radiation) then
      lat%ele(i)%select = .true.
    else
      lat%ele(i)%select = .false.
    endif
  enddo
  call make_hybrid_lat(lat, hybrid, .true.)
  !call twiss_and_track(hybrid,coh)
  !write(*,*) "hybrid tunes: ", hybrid%ele(hybrid%n_ele_track)%a%phi/twopi, hybrid%ele(hybrid%n_ele_track)%b%phi/twopi
  do i=1,hybrid%n_ele_track
    if(associated(hybrid%ele(i)%wake)) then
      wake_ele => hybrid%ele(i)
    endif
  enddo

  if(associated(wake_ele)) then
    do i=1,size(wake_ele%wake%lr_mode)
      write(filename,'(a,i0.4,a)') '~wake_',i,'.out'
      open(2000+i,file=filename)
    enddo
  endif

  call report(beam,0)
  if(associated(wake_ele)) call report_wake(wake_ele,0)

  do j=1,n_turns
    if(mod(j,100)==0) write(*,'(a,i8)') "Tracking turn ", j
    call track_beam(hybrid,beam,err=err)
    if(any(beam%bunch(:)%n_live .ne. n_part)) then
      write(*,'(a)') "bunch lost"
      exit
    endif
    call report(beam,j)
    if(associated(wake_ele)) call report_wake(wake_ele,j)
    !MBF_kick = FIR(, x, a
  enddo

  close(45)
  do i=1,n_bunches
    close(1000+i)
  enddo

  contains
    function FIR(xnew,x,a) result(y)
      real(rp) xnew, x(20), a(20)
      real(rp) y
      
      x(2:20) = x(1:19)
      x(1) = xnew

      y = dot_product(x,a)
    end function

    subroutine report_wake(ele,turn)
      type(ele_struct) ele
      integer i, turn

      do i=1, size(ele%wake%lr_mode)
        write(2000+i,'(i8,4es14.5)') turn, ele%wake%lr_mode(i)%a_sin, ele%wake%lr_mode(i)%a_cos, ele%wake%lr_mode(i)%b_sin, ele%wake%lr_mode(i)%b_cos
      enddo
    end subroutine
      
    subroutine report(beam,turn)
      type(beam_struct) beam
      integer i, turn
      real(rp) centroid(6), moments(6)

      do i=1,size(beam%bunch(:))
        call make_moments(beam%bunch(i)%particle(:), centroid, moments)
        write(45,'(i12,i6,12es14.5)') turn, i, centroid(1:6)
        write(1000+i,'(i12,i6,12es14.5)') turn, i, centroid(1:6)
        !write(45,'(i12,i6,12es14.5)') turn, i, centroid(1:6), moments(1:6)
        !write(1000+i,'(i12,i6,12es14.5)') turn, i, centroid(1:6), moments(1:6)
      enddo
      write(45,*)
      write(45,*)
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






