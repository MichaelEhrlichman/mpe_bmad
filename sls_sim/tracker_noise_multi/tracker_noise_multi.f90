program tracker_noise_multi
  use mpi
  use bmad_parser_mod, only: bp_com
  use bmad
  use sls_lib
  use naff_mod, except_rp=>rp
  use noise_mod
  use beam_mod
 
  implicit none

  integer, parameter :: N_MAX_NOISE = 20

  type(lat_struct) lat, lat_ideal
  type(coord_struct), allocatable :: co(:)
  type(coord_struct), allocatable :: orbit(:)
  type(coord_struct), allocatable :: tbt_coords(:)
  type(coord_struct), allocatable :: norm_coords(:)
  type(ele_noise_struct) ele_noise(N_MAX_NOISE)
  type(beam_init_struct) beam_init
  type(beam_struct) beam
  real(rp), allocatable :: particles(:,:)
  real(rp), allocatable :: particles_alive(:,:)

  integer i, j, m
  integer track_state, n_alive
  integer ix_a, ix_b
  integer n_seg
  integer status

  logical, allocatable :: particle_alive(:)

  ! NAFF, variables for spectral analysis
  real(rp), allocatable :: x_tunes(:,:)
  real(rp), allocatable :: y_tunes(:,:)
  complex(rp), allocatable :: x_amp(:,:)
  complex(rp), allocatable :: y_amp(:,:)
  complex(rp), allocatable :: cdata(:)
  integer, parameter :: max_naff = 10
  complex(rp) x_basis(max_naff,max_naff)
  complex(rp) y_basis(max_naff,max_naff)
  real(rp) x_norm(max_naff)
  real(rp) y_norm(max_naff)

  real(rp) vec(6)
  real(rp) harvest_r
  real(rp) sigma_mat(6,6), normal(3)

  logical ok, radiation

  character(100) lat_file
  character(100) in_file
  character(100) lat_file_override

  integer iargs

  !parameters from .in
  integer ix_inj
  real(rp) init_vec(6)
  integer harvest_i
  integer n_turns
  integer naff_turns
  integer dims
  integer phase_seed
  integer n_part, n_periods

  !mpi housekeeping
  integer myrank, from_id
  integer n_workers, cluster_size
  integer mpierr
  integer mpistatus(MPI_STATUS_SIZE), tag
  logical master
  character*9 mode
  integer n_sent, n_recv
  integer mpi_i, incoming_ix

  namelist / tracker_noise_params / lat_file, ix_inj, init_vec, n_turns, naff_turns, &
                                    dims, radiation, ele_noise, phase_seed, n_part, &
                                    n_periods

  iargs = iargc()
  lat_file_override = ''
  if( iargs == 1 ) then
    call getarg(1,in_file)
  elseif (iargs == 2 ) then
    call getarg(1,in_file)
    call getarg(2,lat_file_override)
  endif

  n_part = 1
  n_periods = 1
  dims = 4
  radiation = .false.
  phase_seed = 12345

  open (unit = 10, file = in_file)
  read (10, nml = tracker_noise_params)
  close (10)

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  call mpi_init(mpierr)                             ! Introduce yourself to the MPI daemon
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,mpierr)  ! Get your rank number, store in myrank.  Master is rank 0.
  if(myrank .eq. 0) then
    master=.true.
    call mpi_comm_size(MPI_COMM_WORLD,cluster_size,mpierr)
    n_workers=cluster_size-1
    write(*,*) "Number of workers: ", n_workers
    if(n_workers .gt. 0) then
      mode = 'multicore'
    else
      mode = 'singlcore'
    endif
  else
    master=.false.
    mode = ''
  endif

  if(master) write(*,*) "Preparing lattice..."

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, lat)
  
  if(dims==4) then
    CALL set_on_off(rfcavity$, lat, off$)
  elseif (dims==6) then
    CALL set_on_off(rfcavity$, lat, on$)
  endif
  bmad_com%radiation_damping_on = radiation
  bmad_com%radiation_fluctuations_on = radiation
  
  call twiss_and_track(lat,co,status)

  allocate(tbt_coords(n_turns))
  allocate(norm_coords(n_turns))
  allocate(cdata(n_turns))
  allocate(orbit(0:lat%n_ele_track))

  if(master) write(*,*) "Tracking..."

  lat_ideal = lat

  call ran_seed_put(phase_seed) !random phase of element properties should be same for all workers
  call init_ele_noise(lat, ele_noise)
  call ran_seed_put(0) !radiation integrals feeds off same random number generator, this resets
                       !generator based on system clock so each particle see different radiation kicks
  call ran_uniform(harvest_r)
  harvest_i = harvest_r*10000
  call ran_seed_put(harvest_i + myrank) !avoid possibility that workers could have same random harvest
                                        !due to polling same clock

  if(master) then
    open(50,file='distribution.dat')
    write(50,'(6A14)') "# x", "px", "y", "py", "z", "pz"
    open(51,file='moments.dat')
    write(51,'(A9, 3A14)') "# ix", "sqrt(xx)", "sqrt(yy)", "sqrt(zz)"
    open(52,file='sig_emit.dat')
    write(52,'(A9, 3A14)') "# ix", "e1", "e2", "e3"
  endif

  if(master) then
    beam_init%distribution_type(1:3) = 'ran_gauss'
    beam_init%n_particle = n_part
    beam_init%a_emit = 1.806E-09
    beam_init%b_emit = 1.0e-12
    beam_init%sig_e = 8.527E-04
    beam_init%sig_z = 0.00571
    beam_init%species = 'electron'
    beam_init%center = co(ix_inj)%vec

    call init_beam_distribution(lat%ele(ix_inj), lat%param, beam_init, beam)

    allocate(particles(n_part,6))
    allocate(particles_alive(n_part,6))
    do i=1,n_part
      particles(i,:) = beam%bunch(1)%particle(i)%vec(:)
    enddo
    i=0
    call write_report()
  endif

  !FOO test
  !n_part = 6
  !allocate(particles(6,6))
  !particles(1,:) = (/0d0,0d0,0d0,0d0,0d0,0d0/)
  !particles(2,:) = (/1d-5,0d0,0d0,0d0,0d0,0d0/)
  !particles(3,:) = (/-1d-5,0d0,0d0,0d0,0d0,0d0/)
  !particles(4,:) = (/0d0,1d-5,0d0,0d0,0d0,0d0/)
  !particles(5,:) = (/0d0,-1d-5,0d0,0d0,0d0,0d0/)
  !particles(6,:) = (/1d-5,-1d-5,0d0,0d0,0d0,0d0/)

  !MPI Tags
  ! 1: vec(6) from master to slave
  ! 2: state from slave to master
  ! 3: vec(6) from slave to master
  ! 4: signals slave to terminate
  if(master) then
    if(mode == 'singlcore') then
      write(*,*) "Non-mpi environment or cluster size is 1.  Running in single core mode."
      do j=1, n_periods
        n_alive = 0
        do i=1,n_part
          call track_particle(particles(i,:),track_state)
          if(track_state == moving_forward$) then
            n_alive = n_alive + 1 
            particles_alive(n_alive,:) = particles(i,:)
          endif
          write(*,'(a,a,i9,a,i9)') char(10), "particle ", i, " of ", n_part
        enddo
        particles(:,:) = particles_alive(:,:)
        n_part = n_alive

        !end of period
        !write coordinates of each particle in distribution
        call write_report()
      enddo
    else
      write(*,*) "MPI environment detected.  Running with ", n_workers, " workers."
      do i=1, n_periods
        write(*,*) "Tracking period ", i
        n_sent = 0
        !initial seeding
        do mpi_i=1, min(n_workers,n_part)
          n_sent = n_sent + 1
          call mpi_send(particles(mpi_i,:), 6, MPI_DOUBLE_PRECISION, mpi_i, 1, MPI_COMM_WORLD, mpierr)
        enddo
        !receive and reseed
        n_recv = 0
        n_alive = 0
        do while(n_recv .lt. n_part)
          call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
          from_id = mpistatus(MPI_SOURCE)
          n_recv = n_recv + 1
          call mpi_recv(track_state, 1, MPI_INTEGER, from_id, 2, MPI_COMM_WORLD, mpistatus, mpierr)
          if(track_state == moving_forward$) then
            n_alive = n_alive + 1
            call mpi_recv(particles(n_alive,:), 6, MPI_DOUBLE_PRECISION, from_id, 3, MPI_COMM_WORLD, mpistatus, mpierr)
          endif
          if(n_sent .lt. n_part) then
            n_sent = n_sent + 1
            call mpi_send(particles(n_sent,:), 6, MPI_DOUBLE_PRECISION, from_id, 1, MPI_COMM_WORLD, mpierr)
          endif
          write(*,'(a,i6,a,i6,a)') "Multi mode: point ", n_recv, " of ", n_part, " complete."
        enddo
        n_part = n_alive

        !end of period
        !write coordinates of each particle in distribution
        call write_report()
      enddo
      do i=1,n_workers
        call mpi_send(-1, 1, MPI_INTEGER, i, 4, MPI_COMM_WORLD, mpierr)
      enddo
    endif
    close(50)
    close(51)
    close(52)
  else
    do while(.true.)
      call mpi_probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
      tag = mpistatus(MPI_TAG)
      if(tag .eq. 4) exit
      call mpi_recv(vec, 6, MPI_DOUBLE_PRECISION, 0, 1, MPI_COMM_WORLD, mpistatus, mpierr)
      call track_particle(vec,track_state)
      call mpi_send(track_state, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, mpierr) 
      if( track_state .eq. moving_forward$) then
        call mpi_send(vec, 6, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, mpierr)
      endif
    enddo
    write(*,*) "Worker ", myrank, " peacefully exited loop."
  endif

  !if(master) then
  !  if( track_state == moving_forward$ ) then
  !    do i = 1, n_turns
  !      call xy_to_action(lat, ix_inj, tbt_coords(i)%vec, norm_coords(i)%vec, ok) !4D + dispersion
  !      if (.not. ok) then
  !        write(*,*) "Error from xy_to_action"
  !        stop
  !      endif
  !    enddo

  !    open(100,file='tbt_Anorm.dat')
  !    write(100,'(A6,4A14)') "# turn", "Jx", "Jx'", "Jy", "Jy'"
  !    write(100,'(a,es14.6)') "# Jx is ", (norm_coords(1)%vec(1)**2 + norm_coords(1)%vec(2)**2)/2.0d0
  !    write(100,'(a,es14.6)') "# Jy is ", (norm_coords(1)%vec(3)**2 + norm_coords(1)%vec(4)**2)/2.0d0
  !    do i=1,n_turns
  !      write(100,'(I6,4ES14.5)') i, norm_coords(i)%vec(1:4)
  !    enddo
  !    close(100)

  !    !subtract centroid
  !    norm_coords(:)%vec(1) = norm_coords(:)%vec(1) - sum(norm_coords(:)%vec(1))/size(norm_coords(:)%vec(1))
  !    norm_coords(:)%vec(2) = norm_coords(:)%vec(2) - sum(norm_coords(:)%vec(2))/size(norm_coords(:)%vec(2))
  !    norm_coords(:)%vec(3) = norm_coords(:)%vec(3) - sum(norm_coords(:)%vec(3))/size(norm_coords(:)%vec(3))
  !    norm_coords(:)%vec(4) = norm_coords(:)%vec(4) - sum(norm_coords(:)%vec(4))/size(norm_coords(:)%vec(4))

  !    write(*,*) "Spectral Analysis ..."

  !    n_seg = n_turns/naff_turns
  !    allocate(x_tunes(n_seg,max_naff))
  !    allocate(y_tunes(n_seg,max_naff))
  !    allocate(x_amp(n_seg,max_naff))
  !    allocate(y_amp(n_seg,max_naff))
  !    x_tunes = 0.0d0
  !    y_tunes = 0.0d0
  !    do i=1, n_seg
  !      ix_a = 1 + (i-1)*naff_turns
  !      ix_b = i*naff_turns
  !      x_basis = 0.0d0
  !      x_norm = 0.0d0
  !      !call naff( (norm_coords(ix_a:ix_b)%vec(1), norm_coords(ix_a:ix_b)%vec(2)), x_tunes(i,:), x_amp(i,:),opt_dump_spectra=40+i)
  !      cdata = cmplx(norm_coords(:)%vec(1),norm_coords(:)%vec(2))
  !      call naff( cdata(ix_a:ix_b), x_tunes(i,:), x_amp(i,:))
  !      !call naff( (norm_coords(ix_a:ix_b)%vec(3), norm_coords(ix_a:ix_b)%vec(4)), y_tunes(i,:), y_amp(i,:),opt_dump_spectra=50+i)
  !      cdata = cmplx(norm_coords(:)%vec(3),norm_coords(:)%vec(4))
  !      call naff( cdata(ix_a:ix_b), y_tunes(i,:), y_amp(i,:))
  !    enddo

  !    open(100,file='tunes_x_slicingB.dat')
  !    open(101,file='tunes_y_slicingB.dat')
  !    write(100,'(4a18)') "# segment", "tune", "amp", "angle (rad)"
  !    write(101,'(4a18)') "# segment", "tune", "amp", "angle (rad)"
  !    do i=1,n_seg
  !      do j=1,max_naff
  !        write(100,'(i18,3es18.8)') i, x_tunes(i,j), abs(x_amp(i,j)), atan2(aimag(x_amp(i,j)),real(x_amp(i,j)))
  !      enddo
  !      do j=1,max_naff
  !        write(101,'(i18,3es18.8)') i, y_tunes(i,j), abs(y_amp(i,j)), atan2(aimag(y_amp(i,j)),real(y_amp(i,j)))
  !      enddo
  !      write(100,*)
  !      write(100,*)
  !      write(101,*)
  !      write(101,*)
  !    enddo
  !    close(100)
  !    close(101)

  !    deallocate(x_tunes)
  !    deallocate(y_tunes)
  !    deallocate(x_amp)
  !    deallocate(y_amp)
  !  else
  !    write(*,*) "Particle was lost.  Not calculating tunes."
  !  endif
  !endif

  call mpi_finalize(mpierr)

  contains
    subroutine track_particle(vec,state)
      implicit none
      real(rp) vec(6)
      integer state
      integer i

      !orbit(ix_inj)%vec = co(ix_inj)%vec + init_vec

      do i = 1, n_turns
        call apply_ele_noise(lat, lat_ideal, ele_noise, i)

        if(mode .eq. 'singlcore') then
          if(mod(i,100).eq.0) write(*,'(a,a,i9,a,i9)',advance='no') achar(13), "turn ", i, " of ", n_turns
        endif
        call init_coord(orbit(ix_inj), vec, ele=lat%ele(ix_inj), element_end=upstream_end$, particle=electron$)
        call track_many(lat, orbit, ix_inj, ix_inj, 1, track_state=state)
        if(state /= moving_forward$) then
          exit
        endif
      enddo
      vec = orbit(ix_inj)%vec
    end subroutine

    subroutine make_sig_mat(particles,sigma_mat)
      implicit none

      real(rp) particles(:,:)
      real(rp) particles_cent(size(particles(:,1)),size(particles(1,:)))
      real(rp) sigma_mat(6,6)

      integer j, k
      integer nparts

      Nparts = size(particles(:,1))

      do j=1,6
        particles_cent(:,j) = particles(:,j) - sum(particles(:,j))/Nparts
      enddo
      
      do j=1,6
        do k=j,6
          sigma_mat(j,k) = sum(particles_cent(:,j)*particles_cent(:,k)) / Nparts
          if( j .ne. k ) then
            sigma_mat(k,j) = sigma_mat(j,k)
          endif
        enddo
      enddo
    end subroutine

    subroutine normal_sig_mat(sigma_mat,normal)
      use bmad
      use eigen_mod

      implicit none

      real(rp) sigma_mat(1:6,1:6)
      real(rp) normal(1:3)
      real(rp) sigmaS(1:6,1:6)
      real(rp) S(1:6,1:6)
      real(rp), parameter :: S2(1:2,1:2) = reshape([0,-1,1,0], [2,2])

      complex(rp) eigen_val(6), eigen_vec(6,6)

      logical error

      S = 0.0d0
      S(1:2,1:2) = S2
      S(3:4,3:4) = S2
      S(5:6,5:6) = S2

      sigmaS = matmul(sigma_mat,S)

      call mat_eigen(sigmaS, eigen_val, eigen_vec, error)

      normal(1) = abs(aimag(eigen_val(1)))
      normal(2) = abs(aimag(eigen_val(3)))
      normal(3) = abs(aimag(eigen_val(5)))
    end subroutine

    subroutine write_report()
      do j=1, n_part
        write(50,'(6es14.5)') particles(j,:)
      enddo
      write(50,*)
      write(50,*)

      call make_sig_mat(particles,sigma_mat)
      write(51,'(i9,3es14.5)') i*n_turns, sqrt(sigma_mat(1,1)), sqrt(sigma_mat(3,3)), sqrt(sigma_mat(5,5))
      call normal_sig_mat(sigma_mat,normal)
      write(52,'(i9,3es14.5)') i*n_turns, normal(1:3)

      call flush(50)
      call flush(51)
      call flush(52)
    end subroutine
end program






