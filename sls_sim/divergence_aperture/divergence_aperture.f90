program divergence_aperture
  use mpi
  use bmad
  use bmad_parser_mod, only: bp_com
  use momentum_aperture_mod
  use touschek_mod
  use sls_lib

  use namelist_general ! lat_file, use_hybrid, use_line, periodicity
  use namelist_touschek ! tracking_method, use_line, rf_bucket, n_turn, dims, current, bunch_length, 
                        ! energy_spread, vertical_emittance, horizontal_emittance, aperture_file, stepping 

  implicit none

  type (lat_struct) ring
  type (coord_struct), allocatable :: co(:)
  type (normal_modes_struct) mode
  type (momentum_aperture_struct), allocatable :: max_delta_xp(:)

  integer i,j
  integer throw_away
  integer iargs
  integer iostat
  integer n_aperture_test
  integer ixer
  integer n_xp_locs

  real(rp) x_init, y_init
  real(rp) Tl
  real(rp) delta_s
  real(rp) high_water
  real(rp), parameter :: dhd = 0.01

  character*60 in_file
  character*100 lat_file_override
  character*200 line

  logical err
  integer status

  !mpi housekeeping
  integer myrank, from_id
  integer n_slave, cluster_size
  integer mpierr
  integer mpistatus(MPI_STATUS_SIZE)
  logical master
  integer n_sent, n_recv, mpi_i
  integer ix, incoming_ix

  !set default values
  use_line = ''
  periodicity = 1
  stepping = 'by_ix'
  aperture_file = 'generate_new'
  rf_bucket = 0.05
  n_xp_locs = 0
  energy_spread = -1
  bunch_length = -1
  vertical_emittance = -1
  horizontal_emittance = -1

  iargs = iargc()
  lat_file_override = ''
  if( iargs == 1 ) then
    call getarg(1,in_file)
  elseif (iargs == 2 ) then
    call getarg(1,in_file)
    call getarg(2,lat_file_override)
  endif

  call mpi_init(mpierr)                             ! Introduce yourself to the MPI daemon
  call mpi_comm_rank(MPI_COMM_WORLD,myrank,mpierr)  ! Get your rank number, store in myrank.  Master is rank 0.
  if(myrank .eq. 0) then
    master=.true.
  else
    master=.false.
  endif

  if(master) then
    !Check that cluster has at least two nodes
    call mpi_comm_size(MPI_COMM_WORLD,cluster_size,mpierr)
    n_slave=cluster_size-1
    if(n_slave .eq. 0) then
      write(*,*) "ERROR: no slaves found in cluster.  At least two nodes"
      write(*,*) "must be available to run this program."
      call mpi_finalize(mpierr)
      error stop
    endif
  endif
    
  open (unit = 10, file = in_file, action='read')
  read (10, nml = general)
  read (10, nml = touschek)
  close (10)

  n_turn = n_turn * periodicity

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, ring, use_line=use_line)
  bmad_com%aperture_limit_on = .true.
  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .false.

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
  if( (1.0*n_aperture_test)/ring%n_ele_track .lt. 0.5 ) then
    write(*,*) "Less than half the elements do not have a x1 physical aperture."
    write(*,*) "Probably something is wrong.  Check that lattice file defines aperture."
    write(*,*) "Aborting"
    stop
  endif

  if (tracking_method .gt. 0) then
    if(tracking_method .eq. 4) then
      write(*,*) "Tracking_method = 4 is linear tracking, which is a bad idea for a momentum aperture calculation."
      call mpi_finalize(mpierr)
      stop
    endif
    do i = 1, ring%n_ele_track
      if(ring%ele(i)%key .ne. wiggler$) then
        ring%ele(i)%tracking_method = tracking_method
      endif
    enddo
  endif

  !call calc_ring(ring,6,co,err)
  call twiss_and_track(ring,co,status)
  call radiation_integrals(ring, co, mode)
  if(energy_spread .gt. 0.0d0) mode%sigE_E = energy_spread
  if(bunch_length .gt. 0.0d0) mode%sig_z = bunch_length
  if(vertical_emittance .gt. 0.0d0) mode%b%emittance = vertical_emittance
  if(horizontal_emittance .gt. 0.0d0) mode%a%emittance = horizontal_emittance

  if(master) then
    write(*,*) "Emittances"
    write(*,*) "Horizontal Emittance: ", mode%a%emittance
    write(*,*) "Vertical Emittance:   ", mode%b%emittance
    write(*,*) "Energy Spread:        ", mode%sigE_E
    write(*,*) "Bunch Length:         ", mode%sig_z
  endif

  if(aperture_file .ne. 'generate_new') then
    open(21,file=aperture_file,action='read')
    n_xp_locs = 0
    do while(.true.)
      read(21,'(a)',iostat=iostat) line 
      if(iostat .ne. 0) exit
      line = trim(adjustl(line))
      if( line(1:1) .ne. '#' ) then
        n_xp_locs = n_xp_locs + 1
      endif
    enddo
    close(21)
    allocate(max_delta_xp(n_xp_locs))
  elseif( trim(stepping) == 'by_ix' ) then
    high_water = 0.0 !needed to deal with negative length elements
    n_xp_locs = 0
    do i=1,ring%n_ele_track
      if(ring%ele(i)%value(l$) .gt. 0.0001) then
        if(ring%ele(i)%s .gt. high_water+0.0001) then
          high_water = ring%ele(i)%s
          n_xp_locs = n_xp_locs + 1
        endif
      endif
    enddo
    allocate(max_delta_xp(n_xp_locs))
    ixer = 0
    high_water = 0.0
    do i=1,ring%n_ele_track
      if(ring%ele(i)%value(l$) .gt. 0.0001) then
        if(ring%ele(i)%s .gt. high_water+0.0001) then
          high_water = ring%ele(i)%s
          ixer = ixer + 1
          max_delta_xp(ixer)%s=ring%ele(i)%s
        endif
      endif
    enddo
    max_delta_xp(n_xp_locs)%s = ring%param%total_length - 0.0001
  elseif(trim(stepping) == 'by_n_steps') then
    if(n_xp_locs .eq. 0) then
      write(*,*) "If stepping is by_n_steps, then n_xp_locs must be greater than zero."
      stop
    endif
    allocate(max_delta_xp(n_xp_locs))
    delta_s = ring%param%total_length / (n_xp_locs-1)
    max_delta_xp(1)%s = 0.0d0
    do i=2, n_xp_locs-1
      max_delta_xp(i)%s = delta_s*(i-1)
    enddo
    max_delta_xp(n_xp_locs)%s = ring%param%total_length * 0.999999
  else
    write(*,*) "Unknown setting for stepping: ", trim(stepping)
    stop
  endif

  if( aperture_file == 'generate_new' ) then
    if(master) then
      n_sent = 0
      !initial seeding
      do mpi_i=1, min(n_slave,n_xp_locs)
        n_sent = n_sent + 1
        call mpi_send(n_sent, 1, MPI_INTEGER, mpi_i, 1, MPI_COMM_WORLD, mpierr)
      enddo
      !receive and reseed
      n_recv = 0
      do while(n_recv .lt. n_xp_locs)
        call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
        from_id = mpistatus(MPI_SOURCE)
        call mpi_recv(incoming_ix, 1, MPI_INTEGER, from_id, 2, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(max_delta_xp(incoming_ix)%pos, 1, MPI_DOUBLE_PRECISION, from_id, 3, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(max_delta_xp(incoming_ix)%neg, 1, MPI_DOUBLE_PRECISION, from_id, 4, MPI_COMM_WORLD, mpistatus, mpierr)
        write(*,'(a,i8,a,i8,a,2es14.5)') "Received ", incoming_ix, " of ", n_xp_locs, ". ",max_delta_xp(incoming_ix)%neg, max_delta_xp(incoming_ix)%pos
        n_recv = n_recv + 1
        if(n_sent .lt. n_xp_locs) then
          n_sent = n_sent + 1
          call mpi_send(n_sent, 1, MPI_INTEGER, from_id, 1, MPI_COMM_WORLD, mpierr)
        else
          call mpi_send(-1, 1, MPI_INTEGER, from_id, 1, MPI_COMM_WORLD, mpierr)
        endif
      enddo
    else
      do while(.true.)
        call mpi_recv(ix, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpistatus, mpierr)
        if(ix .lt. 0) then
          exit
        endif
        call divergence_aperture_one(ring,co,n_turn,6,max_delta_xp(ix),.true.)
        call mpi_send(ix, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, mpierr)
        call mpi_send(max_delta_xp(ix)%pos, 1, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, mpierr)
        call mpi_send(max_delta_xp(ix)%neg, 1, MPI_DOUBLE_PRECISION, 0, 4, MPI_COMM_WORLD, mpierr)
      enddo
    endif

    if(master) then
      !write aperture file
      open(21,file='div_aperture.dat')
      write(21,'(a8,a11,8a14)') "# ix", "s", "neg ap", "pos ap"
      do i=1,n_xp_locs
        write(21,'(i8,f11.4,2es14.5)') i, max_delta_xp(i)%s, max_delta_xp(i)%neg, max_delta_xp(i)%pos
      enddo
      close(21)
    endif
  else !read aperture in from file.
    if(master) then
      open(21,file=aperture_file,action='read')
      i = 0
      do while(.true.)
        read(21,'(a)') line
        line = trim(adjustl(line))
        if( line(1:1) .ne. '#' ) then
          i = i + 1
          read(line,*) throw_away, max_delta_xp(i)%s, max_delta_xp(i)%neg, max_delta_xp(i)%pos
          max_delta_xp(i)%s = min(max_delta_xp(i)%s,ring%ele(ring%n_ele_track)%s-0.0001)
          if( i .eq. n_xp_locs ) exit
        endif
      enddo
      close(21)
    endif
  endif

  ! if(master) then
  !   ring%param%n_part = current / e_charge * (periodicity*ring%param%total_length) / c_light
  !   call touschek_lifetime_with_aperture(mode, Tl, ring, max_delta_xp)

  !   open(21,file='tl-6D.dat')
  !   write(21,*) "# 6D Touschek Lifetime (hours)"
  !   write(21,*) Tl/3600
  !   close(21)
  !   write(*,*) "# 6D Touschek Lifetime (hours)"
  !   write(*,*) Tl/3600
  ! endif

  call mpi_finalize(mpierr)

end program
                                                            







