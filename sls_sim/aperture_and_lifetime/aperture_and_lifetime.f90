program aperture_and_lifetime
  use mpi
  use bmad
  use bmad_parser_mod, only: bp_com
  use momentum_aperture_mod
  use touschek_mod

  implicit none

  type (lat_struct) ring
  type (coord_struct), allocatable :: co(:)
  type (normal_modes_struct) mode
  type (momentum_aperture_struct), allocatable :: max_deltam(:)
  type (ele_pointer_struct), allocatable :: eles(:)

  character(5) mpi_mode
  integer i,j
  integer throw_away
  integer iargs
  integer iostat
  integer n_aperture_test
  integer ixer
  integer n_hd

  real(rp) x_init, y_init
  real(rp) Tl
  real(rp) delta_s
  real(rp) high_water
  real(rp) k1, u0, freq

  character*60 in_file
  character*100 lat_file_override
  character*200 line

  logical err
  integer status

  integer given_ma_locs(20)
  integer given_ma_range(2)

  !mpi housekeeping
  integer myrank, from_id
  integer n_slave, cluster_size
  integer mpierr
  integer mpistatus(MPI_STATUS_SIZE)
  logical master
  integer n_sent, n_recv, mpi_i
  integer ix, incoming_ix

  !namelist touschek
  integer tracking_method
  real(rp) rf_bucket
  integer n_turn
  real(rp) current
  real(rp) bunch_length
  real(rp) energy_spread
  real(rp) vertical_emittance
  real(rp) horizontal_emittance
  character*100 aperture_file
  character*10 stepping
  character*20 use_line
  integer n_ma_locs
  integer periodicity
  character(100) lat_file
  character(100) mask
  logical use_hybrid
  real(rp) MV_BL_rfvoltage
  type(ele_pointer_struct), allocatable :: rfcav_eles(:)
  integer n_cavs
  real(rp) voltage

  namelist / general /    lat_file, &
                          use_hybrid

  namelist / touschek /   tracking_method, &
                          use_line, &
                          periodicity, &
                          rf_bucket, &
                          n_turn, &          !Number of turns particle must survive
                          current, &
                          energy_spread, &
                          bunch_length, &
                          MV_BL_rfvoltage, &
                          vertical_emittance, &
                          horizontal_emittance, &
                          aperture_file, &
                          stepping, &
                          n_ma_locs, &
                          given_ma_locs, &
                          given_ma_range


  !set default values
  given_ma_locs = -1
  given_ma_range = -1
  use_line = ''
  periodicity = 1
  stepping = 'by_ix'
  aperture_file = 'generate_new'
  rf_bucket = 0.05
  n_ma_locs = 0
  energy_spread = -1
  bunch_length = -1
  vertical_emittance = -1
  horizontal_emittance = -1
  MV_BL_rfvoltage = -1

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
    call mpi_comm_size(MPI_COMM_WORLD,cluster_size,mpierr)
    n_slave=cluster_size-1
    write(*,*) "Number of slaves: ", n_slave
    if(n_slave .gt. 0) then
      mpi_mode = 'multi'
    else
      mpi_mode = 'singl'
    endif
  else
    master=.false.
    mpi_mode = ''
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

  call twiss_and_track(ring,co,status)
  call radiation_integrals(ring, co, mode)
  if(energy_spread .gt. 0.0d0) mode%sigE_E = energy_spread
  if(bunch_length .gt. 0.0d0) mode%sig_z = bunch_length
  if(vertical_emittance .gt. 0.0d0) mode%b%emittance = vertical_emittance
  if(horizontal_emittance .gt. 0.0d0) mode%a%emittance = horizontal_emittance

  if(MV_BL_rfvoltage .gt. 0) then
    write(*,*) "Setting bunch length according to MV calculation."

    call lat_ele_locator("rfcav::*",ring,rfcav_eles,n_cavs)
    freq = rfcav_eles(1)%ele%value(rf_frequency$)
    !voltage = rfcav_eles(1)%ele%value(voltage$)
    voltage = MV_BL_rfvoltage
    deallocate(rfcav_eles)

    k1 = twopi*freq/2.99792458e8
    u0 = mode%e_loss
    mode%sig_z = 0.765*(mode%synch_int(1)*ring%ele(0)%value(E_TOT$)*mode%sigE_E**2 / &
                 (k1**3 * voltage * sqrt(1-81.0d0/64.0d0*(u0/voltage)**2)))**0.25
    write(*,*) "MV Bunch Length: ", mode%sig_z
  endif

  if(master) then
    write(*,*) "Emittances"
    write(*,*) "Horizontal Emittance: ", mode%a%emittance
    write(*,*) "Vertical Emittance:   ", mode%b%emittance
    write(*,*) "Energy Spread:        ", mode%sigE_E
    write(*,*) "Bunch Length:         ", mode%sig_z
  endif

  if(aperture_file .ne. 'generate_new') then
    open(21,file=aperture_file,action='read')
    n_ma_locs = 0
    do while(.true.)
      read(21,'(a)',iostat=iostat) line 
      if(iostat .ne. 0) exit
      line = trim(adjustl(line))
      if( line(1:1) .ne. '#' ) then
        n_ma_locs = n_ma_locs + 1
      endif
    enddo
    close(21)
    allocate(max_deltam(n_ma_locs))
  elseif( trim(stepping) == 'given_ix' ) then
    n_ma_locs = count(given_ma_locs .ge. 0)
    allocate(max_deltam(n_ma_locs))
    do i=1,n_ma_locs
      max_deltam(i)%s=ring%ele(given_ma_locs(i))%s
    enddo
  elseif( trim(stepping) == 'mask' ) then
    call lat_ele_locator (mask, ring, eles, n_ma_locs, err)
    allocate(max_deltam(n_ma_locs))
    do i=1,n_ma_locs
      max_deltam(i)%s=eles(i)%ele%s
    enddo
  elseif( trim(stepping) == 'range_ix' ) then
    n_ma_locs = given_ma_range(2)-given_ma_range(1)+1
    allocate(max_deltam(n_ma_locs))
    max_deltam(1:n_ma_locs)%s=ring%ele(given_ma_range(1):given_ma_range(2))%s
    max_deltam(1:n_ma_locs)%s=max(max_deltam(1:n_ma_locs)%s,0.0001)
  elseif( trim(stepping) == 'by_ix' ) then
    high_water = 0.0 !needed to deal with negative length elements
    n_ma_locs = 0
    do i=1,ring%n_ele_track
      if(ring%ele(i)%value(l$) .gt. 0.0001) then
        if(ring%ele(i)%s .gt. high_water+0.0001) then
          high_water = ring%ele(i)%s
          n_ma_locs = n_ma_locs + 1
        endif
      endif
    enddo
    allocate(max_deltam(n_ma_locs))
    ixer = 0
    high_water = 0.0
    do i=1,ring%n_ele_track
      if(ring%ele(i)%value(l$) .gt. 0.0001) then
        if(ring%ele(i)%s .gt. high_water+0.0001) then
          high_water = ring%ele(i)%s
          ixer = ixer + 1
          max_deltam(ixer)%s=ring%ele(i)%s
        endif
      endif
    enddo
    max_deltam(n_ma_locs)%s = min(max_deltam(n_ma_locs)%s, ring%param%total_length - 0.0001)
  elseif(trim(stepping) == 'by_n_steps') then
    if(n_ma_locs .eq. 0) then
      write(*,*) "If stepping is by_n_steps, then n_ma_locs must be greater than zero."
      stop
    endif
    allocate(max_deltam(n_ma_locs))
    delta_s = ring%param%total_length / (n_ma_locs-1)
    max_deltam(1)%s = 0.0d0
    do i=2, n_ma_locs-1
      max_deltam(i)%s = delta_s*(i-1)
    enddo
    max_deltam(n_ma_locs)%s = ring%param%total_length * 0.999999
  else
    write(*,*) "Unknown setting for stepping: ", trim(stepping)
    stop
  endif

  if( aperture_file == 'generate_new' ) then
    if(master) then
      if(mpi_mode == 'singl') then
        do ix = 1, n_ma_locs
          call momentum_aperture_one(ring,co,n_turn,6,max_deltam(ix),.true.)
          write(*,'(a,i6,a,i6,a)') "Single mode: point ", ix, " of ", n_ma_locs, " complete."
        enddo
      elseif(mpi_mode == 'multi') then
        n_sent = 0
        !initial seeding
        do mpi_i=1, min(n_slave,n_ma_locs)
          n_sent = n_sent + 1
          call mpi_send(n_sent, 1, MPI_INTEGER, mpi_i, 1, MPI_COMM_WORLD, mpierr)
        enddo
        !receive and reseed
        n_recv = 0
        do while(n_recv .lt. n_ma_locs)
          call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
          from_id = mpistatus(MPI_SOURCE)
          call mpi_recv(incoming_ix, 1, MPI_INTEGER, from_id, 2, MPI_COMM_WORLD, mpistatus, mpierr)
          write(*,'(a,i8,a,i8)') "Received ", incoming_ix, " of ", n_ma_locs
          call mpi_recv(max_deltam(incoming_ix)%pos, 1, MPI_DOUBLE_PRECISION, from_id, 3, MPI_COMM_WORLD, mpistatus, mpierr)
          call mpi_recv(max_deltam(incoming_ix)%neg, 1, MPI_DOUBLE_PRECISION, from_id, 4, MPI_COMM_WORLD, mpistatus, mpierr)
          n_recv = n_recv + 1
          if(n_sent .lt. n_ma_locs) then
            n_sent = n_sent + 1
            call mpi_send(n_sent, 1, MPI_INTEGER, from_id, 1, MPI_COMM_WORLD, mpierr)
          else
            call mpi_send(-1, 1, MPI_INTEGER, from_id, 1, MPI_COMM_WORLD, mpierr)
          endif
        enddo
      endif
    else
      do while(.true.)
        call mpi_recv(ix, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpistatus, mpierr)
        if(ix .lt. 0) then
          exit
        endif
        call momentum_aperture_one(ring,co,n_turn,6,max_deltam(ix),.true.)
        call mpi_send(ix, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, mpierr)
        call mpi_send(max_deltam(ix)%pos, 1, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, mpierr)
        call mpi_send(max_deltam(ix)%neg, 1, MPI_DOUBLE_PRECISION, 0, 4, MPI_COMM_WORLD, mpierr)
      enddo
    endif

    if(master) then
      !write aperture file
      open(21,file='aperture.dat')
      write(21,'(a8,a11,8a14)') "# ix", "s", "neg ap", "pos ap"
      do i=1,n_ma_locs
        write(21,'(i8,f11.4,2es14.5)') i, max_deltam(i)%s, max_deltam(i)%neg, max_deltam(i)%pos
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
          read(line,*) throw_away, max_deltam(i)%s, max_deltam(i)%neg, max_deltam(i)%pos
          max_deltam(i)%s = min(max_deltam(i)%s,ring%ele(ring%n_ele_track)%s-0.0001)
          if( i .eq. n_ma_locs ) exit
        endif
      enddo
      close(21)
    endif
  endif

  if(master) then
    do i=1,n_ma_locs
      max_deltam(i)%neg = max(max_deltam(i)%neg,-1*rf_bucket)
      max_deltam(i)%pos = min(max_deltam(i)%pos,   rf_bucket)
    enddo
    ring%param%n_part = current / e_charge * (periodicity*ring%param%total_length) / c_light
    call touschek_lifetime_with_aperture(mode, Tl, ring, max_deltam)

    open(21,file='tl-6D.dat')
    write(21,*) "# 6D Touschek Lifetime (hours)"
    write(21,*) Tl/3600
    close(21)
    write(*,*) "# 6D Touschek Lifetime (hours)"
    write(*,*) Tl/3600
  endif

  call mpi_finalize(mpierr)

end program
                                                            







