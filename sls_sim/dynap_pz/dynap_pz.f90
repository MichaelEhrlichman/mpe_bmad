program dynap_pz
  use mpi
  use bmad
  use bmad_parser_mod, only: bp_com
  use momentum_aperture_mod
  use sls_lib
  use naff_mod
  use linear_aperture_mod
  use custom_dynamic_aperture_mod

  implicit none

  type raster_points_struct
    integer pz_k
    real(rp) x,pz
    real(rp) alive
    real(rp) nu_x, nu_y
    real(rp) Dfma
  end type

  character(100) lat_file
  character(1) ix_str
  character(60) in_file
  character(100) lat_file_override

  logical use_hybrid
  logical radiation
  logical stable

  real(rp) x_min, x_max
  real(rp) pz_min, pz_max

  type (lat_struct) ring
  type (lat_struct) ring0
  type (coord_struct), allocatable :: co(:)
  type (coord_struct), allocatable :: orb(:)
  type (coord_struct), allocatable :: naff_turns(:)
  type (raster_points_struct), allocatable :: raster_points(:)
  type (raster_points_struct), allocatable :: simple_grid(:,:)
  type (custom_aperture_scan_struct) da_config

  complex(rp), allocatable :: cdata(:)

  integer i,j,k,m
  integer ix, nix
  integer last_k
  integer track_state
  integer iargs
  integer n_aperture_test
  integer toggle
  integer iostat
  integer n_points
  integer status
  integer nx, n_naff, n_turns
  integer npz

  !mpi housekeeping
  integer myrank, from_id
  integer n_slave, cluster_size
  integer mpierr
  integer mpistatus(MPI_STATUS_SIZE)
  logical master
  character(5) mode
  integer n_sent, n_recv
  integer mpi_i, incoming_ix

  logical err

  real(rp) nu_x(2), nu_y(2)
  real(rp) freqs(2)
  complex(rp) amps(2)

  namelist / general /    lat_file, &
                          use_hybrid

  namelist / xpz_raster /  radiation, &
                           x_min, &
                           x_max, &
                           nx, &
                           n_turns, &
                           n_naff, &
                           pz_min, &
                           pz_max, &
                           npz

  n_turns = 100
  use_hybrid = .false.
  n_naff = 64

  call getarg(1,in_file)
  lat_file_override = ''
  iargs = iargc()
  if (iargs == 2 ) then
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
      mode = 'multi'
    else
      mode = 'singl'
    endif
  else
    master=.false.
    mode = ''
  endif

  if(master) then
    open(unit=21, iostat=iostat, file='raster.key', status='old')
    if (iostat .eq. 0) then
      close(21, status='delete')
    endif
    close(21)
  endif

  radiation = .false.
    
  open (unit = 10, file = in_file, action='read')
  read (10, nml = general)
  read (10, nml = xpz_raster)
  close (10)

  if( (mod(nx,2)==0) .or. (mod(npz,2)==0) ) then
    write(*,*) "Please use only odd nx and npz."
    stop
  endif

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, ring)

  if(radiation) then
    bmad_com%radiation_damping_on = .true.
    bmad_com%radiation_fluctuations_on = .true.
    call set_on_off(rfcavity$, ring, on$)
  else
    bmad_com%radiation_damping_on = .false.
    bmad_com%radiation_fluctuations_on = .false.
    call set_on_off(rfcavity$, ring, off$)
  endif

  allocate(co(0:ring%n_ele_track))
  allocate(orb(0:ring%n_ele_track))
  if(n_naff .gt. 0) allocate(naff_turns(n_naff))
  if(n_naff .gt. 0) allocate(cdata(n_naff))

  bmad_com%aperture_limit_on = .true.

  do i=1,ring%n_ele_track
    if(ring%ele(i)%key == wiggler$) then
      ring%ele(i)%value(x1_limit$) = 1.0
      ring%ele(i)%value(x2_limit$) = 1.0
      ring%ele(i)%value(y1_limit$) = 1.0
      ring%ele(i)%value(y2_limit$) = 1.0
      ring%ele(i)%aperture_type = elliptical$
    endif
  enddo

!  if(master) then
!    n_aperture_test = 0
!    do i=1,ring%n_ele_track
!      if( ring%ele(i)%value(x1_limit$) .gt. 1e-4 ) n_aperture_test = n_aperture_test + 1
!    enddo
!    if( (1.0*n_aperture_test)/ring%n_ele_track .lt. 0.5 ) then
!      write(*,*) "Less than half the elements do not have a x1 physical aperture."
!      write(*,*) "Probably something is wrong.  Check that lattice file defines aperture."
!      write(*,*) "Aborting"
!      stop
!    endif
!  endif

!  if(tracking_method .gt. 0) then
!    do i = 1, ring%n_ele_max
!      if(ring%ele(i)%key .ne. wiggler$) then
!        ring%ele(i)%tracking_method = tracking_method
!      endif
!    enddo
!  endif

  n_points = nx*npz
  allocate(raster_points(n_points))
  allocate(simple_grid(nx,npz))

  ix = 0
  do k=1,npz
    do j=1,nx
      ix = ix + 1
      raster_points(ix)%pz_k = k
      raster_points(ix)%pz = pz_min + (pz_max-pz_min)/(npz-1)*(k-1)
      raster_points(ix)%x = x_min + (x_max-x_min)/(nx-1)*(j-1)
    enddo
  enddo

  if(master) write(*,*) "Starting xpz raster calculation..."
  ring0 = ring
  last_k = -1

  if(master) then
    if(mode == 'singl') then
      do ix = 1, n_points
        call raster_one(ix)
        write(*,'(a,i6,a,i6,a)') "Single mode: point ", ix, " of ", n_points, " complete."
      enddo
    elseif(mode == 'multi') then
      n_sent = 0
      !initial seeding
      do mpi_i=1, min(n_slave,n_points)
        n_sent = n_sent + 1
        call mpi_send(n_sent, 1, MPI_INTEGER, mpi_i, 1, MPI_COMM_WORLD, mpierr)
      enddo
      !receive and reseed
      n_recv = 0
      do while(n_recv .lt. n_points)
        call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
        from_id = mpistatus(MPI_SOURCE)
        call mpi_recv(incoming_ix, 1, MPI_INTEGER, from_id, 2, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(raster_points(incoming_ix)%alive, 1, MPI_DOUBLE_PRECISION, from_id, 3, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(raster_points(incoming_ix)%Dfma, 1, MPI_DOUBLE_PRECISION, from_id, 4, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(raster_points(incoming_ix)%nu_x, 1, MPI_DOUBLE_PRECISION, from_id, 5, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(raster_points(incoming_ix)%nu_y, 1, MPI_DOUBLE_PRECISION, from_id, 6, MPI_COMM_WORLD, mpistatus, mpierr)
        n_recv = n_recv + 1
        if(n_sent .lt. n_points) then
          n_sent = n_sent + 1
          call mpi_send(n_sent, 1, MPI_INTEGER, from_id, 1, MPI_COMM_WORLD, mpierr)
        else
          call mpi_send(-1, 1, MPI_INTEGER, from_id, 1, MPI_COMM_WORLD, mpierr)
        endif
        write(*,'(a,i6,a,i6,a)') "Multi mode: point ", n_recv, " of ", n_points, " complete."
      enddo
    else
      write(*,*) "broken"
      call mpi_finalize(mpierr)
    endif

    !convert data format to simple grid, because keeping track of indices sucks.
    ix = 0
    do k=1,npz
      do j=1,nx
        ix = ix + 1
        simple_grid(j,k) = raster_points(ix)
      enddo
    enddo

    !write accumulated survival data
    open(20,file='raster_xpz.dat')
    open(30,file='raster_xpz_fma.dat')
    open(40,file='raster_xpz_fprint.dat')
    write(20,'(i9,500f9.5)') npz, (simple_grid(1,k)%pz,k=1,npz)
    write(30,'(i9,500f9.5)') npz, (simple_grid(1,k)%pz,k=1,npz)
    write(40,'(a,i6,a,i6)') "# nx= ", nx, " npz=", npz
    do j=1,nx
      write(20,'(f9.5)',advance='no') simple_grid(j,1)%x
      write(30,'(f9.5)',advance='no') simple_grid(j,1)%x
      do k=1,npz
        write(20,'(f9.1)',advance='no') simple_grid(j,k)%alive
        write(30,'(f9.3)',advance='no') simple_grid(j,k)%Dfma
        write(40,'(4f10.6,es14.4)') simple_grid(j,k)%x, simple_grid(j,k)%pz, simple_grid(j,k)%nu_x, simple_grid(j,k)%nu_y, simple_grid(j,k)%Dfma
      enddo
      write(20,*)
      write(30,*)
    enddo
    close(20)
    close(30)
    close(40)
  else
    !slave
    do while(.true.)
      call mpi_recv(ix, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpistatus, mpierr)
      if(ix .lt. 0) then
        exit
      endif
      call raster_one(ix)

      call mpi_send(ix, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, mpierr)
      call mpi_send(raster_points(ix)%alive, 1, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, mpierr)
      call mpi_send(raster_points(ix)%Dfma, 1, MPI_DOUBLE_PRECISION, 0, 4, MPI_COMM_WORLD, mpierr)
      call mpi_send(raster_points(ix)%nu_x, 1, MPI_DOUBLE_PRECISION, 0, 5, MPI_COMM_WORLD, mpierr)
      call mpi_send(raster_points(ix)%nu_y, 1, MPI_DOUBLE_PRECISION, 0, 6, MPI_COMM_WORLD, mpierr)
    enddo
  endif

  call mpi_finalize(mpierr)

  contains
    subroutine raster_one(ix)
      integer ix

      if(raster_points(ix)%pz_k .ne. last_k) then
        ring = ring0
        co(0)%vec = 0.0d0
        co(0)%vec(6) = raster_points(ix)%pz
        call clear_lat_1turn_mats(ring)
        call twiss_and_track(ring,co,status)
        stable = .true.
        if(status .ne. ok$) then
          write(*,*) "Lattice unstable at pz = ", raster_points(ix)%pz
          stable = .false.
        endif
        last_k = raster_points(ix)%pz_k
      endif

      orb(0) = co(0)
      orb(0)%vec(1) = raster_points(ix)%x
      orb(0)%vec(3) = 1.0d-6
      raster_points(ix)%alive = 0.0
      raster_points(ix)%nu_x = -1.0
      raster_points(ix)%nu_y = -1.0
      raster_points(ix)%Dfma = -100.0d0
      if(stable) then
        raster_points(ix)%alive = 1.0
        nix = 0
        toggle = 1
        do m=1,n_turns
          call track_all(ring,orb,track_state=track_state)
          if(track_state == moving_forward$) then
            orb(0)=orb(ring%n_ele_track) 
          else
            raster_points(ix)%alive = 0.0
            exit
          endif
          if(n_naff .gt. 0) then
            if( (m .le. n_naff) .or. (m .gt. n_turns-n_naff)) then
              nix = nix + 1
              naff_turns(nix) = orb(ring%n_ele_track)
            endif
            if(nix .eq. n_naff) then
              !get x tune
              cdata = cmplx(naff_turns(:)%vec(1),naff_turns(:)%vec(2))
              call naff(cdata,freqs,amps,opt_zero_first=.true.)
              nu_x(toggle) = freqs(2)

              !get y tune
              cdata = cmplx(naff_turns(:)%vec(3),naff_turns(:)%vec(4))
              call naff(cdata,freqs,amps,opt_zero_first=.true.)
              nu_y(toggle) = freqs(2)

              if(toggle == 1) then
                raster_points(ix)%nu_x = nu_x(1)
                raster_points(ix)%nu_y = nu_y(1)
              endif

              !reset counters
              toggle = merge(2,1,toggle==1)
              nix = 0
            endif
          endif
        enddo
        if(track_state == moving_forward$) then
          raster_points(ix)%Dfma = log(sqrt((nu_x(1)-nu_x(2))**2+(nu_y(1)-nu_y(2))**2))
        else
          raster_points(ix)%Dfma = -100.0d0
        endif
      endif
    end subroutine
end program
                                                            










