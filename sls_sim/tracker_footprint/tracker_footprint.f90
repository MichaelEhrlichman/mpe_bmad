program tracker_footprint
  use mpi
  use bmad
  use bmad_parser_mod, only: bp_com
  use sls_lib
  use naff_mod
  use sim_utils

  use namelist_general !general: lat_file, use_hybrid
  use namelist_fm !fm: dE, xmax, ymax, n_turns, fft_turns, nx, ny
 
  implicit none

  type(lat_struct) ring
  type(coord_struct), allocatable :: co(:)
  type(coord_struct), allocatable :: orbit(:)
  type(coord_struct), allocatable :: tbt_coords(:)
  type(coord_struct), allocatable :: norm_coords(:)
  complex(rp), allocatable :: cdata(:)

  integer i
  integer ix, iy
  integer track_state, status
  integer num_workers, my_worker_num !coarray housekeeping
  integer, allocatable :: task_list(:,:) 
  integer task_at_hand, n_tasks
  integer id, worker_status

  real(rp) x, y, dx, dy
  real(rp) pz
  real(rp) init_vec(6)
  real(rp), allocatable :: results(:,:,:) ![:]  !coarray declaration
  integer, allocatable :: tasker(:) ![:]

  real(rp) tunes(1)
  complex(rp) amps(1)
  integer freqs_found

  logical ok, err

  character(100) in_file
  character(100) lat_file_override

  integer iargs

  !mpi housekeeping
  integer myrank, from_id
  integer n_slave, cluster_size
  integer mpierr
  integer mpistatus(MPI_STATUS_SIZE)
  logical master
  character*5 mode
  integer n_sent, n_recv
  integer mpi_i, incoming_ix

  iargs = iargc()
  lat_file_override = ''
  if( iargs == 1 ) then
    call getarg(1,in_file)
  elseif (iargs == 2 ) then
    call getarg(1,in_file)
    call getarg(2,lat_file_override)
  endif

  ! mpi stuff
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


  open (unit = 10, file = in_file)
  read (10, nml = general)
  read (10, nml = fm)
  close (10)

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  bp_com%always_parse = .true.
  call bmad_parser(lat_file,ring)
  
  call set_on_off(rfcavity$, ring, off$)
  bmad_com%aperture_limit_on = .true.
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  allocate(co(0:ring%n_ele_track))
  co(0)%vec(6) = dE
  call twiss_and_track(ring,co,status)
 
  allocate(results(nx,ny,2)[*])  !coarray allocation
  allocate(tasker(num_workers)[*])  !coarray allocation
  allocate(tbt_coords(naff_turns))
  allocate(norm_coords(naff_turns))
  allocate(cdata(naff_turns))
  allocate(orbit(0:ring%n_ele_track))

  if( my_worker_num == 1 ) write(*,*) "Starting simulation..."

  n_tasks = ny*nx
  allocate(task_list(ny*nx,2))
  do iy = 1, ny
    do ix = 1, nx
      task_list((iy-1)*nx+ix,:) = (/iy, ix/)
    enddo
  enddo

  tasker(:) = 0

  dy = ymax/ny
  dx = 2.0*xmax/(nx-1)
  if( my_worker_num == 1 ) then
    sync all
    id = 2
    do i = 1, n_tasks
      do while(.true.)
        worker_status = tasker(id)[1]
        if(worker_status == -1) then
          tasker(id)[id] = i
          tasker(id)[1] = 0
          sync images(id)
          exit
        endif
        id = id + 1
        if(id .gt. num_workers) id = 2
        call sleepqq(2)
      enddo
    enddo
    do i=1,num_workers
      tasker(i)[i] = -99
    enddo
    sync images(*)
  else
    sync all
    do while(.true.)
      tasker(my_worker_num)[1] = -1
      sync images(1)
      task_at_hand = tasker(my_worker_num)[my_worker_num]
      if( task_at_hand == -99 ) exit

      iy = task_list(task_at_hand,1)
      ix = task_list(task_at_hand,2)
      write(*,'(A,I4,A,I6,A,I6,A,I6,A,I6)') "Image ", my_worker_num, " Processing ", iy, ", ", ix, " of ", ny, ", ", nx
      y = dy*iy + 1.0e-6
      x = -xmax + dx*(ix-1) + 1.0e-6

      init_vec = (/ x, 0.0d0, y, 0.0d0, 0.0d0, 0.0d0 /)
      orbit(0)%vec = co(0)%vec + init_vec

      do i = 1, pre_turns
        call track_all(ring, orbit, track_state=track_state)
        if(track_state /= moving_forward$) then
          exit
        endif
        orbit(0) = orbit(ring%n_ele_track)
      enddo

      if(track_state == moving_forward$ ) then
        do i = 1, naff_turns
          call track_all(ring, orbit, track_state=track_state)
          if(track_state /= moving_forward$) then
            exit
          endif
          tbt_coords(i)%vec(:) = orbit(0)%vec(:) - co(0)%vec(:)
          orbit(0) = orbit(ring%n_ele_track)
        enddo
      endif

      if( track_state == moving_forward$ ) then
        do i=1,naff_turns
          call xy_to_action(ring, 0, tbt_coords(i)%vec, norm_coords(i)%vec, ok)
          if (.not. ok) then
            write(*,*) "Error from xy_to_action"
            stop
          endif
        enddo

        !naff for strongest frequency component
        cdata = cmplx(norm_coords(:)%vec(1), norm_coords(:)%vec(2))
        call naff( cdata, tunes(:), amps(:), opt_zero_first=.true.)
        results(ix,iy,1)[1] = tunes(2)
        cdata = cmplx(norm_coords(:)%vec(3), norm_coords(:)%vec(4))
        call naff( cdata, tunes(:), amps(:), opt_zero_first=.true.)
        results(ix,iy,2)[1] = tunes(2)
      else
        results(ix,iy,1)[1] = 999.
        results(ix,iy,2)[1] = 999.
      endif
    enddo
  endif

  sync all
  write(*,*) "Image ", my_worker_num, " made it!"
  if( my_worker_num == 1 ) then
    open(23,FILE='tracker_footprint.dat')
    write(23,'(4A16)') "# x", "y", "Qx", "Qy"
    do iy = 1, ny
      do ix = 1, nx
        y = dy*iy
        x = -xmax + dx*(ix-1)
        if ( results(ix,iy,1) .lt. 900. ) then
          write(23, '(2F16.6,2ES16.5)') x, y, results(ix,iy,1), results(ix,iy,2)
        else
          write(23, '(2F16.6,2A16)')    x, y, "NaN", "NaN"
        endif
      enddo
    enddo
    close(23)
  endif
  sync all
end program






