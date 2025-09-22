program tracker_icer
  use mpi
  use bmad_parser_mod, only: bp_com
  use bmad
 
  implicit none

  type(lat_struct) lat
  type(coord_struct), allocatable :: co(:), orbit(:)

  procedure(track1_custom_def) :: track1_custom

  integer i, n_turns, track_state, status

  character(100) lat_file

  !mpi housekeeping
  integer num_workers, my_worker_num
  integer task_at_hand, n_tasks
  integer myrank, from_id
  integer n_slave, cluster_size
  integer mpierr
  integer mpistatus(MPI_STATUS_SIZE)
  integer id, worker_status
  integer n_sent, n_recv
  integer mpi_i, incoming_ix
  integer, allocatable :: task_list(:,:) 
  integer, allocatable :: tasker(:) ![:]

  real(rp), allocatable :: results(:,:,:) ![:]  !coarray declaration

  logical master
  character*5 mode
  !end mpi housekeeping

  track1_custom_ptr => track1_custom

  call getarg(1,lat_file)

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

  n_turns = 20000000

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, lat)

  !call set_on_off(rfcavity$, lat, on$)
  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .true.

  call twiss_and_track(lat,co,status)

  allocate(orbit(0:lat%n_ele_track))

  open(100,file='tracker_simple.dat')
  write(100,'(a8,6a14,a14)') "# turn", "x", "px", "y", "py", "z", "pz", "track state"

  write(*,'(a,6es14.5)') "Closed Orbit:   ", co(0)%vec
  orbit(0)%vec = co(0)%vec + (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  0.00d0 /)
  !orbit(0)%vec = (/-3.82749E-05, -2.10461E-05,  8.11201E-07, -1.08117E-06,  1.55813E-02, -7.74824E-04/)
  write(*,'(a,6es14.5)') "Initial Vector: ", orbit(0)%vec


  if(master) then
    if(mode == 'singl') then
      write(*,*) "Non-mpi environment or cluster size is 1.  Running in single core mode."
      orbit(0)%vec = co(0)%vec + (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0,  0.00d0 /)
      !orbit(0)%vec = (/-3.82749E-05, -2.10461E-05,  8.11201E-07, -1.08117E-06,  1.55813E-02, -7.74824E-04/)
      write(*,'(a,6es14.5)') "Initial Vector: ", orbit(0)%vec
      do i = 1, n_turns
        call track_all(lat, orbit, track_state=track_state)
        write(100,'(i8,6es14.5,i14)') i, orbit(0)%vec(1:6), track_state
        if(track_state /= moving_forward$) then
          write(*,*) "Particle lost at turn ", i
          exit
        endif
        orbit(0) = orbit(lat%n_ele_track)
        if (mod(i,10000) == 0) then
          write(*,'(a,i8,a,i8,a)') "Turn ", i, " of ", n_turns, " complete."
        endif
      enddo
    else
      n_sent = 0
      do mpi_i=1,min(n_workers,n_part)
        n_sent = n_sent + 1
        call mpi_send(mpi_i, 1, MPI_INTEGER, mpi_i, 1, MPI_COMM_WORLD, mpierr)
      enddo
    endif
  else
    do while(.true.)
      call mpi_probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
      tag = mpistatus(MPI_TAG)
      if(tag .eq. 4) exit
      call mpi_recv(part_ix, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpistatus, mpierr)

      do i = 1, n_turns
        call track_all(lat, orbit, track_state=track_state)
        write(100,'(i8,6es14.5,i14)') i, orbit(0)%vec(1:6), track_state
        if(track_state /= moving_forward$) then
          write(*,*) "Particle lost at turn ", i
          exit
        endif
        orbit(0) = orbit(lat%n_ele_track)
        if (mod(i,10000) == 0) then
          write(*,'(a,i8,a,i8,a)') "Turn ", i, " of ", n_turns, " complete."
        endif
      enddo

      call mpi_send(track_state, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, mpierr) 
      if( track_state .eq. moving_forward$) then
        call mpi_send(vec, 6, MPI_DOUBLE_PRECISION, 0, 3, MPI_COMM_WORLD, mpierr)
      endif
    enddo
    write(*,*) "Worker ", myrank, " peacefully exited loop."
  endif

  close(100)
end program






