program tracker_icer
  use mpi
  use bmad_parser_mod, only: bp_com
  use bmad
 
  implicit none

  type(lat_struct) lat
  type(coord_struct), allocatable :: co(:), orbit(:)
  real(rp), allocatable :: particle(:,:)

  procedure(track1_custom_def) :: track1_custom

  integer i, j
  integer n_turns, track_state, status
  integer n_part, part_ix

  character(100) lat_file
  character(3) part_ix_str

  real(rp) z_min, z_max, z0, d0

  !mpi housekeeping
  integer num_workers, my_worker_num
  integer task_at_hand, n_tasks
  integer myrank, from_id
  integer cluster_size
  integer mpierr
  integer mpistatus(MPI_STATUS_SIZE)
  integer id, worker_status
  integer n_sent, n_recv
  integer mpi_i, incoming_ix
  integer tag, source_id

  real(rp), allocatable :: results(:,:,:)

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
    num_workers=cluster_size-1
    write(*,*) "Number of slaves: ", num_workers
    if(num_workers .gt. 0) then
      mode = 'multi'
    else
      mode = 'singl'
    endif
  else
    master=.false.
    mode = ''
  endif

  call ran_seed_put(0,myrank)

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, lat)

  !call set_on_off(rfcavity$, lat, on$)
  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .true.

  call twiss_and_track(lat,co,status)

  allocate(orbit(0:lat%n_ele_track))

  n_turns = 100e3
  z_min = -0.5
  z_max =  0.5
  n_part = 1
  d0 = 0.0
  allocate(particle(n_part,6))
  
  !make distribution
  if (n_part > 1) then
    do i=1,n_part
      z0 = z_min + (i-1)*(z_max-z_min)/(n_part-1)
      particle(i,:) = co(0)%vec + (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, z0, d0 /)
    enddo
  else
    z0 = (z_max + z_min) / 2.0d0
    particle(1,:) = co(0)%vec + (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0, z0, d0 /)
  endif 

  if(master) then
    if(mode == 'singl') then
      write(*,*) "Non-mpi environment or cluster size is 1.  Running in single core mode."
      do part_ix = 1, n_part
        write(part_ix_str,'(i3.3)') part_ix
        open(50,file='particle_'//part_ix_str//'.out')
        write(50,'(a8,6a14,a14)') "# turn", "x", "px", "y", "py", "z", "pz", "track state"

        call init_coord(orbit(0), particle(part_ix,:), lat%ele(0), element_end=upstream_end$)
        do i = 1, n_turns
          write(50,'(i8,6es14.5,i14)') i, orbit(0)%vec(1:6), track_state
          call track_all(lat, orbit, track_state=track_state)
          if(track_state /= moving_forward$) then
            write(*,*) "Particle lost at turn ", i
            exit
          endif
          orbit(0) = orbit(lat%n_ele_track)
          !if (mod(i,10000) == 0) then
          !  write(*,'(a,i8,a,i8,a)') "Turn ", i, " of ", n_turns, " complete."
          !endif
        enddo
        close(50)
      enddo
    else
      n_sent = 0
      do mpi_i=1,min(num_workers,n_part)
        n_sent = n_sent + 1
        call mpi_send(mpi_i, 1, MPI_INTEGER, mpi_i, 1, MPI_COMM_WORLD, mpierr)
      enddo

      do while(n_sent < n_part)
        !call mpi_probe(MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
        call mpi_recv(part_ix, 1, MPI_INTEGER, MPI_ANY_SOURCE, 2, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
        source_id = mpistatus(MPI_SOURCE)
        n_sent = n_sent + 1
        call mpi_send(n_sent, 1, MPI_INTEGER, source_id, 1, MPI_COMM_WORLD, mpierr)
      enddo

      do mpi_i=1,num_workers
        call mpi_send(-1, 1, MPI_INTEGER, mpi_i, 4, MPI_COMM_WORLD, mpierr)
      enddo
      call mpi_finalize(mpierr)
    endif
  else
    do while(.true.)
      call mpi_probe(0, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
      tag = mpistatus(MPI_TAG)
      if(tag .eq. 4) exit
      call mpi_recv(part_ix, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpistatus, mpierr)
      write(*,'(a,i2,a,i3)') "worker ", myrank, " processing ", part_ix

      write(part_ix_str,'(i3.3)') part_ix
      open(50,file='particle_'//part_ix_str//'.out')
      write(50,'(a8,6a14,a14)') "# turn", "x", "px", "y", "py", "z", "pz", "track state"
      call init_coord(orbit(0), particle(part_ix,:), lat%ele(0), element_end=upstream_end$)

      do i = 1, n_turns
        call track_all(lat, orbit, track_state=track_state)
        write(50,'(i8,6es14.5,i14)') i, orbit(0)%vec(1:6), track_state
        orbit(0) = orbit(lat%n_ele_track)
      enddo

      close(50)

      call mpi_send(part_ix, 1, MPI_INTEGER, 0, 2, MPI_COMM_WORLD, mpierr)
    enddo
    call mpi_finalize(mpierr)
    write(*,*) "Worker ", myrank, " peacefully exited loop."
  endif
end program






