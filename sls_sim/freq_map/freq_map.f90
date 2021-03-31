program freq_map

  use bmad
  use bmad_parser_mod, only: bp_com
  use four_mod
  use naff_mod
  use sim_utils
  use ptc_layout_mod
  use madx_ptc_module

  use namelist_general !general: lat_file, use_hybrid
  use namelist_fm !fm: dE, xmax, ymax, n_turns, fft_turns, nx, ny
 
  implicit none

  type(lat_struct) ring
  type (lat_struct) ring_hybrid
  type (lat_struct) ring_use
  type(coord_struct), allocatable :: co(:)
  type(coord_struct), allocatable :: orbit(:)
  type(coord_struct), allocatable :: tbt_coords(:)
  type(coord_struct), allocatable :: norm_coords(:)

  type(taylor_struct) one_turn_map(6)
  type(taylor_struct) A(1:6), A_inverse(1:6), dhdj(1:6)

  integer i
  integer ix, iy
  integer track_state
  integer num_workers, my_worker_num !coarray housekeeping
  integer, allocatable :: task_list(:,:) 
  integer task_at_hand, n_tasks
  integer id, worker_status
  integer status

  real(rp) x, y, dx, dy
  real(rp) pz
  real(rp) init_vec(6)
  real(rp), allocatable :: results(:,:,:,:)[:]  !coarray declaration
  integer, allocatable :: tasker(:)[:]
  real(rp) metric

  real(rp) tunes(5)
  complex(rp) amps(5)
  integer freqs_found

  logical rf_on
  logical ok

  character(100) in_file
  character(100) lat_file_override

  integer iargs

  logical, allocatable :: keep_ele(:)

  use_hybrid = .false.

  iargs = iargc()
  lat_file_override = ''
  if( iargs == 1 ) then
    call getarg(1,in_file)
  elseif (iargs == 2 ) then
    call getarg(1,in_file)
    call getarg(2,lat_file_override)
  endif

  ! coarray stuff
  num_workers = num_images()
  my_worker_num = this_image()

  call getarg(1,in_file)
  open (unit = 10, file = in_file)
  read (10, nml = general)
  read (10, nml = fm)
  close (10)

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  if( my_worker_num == 1) WRITE(*,*) "Preparing lattice..."

  bp_com%always_parse = .true.
  call bmad_parser(lat_file,ring)
  
  call set_on_off(rfcavity$, ring, off$)
  bmad_com%aperture_limit_on = .true.
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  call lat_twiss_and_track(ring,co,status)
  
  allocate(results(nx,ny,2,3)[*])  !coarray allocation
  allocate(tasker(num_workers)[*])  !coarray allocation
  allocate(tbt_coords(n_turns))
  allocate(norm_coords(n_turns))
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

    if(use_hybrid) then
      allocate(keep_ele(ring%n_ele_max))
      keep_ele = .false.
      do i=1,ring%n_ele_track
        if( (ring%ele(i)%key == sbend$) .or. &
            (ring%ele(i)%key == sextupole$) .or. &
            (ring%ele(i)%key == rfcavity$) .or. &
            (ring%ele(i)%key == multipole$) ) then
          call update_hybrid_list(ring, i, keep_ele)
        endif
      enddo
      call make_hybrid_lat(ring, keep_ele, .true., ring_hybrid)
      ring_use = ring_hybrid
    else
      ring_use = ring
    endif

    do while(.true.)
      tasker(my_worker_num)[1] = -1
      sync images(1)
      task_at_hand = tasker(my_worker_num)[my_worker_num]
      if( task_at_hand == -99 ) exit

      iy = task_list(task_at_hand,1)
      ix = task_list(task_at_hand,2)
      write(*,'(A,I4,A,I6,A,I6,A,I6,A,I6)') "Image ", my_worker_num, " Processing ", iy, ", ", ix, " of ", ny, ", ", nx
      y = dy*iy
      x = -xmax + dx*(ix-1)

      init_vec = (/ x, 0.0d0, y, 0.0d0, 0.0d0, dE /)
      orbit(0)%vec = co(0)%vec + init_vec

      do i = 1, n_turns
        call track_all(ring_use, orbit, track_state=track_state)

        if(track_state /= moving_forward$) then
          exit
        endif

        tbt_coords(i)%vec(:) = orbit(0)%vec(:) - co(0)%vec(:)
        orbit(0) = orbit(ring_use%n_ele_track)
      enddo

      if( track_state == moving_forward$ ) then
        do i=1,n_turns
          call xy_to_action(ring_use, 0, tbt_coords(i)%vec, norm_coords(i)%vec, ok)
          if (.not. ok) then
            write(*,*) "Error from xy_to_action"
            stop
          endif
        enddo

        ! !subtract centroid
        ! do i=1,4
        !   norm_coords(1:fft_turns)%vec(i) = norm_coords(1:fft_turns)%vec(i) - sum(norm_coords(1:fft_turns)%vec(i))/size(norm_coords(1:fft_turns)%vec(i))
        !   norm_coords(1+n_turns-fft_turns:n_turns)%vec(i) = norm_coords(1+n_turns-fft_turns:n_turns)%vec(i) - sum(norm_coords(1+n_turns-fft_turns:n_turns)%vec(i))/size(norm_coords(1+n_turns-fft_turns:n_turns)%vec(i))
        ! enddo

        !naff for strongest frequency component
        call naff( norm_coords(1:fft_turns)%vec(1), norm_coords(1:fft_turns)%vec(2), tunes(:), amps(:), freqs_found)
        results(ix,iy,1,1)[1] = tunes(1)
        call naff( norm_coords(1+n_turns-fft_turns:n_turns)%vec(1), norm_coords(1+n_turns-fft_turns:n_turns)%vec(2), tunes(:), amps(:), freqs_found)
        results(ix,iy,1,2)[1] = tunes(1)
        call naff( norm_coords(1:fft_turns)%vec(3), norm_coords(1:fft_turns)%vec(4), tunes(:), amps(:), freqs_found)
        results(ix,iy,2,1)[1] = tunes(1)
        call naff( norm_coords(1+n_turns-fft_turns:n_turns)%vec(3), norm_coords(1+n_turns-fft_turns:n_turns)%vec(4), tunes(:), amps(:), freqs_found)
        results(ix,iy,2,2)[1] = tunes(1)

        metric = log10(sqrt( (results(ix,iy,1,2)[1] - results(ix,iy,1,1)[1] )**2 + (results(ix,iy,2,2)[1] - results(ix,iy,2,1)[1] )**2 ))
        results(ix,iy,1,3)[1] = metric
      else
        results(ix,iy,1,3)[1] = 999.
      endif
    enddo
  endif

  sync all
  write(*,*) "Image ", my_worker_num, " made it!"
  if( my_worker_num == 1 ) then
    open(23,FILE='freq_map.dat')
    write(23,'(7A16)') "# x", "y", "initial Qx", "last Qx", "initial Qy", "last Qy", "Freq Map Metric"
    do iy = 1, ny
      do ix = 1, nx
        y = dy*iy
        x = -xmax + dx*(ix-1)
        if ( results(ix,iy,1,3) .lt. 900. ) then
          write(23, '(2F16.6,5ES16.5)') x, y, results(ix,iy,1,1), results(ix,iy,1,2), results(ix,iy,2,1), results(ix,iy,2,2), results(ix,iy,1,3)
        else
          write(23, '(2F16.6,5A16)')    x, y, "NaN", "NaN", "NaN", "NaN", "NaN"
        endif
      enddo
    enddo
    close(23)
  endif
end program






