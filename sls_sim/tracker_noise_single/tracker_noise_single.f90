program tracker_noise_single
  use bmad
  use sls_lib
  use naff_mod, except_rp=>rp
  use apfft_mod
  use noise_mod
  use fgsl
  use mode3_mod
 
  implicit none

  integer, parameter :: N_MAX_NOISE = 100

  type(lat_struct) lat, lat_ideal
  type(coord_struct), allocatable :: co(:)
  type(coord_struct), allocatable :: orbit(:)
  type(coord_struct), allocatable :: tbt_coords(:)
  type(coord_struct), allocatable :: norm_coords(:)
  type(ele_noise_struct_in) ele_noise_in(N_MAX_NOISE)
  type(ele_noise_struct), allocatable :: ele_noise(:)

  integer i, j, k, m
  integer state
  integer ix_a, ix_b
  integer n_seg
  integer status
  integer idx, data_col, n_file_columns

  character(5) spec_algorithm

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

  real(rp) t6(6,6), N(6,6)
  real(rp) abz_tunes(3)
  logical err

  !check
  complex(rp), allocatable :: check_data(:), xdft(:)
  real(rp), allocatable :: psd(:)
  character(20) check_attribute
  integer check_ele_ix
  integer(fgsl_size_t) check_n
  integer half
  real(rp) frev

  !apfft
  real(rp) x_bounds(2), y_bounds(2)
  real(rp) x_phase, x_amp_, x_freq
  real(rp) y_phase, y_amp_, y_freq

  real(rp) vec(6)
  real(rp) harvest_r
  real(rp) real_bucket
  real(rp) f

  logical ok, radiation_damping, radiation_excitation

  character(100) lat_file
  character(100) in_file
  character(100) lat_file_override

  integer iargs

  !parameters from .in
  integer ix_inj, ix_obs
  real(rp) init_vec(6)
  integer harvest_i
  integer n_turns
  integer block_size
  integer dims
  integer seed

  !fgsl vars
  type(fgsl_fft_complex_wavetable) :: wavetable
  type(fgsl_fft_complex_workspace) :: work

  !for spectrum data file parsing
  integer io
  character(2000) line
  real(rp) line_items(100)

  integer nnoise

  namelist / tracker_noise_params / lat_file, ix_inj, ix_obs, init_vec, n_turns, block_size, &
                                    dims, radiation_damping, radiation_excitation, ele_noise_in, seed, spec_algorithm, &
                                    x_bounds, y_bounds, check_n, check_attribute, check_ele_ix

  check_n = -1
  check_attribute = ''
  check_ele_ix = -1
  iargs = iargc()
  lat_file_override = ''
  if( iargs == 1 ) then
    call getarg(1,in_file)
  elseif (iargs == 2 ) then
    call getarg(1,in_file)
    call getarg(2,lat_file_override)
  endif

  dims = 4
  radiation_damping = .false.
  radiation_excitation = .false.
  seed = 12345
  spec_algorithm = 'undef'
  x_bounds = 0.0d0
  y_bounds = 0.0d0
  ix_obs = 0

  open (unit = 10, file = in_file)
  read (10, nml = tracker_noise_params)
  close (10)

  nnoise = 0
  do i=1, N_MAX_NOISE
    if (len(trim(ele_noise_in(i)%attrib_name)) .eq. 0) exit
    nnoise = nnoise + 1
  enddo
  allocate(ele_noise(nnoise)) 
  do i=1, nnoise
    ele_noise(i)%ele_identifier         = ele_noise_in(i)%ele_identifier
    ele_noise(i)%attrib_name            = ele_noise_in(i)%attrib_name
    ele_noise(i)%dist_type              = ele_noise_in(i)%dist_type
    ele_noise(i)%error_type             = ele_noise_in(i)%error_type
    ele_noise(i)%n_turns_recalc         = ele_noise_in(i)%n_turns_recalc
    ele_noise(i)%nrm                    = ele_noise_in(i)%nrm
    ele_noise(i)%dfreq                  = ele_noise_in(i)%dfreq
    ele_noise(i)%start_freq             = ele_noise_in(i)%start_freq
    ele_noise(i)%random_phase_element   = ele_noise_in(i)%random_phase_element
    ele_noise(i)%phase_gang_element     = ele_noise_in(i)%phase_gang_element
    ele_noise(i)%random_phase_component = ele_noise_in(i)%random_phase_component
    ele_noise(i)%spec_data_file         = ele_noise_in(i)%spec_data_file
    if(ele_noise(i)%spec_data_file == 'none') then
      ele_noise(i)%spec_data(:)         = ele_noise_in(i)%spec_data(:)
    else
      idx = index_nocase(ele_noise(i)%spec_data_file,'::')
      if (idx == 0) then
        data_col = 1
      else
        read(ele_noise(i)%spec_data_file(idx+2:),*) data_col
        ele_noise(i)%spec_data_file = ele_noise(i)%spec_data_file(1:idx-1)
      endif
      !--------------------------------------------
      !count number of columns of data in data file
      !--------------------------------------------
      open(800,file=ele_noise(i)%spec_data_file,action='read')
      read(800,'(a)') line !skip header
      read(800,'(a)',iostat=io) line
      if(io < 0) then
        write(*,*) "Error"
        stop
      endif
      do j=1,100
        read(line,*,iostat=io) real_bucket, line_items(1:j)
        if(io==-1) exit
      enddo
      n_file_columns = j-1
      rewind(800)
      !--------------------------------------------
      !read data file
      !--------------------------------------------
      read(800,'(a)') line !skip header
      j = 0
      do while(.true.)
        line_items = 0.0d0
        read(800,'(a)',iostat=io) line
        if(io < 0) exit
        read(line,*) f, line_items(1:n_file_columns)
        j = j + 1
        ele_noise(i)%spec_data(j)%f = f
        ele_noise(i)%spec_data(j)%a = line_items(data_col)
      enddo
      close(800)
    endif
  enddo
  open(5050,file='noise_sources.info')
  do i=1,nnoise
    write(5050,*) "i: ", i
    write(5050,*) "element identifier: ", ele_noise(i)%ele_identifier
    write(5050,*) "attribute name: ", ele_noise(i)%attrib_name
    write(5050,*) "distribution type: ", ele_noise(i)%dist_type
    write(5050,*) "error type: ", ele_noise(i)%error_type
    write(5050,*) "N turns recalc: ", ele_noise(i)%n_turns_recalc
    write(5050,*) "Normal distribution parameters: ", ele_noise(i)%nrm
    write(5050,*) "starting frequency: ", ele_noise(i)%start_freq
    write(5050,*) "delta frequency: ", ele_noise(i)%dfreq
    write(5050,*) "random phase per element: ", ele_noise(i)%random_phase_element
    write(5050,*) "random phase per component: ", ele_noise(i)%random_phase_component
    write(5050,*) "spectrum data file: ", trim(ele_noise_in(i)%spec_data_file)
    write(5050,'(a,5es14.4)') " first 5 amplitudes: ", ele_noise(i)%spec_data(1:5)%a
    write(5050,*)
  enddo
  close(5050)

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  write(*,*) "Preparing lattice..."

  bmad_com%auto_bookkeeper = .false.
  call bmad_parser(lat_file, lat)
  
  if(dims==4) then
    CALL set_on_off(rfcavity$, lat, off$)
  elseif (dims==6) then
    CALL set_on_off(rfcavity$, lat, on$)
  endif
  bmad_com%radiation_damping_on = radiation_damping
  bmad_com%radiation_fluctuations_on = radiation_excitation
  bmad_com%rel_tol_tracking = 1e-12
  bmad_com%abs_tol_tracking = 1e-14
  
  call twiss_and_track(lat,co,status)

  !call transfer_matrix_calc (lat, t6, ix1=ix_obs, one_turn=.true.)
  !call calc_z_tune(lat)
  !abz_tunes(1) = lat%a%tune
  !abz_tunes(2) = lat%b%tune
  !abz_tunes(3) = lat%z%tune
  !call make_N(t6, N, err, abz_tunes)

  allocate(tbt_coords(n_turns))
  allocate(norm_coords(n_turns))
  allocate(cdata(n_turns))
  allocate(orbit(0:lat%n_ele_track))

  write(*,*) "Tracking..."

  lat_ideal = lat

  if(check_n > 0 ) then
    allocate(check_data(check_n))
  endif

  call ran_seed_put(seed)
  call init_ele_noise(lat, ele_noise)

  ! Tracking
  call init_coord(orbit(ix_inj), co(ix_inj)%vec+init_vec, ele=lat%ele(ix_inj), element_end=upstream_end$, &
                  particle=electron$)
  do i = 1, n_turns
    call apply_ele_noise(lat, lat_ideal, ele_noise, i)
    if(check_ele_ix .gt. 0) then
      if( i .lt. check_n ) then
        check_data(i) = cmplx(value_of_attribute(lat%ele(check_ele_ix), check_attribute, err_print_flag = .true.), 0.0d0)
      endif
    endif

    if(mod(i,100).eq.0) write(*,'(a,a,i9,a,i9)',advance='no') achar(13), "turn ", i, " of ", n_turns

    call track_many(lat, orbit, ix_inj, ix_inj, 1, track_state=state)
    if(state /= moving_forward$) then
      write(*,*) "Particle Lost"
      exit
    endif
    tbt_coords(i)%vec(:) = orbit(ix_obs)%vec(:)
  enddo
  write(*,*) 

  open(100,file='tbt.dat')
  write(100,'(A6,6A14)') "# turn", "x", "px", "y", "py", "z", "pz"
  do i=1,n_turns
    write(100,'(I9,6ES16.7)') i, tbt_coords(i)%vec(1:6)
  enddo
  close(100)

  if(check_ele_ix .gt. 0) then
    ! Compute spectrum of checked element attribute
    check_data = check_data - sum(check_data)/size(check_data)
    write(*,'(2a)') "Calculating spectrum of element ", lat%ele(check_ele_ix)%name
    open(201,file='check.element')
    do i=1,check_n
      write(201,'(i8,es14.6)') i, real(check_data(i))
    enddo
    close(201)
    status = fgsl_fft_complex_radix2_forward(check_data,1_fgsl_size_t,check_n)
    half = (check_n/2)+1
    allocate(xdft(1:half))
    xdft = check_data(1:half)
    frev = c_light/lat%ele(lat%n_ele_track)%s
    open(200,file='check.fft')
    do i=1,check_n
      write(200,'(2es14.6)') (i-1.0d0)/check_n*frev, abs(check_data(i)/sqrt(dble(check_n)))
    enddo
    close(200)
    allocate(psd(size(xdft)))
    psd = abs(xdft)**2 / frev / check_n
    psd(2:half) = 2*psd(2:half)
    open(201,file='check.psd')
    do i=1,half
      write(201,'(2es14.6)') (i-1.0d0)/check_n*frev, psd(i)
    enddo
    close(201)
  endif

  ! Analysis
  if( state == moving_forward$ ) then
    do i = 1, n_turns
      call xy_to_action(lat_ideal, ix_obs, tbt_coords(i)%vec, norm_coords(i)%vec, ok) 
      !norm_coords(i)%vec = matmul(mat_symp_conj(N), tbt_coords(i)%vec)
      ok=.true.
      if (.not. ok) then
        write(*,*) "Error from xy_to_action"
        stop
      endif
    enddo

    open(100,file='tbt_Anorm.dat')
    write(100,'(A6,4A14)') "# turn", "Jx", "Jx'", "Jy", "Jy'"
    write(100,'(a,es14.6)') "# Jx is ", (norm_coords(1)%vec(1)**2 + norm_coords(1)%vec(2)**2)/2.0d0
    write(100,'(a,es14.6)') "# Jy is ", (norm_coords(1)%vec(3)**2 + norm_coords(1)%vec(4)**2)/2.0d0
    do i=1,n_turns
      write(100,'(I9,4ES14.5)') i, norm_coords(i)%vec(1:4)
    enddo
    close(100)

    n_seg = n_turns/block_size
    open(100,file='tbt_x.fft')
    open(101,file='tbt_y.fft')
    write(100,*) '# apfft'
    write(101,*) '# apfft'
    write(100,'(4a18)') "# component", "tune", "amp", "angle (rad)"
    write(101,'(4a18)') "# component", "tune", "amp", "angle (rad)"
    do i=1, n_seg
      ix_a = 1 + (i-1)*block_size
      ix_b = i*block_size
      call apfft_corr(norm_coords(ix_a:ix_b)%vec(1), x_bounds, 'han', x_phase, x_amp_, x_freq)
      call apfft_corr(norm_coords(ix_a:ix_b)%vec(3), y_bounds, 'han', y_phase, y_amp_, y_freq)
      write(100,'(f18.1,3es18.8)') (i-0.5)*block_size, x_freq, x_amp_, x_phase
      write(101,'(f18.1,3es18.8)') (i-0.5)*block_size, y_freq, y_amp_, y_phase
    enddo
    close(100)
    close(101)

    !subtract centroid
    norm_coords(:)%vec(1) = norm_coords(:)%vec(1) - sum(norm_coords(:)%vec(1))/size(norm_coords(:)%vec(1))
    norm_coords(:)%vec(2) = norm_coords(:)%vec(2) - sum(norm_coords(:)%vec(2))/size(norm_coords(:)%vec(2))
    norm_coords(:)%vec(3) = norm_coords(:)%vec(3) - sum(norm_coords(:)%vec(3))/size(norm_coords(:)%vec(3))
    norm_coords(:)%vec(4) = norm_coords(:)%vec(4) - sum(norm_coords(:)%vec(4))/size(norm_coords(:)%vec(4))

    write(*,*) "Spectral Analysis ..."

    ! Spectral Analysis of Beam Motion
    if(match_reg('naff',spec_algorithm)) then
      n_seg = n_turns/block_size
      allocate(x_tunes(n_seg,max_naff))
      allocate(y_tunes(n_seg,max_naff))
      allocate(x_amp(n_seg,max_naff))
      allocate(y_amp(n_seg,max_naff))
      x_tunes = 0.0d0
      y_tunes = 0.0d0
      do i=1, n_seg
        ix_a = 1 + (i-1)*block_size
        ix_b = i*block_size
        x_basis = 0.0d0
        x_norm = 0.0d0
        cdata = cmplx(norm_coords(:)%vec(1),norm_coords(:)%vec(2))
        call naff( cdata(ix_a:ix_b), x_tunes(i,:), x_amp(i,:))
        !call naff( cdata(ix_a:ix_b), x_tunes(i,:), x_amp(i,:), opt_dump_spectra=700)
        cdata = cmplx(norm_coords(:)%vec(3),norm_coords(:)%vec(4))
        call naff( cdata(ix_a:ix_b), y_tunes(i,:), y_amp(i,:))
        !call naff( cdata(ix_a:ix_b), y_tunes(i,:), y_amp(i,:), opt_dump_spectra=800)
      enddo

      open(100,file='tunes_x_by_turn.dat')
      open(101,file='tunes_y_by_turn.dat')
      write(100,*) '# naff'
      write(101,*) '# naff'
      write(100,'(4a18)') "# segment", "tune", "amp", "angle (rad)"
      write(101,'(4a18)') "# segment", "tune", "amp", "angle (rad)"
      do i=1,n_seg
        do j=1,max_naff
          write(100,'(f18.1,3es18.8)') (i-0.5)*block_size, x_tunes(i,j), abs(x_amp(i,j)), atan2(aimag(x_amp(i,j)),real(x_amp(i,j)))
          write(101,'(f18.1,3es18.8)') (i-0.5)*block_size, y_tunes(i,j), abs(y_amp(i,j)), atan2(aimag(y_amp(i,j)),real(y_amp(i,j)))
        enddo
        write(100,*)
        write(100,*)
        write(101,*)
        write(101,*)
      enddo
      close(100)
      close(101)

      open(100,file='tunes_x_by_component.dat')
      open(101,file='tunes_y_by_component.dat')
      write(100,*) '# naff'
      write(101,*) '# naff'
      write(100,'(4a18)') "# component", "tune", "amp", "angle (rad)"
      write(101,'(4a18)') "# component", "tune", "amp", "angle (rad)"
      do j=1,max_naff
        do i=1,n_seg
          write(100,'(f18.1,3es18.8)') (i-0.5)*block_size, x_tunes(i,j), abs(x_amp(i,j)), atan2(aimag(x_amp(i,j)),real(x_amp(i,j)))
          write(101,'(f18.1,3es18.8)') (i-0.5)*block_size, y_tunes(i,j), abs(y_amp(i,j)), atan2(aimag(y_amp(i,j)),real(y_amp(i,j)))
        enddo
        write(100,*)
        write(100,*)
        write(101,*)
        write(101,*)
      enddo
      close(100)
      close(101)

      deallocate(x_tunes)
      deallocate(y_tunes)
      deallocate(x_amp)
      deallocate(y_amp)
    elseif(match_reg('apfft',spec_algorithm)) then
      n_seg = n_turns/block_size
      open(100,file='tunes_x_by_turn.dat')
      open(101,file='tunes_y_by_turn.dat')
      write(100,*) '# apfft'
      write(101,*) '# apfft'
      write(100,'(4a18)') "# turn", "tune", "amp", "angle (rad)"
      write(101,'(4a18)') "# turn", "tune", "amp", "angle (rad)"
      do i=1, n_seg
        ix_a = 1 + (i-1)*block_size
        ix_b = i*block_size
        call apfft_corr(norm_coords(ix_a:ix_b)%vec(1), x_bounds, 'han', x_phase, x_amp_, x_freq)
        call apfft_corr(norm_coords(ix_a:ix_b)%vec(3), y_bounds, 'han', y_phase, y_amp_, y_freq)
        write(100,'(f18.1,3es18.8)') (i-0.5)*block_size, x_freq, x_amp_, x_phase
        write(101,'(f18.1,3es18.8)') (i-0.5)*block_size, y_freq, y_amp_, y_phase
      enddo
      close(100)
      close(101)
    else
      write(*,*) "Unknown spec_algorithm: ", spec_algorithm
      write(*,*) "No spectral analysis done."
    endif
  else
    write(*,*) "Particle was lost.  Not calculating tunes."
  endif

end program






