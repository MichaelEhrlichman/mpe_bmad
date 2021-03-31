program moga
  use mpi
  use bmad
  use bmad_parser_mod, only: bp_com
  use custom_dynamic_aperture_mod
  use pisa_mod
  use dynap_mod
  use linear_aperture_mod
 
  implicit none

  !namelist_general
  character*100 lat_file
  logical use_hybrid
  namelist / general /    lat_file, use_hybrid

  !namelist_da
  integer, parameter :: max_dE = 11
  integer tracking_method
  integer n_adts
  integer n_turn
  integer n_angle
  integer track_dims
  real(rp) dE(max_dE)
  real(rp) init_len
  real(rp) adts_x_min, adts_x_max
  real(rp) adts_y_min, adts_y_max
  
  namelist / da /     tracking_method, &
                      n_adts, &
                      n_turn, &          !Number of turns particle must survive
                      n_angle, &
                      track_dims, &            !either 4 or 6
                      dE, &
                      adts_x_min, &
                      adts_x_max, &
                      adts_y_min, &
                      adts_y_max, &
                      init_len

  !namelist_moga
  integer, parameter :: max_mags = 200
  
  character*100 moga_output_file
  real(rp) set_chrom_x
  real(rp) set_chrom_y
  character*100 initial_pop
  integer seed
  integer generate_feasible_seeds_only
  type(breeder_params_struct) breeder_params
  integer max_gen
  real(rp) co_limit
  real(rp) linear_vec_cutoff
  real(rp) x_fp_min, x_fp_max
  real(rp) y_fp_min, y_fp_max
  real(rp) fp_dE_neg, fp_dE_pos
  real(rp) work_pt_x_min, work_pt_x_max
  real(rp) work_pt_y_min, work_pt_y_max
  integer n_fp_steps
  type(mag_struct) mags_in(max_mags)
  type(mag_struct) chrom_mags(2)
  real(rp) emittance_target
  real(rp) beta_x_max, beta_y_max, eta_x_abs_max, global_beta_x_max
  real(rp) tune_delta, coupling_delta_max
  real(rp) chrom_x_min, chrom_x_max, chrom_y_min, chrom_y_max
  logical coupling_con, use_chrom_mags
  
  Namelist / nl_moga /    moga_output_file, &
                          generate_feasible_seeds_only, &
                          set_chrom_x, &
                          set_chrom_y, &
                          emittance_target, &
                          initial_pop, &
                          seed, &
                          breeder_params, &
                          max_gen, &
                          co_limit, &
                          linear_vec_cutoff, &
                          work_pt_x_min, &
                          work_pt_x_max, &
                          work_pt_y_min, &
                          work_pt_y_max, &
                          x_fp_min, &
                          x_fp_max, &
                          y_fp_min, &
                          y_fp_max, &
                          fp_dE_neg, &
                          fp_dE_pos, &
                          n_fp_steps, &
                          mags_in, &
                          chrom_mags, &
                          global_beta_x_max, &
                          beta_x_max, &
                          beta_y_max, &
                          eta_x_abs_max, &
                          coupling_delta_max, &
                          coupling_con, &
                          use_chrom_mags, &
                          chrom_x_min, chrom_x_max, chrom_y_min, chrom_y_max

  type (lat_struct) ring
  type (lat_struct) ring0
  type (lat_struct), save :: ring_on_energy
  type (lat_struct), save :: ring_use
  type (ele_pointer_struct), allocatable :: eles(:)
  type (custom_aperture_scan_struct) da_config
  type (custom_aperture_scan_struct) da_block_linear
  type (coord_struct), allocatable :: co(:)
  type (coord_struct), allocatable :: co_on_energy(:)
  type (coord_struct), allocatable :: orb(:)
  type(normal_modes_struct) mode

  integer i,j,k
  integer n_dE
  integer n_feasible
  integer time_stamp(8)
  integer n_aperture_test
  integer rsize
  integer ix_cache
  integer coupling_con_ix, chrom_cons_x_ix, chrom_cons_y_ix
  integer, allocatable :: seed_arr(:)

  real(rp) metric
  real(rp) chrom_x, chrom_y
  real(rp) tr_x, tr_y
  real(rp) delta, linear_vec, da_vec
  real(rp) nu_x, nu_y, dist_to_coup
  real(rp) nux_cons, nuy_cons
  real(rp) str_cons
  real(rp) min_dist_to_coup
  real(rp) delta_e
  real(rp) c_x, c_y, ccc

  logical linear_ok
  logical feasible
  logical err_flag, mat_err
  logical first_loop
  logical mat_ok
  logical good_wp
  logical sane_chrom

  character*60 in_file
  character*100 thin_file
  character*18 var_str
  character*30 set_str
  integer iostat
  type(mag_struct), pointer :: var_mags(:)

  real(rp), allocatable :: etax_base(:)
  real(rp), allocatable :: etay_base(:)
  real(rp) co_screen, co_screen_x, co_screen_y
  integer status
  
  !pisa vars
  character(20) prefix
  character(10) poll_str
  real poll
  integer polli
  integer gen_num
  integer alpha
  integer curname
  integer lambda
  integer mu
  integer sta
  integer nsel, narc
  integer ix
  integer dims
  integer con_set, con_read
  logical dead
  integer pool_gap

  !pisa statistics
  integer stats_surviving
  integer stats_feasible

  !moga vars
  real(rp), allocatable :: K2(:)
  type(pool_struct), allocatable :: pool(:)
  type(pop_struct), allocatable :: pop(:)
  integer, allocatable :: arc(:)
  integer, allocatable :: last_arc(:)
  integer, allocatable :: sel(:)
  real(rp) omega_bound_lir
  real(rp) omega_bound_uir
  real(rp) r

  !mpi housekeeping
  integer myrank, from_id
  integer n_slave, cluster_size
  integer mpierr
  integer mpistatus(MPI_STATUS_SIZE)
  logical master

  real(rp), allocatable :: vars(:)
  real(rp), allocatable :: objs(:)
  real(rp), allocatable :: cons(:)
  integer worker_id, worker_status
  integer n_vars, n_loc, n_dep
  integer name, pool_ptr, pool_ptr_b, lambda_recv
  integer n_recv, slot_num
  integer, allocatable :: lambda_vec(:)

  !chromaticity setting
  real(rp) A(2,2), Ainv(2,2)
  real(rp) chrom_vars(2), chrom_vec(2)
  real(rp) chrom_rep_x, chrom_rep_y, init_chrom_x, init_chrom_y

  logical rf_on
  logical ok, err
  logical slaves_done

  ! reduce number of error messages
  call output_direct(print_and_capture=.false., min_level=s_blank$, max_level=s_error$)

  ! read command line arguments
  call getarg(1,in_file)
  call getarg(2,prefix)
  call getarg(3,poll_str)
  read(poll_str,*) poll
  polli = floor(poll*1000)

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
    pool_gap = cluster_size-1
    if(n_slave .eq. 0) then
      write(*,*) "ERROR: no slaves found in cluster.  At least two nodes"
      write(*,*) "must be available to run this program."
      call mpi_finalize(mpierr)
      error stop
    endif
    !Clear PISA state file
    call write_state(prefix,0)
  endif

  !Set default for parameters not likely to be set if user is not optimizing linear optics.
  work_pt_x_min = 0.0d0
  work_pt_x_max = 999.9d0
  work_pt_y_min = 0.0d0
  work_pt_y_max = 999.9d0

  eta_x_abs_max = 999.9d0

  coupling_con = .false.

  ! parse parameters file and check for necessary initializations.
  use_hybrid = .false.  !default
  generate_feasible_seeds_only = -1
  call set_params_to_bomb()
  open (unit = 10, file = in_file, action='read')
  read (10, nml = general)
  read (10, nml = da)
  read (10, nml = nl_moga)
  close (10)
  if(master) call check_params_bomb()

  ! process parameters
  i=0
  n_vars = 0
  do while(.true.)
    i=i+1
    if( mags_in(i)%name == '' ) exit
    n_vars = n_vars + 1
  enddo
  allocate(var_mags(n_vars))

  if(use_chrom_mags) then
    n_dep = 2 !two dependent variables for chromaticity
  else
    n_dep = 0
  endif

  !- this section of code assumes that the magnet types in mags_in are ordered.
  var_mags(1:n_vars) = mags_in(1:n_vars)

  do i=1,max_de
    if( de(i) .lt. -998. ) then
      exit
    endif
  enddo
  n_de = i-1

  ! parse shared config params
  call pisa_cfg_parser(prefix, alpha, mu, lambda, dims, con_read)
  con_set = 11 !add beta_x, beta_y constraint
  if(coupling_con) then
    con_set = con_set + 1
    coupling_con_ix = con_set
  endif
  if(.not. use_chrom_mags) then
    con_set = con_set + 1
    chrom_cons_x_ix = con_set
    con_set = con_set + 1
    chrom_cons_y_ix = con_set
  endif

  ! check shared parameters meet program limitations
  if( con_read /= con_set ) then
    write(*,*) "PISA_cfg constraints must be ", con_set, ".  Terminating."
    call mpi_finalize(mpierr)
    error stop
  endif
  if ( mu .ne. lambda ) then
    write(*,*) "this program can handle only mu == lambda.", mu, lambda
    call mpi_finalize(mpierr)
    error stop
  endif
  if ( mod(mu,2) .ne. 0 ) then
    write(*,*) "this program can handle only even mu."
    call mpi_finalize(mpierr)
    error stop
  endif
  if ( dims .eq. 3) then
    write(*,*) "DA-objectives only"
  elseif( dims .eq. 4) then
    write(*,*) "DA and emittance objectives"
  endif


  if( master ) write(*,*) "preparing lattice..."


  bp_com%always_parse = .true.
  call bmad_parser(lat_file,ring)

  allocate(co(0:ring%n_ele_track))
  allocate(co_on_energy(0:ring%n_ele_track))
  allocate(orb(0:ring%n_ele_track))

  call set_on_off(rfcavity$, ring, off$)
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.
  bmad_com%aperture_limit_on = .true.
  rf_on = .false.

  do i=1,ring%n_ele_track
    if(ring%ele(i)%key == wiggler$) then
      ring%ele(i)%value(x1_limit$) = 1.0
      ring%ele(i)%value(x2_limit$) = 1.0
      ring%ele(i)%value(y1_limit$) = 1.0
      ring%ele(i)%value(y2_limit$) = 1.0
      ring%ele(i)%aperture_type = elliptical$
    endif
  enddo

  !do i=1,ring%n_ele_track
  !  if( ring%ele(i)%value(x1_limit$) .lt. 1e-4 ) then
  !    ring%ele(i)%value(x1_limit$) = 0.05 
  !    ring%ele(i)%value(x2_limit$) = 0.05 
  !    ring%ele(i)%value(y1_limit$) = 0.05 
  !    ring%ele(i)%value(y2_limit$) = 0.05 
  !  endif
  !enddo

  n_aperture_test = 0
  do i=1,ring%n_ele_track
    if( ring%ele(i)%value(x1_limit$) .gt. 1e-4 ) n_aperture_test = n_aperture_test + 1
  enddo
  if( (1.0*n_aperture_test)/ring%n_ele_track .lt. 0.5 ) then
    write(*,*) "Less than half the elements do not have a x1 physical aperture."
    write(*,*) "Probably something is wrong.  Check that lattice file defines aperture."
    write(*,*) "Aborting"
    call mpi_finalize(mpierr)
    error stop
  endif

  if (tracking_method .gt. 0) then
    do i = 1, ring%n_ele_max
      if(ring%ele(i)%key .ne. wiggler$) then
        ring%ele(i)%tracking_method = tracking_method
      endif
    enddo
  endif
  call twiss_and_track(ring,co,status)
  if(status .ne. ok$) then
    write(*,'(a,i6)') "FATAL: Optics calculation for stock lattice failed. status returned is ", status
    call mpi_finalize(mpierr)
    error stop
  endif
    
  !+ Allocate space in master for storing gene pool
  if( master ) then
    ! allocate memory for storing vector pool
    allocate(pool(alpha+pool_gap))
    pool(:)%name = -1
    do i=1,size(pool)
      allocate(pool(i)%x(n_vars))
    enddo
  endif

  ! mutator is simply copied over
  do i=1,n_vars
    breeder_params%mutate_delta(i) = var_mags(i)%mutate_delta
  enddo

  !-
  allocate(etax_base(1:ring%n_ele_track))
  allocate(etay_base(1:ring%n_ele_track))

  !- set da configuration structures
  da_config%param%n_turn = n_turn
  da_config%param%accuracy = 0.00001d0
  da_config%param%init_len = init_len
  da_config%param%step_len = 0.0005d0
  da_config%min_angle = 0.2d0
  da_config%max_angle = pi-0.2d0
  da_config%n_angle = n_angle
  da_config%param%adts_x_min = adts_x_min
  da_config%param%adts_x_max = adts_x_max
  da_config%param%adts_y_min = adts_y_min
  da_config%param%adts_y_max = adts_y_max
  allocate(da_config%aperture(1:n_angle))

  da_block_linear%min_angle = da_config%min_angle
  da_block_linear%max_angle = da_config%max_angle
  da_block_linear%n_angle = da_config%n_angle
  allocate(da_block_linear%aperture(1:da_block_linear%n_angle))

  !-
  allocate(vars(n_vars))
  allocate(objs(dims))
  allocate(cons(con_set))  !zero length arrays are ok in fortran
  allocate(lambda_vec(mu))

  ! alpha is population size
  ! mu is number of parents picked out by selector
  ! lambda is number of children created by breeder

  if( master ) then
    ! manager
    write(*,*) "starting simulation..."
    if(seed .gt. 0) then
      call random_seed(size=rsize)
      allocate(seed_arr(rsize))
      do i=1,rsize
        seed_arr(i) = seed*i
      enddo 
      call random_seed(put=seed_arr)
      deallocate(seed_arr)
    else
      call random_seed()
    endif

    ! allocate memory for storing population
    allocate(pop(alpha+lambda))
    pop(:)%name = -1
    do i=1,size(pop)
      allocate(pop(i)%o(dims))
      allocate(pop(i)%x(n_vars))
      allocate(pop(i)%x_dep(n_dep))
      allocate(pop(i)%c(con_set))
    enddo
    allocate(arc(alpha))
    allocate(last_arc(alpha))
    last_arc = 0
    allocate(sel(lambda))

    ! generate or read initial population
    if ( trim(initial_pop) == 'random' ) then
      do i=1,size(pool)
        ! bounds on magnets are simply copied over
        do j=1, n_vars 
          call random_number(r)
          pool(i)%x(j) = r*(var_mags(j)%uir-var_mags(j)%lir) + var_mags(j)%lir
        enddo
        pool(i)%name = i
      enddo
    elseif ( trim(initial_pop) .ne. '' ) then
      call read_initial_population(pool, n_dep, alpha, n_vars, initial_pop, err_flag)
      if(err_flag) then
        write(*,*) "Error reading initial population.  aborting."
        call mpi_finalize(mpierr)
        stop
      endif
      do i=1,pool_gap
        pool(alpha+i)%x = 0.0d0
        pool(alpha+i)%name = pool(i)%name + alpha
      enddo
    endif
    pool_ptr_b = alpha+pool_gap

    ! seed each worker with a trial vector
    write(*,*) "seeding generation 1"
    pool_ptr = 0
    do worker_id=1, min(n_slave,alpha)
      call increment_ptr(pool_ptr,size(pool))
      call mpi_send(pool(pool_ptr)%name, 1, MPI_INTEGER, worker_id, 1, MPI_COMM_WORLD, mpierr)
      call mpi_send(pool(pool_ptr)%x(:), n_vars, MPI_DOUBLE_PRECISION, worker_id, 2, MPI_COMM_WORLD, mpierr)
    enddo

    ! receive objectives from workers
    ! refresh worker with new trial vector
    n_recv = 0
    do while (n_recv .lt. alpha)
      call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
      from_id = mpistatus(MPI_SOURCE)
      call mpi_recv(name, 1, MPI_INTEGER, from_id, 3, MPI_COMM_WORLD, mpistatus, mpierr)
      call find_empty_pop_slot(pop(:),slot_num)
      pop(slot_num)%name = name
      call mpi_recv(pop(slot_num)%x(:), n_vars, MPI_DOUBLE_PRECISION, from_id, 4, MPI_COMM_WORLD, mpistatus, mpierr)
      call mpi_recv(pop(slot_num)%x_dep(:), n_dep, MPI_DOUBLE_PRECISION, from_id, 5, MPI_COMM_WORLD, mpistatus, mpierr)
      call mpi_recv(pop(slot_num)%o(:), dims, MPI_DOUBLE_PRECISION, from_id, 6, MPI_COMM_WORLD, mpistatus, mpierr)
      call mpi_recv(pop(slot_num)%c(:), con_set, MPI_DOUBLE_PRECISION, from_id, 7, MPI_COMM_WORLD, mpistatus, mpierr)

      n_recv = n_recv + 1

      !refresh worker with new trial vector
      call increment_ptr(pool_ptr,size(pool))
      call mpi_send(pool(pool_ptr)%name, 1, MPI_INTEGER, from_id, 1, MPI_COMM_WORLD, mpierr)
      call mpi_send(pool(pool_ptr)%x(:), n_vars, MPI_DOUBLE_PRECISION, from_id, 2, MPI_COMM_WORLD, mpierr)
    enddo

    call write_pop_pisa(pop(1:alpha),trim(prefix)//'ini')
    call write_state(prefix,1)
    call date_and_time(values=time_stamp)
    write(*,'(a,i4,a,3i5)') "generation ", 1, " complete at ", time_stamp(5:7)

    !delete old output files if present
    call file_suffixer(moga_output_file,thin_file,'.thin',.true.)
    open(unit=21, iostat=iostat, file=moga_output_file, status='old')
    if (iostat .eq. 0) then
      write(*,*) "deleting old output file."
      close(21, status='delete')
    endif
    close(21)
    open(unit=21, iostat=iostat, file=thin_file, status='old')
    if (iostat .eq. 0) then
      write(*,*) "deleting old thin output file."
      close(21, status='delete')
    endif
    close(21)
    open(unit=22, iostat=iostat, file='constraint_report.out', status='old')
    if (iostat .eq. 0) then
      write(*,*) "deleting old constraint file."
      close(22, status='delete')
    endif
    close(22)
    open(unit=22, iostat=iostat, file='objective_report.out', status='old')
    if (iostat .eq. 0) then
      write(*,*) "deleting old objective report file."
      close(22, status='delete')
    endif
    close(22)
    open(unit=22, iostat=iostat, file='offspring_report.out', status='old')
    if (iostat .eq. 0) then
      write(*,*) "deleting old offspring_report.out file."
      close(22, status='delete')
    endif
    close(22)

    !make new output files, write header
    open(unit=21, iostat=iostat, file=moga_output_file, access='append')
    write(21,'(a8,50a19)') '# id', (trim(chrom_mags(i)%name)//'['//trim(chrom_mags(i)%property)//']',i=1,2), &
                                   (trim(var_mags(i)%name)//'['//trim(var_mags(i)%property)//']',i=1,n_vars), "o1", "o2", "o3", "o4", "feasible"
    close(21)
    open(unit=21, iostat=iostat, file=thin_file, access='append')
    write(21,'(a8,50a11)') '# id', (trim(chrom_mags(i)%name)//'['//trim(chrom_mags(i)%property)//']',i=1,2), &
                                   (trim(var_mags(i)%name)//'['//trim(var_mags(i)%property)//']',i=1,n_vars), "o1", "o2", "o3", "o4", "feasible"
    close(21)
    open(unit=22, iostat=iostat, file='constraint_report.out', access='append')
    write(22,'(a6,30a18)') '# id', 'max|k2|', 'mats@+de', 'mats@-de', 'co@-de', 'co@+de', 'wk_pt_x', 'wk_pt_y', 'beta_x', 'beta_y', 'eta_x'
    close(22)
    open(unit=22, iostat=iostat, file='objective_report.out', access='append')
    write(22,'(a6,10a18)') '# id', 'o1', 'o2', 'o3', 'o4'
    close(22)
    open(unit=22, iostat=iostat, file='offspring_report.out', access='append')
    write(22,'(a6,30a18)') '# gen', '% surviving', '% feasible'
    close(22)
    if( generate_feasible_seeds_only .gt. 0 ) then
      open(44,file='feasible.log')
    endif

    do gen_num = 2, max_gen
      call block_on_pisa_status(polli,prefix)
      call read_pisa_indexes(prefix,'sel',nsel,sel)
      last_arc = arc
      call read_pisa_indexes(prefix,'arc',narc,arc)
      call delete_the_dead(pop(:),arc,narc)
      call write_population(pop, generate_feasible_seeds_only, gen_num-1, moga_output_file)
      call write_population(pop, generate_feasible_seeds_only, gen_num-1, thin_file, prec=3)
      call write_constraint_report(pop(:),gen_num-1,'constraint_report.out')
      call write_objective_report(pop(:),gen_num-1,'objective_report.out')
      if( generate_feasible_seeds_only .gt. 0 ) then
        call count_feasible_in_pop(pop, n_feasible)
        write(*,'(a,i6,a)') "population contains ", n_feasible, " feasible seeds."
        write(44,'(a,i6,a)') "population contains ", n_feasible, " feasible seeds."
        if(n_feasible .ge. generate_feasible_seeds_only) then
          exit
          !call write_state(prefix,5) !tell pisa selector to shut down
          !call mpi_finalize(mpierr)
          !stop
        endif
      endif

      call kangal_breeder(pop(:), sel, pool, pool_ptr_b, breeder_params)

      ! receive objectives from workers
      ! refresh worker with new trial vector
      n_recv = 0
      stats_feasible = 0
      do while (n_recv .lt. mu)
        call mpi_probe(MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
        from_id = mpistatus(MPI_SOURCE)
        call find_empty_pop_slot(pop(:),slot_num)
        call mpi_recv(pop(slot_num)%name, 1, MPI_INTEGER, from_id, 3, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(pop(slot_num)%x(:), n_vars, MPI_DOUBLE_PRECISION, from_id, 4, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(pop(slot_num)%x_dep(:), n_dep, MPI_DOUBLE_PRECISION, from_id, 5, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(pop(slot_num)%o(:), dims, MPI_DOUBLE_PRECISION, from_id, 6, MPI_COMM_WORLD, mpistatus, mpierr)
        call mpi_recv(pop(slot_num)%c(:), con_set, MPI_DOUBLE_PRECISION, from_id, 7, MPI_COMM_WORLD, mpistatus, mpierr)


        feasible = all( pop(slot_num)%c(:) .ge. 0.0d0 )
        if(feasible) then
          stats_feasible = stats_feasible + 1
        endif
        n_recv = n_recv + 1
        lambda_vec(n_recv) = slot_num

        !refresh worker with new trial vector
        call increment_ptr(pool_ptr,size(pool))
        call mpi_send(pool(pool_ptr)%name, 1, MPI_INTEGER, from_id, 1, MPI_COMM_WORLD, mpierr)
        call mpi_send(pool(pool_ptr)%x(:), n_vars, MPI_DOUBLE_PRECISION, from_id, 2, MPI_COMM_WORLD, mpierr)
      enddo


      stats_surviving = alpha
      do i=1,alpha
        do j=1,alpha
          if ( last_arc(i) .eq. arc(j) ) then
            stats_surviving = stats_surviving - 1 
            exit
          endif
        enddo
      enddo
      open(unit=22, iostat=iostat, file='offspring_report.out', access='append')
      write(22,'(i6,2f18.3)') gen_num, (100.0*stats_surviving)/mu, (100.0*stats_feasible)/mu
      close(22)

      call write_pop_pisa(pop(lambda_vec),trim(prefix)//'var')
      call write_state(prefix,3)
      call date_and_time(values=time_stamp)
      write(*,*) "************************************************************"
      write(*,*) "*"
      write(*,*) "*"
      write(*,*) "*"
      write(*,'(a,i4,a,3i5)') "generation ", gen_num, " complete at ", time_stamp(5:7)
      write(*,*) "*"
      write(*,*) "*"
      write(*,*) "*"
      write(*,*) "************************************************************"
    enddo

    !read final population archive from selector
    call block_on_pisa_status(polli,prefix)
    call read_pisa_indexes(prefix,'arc',narc,arc)
    call delete_the_dead(pop(:),arc,narc)
    call write_population(pop, generate_feasible_seeds_only, max_gen, moga_output_file) !write final population to log file
    call write_population(pop, generate_feasible_seeds_only, max_gen, thin_file, prec=3) !write final population to log file
    call write_constraint_report(pop(:),max_gen,'constraint_report.out')
    call write_objective_report(pop(:),max_gen,'objective_report.out')
    call write_state(prefix,5) !tell pisa selector to shut down

    !tell workers to shut down
    do i=1,n_slave
      call mpi_send(0, 1, MPI_INTEGER, i, 1, MPI_COMM_WORLD, mpierr)
    enddo
    close(23)
    call mpi_finalize(mpierr)
  else
    ! worker
    if(use_hybrid) then
      do i=1,ring%n_ele_track
        if( (ring%ele(i)%key == sbend$) .or. &
            (ring%ele(i)%key == sextupole$) .or. &
            (ring%ele(i)%key == rfcavity$) .or. &
            (ring%ele(i)%key == wiggler$) .or. &
            (ring%ele(i)%key == multipole$) ) then
          ring%ele(i)%select = .true.
        else
          ring%ele(i)%select = .false.
        endif
      enddo
    endif

    ring0 = ring  !stash original, unaltered lattice
    first_loop = .true.
    do while(.true.)
      if( .not. first_loop ) then
        write(*,'(a,i5,a,i8)') "worker ", myrank, " processed name ", name
        call mpi_send(name, 1, MPI_INTEGER, 0, 3, MPI_COMM_WORLD, mpierr)
        call mpi_send(vars, n_vars, MPI_DOUBLE_PRECISION, 0, 4, MPI_COMM_WORLD, mpierr)
        call mpi_send(chrom_vars, n_dep, MPI_DOUBLE_PRECISION, 0, 5, MPI_COMM_WORLD, mpierr)
        call mpi_send(objs, dims, MPI_DOUBLE_PRECISION, 0, 6, MPI_COMM_WORLD, mpierr)
        call mpi_send(cons, con_set, MPI_DOUBLE_PRECISION, 0, 7, MPI_COMM_WORLD, mpierr)
      else
        first_loop = .false.
      endif

      ring = ring0

      ! receive magnet strengths from master
      call mpi_recv(name, 1, MPI_INTEGER, 0, 1, MPI_COMM_WORLD, mpistatus, mpierr)
      if( name .eq. 0 ) then
        call mpi_finalize(mpierr)
        exit
      endif
      call mpi_recv(vars, n_vars, MPI_DOUBLE_PRECISION, 0, 2, MPI_COMM_WORLD, mpistatus, mpierr)

      objs(:) = 1.0d0

      ! Variable bound constraints
      cons(1)  = -100.0    !magnet strengths
      ! Constraint defaults below remain if on-energy lattice us unstable.
      ! Off-energy constraints
      cons(2)  = -200.0    ! 1-turn mats, traces, tunes along fp_dE_neg
      cons(3)  = -200.0    ! 1-turn mats, traces, tunes along fp_dE_neg
      cons(4)  = -100.0    ! nonlinear dispersion at de(1)
      cons(5)  = -100.0    ! nonlinear dispersion at de(2)
      ! On-energy constraints
      cons(6)  = -100.0    ! x working point bounds
      cons(7)  = -100.0    ! y working point bounds
      cons(8)  = -100.0   ! beta_x at element 0
      cons(9)  = -100.0   ! beta_y at element 0
      cons(10)  = -100.0   ! eta_x at element 0
      cons(11)  = -100.0   ! global beta_x max
      if(coupling_con) then
        cons(coupling_con_ix) = -100.0  !nu_x - nu_y
      endif
      if(.not. use_chrom_mags) then
        cons(chrom_cons_x_ix) = -90000.0
        cons(chrom_cons_y_ix) = -90000.0
      endif
      sane_chrom = .false.

      !- calculate magnet strength constraint
      str_cons = 0.00000001d0
      do i=1, n_vars
        if(vars(i) .lt. var_mags(i)%lb) then
          str_cons = str_cons + (vars(i) - var_mags(i)%lb)/(abs(var_mags(i)%ub)+abs(var_mags(i)%lb))
        elseif(vars(i) .gt. var_mags(i)%ub) then
          str_cons = str_cons + (var_mags(i)%ub - vars(i))/(abs(var_mags(i)%ub)+abs(var_mags(i)%lb))
        endif
      enddo
      cons(1) = str_cons

      call set_magnet_strengths(var_mags,ring,vars(1:n_vars))
      if(use_chrom_mags) call set_magnet_strengths(chrom_mags,ring,(/0.0d0,0.0d0/))
      call clear_lat_1turn_mats(ring)
      co(0)%vec(:) = 0.0d0
      call twiss_and_track(ring,co,status)
      co_on_energy = co  !stash for later emittance calculation
      ring_on_energy = ring  !stash for later emittance calculation
      chrom_vars = 0.0d0
      
      !- Calculate working point tune constraints
      tr_x = ring%param%t1_no_RF(1,1)+ring%param%t1_no_RF(2,2)
      nu_x = ring%ele(ring%n_ele_track)%a%phi/twopi
      cons(6) = trace_tunes_con_a(tr_x, nu_x, work_pt_x_min, work_pt_x_max)
      tr_y = ring%param%t1_no_RF(3,3)+ring%param%t1_no_RF(4,4)
      nu_y = ring%ele(ring%n_ele_track)%b%phi/twopi
      cons(7) = trace_tunes_con_a(tr_y, nu_y, work_pt_y_min, work_pt_y_max)

      if(coupling_con) then
        cons(coupling_con_ix) = -100.0
        if( abs(tr_x) .lt. 2.0d0 ) then
          if( abs(tr_y) .lt. 2.0d0 ) then
            tune_delta = abs( (nu_x-int(nu_x)) - (nu_y-int(nu_y)) )
            cons(coupling_con_ix) = -(tune_delta-coupling_delta_max)
          endif
        endif
      endif

      if( status == ok$ ) then
        cons(8) = beta_x_max - ring%ele(0)%a%beta
        cons(9) = beta_y_max - ring%ele(0)%b%beta
        cons(10) = eta_x_abs_max - abs(ring%ele(0)%a%eta)
        cons(11) = global_beta_x_max - maxval(ring%ele(0:ring%n_ele_track)%a%beta)
        !- Stash dispersion for nonlinear dispersion calculation
        etax_base(:) = ring%ele(1:ring%n_ele_track)%a%eta
        etay_base(:) = ring%ele(1:ring%n_ele_track)%b%eta

        call chrom_calc(ring, 1.0d-6, init_chrom_x, init_chrom_y, err_flag)
        if(use_chrom_mags) then
          !- set chromaticity
          call set_magnet_strengths(chrom_mags,ring,(/1.0d0,0.0d0/))
          call twiss_and_track(ring,co,status)
          call chrom_calc(ring, 1.0d-6, chrom_rep_x, chrom_rep_y, err_flag)
          A(1,1) = (chrom_rep_x-init_chrom_x)/1.0d0
          A(2,1) = (chrom_rep_y-init_chrom_y)/1.0d0

          call set_magnet_strengths(chrom_mags,ring,(/0.0d0,1.0d0/))
          call twiss_and_track(ring,co,status)
          call chrom_calc(ring, 1.0d-6, chrom_rep_x, chrom_rep_y, err_flag)
          A(1,2) = (chrom_rep_x-init_chrom_x)/1.0d0
          A(2,2) = (chrom_rep_y-init_chrom_y)/1.0d0

          call mat_inverse(A,Ainv,ok,.true.)
          if(.not. ok) then
            write(*,*) "Chromaticity response matrix inversion failed."
          endif

          chrom_vec(1) = set_chrom_x - init_chrom_x
          chrom_vec(2) = set_chrom_y - init_chrom_y
          chrom_vars = matmul(Ainv,chrom_vec)

          !Update magnet strength constraint for the two chrom_mags
          do i=1, 2
            if(chrom_vars(i) .lt. chrom_mags(i)%lb) then
              cons(1) = cons(1) + (chrom_vars(i) - chrom_mags(i)%lb)/(abs(chrom_mags(i)%ub)+abs(chrom_mags(i)%lb))
            elseif(chrom_vars(i) .gt. chrom_mags(i)%ub) then
              cons(1) = cons(1) + (chrom_mags(i)%ub - chrom_vars(i))/(abs(chrom_mags(i)%ub)+abs(chrom_mags(i)%lb))
            endif
          enddo

          if( (abs(chrom_vars(1)) .lt. 1.0e5) .and. (abs(chrom_vars(2)) .lt. 1.0e5)) then
            sane_chrom = .true.
            call set_magnet_strengths(chrom_mags,ring,chrom_vars)
            call clear_lat_1turn_mats(ring)
            co(0)%vec(:) = 0.0d0
            call twiss_and_track(ring,co,status)
            ring_on_energy = ring  !update on_energy ring for later emittance calculation
            co_on_energy = co
          endif
        else
          if(.not. err_flag) then
            cons(chrom_cons_x_ix) = 1.0d0
            if(init_chrom_x .lt. chrom_x_min) then
              cons(chrom_cons_x_ix) = init_chrom_x - chrom_x_min
            elseif(init_chrom_x .gt. chrom_x_max) then
              cons(chrom_cons_x_ix) = chrom_x_max - init_chrom_x
            endif
            cons(chrom_cons_y_ix) = 1.0d0
            if(init_chrom_y .lt. chrom_y_min) then
              cons(chrom_cons_y_ix) = init_chrom_y - chrom_y_min
            elseif(init_chrom_y .gt. chrom_y_max) then
              cons(chrom_cons_y_ix) = chrom_y_max - init_chrom_y
            endif
            if( cons(chrom_cons_x_ix) > -2.0 .and. cons(chrom_cons_y_ix) > -2.0 ) then
              !chrom is sane if within 2 of the boundaries
              sane_chrom = .true.
            else
              sane_chrom = .false.
            endif
          else
            sane_chrom = .false.
          endif
        endif
      endif

      if( sane_chrom ) then
        !calculate closed orbit amplitude constraints.
        do i=2,3  !assume negative and positive dE at de(2) and de(3)
          co(0)%vec = 0.0d0
          co(0)%vec(6) = de(i)
          call clear_lat_1turn_mats(ring)
          call twiss_and_track(ring,co,status)
          if(status == ok$) then
            co_screen_x = abs(co(1)%vec(1)-etax_base(1)*de(i))
            co_screen_y = abs(co(1)%vec(3)-etay_base(1)*de(i))
            do k=2,ring%n_ele_track
              if( .not. any(ring%ele(k)%key == (/wiggler$, marker$/)) ) then
                co_screen_x = max(co_screen_x,abs(co(k)%vec(1)-etax_base(k)*de(i)))
                co_screen_y = max(co_screen_y,abs(co(k)%vec(3)-etay_base(k)*de(i)))
              endif
            enddo
            co_screen = (co_limit - max(co_screen_x,co_screen_y)) / co_limit
            if(i==2) then
              cons(4) = co_screen
            elseif(i==3) then
              cons(5) = co_screen
            endif
          endif
        enddo

        do i=1,2  ! i==1: fp_dE_neg.  i==2: fp_dE_pos
          call clear_lat_1turn_mats(ring)
          ccc = 0.0d0
          do j=1,n_fp_steps
            !calculate closed orbit and optics at energy offset
            co(0)%vec = 0.0d0
            if(i==1) then !negative, fp_de_neg assumed to be a negative number
              delta_e = (fp_de_neg/n_fp_steps)*j
            elseif(i==2) then !positive
              delta_e = (fp_de_pos/n_fp_steps)*j
            endif
            co(0)%vec(:) = 0.0d0
            co(0)%vec(6) = delta_e
            call twiss_and_track(ring,co,status)

            !calculate mat OK constraints.
            !i.e. how far out along fp_dE_neg and fp_dE_pos can the energy go before the 
            !1-turn mat becomes unstable.
            if( any(status==(/no_closed_orbit$, xfer_mat_calc_failure$/)) .or. ( (status /= ok$).and.(status .gt. 0) ) ) then
              ccc = ccc + (-200.0d0)
            else
              ! (c_x + c_y)/2.0 saturates at -120.0 for big x and y traces 
              tr_x = ring%param%t1_no_RF(1,1)+ring%param%t1_no_RF(2,2)
              nu_x = ring%ele(ring%n_ele_track)%a%phi/twopi
              c_x = trace_tunes_con_b(tr_x, nu_x, x_fp_min, x_fp_max)
              tr_y = ring%param%t1_no_RF(3,3)+ring%param%t1_no_RF(4,4)
              nu_y = ring%ele(ring%n_ele_track)%b%phi/twopi
              c_y = trace_tunes_con_b(tr_y, nu_y, y_fp_min, y_fp_max)

              ccc = ccc + (c_x + c_y)/2.0d0 
            endif
          enddo
          if(i==1) then !negative chromatic footprint constraint
            cons(2) = ccc / n_fp_steps
          elseif(i==2) then !positive chromatic footprint constraint
            cons(3) = ccc / n_fp_steps
          endif
        enddo

        ! feasible if all constraints met, otherwise infeasible
        feasible = all( cons(:) .ge. 0.0d0 )

        ! tracking study begins here
        ! if lattice is infeasible, then the optimizer ignores the objectives, so no point in tracking if infeasible.
        if( feasible .and. (generate_feasible_seeds_only .lt. 0) ) then
          if( dims .eq. 4) then
            ix_cache = -1
            call radiation_integrals(ring_on_energy, co_on_energy, mode, ix_cache)
            objs(4) = (mode%a%emittance - emittance_target) / emittance_target
          endif

          do i=1,n_de
            ! calculate linear aperture (dynamic aperture assuming linear optics)
            co(0)%vec = 0.0d0
            co(0)%vec(6) = de(i)
            call clear_lat_1turn_mats(ring)
            call twiss_and_track(ring,co,status)

            if(status == ok$) then
              da_block_linear%param%closed_orbit = co(0)
              if(use_hybrid) then
                call make_hybrid_lat(ring, ring_use)
              else
                ring_use = ring
              endif
              call linear_aperture(ring_use,da_block_linear)

              ! screen for small linear dynamic aperture dimension
              linear_ok = .true.
              do j=1,n_angle
                linear_vec = sqrt( (da_block_linear%aperture(j)%x-co(0)%vec(1))**2 + &
                                   (da_block_linear%aperture(j)%y-co(0)%vec(3))**2)
                if(linear_vec .lt. linear_vec_cutoff) then
                  write(*,*) "linear aperture too small: ", linear_vec
                  linear_ok = .false.
                endif
              enddo
            else
              linear_ok = .false.
            endif

            if( linear_ok ) then
              !calculate dynamic aperture

              ! only the on-energy adts is constrained
              da_config%param%n_adts = -1
              if(i==1) da_config%param%n_adts = n_adts

              da_config%param%closed_orbit%vec = co(0)%vec
              call custom_dynamic_aperture_scan(ring_use, da_config)
            endif

            ! total up objective values
            if( .not. linear_ok  ) then
              write(*,'(a,f9.4,a)') "linear lattice is not OK for energy offset ", de(i)*100.0, " %"
              metric = 1.0d0
            else
              metric = 0.0d0
              do j=1,n_angle
                linear_vec = sqrt( (da_block_linear%aperture(j)%x-co(0)%vec(1))**2 + &
                                   (da_block_linear%aperture(j)%y-co(0)%vec(3))**2 )
                da_vec = sqrt( (da_config%aperture(j)%x-co(0)%vec(1))**2 + &
                               (da_config%aperture(j)%y-co(0)%vec(3))**2 )
                delta = (linear_vec - da_vec)/linear_vec/sqrt(1.0d0*n_angle)
                if(delta .lt. 0) then
                  delta = 0.0d0  !no contribution from exceeding physical aperture
                endif
                metric = metric + delta**2
              enddo
            endif

            objs(i) = metric
          enddo
        endif
      endif
    enddo
  endif

  contains

    function trace_tunes_con_a(trace,tune,tune_min,tune_max) result(c)
      implicit none
      real(rp) trace, tune, tune_min, tune_max, c

      if( abs(trace) .gt. 2.0d0 ) then
        !unstable.  con gets log of trace as measure of instability
        c = -log(abs(trace)/2.0d0) / log(1.1d0) - 20.0  ! -inf for big trace
      elseif( tune .gt. tune_max ) then
        c = -20.0 * tanh(tune-tune_max)
      elseif( tune .lt. tune_min ) then
        c =  20.0 * tanh(tune-tune_min)
      else
        c = 0.00001d0
      endif
    end function

    function trace_tunes_con_b(trace,tune,tune_min,tune_max) result(c)
      implicit none
      real(rp) trace, tune, tune_min, tune_max, c

      if( abs(trace) .gt. 2.0d0 ) then
        c = -20.0 * (1.0 + 5.0*tanh((abs(trace)-2.0)/1000.0))  !saturates at -120 for big trace
      elseif( tune .gt. tune_max ) then
        c = -20.0 * tanh(tune-tune_max)
      elseif( tune .lt. tune_min ) then
        c =  20.0 * tanh(tune-tune_min)
      else
        c = 0.00001d0
      endif
    end function

    subroutine set_params_to_bomb()
      implicit none
      de(:) = -999.
      set_chrom_x = -999.0
      set_chrom_y = -999.0
      co_limit = -999.0
      linear_vec_cutoff = -999.0
      x_fp_min = -999.0
      x_fp_max = -999.0
      y_fp_min = -999.0
      y_fp_max = -999.0
      adts_x_min = -999.0
      adts_x_max = -999.0
      adts_y_min = -999.0
      adts_y_max = -999.0
      init_len = -999.0
      fp_de_neg = -999.0
      fp_de_pos = -999.0
      n_fp_steps = -999
      seed = -999
      max_gen = -999
      n_turn = -999
      n_angle = -999
      tracking_method = -999
      track_dims = -999
      lat_file = 'bomb'
      moga_output_file = 'bomb'
      initial_pop = 'bomb'
      mags_in(:)%name = ''
    end subroutine

    subroutine check_params_bomb()
      implicit none
      logical fail
      fail = .false.
      fail = fail .or. check_bomb(set_chrom_x,'set_chrom_x')
      fail = fail .or. check_bomb(set_chrom_y,'set_chrom_y')
      fail = fail .or. check_bomb(x_fp_min,'x_fp_min')
      fail = fail .or. check_bomb(x_fp_max,'x_fp_max')
      fail = fail .or. check_bomb(y_fp_min,'y_fp_min')
      fail = fail .or. check_bomb(y_fp_max,'y_fp_max')
      fail = fail .or. check_bomb(adts_x_min,'adts_x_min')
      fail = fail .or. check_bomb(adts_x_max,'adts_x_max')
      fail = fail .or. check_bomb(adts_y_min,'adts_y_min')
      fail = fail .or. check_bomb(adts_y_max,'adts_y_max')
      fail = fail .or. check_bomb(co_limit,'co_limit')
      fail = fail .or. check_bomb(linear_vec_cutoff,'linear_vec_cutoff')
      fail = fail .or. check_bomb(fp_de_neg,'fp_de_neg')
      fail = fail .or. check_bomb(fp_de_pos,'fp_de_pos')
      fail = fail .or. check_bomb(n_fp_steps,'n_fp_steps')
      fail = fail .or. check_bomb(de(1),'de(1)')
      fail = fail .or. check_bomb(init_len,'init_len')
      fail = fail .or. check_bomb(max_gen,'max_gen')
      fail = fail .or. check_bomb(track_dims,'track_dims')
      fail = fail .or. check_bomb(n_turn,'n_turn')
      fail = fail .or. check_bomb(n_angle,'n_angle')
      fail = fail .or. check_bomb(seed,'seed')
      fail = fail .or. check_bomb(tracking_method,'tracking_method')
      fail = fail .or. check_bomb(lat_file,'lat_file')
      fail = fail .or. check_bomb(moga_output_file,'moga_output_file')
      fail = fail .or. check_bomb(initial_pop,'initial_pop')
      fail = fail .or. check_bomb(mags_in(1)%name,'mags_in(1)%name')
      if(fail) then
        write(*,*) "parameters file does not contain necessary settings."
        write(*,*) "terminating"
        call mpi_finalize(mpierr)
        error stop
      endif
    end subroutine

    function check_bomb(var,var_str)
      implicit none
      class(*) :: var
      character(*) var_str
      logical check_bomb

      check_bomb = .false. 
      select type(var)
        type is (real(rp))
          if(var .lt. -900.0) then
            write(*,*) "parameter ", trim(var_str), " is uninitialized."
            check_bomb = .true. 
          endif
        type is (integer)
          if(var .lt. -900) then
            write(*,*) "parameter ", trim(var_str), " is uninitialized."
            check_bomb = .true.
          endif
        type is (character(*))
          if(trim(var) .eq. 'bomb') then
            write(*,*) "parameter ", trim(var_str), " is uninitialized."
            check_bomb = .true.
          endif
          if(trim(var) .eq. '') then
            write(*,*) "parameter ", trim(var_str), " is uninitialized."
            check_bomb = .true.
          endif
      end select
    end function

end program






