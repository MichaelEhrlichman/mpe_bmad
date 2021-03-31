program scan_dp
  use bmad
  use sls_lib
  use namelist_general !general: lat_file, use_hybrid

  implicit none

  type (lat_struct) ring
  type (lat_struct) ring_use
  type (coord_struct), allocatable :: co(:)
  type (coord_struct), allocatable :: co_use(:)
  type (coord_struct), allocatable :: orb(:)

  integer i,j,k
  integer dims, tracking_method
  integer track_state
  integer n_turns
  integer ix_inj, n_dp
  integer ix_inj_use
  integer status

  real(rp) dp_min, dp_max, delta_dp
  real(rp) alive

  logical err

  character*100 parameter_file

  namelist / scan_dp_params / dp_min, &
                       dp_max, &
                       n_dp, &
                       ix_inj, &
                       dims, &
                       n_turns, &
                       tracking_method

  ! set defaults
  use_hybrid = .false.
  dims = 4
  tracking_method = -1

  call getarg(1,parameter_file)

  open (unit = 10, file = parameter_file)
  read (10, nml = general)
  read (10, nml = scan_dp_params)
  close (10)

  !force no hybrid
  ! use_hybrid = .false.
    
  call bmad_parser(lat_file, ring)
  bmad_com%aperture_limit_on = .true.
  if(dims==4) then
    call set_on_off(rfcavity$, ring, off$)
    bmad_com%radiation_damping_on = .false.
  elseif(dims==6) then
    use_hybrid = .false.
    call set_on_off(rfcavity$, ring, on$)
    bmad_com%radiation_damping_on = .true.
  endif
  bmad_com%radiation_fluctuations_on = .false.

  if(tracking_method .gt. -1) then
    do i=1,ring%n_ele_track
      if( ring%ele(i)%key .ne. rfcavity$) then
        ring%ele(i)%tracking_method = tracking_method
      endif
    enddo
  endif

  if(use_hybrid) then
    do i=1,ring%n_ele_track
      if( (ring%ele(i)%key == sbend$) .or. &
          (ring%ele(i)%key == sextupole$) .or. &
          (ring%ele(i)%key == rfcavity$) .or. &
          (ring%ele(i)%key == multipole$) ) then
        ring%ele(i)%select = .true.
      else
        ring%ele(i)%select = .false.
      endif
    enddo
  endif

  allocate(orb(0:ring%n_ele_track))

  if(ix_inj_use == 0) then
    write(*,*) "Injection element has been concatenated.  Try again."
    stop
  endif

  ix_inj_use = 0

  open(101,file='scan_dp.dat')
  write(101,'(2a10)') "delta dp", "alive?"
  do j=1,n_dp
    if(n_dp .gt. 1) then
      delta_dp = (dp_max-dp_min)/(n_dp-1)*(j-1) + dp_min
    else
      delta_dp = dp_max
    endif
    co(0)%vec = 0.0d0
    co(0)%vec(6) = delta_dp
    call twiss_and_track(ring,co,status)
    if(use_hybrid) then
      call make_hybrid_lat(ring, ring_use)
    else
      ring_use = ring
    endif
    orb(ix_inj_use) = co(ix_inj)
    orb(ix_inj_use)%vec(6) = orb(ix_inj_use)%vec(6) + delta_dp
    alive = 1.0
    do i=1,n_turns
      call track_many(ring_use,orb,ix_inj_use,ix_inj_use,1,0,track_state)

      if( abs(orb(ix_inj_use)%vec(5)) .gt. c_light/(500.0d6)/2.0 ) then
        write(*,*) "Out of Bucket!"
        track_state = lost$
      endif

      if (track_state .ne. moving_forward$) then
        alive = 0.0
        exit
      endif
    enddo
    write(*,'(a,i4,a,i4,a,f13.9,a,f6.1,a,i6)') "Processed ", j, " of ", n_dp, ". delta_dp = ", delta_dp, " Alive? ", alive, "  ",  i
    write(101,'(f12.5,f6.1)') delta_dp, alive
  enddo
  close(101)
end program





