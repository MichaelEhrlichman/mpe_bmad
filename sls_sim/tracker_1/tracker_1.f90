program tracker_1
  use bmad
  use mode3_mod
  use sls_lib
  use naff_mod, except_rp=>rp
 
  implicit none

  type(lat_struct) lat
  type(coord_struct), allocatable :: co(:)
  type(coord_struct), allocatable :: orbit(:)
  type(coord_struct), allocatable :: tbt_coords(:)
  type(coord_struct), allocatable :: norm_coords(:)

  integer i, j, m
  integer track_state
  integer ix_a, ix_b
  integer n_seg

  real(rp), allocatable :: x_tunes(:,:)
  real(rp), allocatable :: y_tunes(:,:)
  real(rp), allocatable :: z_tunes(:,:)
  complex(rp), allocatable :: x_amp(:,:)
  complex(rp), allocatable :: y_amp(:,:)
  complex(rp), allocatable :: z_amp(:,:)
  complex(rp), allocatable :: cdata(:)

  integer, parameter :: max_naff = 3
  complex(rp) x_basis(max_naff,max_naff)
  complex(rp) y_basis(max_naff,max_naff)
  real(rp) x_norm(max_naff)
  real(rp) y_norm(max_naff)

  real(rp) t

  logical err_flag

  character(100) lat_file
  character(100) in_file
  character(100) lat_file_override

  integer iargs

  !parameters from .in
  integer ix_inj
  real(rp) init_vec(6)
  integer n_turns
  integer fft_turns
  integer dims

  namelist / tracker_1_params / lat_file, ix_inj, init_vec, n_turns, fft_turns, dims

  iargs = iargc()
  lat_file_override = ''
  if( iargs == 1 ) then
    call getarg(1,in_file)
  elseif (iargs == 2 ) then
    call getarg(1,in_file)
    call getarg(2,lat_file_override)
  endif

  dims = 4

  open (unit = 10, file = in_file)
  read (10, nml = tracker_1_params)
  close (10)

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  write(*,*) "Preparing lattice..."

  call bmad_parser(lat_file, lat)
  
  if(dims==4) then
    CALL set_on_off(rfcavity$, lat, off$)
  elseif (dims==6) then
    CALL set_on_off(rfcavity$, lat, on$)
  endif
  bmad_com%radiation_damping_on = .true.
  bmad_com%radiation_fluctuations_on = .true.

  call twiss_and_track(lat,co)

  allocate(tbt_coords(n_turns))
  allocate(norm_coords(n_turns))
  allocate(cdata(n_turns))
  allocate(orbit(0:lat%n_ele_track))

  write(*,*) "Tracking..."

  orbit(ix_inj)%vec = co(ix_inj)%vec + init_vec
  write(*,'(a,6es14.5)') "Initial Vector: co(ix_inj) + init_vec = ", orbit(ix_inj)%vec

  tbt_coords(1)%vec(:) = orbit(ix_inj)%vec(:)
  do i = 2, n_turns-1
    call track_many(lat, orbit, ix_inj, ix_inj, 1, track_state=track_state)
    if(track_state /= moving_forward$) then
      write(*,*) "Particle lost.  turn, track_state: ", i,  track_state
      exit
    endif
    tbt_coords(i)%vec(:) = orbit(ix_inj)%vec(:)
  enddo

  open(100,file='tbt.dat')
  write(100,'(A6,4A14)') "# turn", "x", "px", "y", "py"
  do i=1,n_turns
    write(100,'(I6,6ES14.5)') i, tbt_coords(i)%vec(1:6)
  enddo
  close(100)

  if( track_state == moving_forward$ ) then
    do i = 1, n_turns
      call xyz_to_action(lat, ix_inj, tbt_coords(i)%vec-co(ix_inj)%vec, norm_coords(i)%vec, err_flag)
      if (err_flag) then
        write(*,*) "Error from xy_to_action"
        stop
      endif
    enddo

    open(100,file='tbt_Anorm.dat')
    write(100,'(A6,6A14)') "# turn", "Jx", "Jx'", "Jy", "Jy'", "Jz", "Jz'"
    write(100,'(a,es14.6)') "# Jx is ", (norm_coords(1)%vec(1)**2 + norm_coords(1)%vec(2)**2)/2.0d0
    write(100,'(a,es14.6)') "# Jy is ", (norm_coords(1)%vec(3)**2 + norm_coords(1)%vec(4)**2)/2.0d0
    write(100,'(a,es14.6)') "# Jz is ", (norm_coords(1)%vec(5)**2 + norm_coords(1)%vec(6)**2)/2.0d0
    do i=1,n_turns
      write(100,'(I6,6ES14.5)') i, norm_coords(i)%vec(1:6)
    enddo
    close(100)

    !subtract centroid
    norm_coords(:)%vec(1) = norm_coords(:)%vec(1) - sum(norm_coords(:)%vec(1))/size(norm_coords(:)%vec(1))
    norm_coords(:)%vec(2) = norm_coords(:)%vec(2) - sum(norm_coords(:)%vec(2))/size(norm_coords(:)%vec(2))
    norm_coords(:)%vec(3) = norm_coords(:)%vec(3) - sum(norm_coords(:)%vec(3))/size(norm_coords(:)%vec(3))
    norm_coords(:)%vec(4) = norm_coords(:)%vec(4) - sum(norm_coords(:)%vec(4))/size(norm_coords(:)%vec(4))
    norm_coords(:)%vec(5) = norm_coords(:)%vec(5) - sum(norm_coords(:)%vec(5))/size(norm_coords(:)%vec(5))
    norm_coords(:)%vec(6) = norm_coords(:)%vec(6) - sum(norm_coords(:)%vec(6))/size(norm_coords(:)%vec(6))

    write(*,*) "Processing..."

    n_seg = n_turns/fft_turns
    allocate(x_tunes(n_seg,max_naff))
    allocate(y_tunes(n_seg,max_naff))
    allocate(z_tunes(n_seg,max_naff))
    allocate(x_amp(n_seg,max_naff))
    allocate(y_amp(n_seg,max_naff))
    allocate(z_amp(n_seg,max_naff))
    x_tunes = 0.0d0
    y_tunes = 0.0d0
    z_tunes = 0.0d0
    do i=1, n_seg
      ix_a = 1 + (i-1)*fft_turns
      ix_b = i*fft_turns
      x_basis = 0.0d0
      x_norm = 0.0d0
      cdata = cmplx(norm_coords(:)%vec(1),norm_coords(:)%vec(2))
      call naff( cdata(ix_a:ix_b), x_tunes(i,:), x_amp(i,:),300)
      cdata = cmplx(norm_coords(:)%vec(3),norm_coords(:)%vec(4))
      call naff( cdata(ix_a:ix_b), y_tunes(i,:), y_amp(i,:),400)
      cdata = cmplx(norm_coords(:)%vec(5),norm_coords(:)%vec(6))
      call naff( cdata(ix_a:ix_b), z_tunes(i,:), z_amp(i,:),500)
    enddo

    open(100,file='tunes_x_slicingB.dat')
    open(101,file='tunes_y_slicingB.dat')
    open(102,file='tunes_z_slicingB.dat')
    write(100,'(4a18)') "# segment", "tune", "amp", "angle (rad)"
    write(101,'(4a18)') "# segment", "tune", "amp", "angle (rad)"
    write(102,'(4a18)') "# segment", "tune", "amp", "angle (rad)"
    do i=1,n_seg
      do j=1,max_naff
        write(100,'(i18,3es18.8)') i, x_tunes(i,j), abs(x_amp(i,j)), atan2(aimag(x_amp(i,j)),real(x_amp(i,j)))
      enddo
      do j=1,max_naff
        write(101,'(i18,3es18.8)') i, y_tunes(i,j), abs(y_amp(i,j)), atan2(aimag(y_amp(i,j)),real(y_amp(i,j)))
      enddo
      do j=1,max_naff
        write(102,'(i18,3es18.8)') i, z_tunes(i,j), abs(z_amp(i,j)), atan2(aimag(z_amp(i,j)),real(z_amp(i,j)))
      enddo
      write(100,*)
      write(100,*)
      write(101,*)
      write(101,*)
      write(102,*)
      write(102,*)
    enddo
    close(100)
    close(101)
    close(102)

    deallocate(x_tunes)
    deallocate(y_tunes)
    deallocate(z_tunes)
    deallocate(x_amp)
    deallocate(y_amp)
    deallocate(z_amp)
  else
    write(*,*) "Particle was lost.  Not calculating tunes."
  endif

end program






