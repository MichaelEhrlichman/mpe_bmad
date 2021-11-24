program orm
  use bmad
  use bmad_parser_mod, only: bp_com
  use sls_lib
  use analytic_rm_mod
  use f95_lapack

  implicit none

  type(lat_struct) ring
  type(coord_struct), allocatable :: co(:)
  type(normal_modes_struct) mode

  character(100) in_file
  character(200) lat_file
  character(200) lat_file_override

  integer i,j,k

  character(20) h_corrector_mask
  character(20) v_corrector_mask
  character(20) bpm_mask
  type(ele_pointer_struct), allocatable :: h_correctors(:)
  type(ele_pointer_struct), allocatable :: v_correctors(:)
  type(ele_pointer_struct), allocatable :: bpms(:)
  integer n_h_correctors
  integer n_v_correctors
  integer n_bpms
  integer n_loops
  integer status

  real(rp), allocatable :: dx(:), dy(:)
  real(rp), allocatable :: tx(:), ty(:)
  real(rp), allocatable :: a(:,:)
  real(rp), allocatable :: ap(:,:)

  real(rp) I1
  real(rp) rms_x, rms_y

  character(6) ix_str
  character(17) str

  logical err

  namelist /parameters/  lat_file,&
                         n_loops,&
                         h_corrector_mask,&
                         v_corrector_mask,&
                         bpm_mask

  call getarg(1,in_file)
  open(10,file=in_file,status='old')
  read(10,nml=parameters)
  close(10)

  lat_file_override = ''
  call getarg(2,lat_file_override)
  if(lat_file_override .ne. '') then
    lat_file = lat_file_override
  endif

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, ring)
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  !call closed_orbit_calc (ring, co, 6, 1, err_flag = err)
  call closed_orbit_from_tracking(ring, co, 6)
  write(*,'(a,6es12.3)') "Found closed orbit: ", co(0)%vec(1:6)
  call lat_make_mat6 (ring, -1, co, err_flag = err)
  call twiss_at_start (ring, status)
  call twiss_propagate_all (ring, err_flag=err)

  !call twiss_and_track(ring,co,status)
  call radiation_integrals(ring, co, mode)
  I1 = mode%synch_int(1)

  write(*,*) "Beam distribution parameters before correction:"
  write(*,*) "   emit_a      : ", mode%a%emittance
  write(*,*) "   emit_b      : ", mode%b%emittance
  write(*,*) "   sigmaE_E    : ", mode%sigE_E
  write(*,*) "   sigma_z     : ", mode%sig_z

  call lat_ele_locator(h_corrector_mask, ring, h_correctors, n_h_correctors)
  call lat_ele_locator(v_corrector_mask, ring, v_correctors, n_v_correctors)
  call lat_ele_locator(bpm_mask, ring, bpms, n_bpms)

  write(*,'(a,i6)') "Number of horizontal correctors found: ", n_h_correctors
  write(*,'(a,i6)') "Number of vertical correctors found: ", n_v_correctors
  write(*,'(a,i6)') "Number of bpms found: ", n_bpms

  allocate(tx(n_h_correctors))
  allocate(ty(n_v_correctors))
  allocate(dx(n_bpms))
  allocate(dy(n_bpms))

  !Write out uncorrected closed orbit
  open(10,file='uncorrected_orbit.out')
  write(10,'(6a12)') 'ix', 's', 'x', 'px', 'y', 'py'
  do i=1,n_bpms
    write(10,'(i12,f12.4,4es15.4)') i, bpms(i)%ele%s, co(bpms(i)%ele%ix_ele)%vec(1:4)
  enddo
  close(10)


  do j=1,n_loops
    !Make measurement vector
    do i=1,n_bpms
      dx(i) = co(bpms(i)%ele%ix_ele)%vec(1)
      dy(i) = co(bpms(i)%ele%ix_ele)%vec(3)
      rms_x = sqrt(sum(dx(:)*dx(:))/size(dx))
      rms_y = sqrt(sum(dy(:)*dy(:))/size(dy))
    enddo
    write(*,'(a,2es14.5)') "rms_x rms_y: ", rms_x, rms_y

    !Correct horizontal orbit
    allocate(A(n_bpms,n_h_correctors))
    allocate(Ap(n_h_correctors,n_bpms))
    call analytic_orm(ring,I1,'x',A,bpms,h_correctors)
    call mat_pseudoinverse(A,Ap)
    !Make correction vector
    tx = -matmul(Ap,dx)
    !Apply correction vector
    do i=1,n_h_correctors
      h_correctors(i)%ele%value(hkick$) = h_correctors(i)%ele%value(hkick$) + tx(i)/10.0
    enddo
    deallocate(A)
    deallocate(Ap)

    !Correct vertical orbit
    allocate(A(n_bpms,n_v_correctors))
    allocate(Ap(n_v_correctors,n_bpms))
    call analytic_orm(ring,I1,'y',A,bpms,v_correctors)
    call mat_pseudoinverse(A,Ap)
    !Make correction vector
    ty = -matmul(Ap,dy)
    !Apply correction vector
    do i=1,n_v_correctors
      v_correctors(i)%ele%value(vkick$) = v_correctors(i)%ele%value(vkick$) + ty(i)/10.0
    enddo
    deallocate(A)
    deallocate(Ap)

    !Calculate new closed orbit

    call twiss_and_track(ring,co,status)
  enddo
  do i=1,n_bpms
    dx(i) = co(bpms(i)%ele%ix_ele)%vec(1)
    dy(i) = co(bpms(i)%ele%ix_ele)%vec(3)
    rms_x = sqrt(sum(dx(:)*dx(:))/size(dx))
    rms_y = sqrt(sum(dy(:)*dy(:))/size(dy))
  enddo
  write(*,'(a,2es14.5)') "rms_x rms_y: ", rms_x, rms_y


  call radiation_integrals(ring, co, mode)

  WRITE(*,*) "Beam distribution parameters after correction:"
  WRITE(*,*) "   emit_a      : ", mode%a%emittance
  WRITE(*,*) "   emit_b      : ", mode%b%emittance
  WRITE(*,*) "   sigmaE_E    : ", mode%sigE_E
  WRITE(*,*) "   sigma_z     : ", mode%sig_z

  !Write out corrected closed orbit
  open(11,file='corrected_orbit.out')
  write(11,'(6a12)') 'ix', 's', 'x', 'px', 'y', 'py'
  do i=1,n_bpms
    write(11,'(i12,f12.4,4es15.4)') i, bpms(i)%ele%s, co(bpms(i)%ele%ix_ele)%vec(1:4)
  enddo
  close(11)

  !Write out corrections
  open(12,file='orbit_correction.lat')
  write(12,*) "! Horizontal and vertical orbit correction determined by ORM"
  write(12,*) 'call, file = '//lat_file
  do i=1,n_h_correctors
    write(str,'(es17.7)') h_correctors(i)%ele%value(hkick$)
    write(ix_str,'(i6.6)') h_correctors(i)%ele%ix_ele
    write(12,'(a)') trim(ix_str)//'[hkick] ='//trim(str)//'   ! '//trim(h_correctors(i)%ele%name)
  enddo
  write(12,*)
  do i=1,n_v_correctors
    write(str,'(es17.7)') v_correctors(i)%ele%value(vkick$)
    write(ix_str,'(i6.6)') v_correctors(i)%ele%ix_ele
    write(12,'(a)') trim(ix_str)//'[vkick] ='//trim(str)//'   ! '//trim(v_correctors(i)%ele%name)
  enddo
  close(12)

  deallocate(bpms,h_correctors,v_correctors,dx,dy,tx,ty)

end program orm






