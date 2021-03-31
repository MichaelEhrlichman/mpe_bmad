program prm
  use bmad
  use bmad_parser_mod, only: bp_com
  use sls_lib
  use analytic_rm_mod
  use f95_lapack
  use make_pseudoinverse_mod  !from util_programs

  implicit none

  type(lat_struct) ideal_ring, ma_ring
  type(coord_struct), allocatable :: co(:)
  type(normal_modes_struct) mode

  character(100) in_file
  character(200) ideal_lat_file
  character(200) ma_lat_file
  character(200) ma_lat_file_override

  integer i,j,k
  integer status

  logical err

  character(200) quad_corrector_mask
  character(20) bpm_mask
  type(ele_pointer_struct), allocatable :: quad_correctors(:)
  type(ele_pointer_struct), allocatable :: bpms(:)
  integer n_quad_correctors
  integer n_bpms
  integer n_loops

  real(rp), allocatable :: res(:)
  real(rp), allocatable :: trim_k1(:)
  real(rp), allocatable :: sum_trim_k1(:)
  real(rp), allocatable :: A(:,:)
  real(rp), allocatable :: Ap(:,:)

  character(6) ix_str
  character(17) str

  namelist /parameters/  ideal_lat_file,&
                         ma_lat_file, &
                         n_loops,&
                         quad_corrector_mask,&
                         bpm_mask

  call getarg(1,in_file)
  open(10,file=in_file,status='old')
  read(10,nml=parameters)
  close(10)

  ma_lat_file_override = ''
  call getarg(2,ma_lat_file_override)
  if(ma_lat_file_override .ne. '') then
    ma_lat_file = ma_lat_file_override
  endif

  bmad_com%auto_bookkeeper = .false.

  bp_com%always_parse = .true.
  call bmad_parser(ma_lat_file, ma_ring)
  call bmad_parser(ideal_lat_file, ideal_ring)
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  call twiss_and_track(ideal_ring,co,status)
  deallocate(co)

  call twiss_and_track(ma_ring,co,status)
  call radiation_integrals(ma_ring, co, mode)

  write(*,*) "Beam distribution parameters before correction:"
  write(*,*) "   emit_a      : ", mode%a%emittance
  write(*,*) "   emit_b      : ", mode%b%emittance
  write(*,*) "   sigmaE_E    : ", mode%sigE_E
  write(*,*) "   sigma_z     : ", mode%sig_z

  call lat_ele_locator(quad_corrector_mask, ma_ring, quad_correctors, n_quad_correctors)
  call lat_ele_locator(bpm_mask, ma_ring, bpms, n_bpms)

  write(*,'(a,i6)') "Number of quadrupole correctors found: ", n_quad_correctors
  write(*,'(a,i6)') "Number of bpms found: ", n_bpms

  allocate(trim_k1(n_quad_correctors))
  allocate(sum_trim_k1(n_quad_correctors))
  allocate(res(3*n_bpms))

  !Write out uncorrected phase advance and dispersion residual
  open(10,file='uncorrected_phi_and_eta_residuals.out')
  write(10,'(4a12)') 'ix', 's', 'dphi_x', 'dphi_y'
  do i=1,n_bpms
    write(10,'(i12,f12.4,4es15.4)') i, bpms(i)%ele%s, ideal_ring%ele(bpms(i)%ele%ix_ele)%a%phi-ma_ring%ele(bpms(i)%ele%ix_ele)%a%phi, &
                                                      ideal_ring%ele(bpms(i)%ele%ix_ele)%b%phi-ma_ring%ele(bpms(i)%ele%ix_ele)%b%phi, &
                                                      ideal_ring%ele(bpms(i)%ele%ix_ele)%a%eta-ma_ring%ele(bpms(i)%ele%ix_ele)%a%eta, &
                                                      ideal_ring%ele(bpms(i)%ele%ix_ele)%b%eta-ma_ring%ele(bpms(i)%ele%ix_ele)%b%eta
  enddo
  close(10)

  sum_trim_k1(:) = 0.0d0
  do j=1,n_loops
    !Make measurement vector
    do i=1,n_bpms
      res(i)          = ma_ring%ele(bpms(i)%ele%ix_ele)%a%phi-ideal_ring%ele(bpms(i)%ele%ix_ele)%a%phi
      res(i+n_bpms)   = -(ma_ring%ele(bpms(i)%ele%ix_ele)%b%phi-ideal_ring%ele(bpms(i)%ele%ix_ele)%b%phi)
      res(i+2*n_bpms) = (ma_ring%ele(bpms(i)%ele%ix_ele)%a%eta-ideal_ring%ele(bpms(i)%ele%ix_ele)%a%eta)
      !res(i+3*n_bpms) = -(ma_ring%ele(bpms(i)%ele%ix_ele)%b%eta-ideal_ring%ele(bpms(i)%ele%ix_ele)%b%eta)
    enddo

    allocate(A(3*n_bpms,n_quad_correctors))
    allocate(Ap(n_quad_correctors,3*n_bpms))
    call analytic_prm(ma_ring,A(1:2*n_bpms,:),bpms,quad_correctors)
    call numerical_drm(ma_ring,A(2*n_bpms+1:3*n_bpms,:),bpms,quad_correctors)
    call make_pseudoinverse(A,Ap,0.0001d0)
    !Make correction vector
    trim_k1 = -matmul(Ap,res)
    sum_trim_k1(:) = sum_trim_k1(:) + trim_k1(:)
    !Apply correction vector
    do i=1,n_quad_correctors
      quad_correctors(i)%ele%value(k1$) = quad_correctors(i)%ele%value(k1$) + trim_k1(i)
      call set_flags_for_changed_attribute(quad_correctors(i)%ele, quad_correctors(i)%ele%value(k1$))
    enddo
    deallocate(A)
    deallocate(Ap)

    call lattice_bookkeeper(ma_ring)
    call twiss_and_track(ma_ring,co,status)

    write(*,'(a,i3,a,f15.7)') "Iteration ", j, " complete!", sum(abs(trim_k1))/size(trim_k1)
  enddo

  call radiation_integrals(ma_ring, co, mode)

  WRITE(*,*) "Beam distribution parameters after correction:"
  WRITE(*,*) "   emit_a      : ", mode%a%emittance
  WRITE(*,*) "   emit_b      : ", mode%b%emittance
  WRITE(*,*) "   sigmaE_E    : ", mode%sigE_E
  WRITE(*,*) "   sigma_z     : ", mode%sig_z

  !Write out corrected closed orbit and residual dispersion
  open(11,file='corrected_phi_and_eta_residual.out')
  write(11,'(4a12)') 'ix', 's', 'dphi_x', 'dphi_y'
  do i=1,n_bpms
    write(11,'(i12,f12.4,4es15.4)') i, bpms(i)%ele%s, ideal_ring%ele(bpms(i)%ele%ix_ele)%a%phi-ma_ring%ele(bpms(i)%ele%ix_ele)%a%phi, &
                                                      ideal_ring%ele(bpms(i)%ele%ix_ele)%b%phi-ma_ring%ele(bpms(i)%ele%ix_ele)%b%phi, &
                                                      ideal_ring%ele(bpms(i)%ele%ix_ele)%a%eta-ma_ring%ele(bpms(i)%ele%ix_ele)%a%eta, &
                                                      ideal_ring%ele(bpms(i)%ele%ix_ele)%b%eta-ma_ring%ele(bpms(i)%ele%ix_ele)%b%eta
  enddo
  close(11)

  !Write out corrections
  open(12,file='phase_correction.lat')
  write(12,*) "! Horizontal and vertical phase correction determined by PRM"
  write(12,*) 'call, file = '//ma_lat_file
  do i=1,n_quad_correctors
    write(str,'(es17.7)') quad_correctors(i)%ele%value(k1$)
    write(ix_str,'(i6.6)') quad_correctors(i)%ele%ix_ele
    write(12,'(a,f14.5)') trim(ix_str)//'[k1] ='//trim(str)//'   ! '//trim(quad_correctors(i)%ele%name)//' total trim: ',sum_trim_k1(i)
  enddo
  close(12)

  deallocate(bpms,quad_correctors,res,trim_k1)

end program






