program orm_print
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

  character(18) ix_str
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

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, ring)
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  call closed_orbit_from_tracking(ring, co, 6)
  call lat_make_mat6 (ring, -1, co, err_flag = err)
  call twiss_at_start (ring, status)
  call twiss_propagate_all (ring, err_flag=err)

  !call twiss_and_track(ring,co,status)
  call radiation_integrals(ring, co, mode)
  I1 = mode%synch_int(1)

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

  !  !Make measurement vector
  !  do i=1,n_bpms
  !    dx(i) = co(bpms(i)%ele%ix_ele)%vec(1)
  !    dy(i) = co(bpms(i)%ele%ix_ele)%vec(3)
  !  enddo

  !Correct horizontal orbit
  allocate(A(n_bpms,n_h_correctors))
  allocate(Ap(n_h_correctors,n_bpms))
  call analytic_orm(ring,I1,'x',A,bpms,h_correctors)
  !call make_pseudoinverse(A,Ap)
  call mat_inverse(A,Ap)
  open(45,file='orm_Ap_h.dat')
  !write(45,'(a14)',advance='no') "#             "
  !do i=1,n_bpms
  !  write(ix_str,*) (i-1)/4+1
  !  write(45,'(a14)',advance='no') trim(adjustl(ix_str))//'-'//trim(bpms(i)%ele%name)
  !enddo
  !write(45,*)
  do i=1,size(Ap(:,1))
    !write(45,'(a14)', advance='no') h_correctors(i)%ele%name
    do j=1,size(Ap(1,:))
      write(45,'(es14.4)',advance='no') Ap(i,j)
    enddo
    write(45,*)
  enddo
  close(45)
  open(45,file='orm_A_h.dat')
  do i=1,size(A(:,1))
    do j=1,size(A(1,:))
      write(45,'(es14.4)',advance='no') A(i,j)
    enddo
    write(45,*)
  enddo
  close(45)
  !   !Make correction vector
  !   tx = -matmul(Ap,dx)
  !   !Apply correction vector
  !   do i=1,n_h_correctors
  !     h_correctors(i)%ele%value(hkick$) = h_correctors(i)%ele%value(hkick$) + tx(i)/10.0
  !   enddo
  deallocate(A)
  deallocate(Ap)

  !Correct vertical orbit
  allocate(A(n_bpms,n_v_correctors))
  allocate(Ap(n_v_correctors,n_bpms))
  call analytic_orm(ring,I1,'y',A,bpms,v_correctors)
  !call make_pseudoinverse(A,Ap)
  call mat_inverse(A,Ap)
  open(45,file='orm_Ap_v.dat')
  do i=1,size(Ap(:,1))
    do j=1,size(Ap(1,:))
      write(45,'(es14.4)',advance='no') Ap(i,j)
    enddo
    write(45,*)
  enddo
  close(45)
  open(45,file='orm_A_v.dat')
  do i=1,size(A(:,1))
    do j=1,size(A(1,:))
      write(45,'(es14.4)',advance='no') A(i,j)
    enddo
    write(45,*)
  enddo
  close(45)
  !   !Make correction vector
  !   ty = -matmul(Ap,dy)
  !   !Apply correction vector
  !   do i=1,n_v_correctors
  !     v_correctors(i)%ele%value(vkick$) = v_correctors(i)%ele%value(vkick$) + ty(i)/10.0
  !   enddo
  deallocate(A)
  deallocate(Ap)
end program






