program test_lsq_h

USE pointer_lattice
USE orm_mod

implicit none

type :: terms_struct
  integer stage
  integer F_index  !used to index complex vector field.  Only for dhdj
  character(5) name
  real(dp) target
  real(dp) weight
end type
type(terms_struct), allocatable :: termdefs(:)

integer, parameter :: max_fam = 500

character(100) lat_file
character(nlp) magname(max_fam)
character(50) jac_file_name
character(30) format_spec
character(2) nvar_str
character(5) c5
character(1) F_index_str
character(7) current_terms_str(100)

type(layout), pointer :: r
TYPE(FIBRE),POINTER :: p

!PTC datatypes - Not allocatable
type(pol_block) polyb
type(probe) ray
type(internal_state), target :: state

!PTC datatypes - Must be allocated
type(c_damap) cmap, L
type(c_vector_field) F
type(c_normal_form) cnorm
type(damapspin) maps
type(damapspin) ids   ! Initialized to identity map
type(c_taylor) t
type(c_taylor) h
type(probe_8) ray8
type(normalform) normal
type(damap) rdamap

integer n_terms
integer repeet, nfam, n_stages
integer magtype(max_fam)
integer i, j, k, ks
integer fid
integer report_lun
integer nvars
integer iterations
integer selection

logical menu
logical is_new
logical, allocatable :: printed_update(:)

real(dp) tx_limiter, limiter_norm
real(dp) x(6)
real(dp) prec
real(dp) kmag(max_fam)
real(dp) der
real(dp), allocatable :: jac_p(:,:)
real(dp), allocatable :: terms(:)
real(dp), allocatable :: ve(:)
real(dp), allocatable :: merits(:)
real(dp), allocatable :: jac_buf(:)

! coarray (shared) variables
real(dp), allocatable :: tx(:)[:]
real(dp), allocatable :: jac(:,:)[:]
logical loop_optimizer[*]
logical loop_stages[*]
integer stage[*]

! coarray bookkeeping stuff
integer my_worker_num
integer num_workers
integer vec_num
integer counter
integer, allocatable :: work_schedule(:)

! coarray stuff
num_workers = num_images()
my_worker_num = this_image()

! Set global parameters
use_info=.true.
my_estate=>state

! Set local parameters
prec=1.d-8

!----------------------------------------------Initialize PTC and read in universe
call ptc_ini_no_append
call getarg(1,lat_file)
call READ_AND_APPEND_VIRGIN_general(M_U,lat_file)

c_verbose=.false.

r=>m_u%start

! state=delta0
state = default + nocavity0

use_complex_in_ptc=my_true

call init(state,5,1)

call alloc(ids)
call alloc(maps)
call alloc(ray8)
call alloc(cmap)
call alloc(cnorm)
call alloc(L)
call alloc(F)
call alloc(t)
call alloc(h)
call alloc(rdamap)
call alloc(normal)

n_stages = 5
n_terms = 29
tx_limiter = 0.1  !normalize strength adjustment so that no magnet changes by more than this relative amount
allocate(termdefs(n_terms))

           !terms_struct(stage, name, target, weight)
termdefs(1)  = terms_struct(3, 0, '21000', 0.0, 0.1)         ! Qx
termdefs(2)  = terms_struct(3, 0, '30000', 0.0, 0.1)         ! 3Qx
termdefs(3)  = terms_struct(3, 0, '10110', 0.0, 0.1)         ! Qx
termdefs(4)  = terms_struct(3, 0, '10020', 0.0, 0.1)         ! Qx-2Qy
termdefs(5)  = terms_struct(3, 0, '10200', 0.0, 0.1)         ! Qx+2Qy
termdefs(6)  = terms_struct(3, 0, '20001', 0.0, 0.1)         ! Synchro-betatron resonances ???
termdefs(7)  = terms_struct(3, 0, '00201', 0.0, 0.1)         ! Momentum-dependence of beta functions ???
termdefs(8)  = terms_struct(3, 0, '10002', 0.0, 0.1)         ! 2nd order dispersion

termdefs(9)  = terms_struct(3, 0, '31000', 0.0, 0.0001)         ! 2Qx
termdefs(10) = terms_struct(3, 0, '00310', 0.0, 0.0001)         ! 2Qy
termdefs(11) = terms_struct(3, 0, '11200', 0.0, 0.0001)         ! 2Qy
termdefs(12) = terms_struct(3, 0, '20110', 0.0, 0.0001)         ! 2Qx

termdefs(13) = terms_struct(5, 0, '22000', 0.0, 0.01)         ! dQx/dJx
termdefs(14) = terms_struct(5, 0, '00220', 0.0, 0.01)         ! dQy/dJy
termdefs(15) = terms_struct(5, 0, '11110', 0.0, 0.01)         ! dQx/dJy

termdefs(16) = terms_struct(4, 0, '20200', 0.0, 0.0001)         ! 2Qx+2Qy
termdefs(17) = terms_struct(4, 0, '20020', 0.0, 0.0001)         ! 2Qx-2Qy
termdefs(18) = terms_struct(4, 0, '40000', 0.0, 0.0001)         ! 4Qx
termdefs(19) = terms_struct(4, 0, '00400', 0.0, 0.0001)         ! 4Qy

termdefs(20)  = terms_struct(1, 1, '00001',  1.0, 10.)     ! chrom x
termdefs(21)  = terms_struct(1, 1, '00002',  0.0, 1.0)     ! 2nd order chrom x
termdefs(22)  = terms_struct(1, 1, '00003',  0.0, 0.1)     ! 3rd order chrom x
termdefs(23)  = terms_struct(1, 2, '00001',  1.0, 10.)     ! chrom y
termdefs(24)  = terms_struct(1, 2, '00002',  0.0, 1.0)     ! 2nd order chrom y
termdefs(25)  = terms_struct(1, 2, '00003',  0.0, 0.1)     ! 3rd order chrom y
termdefs(26)  = terms_struct(2, 6, '00001', -0.1, 0.0)     ! dz/dp: first order momentum compaction (related by negative sign)
termdefs(27)  = terms_struct(2, 6, '00002',  0.0, 0.0)     ! d2x/dp2: second order momentum compaction (related by negative sign)
termdefs(28)  = terms_struct(1, 6, '00003', -1.0, 10.0)     ! d3x/dp3: third order momentum compaction (related by negative sign)
termdefs(29)  = terms_struct(2, 6, '00004', -5.0, 0.01)     ! d3x/dp3: third order momentum compaction (related by negative sign)


!--------------------------------------------- Locate multipole families
k=0
p=>r%start
do i=1,r%n
  if(p%mag%p%nmul > 2) then
    if(p%mag%bn(3) .ne. 0.d0) then
      is_new = .true.
      do j=1,k
        if(p%mag%name == magname(j)) then  !already found this family
          is_new = .false.
          exit
        endif
      enddo
      if(is_new) then 
        if( my_worker_num == 1) write(6,*) i,p%mag%name,p%mag%bn(3)
        k=k+1
        magname(k)=p%mag%name
        magtype(k)=3
      endif
      !call add(p,3,0,0.d0)   ! sets multipole field to zero
    elseif(p%mag%bn(4) .ne. 0.d0) then
      is_new = .true.
      do j=1,k
        if(p%mag%name == magname(j)) then  !already found this family
          is_new = .false.
          exit
        endif
      enddo
      if(is_new) then 
        if( my_worker_num == 1) write(6,*) i,p%mag%name,p%mag%bn(4)
        k=k+1
        magname(k)=p%mag%name
        magtype(k)=4
      endif
      !call add(p,3,0,0.d0)   ! sets multipole field to zero
    endif
  endif
  p=>p%next
enddo
nfam=k ! Total number of multipoles
if( my_worker_num == 1) write(6,*) "Number of multitupole families found: ", nfam

allocate(printed_update(nfam))
allocate(terms(n_terms))
allocate(ve(n_terms))
allocate(jac_p(nfam,n_terms))
allocate(jac_buf(n_terms))
allocate(merits(n_stages))
! coarray stuff
allocate(jac(n_terms,nfam)[*])
allocate(tx(nfam)[*])

! coarray stuff
IF(my_worker_num == 1) write(*,*) "Number of Images: ", num_workers
allocate(work_schedule(num_workers))
work_schedule = 0
counter = 1
do i=1, nfam
  work_schedule(counter) = work_schedule(counter) + 1
  counter = counter + 1
  IF(counter .gt. num_workers) counter=1
enddo

!--------------------------------------------------------------------
! Calculate initial terms
!--------------------------------------------------------------------
if(my_worker_num == 1) THEN
  ids = 1
  x=0.d0 
  ray=0  
  ray=x  
  ray8=ray+ids  ! ray8 is a probe_8
  call track_probe(r,ray8,state,fibre1=1)

  ! Calculate nonlinear driving terms
  maps=ray8   ! damapspin <- probe_8
  cmap=maps   ! c_damap <- damapspin
  call c_normal(cmap,cnorm,dospin=my_false,no_used=1)  ! c_normal_form <- c_damap
  cmap=cnorm%n
  cmap=from_phasor(-1)*cmap*from_phasor(1)
  call c_factor_map(cmap,L,F,1)
  call c_clean_vector_field(F,F,prec)

  h=cgetpb(F,-1)   ! cgetpb includes 1/(-2*i_) 

  rdamap = ray8    ! damap <- probe_8
  normal = rdamap  ! normalform <- damap

  do i=1, n_terms
    call process_term(h, normal%dhdj%v, termdefs(i)%F_index, termdefs(i)%name, terms(i))
  enddo

  merits = 0
  do i=1,n_terms
    merits(termdefs(i)%stage) = merits(termdefs(i)%stage) + termdefs(i)%weight*abs(terms(i) - termdefs(i)%target)
  enddo

  call kanalnummer(report_lun,'report.out')
  printed_update = .false.
  write(report_lun,*) "***********************************************"
  p=>r%start
  do i=1,r%n
    do j=1,nfam
      if(p%mag%name == magname(j)) then
        if(.not. printed_update(j)) then
          printed_update(j) = .true.
          write(report_lun,'(A,ES17.7)') magname(j), p%mag%bn(magtype(j))
          exit
        endif
      endif
    enddo
    p=>p%next
  enddo
endif
!--------------------------------------------------------------------

!--------------------------------------------------------------------
! Staged Optimization Loop
!--------------------------------------------------------------------
stage = 1
loop_stages = .true.
do while( loop_stages )
  loop_optimizer = .true.
  iterations = 1
  do while( loop_optimizer )
    jac=0.0_dp
    jac_buf = 0.0_dp
    jac_p=0.0_dp
    terms=0.0_dp
    ve=0.0_dp
    tx=0.0_dp

    !------------------------------------- 
    ! For each multipole family i, get dh_xxxxx/dSi, for each driving term h_xxxxx.
    !------------------------------------- 
    if(my_worker_num==1) write(*,*)
    do i=1,work_schedule(my_worker_num)
      ks = my_worker_num + (i-1)*num_workers
      polyb=0                     ! pol_block:  polymorphic_block.  sets a knob.
      polyb%ibn(magtype(ks))=1    ! set normal multipole strength as first parameter
      polyb%name=magname(ks)      ! set name of family
      r=polyb                     ! Overloaded to scan_for_polymorphs(layout,pol_block).  
                                  ! For all magnets in r with name == polyb%name, bn(magtype(ks)) becomes
                                  ! a knob.

      ids = 1
      x=0.d0 
      ray=0  
      ray=x  
      ray8=ray+ids  ! ray8 is a probe_8

      call track_probe(r,ray8,+state,fibre1=1)  ! get the 1-turn map by mapping the identity map.
                                                ! + symbol activates the knob
      call kill_para(r)  ! removes all polymorphic blocks (pol_block) from a layout

      ! Process 1-turn map into dhdj, which contains the amplitude-dependent tune shifts.
      rdamap = ray8    ! damap <- probe_8
      normal = rdamap  ! normalform <- damap

      ! call print(normal%dhdj,6)

      ! Process 1-turn map into h, which contains the h_xxxxx driving terms.
      maps=ray8   ! damapspin <- probe_8
      cmap=maps   ! c_damap <- damapspin
      call c_normal(cmap,cnorm,dospin=my_false,no_used=1)  ! c_normal_form <- c_damap
      cmap=cnorm%n
      cmap=from_phasor(-1)*cmap*from_phasor(1)
      call c_factor_map(cmap,L,F,1)
      call c_clean_vector_field(F,F,prec)
      h=cgetpb(F,-1)   ! cgetpb includes 1/(-2*i_)   ! TODO: figure out crash if cgetpb replaced by getpb

      nvars = 0
      do j=1, n_terms
        call process_term(h, normal%dhdj%v, termdefs(j)%F_index, termdefs(j)%name, terms(j), der)

        if(termdefs(j)%stage .le. stage) then
          nvars = nvars + 1
          ve(nvars) = termdefs(j)%weight * (terms(j) - termdefs(j)%target)
          jac_buf(nvars) = der
          write(F_index_str,'(I1)') termdefs(j)%F_index !FOO
          current_terms_str(nvars) = F_index_str//'_'//termdefs(j)%name
        endif
      enddo

      write(*,*) "Completed Family Number ",ks, " of ", nfam
      jac(1:nvars,ks)[1] = jac_buf(1:nvars)
    enddo

    sync all  !coarray: wait here for all processes to finish calculating terms
    if(my_worker_num == 1) then
      merits = 0
      do i=1,n_terms
        merits(termdefs(i)%stage) = merits(termdefs(i)%stage) + termdefs(i)%weight*abs(terms(i) - termdefs(i)%target)
      enddo
      
      write(*,*) "Norms:"
      do i=1,n_stages
        write(*,'(A,I3,A,F20.5)') "Stage ", i, " norm: ", merits(i)
      enddo
      write(*,*)

      write(*,*) "Terms:"
      do i=1,n_terms
        write(*,'(A,I3,A,A,F25.10)') "     ", termdefs(i)%F_index, " ", termdefs(i)%name,  terms(i)
        write(report_lun,'(A,I3,A,A,F25.10)') "     ", termdefs(i)%F_index, " ", termdefs(i)%name,  terms(i)
      enddo
      write(*,*)

      !--------------------------------- Use pseudoinverse to get least squares solution
      call make_pseudoinverse(jac(1:nvars,1:nfam),jac_p(1:nfam,1:nvars))
      tx(1:nfam) = -matmul(jac_p(1:nfam,1:nvars),ve(1:nvars))
      limiter_norm = 1.0_dp
      p=>r%start
      do i=1,r%n
        do j=1,nfam
          if(p%mag%name == magname(j)) then
            if( abs(tx(j)/p%mag%bn(magtype(j))) .gt. tx_limiter ) then 
              if ( abs(tx(j)/p%mag%bn(magtype(j)))/tx_limiter .gt. limiter_norm ) then
                limiter_norm = abs(tx(j)/p%mag%bn(magtype(j)))/tx_limiter
                write(*,'(A,F20.10,A,A)') "Step limiter adjusted to: ", limiter_norm, " because of magnet family ", magname(j)
              endif
            endif
            exit
          endif
        enddo
        p=>p%next
      enddo
      tx = tx/limiter_norm
      do i=2,num_images()
        tx(1:nfam)[i] = tx(1:nfam)  !distribute multipole strength adjustments to all other processes 
      enddo
    endif

    !--------------------------------- Adjust multitupoles according to least squares solution
    sync all  !coarray: wait here for worker 1 to finish calculating and sending updated multitupole strengths
    printed_update = .false.
    if(my_worker_num == 1) write(report_lun,*) "***********************************************"
    p=>r%start
    do i=1,r%n
      do j=1,nfam
        if(p%mag%name == magname(j)) then
          call add(p,magtype(j),1,tx(j))
          kmag(j)=p%mag%bn(magtype(j))
          if(my_worker_num == 1) then
            if(.not. printed_update(j)) then
              printed_update(j) = .true.
              write(*,'(A,I3,A,A,A,ES17.7,A,ES17.7,A)') "Adjusted family ", j, "   ", p%mag%name, " by ", tx(j), " (new strength: ", p%mag%bn(magtype(j)), " )"
              write(report_lun,'(A,ES17.7)') magname(j), kmag(j)
            endif
          endif
          exit
        endif
      enddo
      p=>p%next
    enddo
      
    if(my_worker_num == 1) then
      !---------------------------------- Prepare for next loop
      iterations = iterations - 1
      if(iterations .lt. 1) then
        write(6,*) "Finished stage = ", stage
        menu = .true.
        do while(menu)
          write(6,*) "Make a selection:"
          write(6,*) "  1) Iterate at current stage"
          write(6,*) "  2) Move on to next stage"
          write(6,*) "  3) Write Jacobian to file"
          write(6,*) "  4) Terminate program"
          write(6,'(A)',advance='no') "Selection: "
          read(5,*) selection
          if( selection == 1 ) then
            write(6,'(A)',advance='no') "Enter number of iterations to repeat at this stage: "
            read(5,*) iterations
            loop_optimizer = .true.
            menu = .false.
          elseif( selection == 2 ) then
            iterations = 0
            stage = stage + 1
            loop_optimizer = .false.
            menu = .false.
          elseif( selection == 3 ) then
            write(6,'(A)',advance='no') "Enter name for Jacobian file: "
            read(5,*) jac_file_name
            write(nvar_str,'(I2)') nfam
            open(111,file=trim(jac_file_name))
            format_spec = '(A,'//nvar_str//'A14)'
            write(111,format_spec) '              ', magname(1:nfam)
            format_spec = '(A,'//nvar_str//'ES14.4)'
            do i=1, nvars
              write(111,format_spec) current_terms_str(i), jac(i,1:nfam)
            enddo
            close(111)
            menu = .true.
          elseif( selection == 4 ) then
            iterations = 0
            menu = .false.
            loop_optimizer = .false.
            loop_stages = .false.
          endif
        enddo
      endif
      do i=2,num_images()
        loop_optimizer[i] = loop_optimizer  !distribute loop_optimizer to all other processes 
        loop_stages[i] = loop_stages 
        stage[i] = stage
      enddo
    endif
    sync all  !coarray stuff
  enddo
enddo

deallocate(jac_buf)
deallocate(jac)
deallocate(jac_p)
deallocate(terms)
deallocate(termdefs)
deallocate(ve)
deallocate(tx)

if(my_worker_num == 1) then
  close(report_lun)
  CALL print_COMPLEX_SINGLE_STRUCTURE(r,"ptc_new.lat")
endif

call ptc_end

contains

  subroutine process_term(h_local, v_local, F_index, lie_index, mag, der)
  
    implicit none

    type(c_taylor) h_local
    type(taylor) v_local(ndim2)
    integer F_index
    character(5) lie_index
    real(dp) mag
    real(dp), optional :: der

    !working space
    type(taylor) t_local
    type(c_taylor) c_t_local
    real(dp) hr, hi, dhrdS, dhidS


    if( F_index .gt. 0 ) then
      call alloc(t_local)
      t_local = v_local(F_index)
      t_local = t_local.par.lie_index
      t_local=t_local<=6  ! shift indexes down by 6, which moves the knob index 7 into
                          ! the first position.

      mag = t_local.sub.'0'
      if( present(der) ) then
        der = t_local.sub.'1'
      endif
      call kill(t_local)
    else
      call alloc(c_t_local)
      c_t_local=h_local.par.lie_index  ! pull desired polynomial coefficient from h.
      c_t_local=c_t_local<=6  ! shift indexes down by 5, which moves the knob index 6 into
                          ! the first position.

      dhrdS = real(c_t_local).sub.'1'
      dhidS = aimag(c_t_local).sub.'1'
      hr = real(c_t_local).sub.'0'
      hi = aimag(c_t_local).sub.'0'
      mag = sqrt(hr**2 + hi**2)
      if( present(der) ) then
        der = (hr*dhrdS+hi*dhidS) / mag
      endif
      call kill(c_t_local)
    endif

  end subroutine process_term

end program test_lsq_h








