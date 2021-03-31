program ptc_response_sext

USE pointer_lattice

implicit none

integer, parameter :: max_fam = 100
integer, parameter :: max_per_fam = 100

character(nlp) sexname(max_fam)
integer possex(max_fam,max_per_fam)
real(dp) ksex(max_fam),ksex0(max_fam)
integer nfam, fam_num
integer nper(max_fam)

character(5) c5

type(layout), pointer :: r
TYPE(FIBRE),POINTER :: p

type(pol_block) sexb
type(internal_state), target :: state
type(c_damap) cmap, L
type(c_vector_field) F
type(c_normal_form) cnorm
type(damapspin) maps
type(damapspin) ids   ! Initialized to identity map
type(c_taylor) t
type(probe_8) ray8
type(probe) ray
type(normalform) normal
type(damap) rdamap

integer neq,repeet,run_number
integer i,j,k,mf,ks

logical run_lsq, is_new

real(dp) xv
real(dp) x(6)
real(dp) prec
real(dp), allocatable :: ma(:,:),ve(:),gt(:,:)

! Set global parameters
lmax=100.0d0
use_info=.true. 
my_estate=>state

! Set local parameters
prec=1.d-8

!! This program reads a flat file
!! and produce by PTC via BMAD

! It runs the iterations in four groups

! Runs=1 Fits in one iteration the tunes and all sextupole geometric resonances
! Runs=2 Fits first most deadly quartic resonances : interactive
! Runs=3 Fits additional quartic resonances
! Runs=4 Fits all the rest as well (difference resonances and tune shifts)

! The file "ksex.txt" contains the Lie polynomial in phasors for easier checking that all is indeed zero
! It also contains the old and new sextupole strengths in PTC units

!----------------------------------------------Initialize PTC and read in universe
call ptc_ini_no_append
call READ_AND_APPEND_VIRGIN_general(M_U,'ptc.lat')

c_verbose=.false.
call kanalnummer(mf,"eq.txt")

r=>m_u%start

state=delta0

use_complex_in_ptc=my_true

call init(state,4,1)   ! state, order of polynomials, number of parameters
                       !init is overloaded so S_init, resulting call: init(NO1=4, ND1=2, NP1+NDEL=1+1=2, NDPT1=0, PACKAGE is not present)

call alloc(ids)     ! damapspin:  set to identity map
call alloc(maps)    ! damapspin
call alloc(ray8)    ! probe_8
call alloc(cmap)    ! c_damap
call alloc(cnorm)   ! c_normal_form
call alloc(L)       ! c_damap
call alloc(F)       ! c_vector_vield
call alloc(t)       ! c_taylor

!--------------------------------------------- Locate sextupoles families
possex(:,:) = 0  ! 2d array of sextupole locations.
nfam=0  ! number of unique sextupole names in lattice
nper(:)=0
p=>r%start
do i=1,r%n
  if(p%mag%p%nmul>2) then
    if(p%mag%bn(3)/=0.d0) then
      is_new = .true.
      do j=1,nfam
        if(p%mag%name == sexname(j)) then  !already found this family
          is_new = .false.
          fam_num = j
          exit
        endif
      enddo
      if( is_new ) then
        write(6,*) i,p%mag%name,p%mag%bn(3)
        nfam=nfam+1
        sexname(nfam)=p%mag%name
        ksex0(nfam)=p%mag%bn(3)
        fam_num = nfam
        nper(fam_num) = 1
      else
        nper(fam_num) = nper(fam_num) + 1
      endif
      possex(fam_num,nper(fam_num))=i
      call add(p,3,0,0.d0)
    endif
  endif
  p=>p%next
enddo
write(6,*) "Number of sextupoles found: ", nfam

do run_number=1,4
  neq= 12
  if(run_number==2) then
    neq=18
  endif
  if(run_number==3) then
    neq=24
  endif
  if(run_number==4) then
    neq=31
  endif

  allocate(ma(neq,nfam))
  allocate(ve(neq+nfam))
  allocate(gt(neq+nfam,neq+nfam))

  run_lsq = .true.
  do while( run_lsq )
    ma=0.0_dp
    ve=0.0_dp
    gt=0.0_dp

    !------------------------------------- 
    ! For each sextupole: 
    do ks=1,nfam
      sexb=0
      sexb%ibn(3)=1
      sexb%name=sexname(ks)
      r=sexb  ! nested overloaded assignment operators.  All elements in r with name matching sexb%name, have
              ! their ibn(3) overwritten with 1.  ???

      x(:)=0.d0  !real
      ray=0      !probe
      ray=x
      ids=1      !damapspin, initialized to identity
      ray8=ray+ids

      ! Obtain 1-turn taylor map by tracking the identity map around the ring
      call track_probe(r,ray8,+state,fibre1=1)

      maps=ray8   ! damapspin <- probe_8
      cmap=maps   ! c_damap <- damapspin
      call c_normal(cmap,cnorm,dospin=my_false,no_used=1)  ! c_normal_form <- c_damap
      cmap=cnorm%a_t**(-1)*cmap*cnorm%a_t
      cmap=from_phasor(-1)*cmap*from_phasor(1)
      call c_factor_map(cmap,L,F,1)
      call c_clean_vector_field(F,F,prec)

!      rdamap = ray8    ! damap <- probe_8
!      normal = rdamap  ! normalform <- damap
!      normal%dhdj%v.par.

      c5='10001'
      t=F%v(1).par.c5
      t=t<=5
      if(c_verbose) call print(t,mf)
      ma(1,ks)=aimag(t).sub.'1'
      ve(1)=aimag(t).sub.'0'

      c5='00101'
      t=F%v(3).par.c5
      t=t<=5
      if(c_verbose) call print(t,mf)
      ma(2,ks)=aimag(t).sub.'1'
      ve(2)=aimag(t).sub.'0'

      c5='20000' ! 1 0 J_x
      t=F%v(1).par.c5
      t=t<=5
      if(c_verbose) call print(t,mf)
      ma(3,ks)=real(t).sub.'1'
      ma(4,ks)=aimag(t).sub.'1'
      ve(3)=real(t).sub.'0'
      ve(4)=aimag(t).sub.'0'

      c5='02000' ! 3 0
      t=F%v(1).par.c5
      t=t<=5
      if(c_verbose) call print(t,mf)
      ma(5,ks)=real(t).sub.'1'
      ma(6,ks)=aimag(t).sub.'1'
      ve(5)=real(t).sub.'0'
      ve(6)=aimag(t).sub.'0'

      c5='00200' ! -1 2
      t=F%v(1).par.c5
      t=t<=5
      if(c_verbose) call print(t,mf)
      ma(7,ks)=real(t).sub.'1'
      ma(8,ks)=aimag(t).sub.'1'
      ve(7)=real(t).sub.'0'
      ve(8)=aimag(t).sub.'0'

      c5='00110' ! 1 0 J_y
      t=F%v(1).par.c5
      t=t<=5
      if(c_verbose) call print(t,mf)
      ma(9,ks)=real(t).sub.'1'
      ma(10,ks)=aimag(t).sub.'1'
      ve(9)=real(t).sub.'0'
      ve(10)=aimag(t).sub.'0'

      c5='00020' ! 1 2
      t=F%v(1).par.c5
      t=t<=5
      if(c_verbose) call print(t,mf)
      ma(11,ks)=real(t).sub.'1'
      ma(12,ks)=aimag(t).sub.'1'
      ve(11)=real(t).sub.'0'
      ve(12)=aimag(t).sub.'0'

      if( neq==18 .or. neq==24 .or. neq==31 ) then
        c5='03000' ! 4 0
        t=F%v(1).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)

        ma(13,ks)=real(t).sub.'1'
        ma(14,ks)=aimag(t).sub.'1'
        ve(13)=real(t).sub.'0'
        ve(14)=aimag(t).sub.'0'

        c5='01020' ! 2 2
        t=F%v(1).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)

        ma(15,ks)=real(t).sub.'1'
        ma(16,ks)=aimag(t).sub.'1'
        ve(15)=real(t).sub.'0'
        ve(16)=aimag(t).sub.'0'

        c5='00030' ! 0 4
        t=F%v(3).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(17,ks)=real(t).sub.'1'
        ma(18,ks)=aimag(t).sub.'1'
        ve(17)=real(t).sub.'0'
        ve(18)=aimag(t).sub.'0'
      endif

      if( neq==24 .or. neq==31 ) then
        c5='12000' ! 2 0 J_x
        t=F%v(1).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(19,ks)=real(t).sub.'1'
        ma(20,ks)=aimag(t).sub.'1'
        ve(19)=real(t).sub.'0'
        ve(20)=aimag(t).sub.'0'


        c5='10020' ! 0 2 J_x
        t=F%v(1).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(21,ks)=real(t).sub.'1'
        ma(22,ks)=aimag(t).sub.'1'
        ve(21)=real(t).sub.'0'
        ve(22)=aimag(t).sub.'0'


        c5='01110' ! 2 0 J_y
        t=F%v(1).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(23,ks)=real(t).sub.'1'
        ma(24,ks)=aimag(t).sub.'1'
        ve(23)=real(t).sub.'0'
        ve(24)=aimag(t).sub.'0'
      endif

      if( neq==31 ) then
        c5='21000' ! J_x**2 tune
        t=F%v(1).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(25,ks)=aimag(t).sub.'1'
        ve(25)=aimag(t).sub.'0'

        c5='10110' ! J_x*J_y tune
        t=F%v(1).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(26,ks)=aimag(t).sub.'1'
        ve(26)=aimag(t).sub.'0'

        c5='00210' ! J_y**2 tune
        t=F%v(3).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(27,ks)=aimag(t).sub.'1'
        ve(27)=aimag(t).sub.'0'

        c5='01200' ! 2 -2
        t=F%v(1).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(28,ks)=real(t).sub.'1'
        ma(29,ks)=aimag(t).sub.'1'
        ve(28)=real(t).sub.'0'
        ve(29)=aimag(t).sub.'0'

        c5='00120' ! 0 2 J_y
        t=F%v(3).par.c5
        t=t<=5
        if(c_verbose) call print(t,mf)
        ma(30,ks)=real(t).sub.'1'
        ma(31,ks)=aimag(t).sub.'1'
        ve(30)=real(t).sub.'0'
        ve(31)=aimag(t).sub.'0'
      endif

      write(6,*) "Tunes: ",  cnorm%tune(1:2)
      write(6,*) "ks ",ks

      call kill_para(r)  ! resets all parameters in layout
    enddo
    close(mf)

    !--------------------------------- Make response matrix
    do i=1,nfam
      gt(neq+i,neq+i)=1.0_dp
    enddo

    do i=1,neq
      do ks=1,nfam
        gt(i,neq+ks)=ma(i,ks)
        gt(neq+ks,i)=ma(i,ks)
      enddo
    enddo
    ks=neq+nfam
    call matinv(gt,gt,ks,ks,i)
    write(6,*) " matinv error status (should be zero): ",i

    !--------------------------------- Calculate Norms
    xv=0
    do i=1,12
      xv=xv+abs(ve(i))
    enddo
    write(6,*) " Norm of sextupole driving terms: ",xv

    xv=0
    do i=13,neq
      xv=xv+abs(ve(i))
    enddo
    write(6,*) " Norm of octupole driving terms: ",xv

    !--------------------------------- Adjust sextupoles according to response matrix
    ve=-ve
    ve=matmul(gt,ve)
    k=1
    p=>r%start
    do i=1,r%n
      IF(possex(k)==i) then
        call add(p,3,1,VE(neq+k))
        ksex(k)=p%mag%bn(3)
        k=k+1
      endif

      p=>p%next
    enddo

    !---------------------------------- Prepare for next loop
    if(run_number>1) then
      write(6,*) "Finished iteration at neq = ",neq
      write(6,*) "To iterate again at same neq, enter 1"
      read(5,*) repeet
      if(repeet==1) then
        run_lsq = .true.
      else
        run_lsq = .false.
        deallocate(ma)
        deallocate(ve)
        deallocate(gt)
      endif
    else
      run_lsq = .false.
      deallocate(ma)
      deallocate(ve)
      deallocate(gt)
    endif
  enddo
enddo

CALL print_COMPLEX_SINGLE_STRUCTURE(r,"ptc_new.lat")

!-------------------------------------- Print resulting driving terms in Lie Polynomial form
call kanalnummer(mf,"ksex.txt")

t=getpb(F,-1)
t=t/(-2*i_)
write(mf,*) " %%%%%%%%%%%%%%%%%%%% Lie polynomial %%%%%%%%%%%%%%%%%%%%% "
call print(t,mf)
write(mf,*) " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% "
do i=1,nfam
  write(MF,*) i,ksex0(i),ksex(i)
enddo

close(mf)

call ptc_end

end program 


