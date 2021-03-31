program ptc_et

use pointer_lattice

implicit none

type :: terms_struct
  integer F_index  !used to index complex vector field.  Only for dhdj
  character(5) name
  integer, allocatable :: j(:)
  complex(dp) c_val
  character(100) desc
end type
type(terms_struct), allocatable :: terms(:)

character(100) lat_file

type(layout), pointer :: r

!PTC datatypes - Not allocatable
type(probe) ray
type(internal_state), target :: state

!PTC datatypes - Must be allocated
type(c_damap) cmap, L
type(c_damap) cmap_test
type(c_vector_field) F
type(c_normal_form) cnorm
type(c_damap) maps
type(c_damap) ids   ! Initialized to identity map
type(c_taylor) h
type(c_taylor) test_t
type(c_taylor) tunes(3)
type(probe_8) ray8
type(normalform) normal
type(c_damap) a0, a1, a2, acs
type(fibre), pointer:: p
complex(dp) test_c
logical ch
real(dp)norm,nm
integer nj
integer, allocatable :: test_i(:)
integer nd(11)
integer i, j, k,no,np
integer pos_a, pos_b
integer n_terms, nix

real(dp) x(6)
real(dp) prec,thi

! Set global parameters
use_info=.true.
my_estate=>state

! Set local parameters
prec=1.d-8
lat_file="ptc.lat"
!----------------------------------------------Initialize PTC and read in universe
call ptc_ini_no_append
!call getarg(1,lat_file)
call READ_AND_APPEND_VIRGIN_general(M_U,lat_file)


c_verbose=.false.

r=>m_u%start

thi=-1
thi=1.d0
!split for speed.  CALL THIN_LENS_resplit(r,thi)
!goto 999
state=delta0
!state = default + nocavity0
!state = only_4d
!state = default

!use_complex_in_ptc=my_true

no=5
np=0
call init(state,no,np)

call alloc(tunes)
call alloc(a0,a1,a2,acs)
call alloc(ids)
call alloc(maps)
call alloc(ray8)
call alloc(cmap)
call alloc(cmap_test)
call alloc(cnorm)
call alloc(L)
call alloc(F)
call alloc(h)
call alloc(test_t)
call alloc(normal)

n_terms = 29
!    ndct=iabs(ndpt-ndptb)  ! 1 if coasting, otherwise 0
!    ndc2t=2*ndct  ! 2 if coasting, otherwise 0
!    nd2t=nd2-2*rf-ndc2t   !  size of harmonic oscillators minus modulated clocks
!    ndt=nd2t/2        ! ndt number of harmonic oscillators minus modulated clocks(*)
!    nd2harm=nd2t+2*rf  !!!!  total dimension of harmonic phase space
!    ndharm=ndt+rf  !!!! total number of harmonic planes
! n(1)=NO
! n(2)=ND
! n(3)=ND2
! n(4)=NV
! n(5)=Ndpt
! n(6)=ndptb
! n(7)=np
! n(8)=rf*2
! n(9)=ndc2t
! n(10)=nd2t
! n(11)=nd2harm
nix = c_%npara
  
allocate(terms(n_terms))
terms(1)  = terms_structor(0, '21000', 'Qx')
terms(2)  = terms_structor(0, '30000', '3Qx')
terms(3)  = terms_structor(0, '10110', 'Qx')
terms(4)  = terms_structor(0, '10020', 'Qx-2Qy')
terms(5)  = terms_structor(0, '10200', 'Qx+2Qy')
terms(6)  = terms_structor(0, '20001', 'Synchro-betatron resonances ???')
terms(7)  = terms_structor(0, '00201', 'Momentum-dependence of beta functions ???')
terms(8)  = terms_structor(0, '10002', '2nd order dispersion')

terms(9)  = terms_structor(0, '31000', '2Qx')
terms(10) = terms_structor(0, '00310', '2Qy')
terms(11) = terms_structor(0, '11200', '2Qy')
terms(12) = terms_structor(0, '20110', '2Qx')

terms(13) = terms_structor(0, '22000', 'dQx/dJx')
terms(14) = terms_structor(0, '00220', 'dQy/dJy')
terms(15) = terms_structor(0, '11110', 'dQx/dJy')

terms(16) = terms_structor(0, '20200', '2Qx+2Qy')
terms(17) = terms_structor(0, '20020', '2Qx-2Qy')
terms(18) = terms_structor(0, '40000', '4Qx')
terms(19) = terms_structor(0, '00400', '4Qy')

terms(20) = terms_structor(1, '00001', 'chrom x')
terms(21) = terms_structor(1, '00002', '2nd order chrom x')
terms(22) = terms_structor(1, '00003', '3rd order chrom x')
terms(23) = terms_structor(2, '00001', 'chrom y')
terms(24) = terms_structor(2, '00002', '2nd order chrom y')
terms(25) = terms_structor(2, '00003', '3rd order chrom y')
terms(26) = terms_structor(3, '00001', 'dz/dp: first order momentum compaction (related by negative sign)')
terms(27) = terms_structor(3, '00002', 'd2z/dp2: second order momentum compaction (related by negative sign)')
terms(28) = terms_structor(3, '00003', 'd3z/dp3: third order momentum compaction (related by negative sign)')
terms(29) = terms_structor(3, '00004', 'd3z/dp3: third order momentum compaction (related by negative sign)')

!--------------------------------------------------------------------
  !p => r%start
  !do i=1,100
  !  if(p%mag%p%nmul>=2) then
  !    write(*,*) "Knob member: ", p%mag%name
  !    call make_it_knob(p%magp%bn(2), 1)  !all these magnets are on knob 1
  !  endif
  !  p => p%next
  !enddo 
  ids = 1
  x=0.d0 
  pos_a = 1
  pos_b = 500 ! r%n+pos_a for 1 turn
  call find_orbit_x(r,x,state,1.0d-6,fibre1=pos_a)
  ray=0  
  ray=x  
  ray8=ray+ids  ! ray8 is a probe_8
  !call track_probe(r,ray8,state,fibre1=pos_a)
  call track_probe(r,ray8,+state,fibre1=pos_a)

  cmap=ray8  ! c_damap <- probe_8
  call c_normal(cmap,cnorm,phase=tunes)  ! c_normal_form <- c_damap
  call c_canonise(cnorm%a_t,acs,a0=a0,a1=a1)
  acs = a0*a1
  !ray8 = ray + (cnorm%a_t.sub.1)
  ray8 = ray + acs
  !call track_probe(r,ray8,state,fibre1=pos_a,fibre2=pos_b)
  call track_probe(r,ray8,+state,fibre1=pos_a,fibre2=pos_b)
  acs = ray8
  !cmap = (cmap.sub.1)**(-1) * cmap

  !cmap=(cnorm%a_t.sub.1)**(-1)*cmap*(cnorm%a_t.sub.1)
  !test cmap=from_phasor(-1)*cmap*from_phasor(1)

  call c_full_factorise(acs, a0=a0, a1=a1, a2=a2)
  !call print(a2)

cmap=a2
!cmap=(cmap.sub.1)**(-1)*cmap 
  F=log(cmap)
  h=getpb(F)

h=h*c_phasor()
!call print(h)
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!call print(tunes(1))
write(*,*) "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"

  F=c_phasor()*F

  allocate(test_i(4))
  !allocate(test_i(6))
  test_i = 0 
  test_i(1) = 1
  test_i(4) = 2
  h=h.par.test_i
  call print(h)
  call c_identify_resonance(test_i,nj,test_c)
  test_t = (F%v(nj).par.test_i)*test_c
  !call print(test_t)
stop
  F=-getvf(h)
cmap=exp(f,cmap)
cmap=cmap.cut.c_%no
call print(cmap)
 !call check_map(cmap,ch,norm,nm)
write(6,*) norm,nm
stop
  ray=0  
  ray=x  
  ray8=ray+ids
p=>r%start
do k=1,r%n
!if(p%mag%kind/=55) then
  call track_probe(r,ray8,state,fibre1=k,fibre2=k+1)
!endif
cmap=ray8
  ! ray=ray8
  ! ray8=ray+ids
call check_map(cmap,ch,norm,nm)
if(ch.or.p%mag%kind==55) then
write(6,*) k,p%mag%name,norm,ch
endif
p=>p%next
enddo
stop 100

stop
  cmap=from_phasor(-1)*a2*from_phasor(1)

  !call c_factor_map(cmap,L,F,1)
  F=log(cmap)
  call c_clean_vector_field(F,F,prec)

  h=cgetpb(F,S2=-1)   ! cgetpb includes 1/(-2*i_) 
  call print(h,6)
 
write(6,*) (h.sub.'21001')
write(6,*) (h.sub.'12001')
stop
  !test begin
  allocate(test_i(4))
  test_i = 0 
  test_i(1) = 1
  test_i(4) = 2
  write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
  call print(h,6)
  write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
  h=h.par.test_i
  call print(h,6)
  write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
  call c_identify_resonance(test_i,nj,test_c)
  test_t = (F%v(nj).par.test_i)*test_c
  write(*,*) "FOO nj: ", nj
  write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
  call print(F%v(nj),6)
  write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
  call print(test_t,6)
  write(*,*) "++++++++++++++++++++++++++++++++++++++++++++++++++"
  stop
  !test end

  !          call clean(tunes,tunes,prec)
  !          call print(tunes(1),6)
  cmap_test = 1
  j = 0
  if(c_%ndpt .ne. 0) j=1
  do i=1,c_%nd-j
    cmap_test%v(2*i-1) = cmap_test%v(2*i-1)*2
    cmap_test%v(2*i)   = 1
  enddo
  tunes(1) = tunes(1).o.cmap_test
  !          call print(tunes(1),6)

  write(*,*)
  write(*,*)
  do i=1, size(terms)
    !call process_term(h, normal%dhdj%v, termdefs(i)%F_index, termdefs(i)%name, termdefs(i)%c_val)
    call process_term(h, tunes, terms(i))
    write(*,'(i4,"   ",A,"   (",f14.5,", ",f14.5,")   ",a)') terms(i)%F_index, terms(i)%name, terms(i)%c_val, terms(i)%desc
  enddo
  write(*,*)
  write(*,*)

999 continue

call ptc_end(graphics_maybe=1,flat_file=.false.)

contains
subroutine check_map(m,ch,n,nm)
implicit none
type(c_damap) m
type(c_taylor) h
integer i,j
real(dp) n,nm
logical ch
call alloc(h)
ch=.false.
n=0
nm=0
do i=1,c_%nd2
nm=nm+full_abs(m%v(i))
do j=1,c_%nd2
h= (m%v(i).pb.m%v(j))
h=h.cut.c_%no
 n=n+full_abs(h)
 enddo
enddo
n=abs(n-c_%nd2)
nm=abs(nm-c_%nd2)
nm=n/nm

if(n>1.d-15) ch=.true.
 call kill(h)
end subroutine check_map

  function terms_structor(F_index, term_name, desc)
    implicit none

    integer F_index
    character(*) term_name
    character(*) desc

    integer, allocatable :: j_local(:)

    type(terms_struct) terms_structor

    allocate(terms_structor%j(len(term_name)))

    terms_structor%F_index = F_index
    terms_structor%name = term_name
    do i=1, len(term_name)
      read(term_name(i:i),'(i)') terms_structor%j(i)
    enddo
    terms_structor%c_val = (0.0, 0.0)
    terms_structor%desc = desc

  end function

  subroutine process_term(h_local, tunes, term)
    implicit none

    type(c_taylor) h_local
    type(c_taylor) tunes(:)
    !type(taylor) v_local(ndim2)
    type(terms_struct) term

    !working space
    type(taylor) t_local
    type(c_taylor) c_t_local

    if( term%F_index .gt. 0 ) then
      call alloc(t_local)
      t_local = tunes(term%F_index)
      t_local = t_local.par.term%name
      t_local=t_local<=6  ! shift indexes down by 6, which moves the knob index 7 into
                          ! the first position.
      term%c_val = cmplx(t_local.sub.'0', 0.0)

      call kill(t_local)
    else
      call alloc(c_t_local)
      c_t_local=h_local.par.term%name  ! pull desired polynomial coefficient from h.
      !write(*,*) "FOO name: ", term%name
      !call print(c_t_local,6)
      !write(*,*) "+++++++++++++++++++++++++++++++++++++++++++++++"

      c_t_local=c_t_local<=6  ! shift indexes down by 5, which moves the knob index 6 into
                              ! the first position.

      term%c_val = c_t_local.sub.'0'

      call kill(c_t_local)
    endif

  end subroutine process_term

end program ptc_et








