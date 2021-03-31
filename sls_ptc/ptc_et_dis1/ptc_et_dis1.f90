program ptc_et_dis1

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
type(probe) co_pr
type(internal_state), target :: state

!PTC datatypes - Must be allocated
type(c_damap) cmap
type(c_damap) a_ph
type(c_damap) region_map
type(c_vector_field) F
type(c_normal_form) cnorm
type(c_damap) ids
type(c_taylor) tunes(3), tunes2(3)
type(probe_8) ray8
type(c_ray) x1_ray, x2_ray
type(c_damap) a0, a1, a2
type(c_damap) acs_full, acs_a, acs_b
integer i, j, k
integer map_order,np
integer pos_a, pos_b
integer n_terms
integer flag
character(1) flag_chr

real(dp) nu_x, nu_y
real(dp) co_dp(6)
real(dp) prec,thi,emi(2)

! Set global parameters
use_info=.true.
my_estate=>state

! Set local parameters
prec=1.d-8
call getarg(1,lat_file)

!call getarg(2,flag_chr)
!read(flag_chr,*) flag

!----------------------------------------------Initialize PTC and read in universe
call ptc_ini_no_append
call read_and_append_virgin_general(M_U,lat_file)

c_verbose=.false.

!n_cai = (0.0d0, -1.0d0) !FOO

r=>m_u%start

!thi=-1
!thi=0.01d0
!CALL THIN_LENS_resplit(r,thi)
!state=delta0
state = default + nocavity0
!state = only_4d
!state = default

!use_complex_in_ptc=my_true

map_order=4
np=0 !number of parameters
call init(state,map_order,np)

call alloc(tunes)
call alloc(tunes2)
call alloc(acs_full, acs_a, acs_b)
call alloc(a0, a1, a2)
call alloc(ids)
call alloc(ray8)
call alloc(cmap)
call alloc(cnorm)
call alloc(F)
call alloc(a_ph)

x1_ray = 0 !FOO
x2_ray = 0 !FOO

n_terms = 29
  
allocate(terms(n_terms))
!                             'x px y py delta z'
terms(1)  = terms_structor(0, '210000', 'Qx')
terms(2)  = terms_structor(0, '300000', '3Qx')
terms(3)  = terms_structor(0, '101100', 'Qx')
terms(4)  = terms_structor(0, '100200', 'Qx-2Qy')
terms(5)  = terms_structor(0, '102000', 'Qx+2Qy')
terms(6)  = terms_structor(0, '200010', 'Synchro-betatron resonances ???')
terms(7)  = terms_structor(0, '002010', 'Momentum-dependence of beta functions ???')
terms(8)  = terms_structor(0, '100020', '2nd order dispersion')

terms(9)  = terms_structor(0, '310000', '2Qx')
terms(10) = terms_structor(0, '003100', '2Qy')
terms(11) = terms_structor(0, '112000', '2Qy')
terms(12) = terms_structor(0, '201100', '2Qx')

terms(13) = terms_structor(1, '110000', 'dQx/dJx')
terms(14) = terms_structor(2, '001100', 'dQy/dJy')
terms(15) = terms_structor(1, '001100', 'dQx/dJy')

terms(16) = terms_structor(0, '202000', '2Qx+2Qy')
terms(17) = terms_structor(0, '200200', '2Qx-2Qy')
terms(18) = terms_structor(0, '400000', '4Qx')
terms(19) = terms_structor(0, '004000', '4Qy')

terms(20) = terms_structor(1, '000010', 'chrom x')
terms(21) = terms_structor(1, '000020', '2nd order chrom x')
terms(22) = terms_structor(1, '000030', '3rd order chrom x')
terms(23) = terms_structor(2, '000010', 'chrom y')
terms(24) = terms_structor(2, '000020', '2nd order chrom y')
terms(25) = terms_structor(2, '000030', '3rd order chrom y')
terms(26) = terms_structor(3, '000010', 'dz/dp: first order momentum compaction (related by negative sign)')
terms(27) = terms_structor(3, '000020', 'd2z/dp2: second order momentum compaction (related by negative sign)')
terms(28) = terms_structor(3, '000030', 'd3z/dp3: third order momentum compaction (related by negative sign)')
terms(29) = terms_structor(3, '000040', 'd3z/dp3: third order momentum compaction (related by negative sign)')

!--------------------------------------------------------------------
  ids = 1  !c_damap to identity map.
  co_dp=0.d0   ! real(rp)*6
  pos_a = 1
  pos_b = r%n+pos_a !r%n+pos_a for 1 turn
  call find_orbit_x(r,co_dp,state,1.0d-6,fibre1=pos_a) !find closed orbit at pos_a, store as rp*6 in co_dp.
  co_pr=0   !co_pr is probe
  co_pr=co_dp   !probe <= real(rp)*6.  Put the closed orbit into a probe.
  ray8=co_pr+ids  ! probe_8 <= probe + c_damap
  call track_probe(r,ray8,state,fibre1=pos_a)

  cmap=ray8  ! c_damap <- probe_8 
  !looking at output of 'call print(ray8%x)' and 'call print(cmap)', 
  !cmap contains the same numbers as ray8%x, but as complex numbers (with imag part zero).

  call c_normal(cmap,cnorm,phase=tunes)  ! c_normal_form <- c_damap
                                         ! get normal form for 1-turn map.
  call c_canonise(cnorm%a_t,acs_full,a0=a0,a1=a1)  ! Enforce CS on 1-turn map.

  !For given physical coord
 ! a_ph = ci_phasor()*cnorm%a_t**(-1)
 ! x1_ray%x(1) = 0.001 !mm relative to CO
 ! x2_ray = a_ph.o.x1_ray
 !prec=1.e-8

 ! call clean(tunes(1),tunes(1),prec)
 ! call clean(tunes(2),tunes(2),prec)

  !emi(1)=1.d-6
  !emi(2)=1.d-6
  !region_map=1
  !do i=1,2 
  !  region_map%v(2*i-1)=  sqrt(emi(i)).cmono.2*i-1
  !  region_map%v(2*i)=  sqrt(emi(i)).cmono.2*i
  !enddo
  !region_map%v(5)=(1.d-3).cmono.5
  !call print(tunes(1:2))

 ! tunes(1)=tunes(1).o.region_map
 ! tunes(2)=tunes(2).o.region_map
 ! call clean(tunes(1),tunes(1),prec)
 ! call clean(tunes(2),tunes(2),prec)
 ! call print(tunes(1:2))
 ! 
 ! stop
 
  ! write(*,'(a,8es14.5)') "FOO x2_ray%x(1:4): ", x2_ray%x(1:4)
  ! write(*,*) x2_ray%x(1)*x2_ray%x(2), x2_ray%x(3)*x2_ray%x(4)
  ! nu_x = tunes(1).o.x2_ray
  ! nu_y = tunes(2).o.x2_ray
  ! write(*,*) "FOO nu_x0, nu_y0: ", cnorm%tune(1), cnorm%tune(2)
  ! write(*,*) "FOO nu_x, nu_y:   ", nu_x, nu_y
  ! stop !FOO

  !acs_a contains those elements of acs_full that are linear in the first
  !three indices.  So perhaps acs_a is the linear normalizing transformation, while
  !acs_full is the full nonlinear normalizing transformation.
  acs_a = a0*a1

  !We can track an acs_a already in hand from the location it was generated to some other location.
  ray8 = co_pr + acs_a
  call track_probe(r,ray8,state,fibre1=pos_a,fibre2=pos_b) !Track the CS linearizing A1 from point a to point b.
  acs_b = ray8  !c_damap <- probe_8

  call c_canonise(acs_b,acs_b,phase=tunes2)
  !call print(tunes2(1)) !FOO
  !call print(tunes2(2)) !FOO

  call c_full_factorise(acs_b, a0=a0, a1=a1, a2=a2)

  F=log(a2) ! c_vector_field <- c_damap
  F=c_phasor()*F

  write(*,*)
  write(*,*)
  do i=1, size(terms)
    call process_term(F, tunes, terms(i))
    write(*,'(i4,"   ",A,"   (",f20.10,", ",f20.10,")   ",a)') terms(i)%F_index, terms(i)%name, terms(i)%c_val, terms(i)%desc
  enddo
  write(*,*)
  write(*,*)

contains
  subroutine process_term(F, tunes, term)
    implicit none

    type(c_vector_field) F
    type(c_taylor) tunes(:)
    type(c_taylor) t_local
    type(terms_struct) term

    !working space
    complex(dp) res_c
    integer nj

    if( term%F_index .gt. 0 ) then
      call alloc(t_local)
      t_local = tunes(term%F_index)
      t_local = t_local.par.term%j
      term%c_val = t_local.sub.'0'

      call kill(t_local)
    else
      call c_identify_resonance(term%j,nj,res_c)
      term%c_val = (F%v(nj).par.term%j)*res_c
    endif

  end subroutine process_term

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

end program ptc_et_dis1








