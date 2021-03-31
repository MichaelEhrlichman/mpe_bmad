module ptc_mod

use pointer_lattice

implicit none

type :: terms_struct
  integer F_index  !used to index complex vector field.  Only for dhdj
  integer j(6)
  character(6) name
  character(100) desc
  complex(dp) c_val
end type

type(terms_struct) :: terms(29) = (/ &
                                  terms_struct(0, (/2,1,0,0,0,0/), '210000', 'Qx', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/3,0,0,0,0,0/), '300000', '3Qx', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/1,0,1,1,0,0/), '101100', 'Qx', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/1,0,0,2,0,0/), '100200', 'Qx-2Qy', (0.0d0,0.0d0)),&
                                  terms_struct(0, (/1,0,2,0,0,0/), '102000', 'Qx+2Qy', (0.0d0,0.0d0)),&
                                  terms_struct(0, (/2,0,0,0,1,0/), '200010', 'Synchro-betatron resonances ???', (0.0d0,0.0d0)),&
                                  terms_struct(0, (/0,0,2,0,1,0/), '002010', 'Momentum-dependence of beta functions ???', (0.0d0,0.0d0)),&
                                  terms_struct(0, (/1,0,0,0,2,0/), '100020', '2nd order dispersion', (0.0d0,0.0d0)),&
                                  terms_struct(0, (/3,1,0,0,0,0/), '310000', '2Qx', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/0,0,3,1,0,0/), '003100', '2Qy', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/1,1,2,0,0,0/), '112000', '2Qy', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/2,0,1,1,0,0/), '201100', '2Qx', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/2,2,0,0,0,0/), '220000', 'Horizontal ADTS', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/1,1,1,1,0,0/), '111100', 'Cross ADTS', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/0,0,2,2,0,0/), '002200', 'Vertical ADTS', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/2,0,2,0,0,0/), '202000', '2Qx+2Qy', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/2,0,0,2,0,0/), '200200', '2Qx-2Qy', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/4,0,0,0,0,0/), '400000', '4Qx', (0.0d0,0.0d0)), &
                                  terms_struct(0, (/0,0,4,0,0,0/), '004000', '4Qy', (0.0d0,0.0d0)), &
                                  terms_struct(1, (/0,0,0,0,1,0/), '000010', 'chrom x', (0.0d0,0.0d0)), &
                                  terms_struct(1, (/0,0,0,0,2,0/), '000020', '2nd order chrom x', (0.0d0,0.0d0)), &
                                  terms_struct(1, (/0,0,0,0,3,0/), '000030', '3rd order chrom x', (0.0d0,0.0d0)), &
                                  terms_struct(2, (/0,0,0,0,1,0/), '000010', 'chrom y', (0.0d0,0.0d0)), &
                                  terms_struct(2, (/0,0,0,0,2,0/), '000020', '2nd order chrom y', (0.0d0,0.0d0)), &
                                  terms_struct(2, (/0,0,0,0,3,0/), '000030', '3rd order chrom y', (0.0d0,0.0d0)), &
                                  terms_struct(3, (/0,0,0,0,1,0/), '000010', 'dz/dp: first order momentum compaction (related by negative sign)', (0.0d0,0.0d0)), &
                                  terms_struct(3, (/0,0,0,0,2,0/), '000020', 'd2z/dp2: second order momentum compaction (related by negative sign)', (0.0d0,0.0d0)), &
                                  terms_struct(3, (/0,0,0,0,3,0/), '000030', 'd3z/dp3: third order momentum compaction (related by negative sign)', (0.0d0,0.0d0)), &
                                  terms_struct(3, (/0,0,0,0,4,0/), '000040', 'd3z/dp3: third order momentum compaction (related by negative sign)', (0.0d0,0.0d0)) &
                                  /)

  type(internal_state), target :: state 

contains

subroutine my_ptc_parse(ptc_lat,r,map_order)
  !parses ptc flat file ptc_lat into global data structure m_u from module pointer lattice.

  character(*) ptc_lat
  type(layout), pointer :: r
  integer map_order

  real(dp) thi
  integer np

  ! Set global parameters
  !state = default0
  !state=delta0
  state = default0 + nocavity0
  ndpt_bmad = 1

  use_info=.true.
  my_estate=>state

  !----------------------------------------------Initialize PTC and read in universe
  call ptc_ini_no_append
  call read_and_append_virgin_general(m_u,ptc_lat)

  c_verbose=.false.

  r=>m_u%start

  !thi=-1 !slow
  ! thi=0.01d0
  ! CALL thin_lens_resplit(r,thi)

  !map_order=3
  np=0 !number of parameters
  call init(state,map_order,np)
end subroutine

subroutine ptc_one_turn_taylor(r,pos,one_turn_ray8,co_pr)
  implicit none

  type(layout) r
  integer pos
  type(probe_8) one_turn_ray8
  type(probe) co_pr

  type(c_damap) ids

  real(dp) co_dp(6)

  call alloc(ids)
  call alloc(one_turn_ray8)

  ids = 1  !c_damap to identity map.
  co_dp=0.d0   ! real(rp)*6

  call find_orbit_x(r,co_dp,state,1.0d-6,fibre1=pos) !find closed orbit at pos_a, store as rp*6 in co_dp.
  co_pr=0   !co_pr is probe
  co_pr=co_dp   !probe <= real(rp)*6.  Put the closed orbit into a probe.
  one_turn_ray8=co_pr+ids  ! probe_8 <= probe + c_damap
  call track_probe(r,one_turn_ray8,state,fibre1=pos)  !track the identity map around the closed orbit, to get map about closed orbit.
  !one_turn_ray8 now contains the 1-turn taylor map about co at pos_a

  call kill(ids)
end subroutine

subroutine ptc_populate_terms(playout,one_turn_ray8,co_pr)

  implicit none

  type(layout) playout

  !PTC datatypes - Not allocatable
  type(probe) co_pr

  !PTC datatypes - Must be allocated
  type(c_damap) cmap
  type(c_vector_field) F
  type(c_normal_form) cnorm
  type(c_taylor) tunes_1turn(3), tunes_atob(3)
  type(probe_8) one_turn_ray8
  type(probe_8) ray8
  type(c_damap) a0, a1, a2
  type(c_damap) acs_full, acs1_a, acs1_b

  integer i, j, k
  integer pos_a, pos_b

  real(dp) co_dp(6)

  call alloc(tunes_1turn)
  call alloc(tunes_atob)
  call alloc(a0,a1,a2,acs_full)
  call alloc(ray8)
  call alloc(acs1_a,acs1_b)
  call alloc(cmap)
  call alloc(cnorm)
  call alloc(F)

  !--------------------------------------------------------------------
  pos_a = 1
  pos_b = playout%n+pos_a !playout%n+pos_a for 1-turn

  cmap=one_turn_ray8  ! c_damap <- probe_8 
  !looking at output of 'call print(one_turn_ray8%x)' and 'call print(cmap)', 
  !cmap contains the same numbers as one_turn_ray8%x, but as complex numbers (with imag part zero).
  call c_normal(cmap,cnorm,phase=tunes_1turn,dospin=my_false,no_used=1)  ! c_normal_form <- c_damap
                                               ! get normal form for 1-turn map.
  call c_canonise(cnorm%a_t,acs_full,a0=a0,a1=a1, a2=a2)  ! Enforce CS on normalising map.

  acs1_a = a0*a1  ! acs1_a is first order of normalizing map.

  ray8 = co_pr + acs1_a 
  !call track_probe(playout,ray8,state,fibre1=pos_a,fibre2=pos_b) !Track acs1_a from point a to point b.
  call track_probe(playout,ray8,state,fibre1=pos_a) !Track acs1_a from point a to point b.
  acs1_b = ray8  !c_damap <- probe_8

  call c_canonise(acs1_b,acs1_b, a0=a0, a1=a1, a2=a2, phase=tunes_atob)
  !write(*,*) "tune 1: ", real((tunes_atob(1).par.'000000').sub.'0')
  !write(*,*) "tune 2: ", real((tunes_atob(2).par.'000000').sub.'0')
  !write(*,*) "tune 3: ", real((tunes_atob(3).par.'000000').sub.'0')

  F=log(a2) ! c_vector_field <- c_damap
  F=c_phasor()*F
  ! write(*,*) "********************** FOO F *********************"
  ! call print(F) !FOO
  ! write(*,*) "********************** FOO *********************"

  !call clean(t_local,t_local,1.0d-6)
  !call print(tunes_atob)

  do i=1, size(terms)
    call process_term(F, tunes_atob, terms(i))
    !write(*,*) "FOO: ", i, terms(i)%c_val
  enddo
  !--------------------------------------------------------------------

  call kill(tunes_1turn)
  call kill(tunes_atob)
  call kill(a0,a1,a2,acs_full)
  call kill(acs1_a)
  call kill(acs1_b)
  call kill(cmap)
  call kill(cnorm)
  call kill(F)

end subroutine

!+
!
!-
subroutine ptc_print_terms()
  implicit none

  integer i
  do i=1, size(terms)
    write(*,'(i4,"   ",A,"   (",f20.10,", ",f20.10,")   ",f20.10,"   ",a)') terms(i)%F_index, terms(i)%name, terms(i)%c_val, abs(terms(i)%c_val), terms(i)%desc
  enddo
end subroutine

!+
!
!-
subroutine process_term(F, tunes, term)

  implicit none

  type(c_vector_field) F
  type(c_taylor) tunes(:)
  type(c_taylor) t_local
  type(terms_struct) term

  !working space
  complex(dp) res_c
  integer nj
  integer j_temp(6)

  if( term%F_index .gt. 0 ) then
    call alloc(t_local)
    t_local = tunes(term%F_index)
    t_local = t_local.par.term%j
    term%c_val = t_local.sub.'0'

    call kill(t_local)
  else
    j_temp = term%j
    call c_identify_resonance(j_temp,nj,res_c)
    term%c_val = (F%v(nj).par.j_temp)*res_c
  endif
end subroutine process_term

! !+
! !
! !-
! 
! function terms_structor(F_index, term_name, desc)
!   implicit none
! 
!   integer F_index
!   character(*) term_name
!   character(*) desc
! 
!   integer, allocatable :: j_local(:)
! 
!   type(terms_struct) terms_structor
! 
!   allocate(terms_structor%j(len(term_name)))
! 
!   terms_structor%F_index = F_index
!   terms_structor%name = term_name
!   do i=1, len(term_name)
!     read(term_name(i:i),'(i)') terms_structor%j(i)
!   enddo
!   terms_structor%c_val = (0.0, 0.0)
!   terms_structor%desc = desc
! 
! end function

end module








