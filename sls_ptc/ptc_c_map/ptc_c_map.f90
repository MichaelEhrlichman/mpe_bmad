program ptc_c_map

use pointer_lattice

implicit none

character(100) lat_file

type(layout), pointer :: r

!PTC datatypes - Not allocatable
type(probe) co_pr
type(internal_state), target :: state

!PTC datatypes - Must be allocated
type(c_damap) cmap
type(c_normal_form) cnorm
type(c_damap) ids, N, A, M
type(probe_8) ray8
type(c_damap) a0, a1, a2
integer i, j, k
integer map_order,np
integer pos_a, pos_b
integer flag
character(1) flag_chr

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

r=>m_u%start

!thi=-1
!thi=0.01d0
!CALL THIN_LENS_resplit(r,thi)
!state=delta0
!state = default + nocavity0
state = only_4d
!state = default

use_complex_in_ptc=my_true

map_order=15
np=0 !number of parameters
call init(state,map_order,np)

call alloc(ids)
call alloc(ray8)
call alloc(cmap)
call alloc(cnorm)
call alloc(A)
call alloc(M)
call alloc(N)

!                             'x px y py delta z'
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

  call c_normal(cmap,cnorm)  ! c_normal_form <- c_damap

  A = cnorm%atot.sub.1  !get first order of normalizing map
  M = A**(-1)*cmap*A  !map to floquet variables
  N = from_phasor(-1)*M*from_phasor(1)

  call clean(N,N,1d-8)
  call print(N)

end program








