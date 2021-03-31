!+
!-
module momentum_aperture_mod

use bmad

use touschek_mod, only: momentum_aperture_struct

implicit none

contains

subroutine momentum_aperture_one(ring,co,n_turns,dims,max_deltam,verbose)

implicit none

logical lost
logical, optional :: verbose
logical err

real(rp) delta_m, accuracy
real(rp) target_accuracy, first_guess
real(rp) rf_freq

integer track_state
integer ix, i, k, j, o
integer n_turns, dims, n_cavs

type(momentum_aperture_struct) max_deltam
type(lat_struct) ring
type(coord_struct) co(0:)
type(coord_struct) co_inj, orb1

type(ele_pointer_struct), allocatable :: rfcav_eles(:)

if(dims .eq. 6) then
  call lat_ele_locator("rfcav::*",ring,rfcav_eles,n_cavs)
  if( n_cavs .lt. 1) then
    write(*,*) "FATAL:  dims is 6 but no rf cavities found"
    call err_exit()
  endif
  rf_freq = rfcav_eles(1)%ele%value(rf_frequency$)
  deallocate(rfcav_eles)
endif

!------------------------------Set default parameters
first_guess = 0.02
target_accuracy = 0.00005

!-----------------------------Preparation is complete.  Begin computing momentum apertures.

call twiss_and_track_at_s(ring,max_deltam%s,orb_at_s=co_inj,orb=co)
do k= 1, 2
  !If k=1, then pos aperture.  if k=2, then neg aperture.
  !The following loop conducts a binary search for the largest energy kick that can be given to a particle
  !such that it is lost.
  !This loop first brackets the largest energy kick, then uses a binary search to obtain the
  !desired precision.

  !Initialization
  call binary_search(.false.,0._rp,accuracy,.true.) !calling binary_search with last arg .true. resets its state variables
  delta_m = first_guess

  do while (accuracy .gt. target_accuracy)
    orb1 = co_inj                               ! start from the closed orbit
    orb1%vec(1) = orb1%vec(1) + 1.0e-6
    orb1%vec(3) = orb1%vec(3) + 1.0e-6   
    orb1%vec(6) = orb1%vec(6) + ((-1)**(k+1))*delta_m   ! add the momentum kick

    do j=1,n_turns
      call track_from_s_to_s(ring, max_deltam%s, max_deltam%s, orb1, orb1, track_state=track_state)
      if( dims == 6 ) then
        if( abs(orb1%vec(5)) .gt. c_light/rf_freq/2.0 ) then
          track_state = lost$
        endif
      endif
      if ( track_state .ne. moving_forward$ ) exit
    enddo
    lost = (track_state .ne. moving_forward$)
    call binary_search(lost, delta_m, accuracy, .false.) !updates delta_m according to one iteration of binary search
  enddo

  if(k==1) then
    max_deltam%pos = ((-1)**(k+1))*delta_m
  else
    max_deltam%neg = ((-1)**(k+1))*delta_m
  endif
enddo

end subroutine

!+
!-

subroutine divergence_aperture_one(ring,co,n_turns,dims,max_delta_xp,verbose)

implicit none

logical lost
logical, optional :: verbose
logical err

real(rp) delta_xp
real(rp) accuracy
real(rp) target_accuracy
real(rp) first_guess
real(rp) rf_freq

integer track_state
integer ix, i, k, j, o
integer n_turns
integer dims
integer n_cavs

type(momentum_aperture_struct) max_delta_xp
type(lat_struct) ring
type(coord_struct) co(0:)
type(coord_struct) co_inj
type(coord_struct) orb1

type(ele_pointer_struct), allocatable :: rfcav_eles(:)

if(dims .eq. 6) then
  call lat_ele_locator("rfcav::*",ring,rfcav_eles,n_cavs)
  if( n_cavs .lt. 1) then
    write(*,*) "FATAL:  dims is 6 but no rf cavities found"
    call err_exit()
  endif
  rf_freq = rfcav_eles(1)%ele%value(rf_frequency$)
  deallocate(rfcav_eles)
endif

!------------------------------Set default parameters
first_guess = 0.02
target_accuracy = 0.00005  ! 0.01 is accurate to 1 part in a hundred

!-----------------------------Preparation is complete.  Begin computing momentum apertures.

call twiss_and_track_at_s(ring,max_delta_xp%s,orb_at_s=co_inj,orb=co)
do k= 1, 2
  !If k=1, then pos aperture.  if k=2, then neg aperture.
  !The following loop conducts a binary search for the largest horizontal div shange can be given to a particle
  !such that it is lost.
  !This loop first brackets the largest change, then uses a binary search to obtain the
  !desired precision.

  !Initialization
  call binary_search(.false.,0._rp,accuracy,.true.) !calling binary_search with last arg .true. resets its state variables
  delta_xp = first_guess

  do while (accuracy .gt. target_accuracy)
    orb1 = co_inj                               ! start from the closed orbit
    orb1%vec(2) = orb1%vec(2) + ((-1)**(k+1))*delta_xp   ! add the momentum kick

    do j=1,n_turns
      call track_from_s_to_s(ring, max_delta_xp%s, max_delta_xp%s, orb1, orb1, track_state=track_state)
      if( dims == 6 ) then
        if( abs(orb1%vec(5)) .gt. c_light/rf_freq/2.0 ) then
          track_state = lost$
        endif
      endif
      if ( track_state .ne. moving_forward$ ) exit
    enddo

    lost = (track_state .ne. moving_forward$)

    call binary_search(lost, delta_xp, accuracy, .false.) !updates delta_m according to one iteration of binary search
  enddo

  if(k==1) then
    max_delta_xp%pos = ((-1)**(k+1))*delta_xp
  else
    max_delta_xp%neg = ((-1)**(k+1))*delta_xp
  endif
enddo

end subroutine

end module




