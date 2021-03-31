module dynap_mod

use bmad, only: rp
use moga_struct_mod

implicit none

contains

subroutine get_magnet_strengths(mags,ring,strengths)
  use bmad

  implicit none

  type(mag_struct) mags(:)
  type(lat_struct) ring
  real(rp) strengths(:)

  integer n_mags, n_loc
  type (ele_pointer_struct), allocatable :: eles(:)
  logical err
  integer i

  n_mags = size(mags)

  do i=1, n_mags
    if(mags(i)%name == '') exit
    call lat_ele_locator(mags(i)%name, ring, eles, n_loc, err)
    if(err .or. (n_loc .lt. 1)) then
      write(*,*) "Get ele attribute error.  Terminating."
      call err_exit
    endif
    strengths(i) = eles(1)%ele%value(k1$)
  enddo

  if(allocated(eles)) deallocate(eles)
end subroutine

subroutine set_magnet_strengths(mags,ring,strengths)
  use bmad

  implicit none

  type(mag_struct) mags(:)
  type(lat_struct) ring
  real(rp) strengths(:)

  integer n_mags, n_loc
  type (ele_pointer_struct), allocatable :: eles(:)
  logical err
  character*21 var_str
  character*50 set_str
  integer i,j

  n_mags = size(mags)

  do i=1, n_mags
    if(mags(i)%name == '') exit
    call lat_ele_locator(mags(i)%name, ring, eles, n_loc, err)
    if(n_loc == 0) then
      write(*,*) "Magnet ", mags(i)%name, " not found!"
      call err_exit
    endif
    !write(var_str,'(f21.11)') strengths(i)
    write(var_str,'(es21.11)') strengths(i)
    set_str = trim(adjustl(mags(i)%property))//'='//trim(adjustl(var_str))

    do j=1, n_loc
      call set_ele_attribute (eles(j)%ele, set_str, err)
      if(err) then
        write(*,*) "Set ele attribute error.  Terminating.", set_str, eles(j)%ele%name, strengths(i)
        !error stop
        call err_exit
      endif
    enddo
  enddo

  if(allocated(eles)) deallocate(eles)
  call lattice_bookkeeper(ring)
end subroutine

function is_linear_ele(ele)
  use bmad
  implicit none
  logical is_linear_ele
  type(ele_struct) ele

  if( (ele%key == sbend$) .or. &
      (ele%key == sextupole$) .or. &
      (ele%key == multipole$) ) then
    is_linear_ele = .false.
  else
    is_linear_ele = .true.
  endif
end function

function frac(x)
  use bmad, only: rp
  real(rp) frac, x
  frac = x - int(x,kind(x))
end function


end module
