!+
! Subroutine make_mat6_custom (ele, param, start_orb, end_orb, err_flag)
!
! Prototype routine for custom tracking. 
! To use, see the Bmad manual.
!
! General rule: Your code may NOT modify any argument that is not listed as
! an output agument below."
!
! Input:
!   ele       -- Ele_struct: Element with transfer matrix
!   param     -- lat_param_struct: Parameters are needed for some elements.
!   start_orb -- Coord_struct: Coordinates at the beginning of element. 
!
! Output:
!   ele       -- Ele_struct: Element with transfer matrix.
!     %mat6     -- 6x6 transfer matrix.
!   end_orb   -- Coord_struct: Coordinates at the end of element.
!   err_flag  -- Logical: Set true if there is an error. False otherwise.
!+

subroutine make_mat6_custom (ele, param, start_orb, end_orb, err_flag)

use bmad_interface

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: start_orb, end_orb
type (lat_param_struct)  param
type(all_pointer_struct) v_ptr, s_ptr

logical :: err_flag
character(32) :: r_name = 'make_mat6_custom'

real(rp) voltage_nominal, slope, voltage

integer i

!

call pointer_to_attribute(ele,'VOLTAGE', .true., v_ptr, err_flag)
call pointer_to_attribute(ele,'SLOPE', .true., s_ptr, err_flag)
voltage_nominal = v_ptr%r
slope = s_ptr%r

err_flag = .false.

ele%mat6 = 0.0d0
do i=1,6
  ele%mat6(i,i) = 1.0d0
enddo
ele%vec0 = 0.0d0
ele%vec0(6) = voltage_nominal

end subroutine
