!+
! Subroutine make_mat6_custom (ele, param, start_orb, end_orb, err_flag)
!
! Dummy routine for custom tracking. 
! This routine needs to be replaced for a custom calculation.
! If not replaced, and this routine is called, this routine will generate an error message.
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

use bmad_struct
use bmad_interface, except_dummy => make_mat6_custom

implicit none

type (ele_struct), target :: ele
type (coord_struct) :: start_orb, end_orb
type (lat_param_struct)  param

real(rp) m11, m12, m21, m22
real(rp) L0, K_und, B, vbar, gamma_rel
real(rp) lambda, k1

integer i, j

logical :: err_flag
character(32) :: r_name = 'make_mat6_custom'

L0 = ele%value(n_pole$)*ele%value(l_pole$) !undulator length
lambda = ele%value(l_pole$)*2 !undulator period
B = ele%value(b_max$) !peak field
gamma_rel = sqrt((start_orb%p0c*(1.+start_orb%vec(6)))**2+m_electron**2)/m_electron !relativistic gamma

K_und = e_charge*B*lambda/(twopi*m_electron*e_charge/c_light) !undulatr parameter
vbar = c_light*(1.-1./2./gamma_rel/gamma_rel*(1.+K_und*K_und)) !average longitudinal velocity

k1 = 1./2.*((e_charge*B)/(vbar*gamma_rel*m_electron*e_charge/c_light**2))**2 !effective focusing strength


!transfer matrix - initially just zeros
do i = 1,6
  do j= 1,6
    ele%mat6(i,j)=0
  end do
end do

m11 = cos(sqrt(k1)*L0)
m12 = 1./sqrt(k1)*sin(sqrt(k1)*L0)
m21 = -sqrt(k1)*sin(sqrt(k1)*L0)
m22 = cos(sqrt(k1)*L0)

ele%mat6(1,1) = m11
ele%mat6(1,2) = m12
ele%mat6(2,1) = m21
ele%mat6(2,2) = m22
ele%mat6(3,3) = m11
ele%mat6(3,4) = m12
ele%mat6(4,3) = m21
ele%mat6(4,4) = m22
ele%mat6(5,5) = 1
ele%mat6(6,6) = 1


end_orb = start_orb

end_orb%vec(1) = m11*start_orb%vec(1) + m12*start_orb%vec(2)
end_orb%vec(2) = m21*start_orb%vec(1) + m22*start_orb%vec(2)
end_orb%vec(3) = m11*start_orb%vec(3) + m12*start_orb%vec(4)
end_orb%vec(4) = m21*start_orb%vec(3) + m22*start_orb%vec(4)


end_orb%s = start_orb%s + L0
end_orb%t = start_orb%t + L0/vbar



err_flag = .false.

end subroutine
