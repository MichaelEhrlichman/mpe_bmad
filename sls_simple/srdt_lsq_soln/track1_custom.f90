!+
! Subroutine track1_custom (start_orb, ele, param, end_orb, err_flag, finished, track)
!
! Dummy routine for custom tracking. 
! This routine needs to be replaced for a custom calculation.
! If not replaced and this routine is called, this routine will generate an error message.
!
! Also see:
!   track1_preprocess
!   track1_postprocess
!
! If this routine takes into account radiation damping and/or excitation when bmad_com%radiation_damping_on 
! and/or bmad_com%radiation_fluctuations_on is True, a custom version of track1_preprocess should be 
! constructed to set its radiation_included argument to True.
! If not, the track1 routine will use track1_radiation to include the radiation effects.
! Note: If this routine calles symp_lie_bmad, the symp_lie_bmad routine does take into account radiation effects.
!
! Note: If ele%spin_tracking_method = tracking$, then it is expected that this routine will also handle
! spin tracking. The alternative is when ele%spin_tracking_method = custom$ in which case track1_spin_custom will
! be called after this routine. If doing spin tracking here, bmad_com%spin_tracking_on should be checked
! to see if spin tracking is actually wanted.
!
! General rule: Your code may NOT modify any argument that is not listed as an output agument below.
!
! Input:
!   start_orb  -- coord_struct: Starting position.
!   ele        -- ele_struct: Element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   end_orb     -- coord_struct: End position.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt processing and return to its calling routine.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!-

subroutine track1_custom (start_orb, ele, param, end_orb, err_flag, finished, track)

use bmad, except_dummy => track1_custom

implicit none

type (coord_struct) :: start_orb
type (coord_struct) :: end_orb
type (ele_struct) :: ele
type (lat_param_struct) :: param
type (track_struct), optional :: track

real(rp) m11, m12, m21, m22
real(rp) L0, K_und, B, vbar, gamma_rel
real(rp) lambda, k1

logical err_flag, finished

character(*), parameter :: r_name = 'track1_custom'

L0 = ele%value(n_pole$)*ele%value(l_pole$) !undulator length
lambda = ele%value(l_pole$)*2 !undulator wavelength
B = ele%value(b_max$) !peak field (T)
gamma_rel = sqrt((start_orb%p0c*(1.+start_orb%vec(6)))**2+m_electron**2)/m_electron !relativistic gamma

K_und = e_charge*B*lambda/(twopi*m_electron*e_charge/c_light) !undualtor paramater
vbar = c_light*(1.-1./2./gamma_rel/gamma_rel*(1.+K_und*K_und)) !average longitudinal velocity (m/s)

k1 = 1./2.*((e_charge*B)/(vbar*gamma_rel*m_electron*e_charge/c_light**2))**2 !effective focusing strength

!transfer matrix
m11 = cos(sqrt(k1)*L0)
m12 = 1./sqrt(k1)*sin(sqrt(k1)*L0)
m21 = -sqrt(k1)*sin(sqrt(k1)*L0)
m22 = cos(sqrt(k1)*L0)

end_orb = start_orb

if(L0<1.29 .or. L0>1.31 .or. L0>-1) then
end_orb%vec(1) = m11*start_orb%vec(1) + m12*start_orb%vec(2)
end_orb%vec(2) = m21*start_orb%vec(1) + m22*start_orb%vec(2)
end_orb%vec(3) = m11*start_orb%vec(3) + m12*start_orb%vec(4)
end_orb%vec(4) = m21*start_orb%vec(3) + m22*start_orb%vec(4)
else
end_orb%vec(1) = 0
end_orb%vec(2) = 0
end_orb%vec(3) = 0
end_orb%vec(4) = 0
end_orb%vec(1) = end_orb%vec(1) + &
0.999264*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
1.29968*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
4.39598e-05*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
7.61458e-05*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.00113243*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.999264*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
6.63486e-08*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.000161052*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
4.39026e-05*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
7.61469e-05*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.999264*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
1.29968*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
6.63461e-08*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
7.31893e-05*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00113241*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.999264*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
0.131118*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
2.5975e-06*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
2.79955e-05*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
0.013569*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
2.5975e-06*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) - &
0.00210267*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) - &
0.0135587*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) - &
0.0088139*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
2.79955e-05*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) - &
0.0135587*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
0.393338*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
2.59387e-06*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
0.013569*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) - &
0.0088139*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) + &
2.59387e-06*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) - &
0.000692918*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) + &
7.42156e-05*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
0.131023*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) + &
4.44797e-06*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
2.22201e-05*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
0.131023*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
0.170414*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
2.60544e-05*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
0.0135953*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) + &
4.44797e-06*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
2.60544e-05*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) + &
0.000222627*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
0.393051*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
2.22201e-05*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
0.0135953*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
0.393051*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(2) = end_orb%vec(2) - &
0.511208*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) - &
5.32615e-05*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
0.013549*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) - &
0.13111*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) - &
1.0071e-06*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
0.013549*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
0.0264427*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
2.20253e-06*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) - &
0.000704878*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) - &
0.13111*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
2.20253e-06*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
2.72924e-06*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
0.0135595*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) - &
1.0071e-06*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) - &
0.000704878*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
0.0135595*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(3) = end_orb%vec(3) + &
0.00881498*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
2.21249e-05*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
5.13052e-05*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) - &
7.42029e-05*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
0.131015*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
5.13052e-05*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
0.0136377*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
0.131012*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
0.170396*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) - &
7.42029e-05*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
0.131012*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
4.87227e-05*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
2.60716e-05*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
0.131015*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
0.170396*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) + &
2.60716e-05*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))
end_orb%vec(4) = end_orb%vec(4) - &
0.0135531*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))
end_orb%vec(1) = end_orb%vec(1) - &
0.120333*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0525151*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00364209*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00770228*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0525151*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0341361*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
2.8512e-06*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00295354*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
0.000452653*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
0.00355523*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0458411*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0197219*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00414973*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.000677368*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0201448*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0132175*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0525151*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0340444*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00177617*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00414273*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0341361*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.026793*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00113988*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00266167*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
0.00355523*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
0.00232541*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0203339*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0127581*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.000677368*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
0.000337507*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0130337*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.010112*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00364209*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00177617*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0286512*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0125558*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) + &
2.8512e-06*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00113988*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.012037*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00782406*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0458411*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0203339*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00864469*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00947659*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0201448*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0130337*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00414127*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00475269*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00770228*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00414273*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0125558*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0086514*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00295354*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00266167*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00782406*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00625337*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0197219*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0127581*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00947659*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0059196*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.0132175*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.010112*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00475269*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(1) = end_orb%vec(1) - &
0.00500516*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.211478*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.13747*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
2.89658e-05*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0122762*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.120278*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.103883*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.00181192*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0100606*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
6.45047e-07*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.00364561*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0704947*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.00574904*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0123053*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.00651783*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0458191*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0399059*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.13747*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.103883*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0059138*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0153964*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.103883*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.100854*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.00769639*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0158796*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.00364561*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.00473882*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0973956*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0394368*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.00651783*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.00205012*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0395757*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0382134*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
2.89658e-05*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0059138*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.0352442*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.0229341*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.00181192*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.00769639*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0286391*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0243983*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0704947*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0973956*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
3.56224e-05*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0182115*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0458191*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0395757*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0141028*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0183499*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0122762*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0153964*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.0229341*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0255753*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0100606*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0158796*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0243983*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0236841*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) + &
0.00574904*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0394368*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0182115*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0236846*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0399059*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0382134*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0183499*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(2) = end_orb%vec(2) - &
0.0215272*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0104659*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0094756*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0458415*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0196751*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00413977*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00471081*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0201445*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0129724*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0286492*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0125535*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0072812*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00177301*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0116597*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0075785*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
5.5607e-06*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00112362*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0094756*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00590002*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0197199*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0123287*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00471081*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0049102*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0131558*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0100322*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0125535*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0087736*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0148144*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00410077*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0075785*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.00641275*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00293372*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0026398*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0458415*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0197199*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.00228336*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.00355847*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0201445*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0131558*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00414818*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.000635372*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0072812*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0148144*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.120333*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.052045*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
5.5607e-06*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00293372*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0517603*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0336458*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0196751*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0123287*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.00355847*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.00234194*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0129724*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0100322*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.000635372*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.000350068*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00177301*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00410077*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.052045*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0335542*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.00112362*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) + &
0.0026398*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0336458*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(3) = end_orb%vec(3) - &
0.0261615*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
3.27774e-05*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00545015*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0704908*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0286286*start_orb%vec(1)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00134429*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00531445*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0458175*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0398577*start_orb%vec(1)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0704891*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0458247*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
6.00876e-06*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00226996*start_orb%vec(1)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0286367*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0246809*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00183491*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.000599898*start_orb%vec(1)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00545015*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.0070934*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0630071*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0394326*start_orb%vec(2)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00531445*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00921529*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0400922*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0392549*start_orb%vec(2)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0458247*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0251992*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00867191*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00354713*start_orb%vec(2)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0246809*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0245469*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00176929*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.000513117*start_orb%vec(2)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0704908*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0630071*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
8.39187e-06*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00182462*start_orb%vec(3)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0458175*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0400922*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00318772*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00177997*start_orb%vec(3)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
6.00876e-06*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00867191*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.105737*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.068702*start_orb%vec(3)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00183491*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00176929*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.12027*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.104628*start_orb%vec(3)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0286286*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0394326*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00182462*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00237042*start_orb%vec(4)*start_orb%vec(1)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0398577*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.0392549*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00177997*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.00409324*start_orb%vec(4)*start_orb%vec(2)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00226996*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.00354713*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.068702*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.104629*start_orb%vec(4)*start_orb%vec(3)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.000599898*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(1)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) + &
0.000513117*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(2)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.104628*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(3)/(1.+start_orb%vec(6))**2
end_orb%vec(4) = end_orb%vec(4) - &
0.102314*start_orb%vec(4)*start_orb%vec(4)*start_orb%vec(4)/(1.+start_orb%vec(6))**2
endif


end_orb%s = start_orb%s + L0
end_orb%t = start_orb%t + L0/vbar

err_flag = .false.
finished = .false.


end subroutine
