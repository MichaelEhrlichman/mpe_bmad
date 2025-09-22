!+
! Subroutine track1_custom (orbit, ele, param, err_flag, finished, track)
!
! Prototype routine for custom tracking. 
!
! Also see:
!   track1_preprocess
!   track1_postprocess
!
! If this routine handles radiation damping and/or excitation when bmad_com%radiation_damping_on 
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
!   orbit      -- coord_struct: Starting position.
!   ele        -- ele_struct: Lattice element.
!   param      -- lat_param_struct: Lattice parameters.
!
! Output:
!   orbit       -- coord_struct: End position.
!   err_flag    -- logical: Set true if there is an error. False otherwise.
!   finished    -- logical: When set True, track1 will halt processing and return to its calling routine.
!   track       -- track_struct, optional: Structure holding the track information if the 
!                    tracking method does tracking step-by-step.
!-

subroutine track1_custom (orbit, ele, param, err_flag, finished, track)

use bmad

implicit none

type (coord_struct) :: orbit
type (ele_struct) :: ele
type (ele_struct), pointer :: lord
type (lat_param_struct) :: param
type (track_struct), optional :: track
type(all_pointer_struct) v_ptr, s_ptr, hv_ptr, ht_ptr

logical err_flag, finished

character(*), parameter :: r_name = 'track1_custom'

real(rp) voltage_nominal, slope, voltage
real(rp) hornV, hornT
real(rp) dE, length, p0c, time
real(rp) beta0, T0

integer i, n_slice, Nturns

err_flag = .false.

call pointer_to_attribute(ele,'VOLTAGE', .true., v_ptr, err_flag)
call pointer_to_attribute(ele,'SLOPE', .true., s_ptr, err_flag)
call pointer_to_attribute(ele,'HORNV', .true., hv_ptr, err_flag)
call pointer_to_attribute(ele,'HORNT', .true., ht_ptr, err_flag)
voltage_nominal = v_ptr%r
slope = s_ptr%r
hornV = hv_ptr%r
hornT = ht_ptr%r

!zero-length induction cell for now
n_slice = 1
length = 0

p0c = orbit%p0c

!call convert_pc_to(p0c,orbit%species,beta=beta0)
!T0 = param%total_length/c_light/beta0

call offset_particle (ele, set$, orbit)

do i = 0, n_slice
  ! if (logic_option(.false., make_matrix)) then
  !   factor = voltage_nominal / n_slice
  !   if (i == 0 .or. i == n_slice) factor = factor / 2

  !   dE = factor
  !   pc = (1 + orbit%vec(6)) * p0c 
  !   E = p0c / orbit%beta
  !   call convert_total_energy_to (E + dE, orbit%species, beta = new_beta)

  !   m2(2,1) = 0.0
  !   m2(2,2) = orbit%beta / new_beta 
  !   m2(1,1) = new_beta / orbit%beta
  !   m2(1,2) = orbit%vec(5) * mc2**2 * p0c * (m2(2,2) / ((E+dE)**3 * orbit%beta) - new_beta / (pc**2 * E))
  !   if (orbit%time_dir == -1) call mat_inverse(m2, m2)

  !   mat6(5:6, :) = matmul(m2, mat6(5:6, :))
  ! endif

  !Nturns = nint(orbit%t / T0)
  !time = orbit%t - ele%ref_time - T0*Nturns
  time = particle_rf_time(orbit, ele, .false.)

  if (abs(time) .lt. hornT) then
    !write(*,*) "no horn    ", time, orbit%vec(5)
    voltage = voltage_nominal 
  else if (time .ge. hornT) then
    !write(*,*) "ge horn ***", time, orbit%vec(5)
    voltage = voltage_nominal + hornV
  else if (time .le. hornT) then
    !write(*,*) "le horn ***", time, orbit%vec(5)
    voltage = voltage_nominal - hornV
  endif
  dE = orbit%time_dir * voltage / n_slice
  !write(*,*) "time, voltage: ", time, voltage
  if (i == 0 .or. i == n_slice) dE = dE / 2

  call apply_energy_kick (dE, orbit, [0.0_rp, 0.0_rp])
  
  if (orbit%vec(6) == -1) then
    orbit%state = lost_pz$
    return
  endif

  if (i /= n_slice) then
    call track_a_drift (orbit, length/n_slice)
  endif
enddo

! If ele represents a section of the entire element and the element has internal structure (not uniform
! longitudinally), then it will be important to know the position of ele relative to the entire element.
! The entire element will be the lord:

if (ele%slave_status == slice_slave$ .or. ele%slave_status == super_slave$) then
  lord => pointer_to_super_lord(ele)
  ! Now [lord%s_start, lord%s] is the position of the entire element.
endif

call offset_particle (ele, unset$, orbit)

contains
  function induction_cell_voltage(time,vmax) result(v)
    real(rp) v
    real(rp) time, vmax

    if (abs(time) .gt. .667e-9) then 
    endif    

  end function

end subroutine
