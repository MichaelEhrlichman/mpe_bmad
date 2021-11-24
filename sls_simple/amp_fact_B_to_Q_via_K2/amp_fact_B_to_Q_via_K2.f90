program amp_fact_B_to_Q_via_K2

use bmad

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) outer_ele
type(ele_struct) inner_ele

integer i,j,mag_ix
integer ii, jj
integer status
integer nper

real(rp) angle, k2l
real(rp) s_outer, s_inner
real(rp) outer_delta_s, inner_delta_s
real(rp) Qx_eta, Qy_eta
real(rp) nu_x, nu_y
real(rp) Dphi_x, D_phi_y
real(rp) DQxij, DQyij

character(100) line

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)

nu_x = lat%ele(lat%n_ele_track)%a%phi / 2.0d0 / pi
nu_y = lat%ele(lat%n_ele_track)%b%phi / 2.0d0 / pi

open(45,file='B_sensitivity_coefficients.out')

write(45,'(a14,2a17)') "# loc (s)", "dQx/(dI/I)", "dQy/(dI/I)"

DQyij = 0  !assume no vertical bend fields
nper = 10
do i=1,lat%n_ele_track
!do i=1,111
  if(value_of_attribute(lat%ele(i), 'L') .lt. 1.0e-6) cycle 
  if(has_attribute(lat%ele(i), 'G')) then
    outer_delta_s = value_of_attribute(lat%ele(i), 'L') / nper
    angle = value_of_attribute(lat%ele(i), 'G', err_print_flag=.true.) * outer_delta_s
    do ii=1,nper
      s_outer = lat%ele(i-1)%s + ii*outer_delta_s
      call twiss_and_track_at_s(lat, s_outer, outer_ele)
      DQxij = 0.0d0
      do j=1,lat%n_ele_track
        if(value_of_attribute(lat%ele(j), 'L') .lt. 1.0e-6) cycle 
        if(has_attribute(lat%ele(j), 'K2')) then
          inner_delta_s = value_of_attribute(lat%ele(j), 'L') / nper
          k2l = value_of_attribute(lat%ele(j), 'K2', err_print_flag=.true.) * inner_delta_s
          do jj = 1,nper
            s_inner = lat%ele(j-1)%s + jj*inner_delta_s
            call twiss_and_track_at_s(lat, s_inner, inner_ele)
            Dphi_x = inner_ele%a%phi - outer_ele%a%phi
            if(s_inner .lt. s_outer) Dphi_x = Dphi_x + lat%ele(lat%n_ele_track)%a%phi
            DQxij = DQxij + inner_ele%a%beta * k2l / 4.0d0 / pi * &
                    sqrt(inner_ele%a%beta*outer_ele%a%beta) / 2.0d0 / sin(pi*nu_x) * cos(pi*nu_x + Dphi_x) * angle
          enddo
        endif
      enddo
      write(45,'(2f14.5,2es17.5,i8)') s_outer, outer_delta_s, DQxij, DQyij, i
    enddo
  else
    write(45,'(2f14.5,2a17,i8)') lat%ele(i)%s, -1.0d0, 'NaN', 'NaN', -1
  endif
enddo
close(45)

!open(77,file='B_sensitivity_coefficients.out_ix')
!write(77,'(a14,2a17)') "# loc (s)", "dQx/(dI/I)", "dQy/(dI/I)"
!do i = 1,lat%n_ele_track
!  if(has_attribute(lat%ele(i), 'G')) then
!    ! Sensitivity of tune to relative changes (dI/I) in magnet current.
!    ele_l = value_of_attribute(lat%ele(i), 'L', err_print_flag=.true.)
!    angle = value_of_attribute(lat%ele(i), 'G', err_print_flag=.true.) * ele_l
!    DQxij = 0.0d0
!    do j=1,lat%n_ele_track
!      if(has_attribute(lat%ele(j), 'K2')) then
!        k2 = value_of_attribute(lat%ele(j), 'K2', err_print_flag=.true.)
!        inner_ele_l = value_of_attribute(lat%ele(j), 'L', err_print_flag=.true.)
!        Dphi_x = lat%ele(j)%a%phi - lat%ele(i)%a%phi
!        if(j .lt. i) Dphi_x = Dphi_x + lat%ele(lat%n_ele_track)%a%phi
!        write(770,'(3es14.4)')  lat%ele(j)%a%beta , k2 , inner_ele_l
!        write(771,'(5es14.4)')  sqrt(lat%ele(j)%a%beta*lat%ele(i)%a%beta) , nu_x, Dphi_x, angle
!        DQxij = DQxij + lat%ele(j)%a%beta * k2 * inner_ele_l / 4.0d0 / pi * &
!                sqrt(lat%ele(j)%a%beta*lat%ele(i)%a%beta) / 2.0d0 / sin(pi*nu_x) * cos(pi*nu_x + Dphi_x) * angle
!      endif
!    enddo
!    write(77,'(f14.5,2es17.5,i8)') lat%ele(i)%s, DQxij, DQyij, i
!  else
!    write(77,'(f14.5,2a17,i8)') lat%ele(i)%s, 'NaN', 'NaN', -1
!  endif
!enddo
!close(77)

end program












