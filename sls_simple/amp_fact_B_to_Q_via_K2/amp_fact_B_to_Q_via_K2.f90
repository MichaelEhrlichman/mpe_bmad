program amp_fact_B_to_Q_via_K2

use bmad

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s
type(ele_struct) inner_ele_at_s

integer i,j,mag_ix
integer status

real(rp) g, k2, ele_l
real(rp) s, s_inner, delta_s
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
nu_y = lat%ele(lat%n_ele_track)%a%phi / 2.0d0 / pi

open(45,file='B_sensitivity_coefficients.out')

write(45,'(a14,2a17)') "# loc (s)", "dQx/(dI/I)", "dQy/(dI/I)"

s=0.0
delta_s = 0.01 
do while (s .lt. lat%param%total_length)
  call twiss_and_track_at_s(lat, s, ele_at_s)
  if(has_attribute(ele_at_s, 'G')) then
    ! Sensitivity of tune to relative changes (dI/I) in magnet current.
    g = value_of_attribute(ele_at_s, 'G', err_print_flag=.true.)
    ele_l = delta_s
    s_inner = 0.0
    DQxij = 0.0d0
    do while (s_inner .lt. lat%param%total_length)
      call twiss_and_track_at_s(lat, s_inner, inner_ele_at_s)
      if(has_attribute(inner_ele_at_s, 'K2')) then
        k2 = value_of_attribute(inner_ele_at_s, 'K2', err_print_flag=.true.)
        Dphi_x = inner_ele_at_s%a%phi - ele_at_s%a%phi
        DQxij = DQxij + inner_ele_at_s%a%beta * k2 * ele_l / 4.0d0 / pi * &
                sqrt(inner_ele_at_s%a%beta*ele_at_s%a%beta) / 2.0d0 / sin(pi*nu_x) * cos(pi*nu_x - Dphi_x)
      endif
      s_inner = s_inner + delta_s
    enddo
    write(45,'(f14.5,2es17.5,i8)') s, DQxij, DQyij, element_at_s(lat,s,.true.)
  else
    !write(45,'(f14.5,2es14.5)') s, 0.0d0, 0.0d0
    write(45,'(f14.5,2a17,i8)') s, 'NaN', 'NaN', -1
  endif
  s = s + delta_s
enddo
close(45)

end program












