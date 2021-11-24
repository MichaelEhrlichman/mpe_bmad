program amp_fact_B_at_K2

use bmad

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) outer_ele
type(ele_struct) inner_ele
type(ele_pointer_struct), allocatable :: sextupoles(:)

integer i,j,mag_ix
integer ii, jj
integer status
integer nper, nloc

real(rp) angle, k2l
real(rp) s_outer, s_inner
real(rp) outer_delta_s, inner_delta_s
real(rp) Qx_eta, Qy_eta
real(rp) nu_x, nu_y
real(rp) Dphi_x, D_phi_y
real(rp) DQxij, DQyij
real(rp) sext_midpoint
real(rp), allocatable :: one_bend_contrib(:)
real(rp), allocatable :: ripple_at_sextupoles(:)

character(100) line

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)

nu_x = lat%ele(lat%n_ele_track)%a%phi / 2.0d0 / pi
nu_y = lat%ele(lat%n_ele_track)%b%phi / 2.0d0 / pi

call lat_ele_locator("sextupole::*", lat, sextupoles, nloc)
write(*,*) "num sextupoles: ", nloc
allocate(ripple_at_sextupoles(nloc))
allocate(one_bend_contrib(nloc))
ripple_at_sextupoles = 0.0d0

nper = 10
do i=1,lat%n_ele_track
  if(value_of_attribute(lat%ele(i), 'L') .lt. 1.0e-6) cycle 
  if(.not. has_attribute(lat%ele(i), 'G')) cycle

  one_bend_contrib(:) = 0.0d0
  outer_delta_s = value_of_attribute(lat%ele(i), 'L') / nper
  angle = value_of_attribute(lat%ele(i), 'G', err_print_flag=.true.) * outer_delta_s

  do ii=1,nper
    s_outer = lat%ele(i-1)%s + ii*outer_delta_s
    call twiss_and_track_at_s(lat, s_outer, outer_ele)
    do j=1,nloc
      if(value_of_attribute(sextupoles(j)%ele, 'L') .lt. 1.0e-6) cycle 
      sext_midpoint = sextupoles(j)%ele%s-value_of_attribute(sextupoles(j)%ele,'L')/2.0d0
      call twiss_and_track_at_s(lat, sext_midpoint, inner_ele)
      Dphi_x = inner_ele%a%phi - outer_ele%a%phi
      if(s_inner .lt. s_outer) Dphi_x = Dphi_x + lat%ele(lat%n_ele_track)%a%phi
      one_bend_contrib(j) = one_bend_contrib(j) + sqrt(inner_ele%a%beta*outer_ele%a%beta) /2.0d0 / sin(pi*nu_x) * cos(pi*nu_x + Dphi_x) * angle
    enddo
  enddo
  do j=1,nloc
    ripple_at_sextupoles(j) = ripple_at_sextupoles(j) + one_bend_contrib(j)**2
  enddo
enddo

open(45,file='Jitter_at_K2_sensitivity_coefficients.out')
write(45,'(a14,2a17)') "# loc (s)", "dx/(dI/I)", "dy/(dI/I)"
do j=1,nloc
  sext_midpoint = sextupoles(j)%ele%s-value_of_attribute(sextupoles(j)%ele,'L')/2.0d0
  write(45,'(f14.5,2es17.5,i8)') sext_midpoint, sqrt(ripple_at_sextupoles(j)), 0.0d0, j
enddo
close(45)

end program












