program amp_fact_Q

use bmad

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s
type(ele_struct), pointer :: ele_at_kick

integer i,ii
integer status
integer nper

real(rp) k1l
real(rp) s, delta_s
real(rp) Qx_eta, Qy_eta

character(100) line

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)

open(45,file='sensitivity_coefficients.out')

write(45,'(a14,4a17)') "# loc (s)", "dQx/(dI/I)", "dQy/(dI/I)", "dQx/(dI/I)*Inom", "dQy/(dI/I)*Inom"

nper = 10
do i=1,lat%n_ele_track
  if(value_of_attribute(lat%ele(i), 'L') .lt. 1.0e-6) cycle 
  if(has_attribute(lat%ele(i), 'K1')) then
    delta_s = value_of_attribute(lat%ele(i), 'L') / nper
    k1l = value_of_attribute(lat%ele(i), 'K1', err_print_flag=.true.) * delta_s
    do ii=1,nper
      s = lat%ele(i-1)%s + ii*delta_s
      call twiss_and_track_at_s(lat, s, ele_at_s)
      Qx_eta = ele_at_s%a%beta * k1l / 4.0d0 / pi
      Qy_eta = ele_at_s%b%beta * k1l / 4.0d0 / pi
      write(45,'(f14.5,2es17.5,i8)') s, Qx_eta, Qy_eta, element_at_s(lat,s,.true.)
    enddo
  else
    !write(45,'(f14.5,2es14.5)') s, 0.0d0, 0.0d0
    write(45,'(f14.5,2a17,i8)') s, 'NaN', 'NaN', -1
  endif
enddo
close(45)

end program












