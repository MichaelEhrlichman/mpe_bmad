program amp_fact_Q

use bmad

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s
type(ele_struct), pointer :: ele_at_kick

type mag_scaling_struct
  character(10) name
  real(rp) nom_I
  real(rp) lat_K1
  real(rp) nom_K1
end type
type(mag_scaling_struct) mag_scaling_data(18)

integer i,j,mag_ix
integer status

real(rp) k1, ele_l
real(rp) s, delta_s
real(rp) Qx_eta, Qy_eta
real(rp) Qx_eta_scaled, Qy_eta_scaled

character(100) line

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,status)

open(40,file='mag_scaling.dat')
read(40,*) line  !header
do i=1,18
  read(40,*)  mag_scaling_data(i)%name, mag_scaling_data(i)%nom_I, mag_scaling_data(i)%lat_K1, mag_scaling_data(i)%nom_K1
enddo
close(40)

open(45,file='sensitivity_coefficients.out')

write(45,'(a14,4a17)') "# loc (s)", "dQx/(dI/I)", "dQy/(dI/I)", "dQx/(dI/I)*Inom", "dQy/(dI/I)*Inom"

s=0.0
delta_s = 0.002  !0.002
do while (s .lt. lat%param%total_length)
  call twiss_and_track_at_s(lat, s, ele_at_s)
  if(has_attribute(ele_at_s, 'K1')) then
    ! Sensitivity of tune to relative changes (dI/I) in magnet current.
    k1 = value_of_attribute(ele_at_s, 'K1', err_print_flag=.true.)
    do mag_ix=1,18
      if( trim(adjustl(upcase(ele_at_s%name))) == trim(adjustl(mag_scaling_data(mag_ix)%name))) then
        exit
      endif
      if(mag_ix==18) then
        write(*,*) "FOO: ERROR NAME NOT MATCHED", ele_at_s%name
      endif
    enddo
    ele_l = delta_s
    Qx_eta = ele_at_s%a%beta * k1 * ele_l / 4.0d0 / pi
    Qy_eta = ele_at_s%b%beta * k1 * ele_l / 4.0d0 / pi
    Qx_eta_scaled = Qx_eta / abs(mag_scaling_data(mag_ix)%lat_K1) * mag_scaling_data(mag_ix)%nom_I
    Qy_eta_scaled = Qy_eta / abs(mag_scaling_data(mag_ix)%lat_K1) * mag_scaling_data(mag_ix)%nom_I
    write(45,'(f14.5,4es17.5,i8)') s, Qx_eta, Qy_eta, Qx_eta_scaled, Qy_eta_scaled, element_at_s(lat,s,.true.)
  else
    !write(45,'(f14.5,2es14.5)') s, 0.0d0, 0.0d0
    write(45,'(f14.5,4a17,i8)') s, 'NaN', 'NaN', 'NaN', 'NaN', -1
  endif
  s = s + delta_s
enddo
close(45)

end program












