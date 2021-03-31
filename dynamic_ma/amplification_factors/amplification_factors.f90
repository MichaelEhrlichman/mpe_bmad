program amplification_factors

use bmad

implicit none

character(200) lat_file
character(5) ele_ix_str

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
!type(ele_struct), target :: ele_at_eval
!type(ele_struct), pointer :: ele_at_eval
!type(ele_struct), pointer :: ele_at_source
type(ele_struct) :: ele_at_eval
type(ele_struct) :: ele_at_source

integer i,j
integer stat
integer ele_ix

real(rp) eval_s
real(rp) s, delta_s
real(rp) k1l_at_source, angle_at_source
real(rp) nu_x, dphi_x
real(rp) nu_y, dphi_y
real(rp) amp_factor_x
real(rp) amp_factor_y

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call getarg(2, ele_ix_str)

call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,stat)
nu_x = lat%ele(lat%n_ele_track)%a%phi / twopi
nu_y = lat%ele(lat%n_ele_track)%b%phi / twopi

!eval_s = 13.06  !0.4 m past end of SHF
!call twiss_and_track_at_s(lat, eval_s, ele_at_eval)
!ele_ix = 112
read(ele_ix_str,*) ele_ix
!ele_at_eval => lat%ele(ele_ix)
call twiss_and_track_at_s(lat, lat%ele(ele_ix)%s-lat%ele(ele_ix)%value(l$)/2.0d0, ele_at_eval)

open(40,FILE='amp_factors_ele'//trim(ele_ix_str)//'.dat')
write(40,'(a,f12.4)') "# amplification factors from source to s= ", eval_s
write(40,'(a)') "# source index, amp_factor_x, amp_factor_y, source name"
do i=1,lat%n_ele_track
  !ele_at_source => lat%ele(i)
  call twiss_and_track_at_s(lat,lat%ele(i)%s-lat%ele(i)%value(l$)/2.0d0,ele_at_source)
  if(has_attribute(ele_at_source, 'K1')) then
    dphi_x = ele_at_eval%a%phi - ele_at_source%a%phi
    dphi_y = ele_at_eval%b%phi - ele_at_source%b%phi
    k1l_at_source = ele_at_source%value(l$) * value_of_attribute(ele_at_source, 'K1', err_print_flag = .true.)
    amp_factor_x = sqrt(ele_at_source%a%beta*ele_at_eval%a%beta)*cos(pi*nu_x-abs(dphi_x))/2.0d0/sin(pi*nu_x) * k1l_at_source
    amp_factor_y = sqrt(ele_at_source%b%beta*ele_at_eval%b%beta)*cos(pi*nu_y-abs(dphi_y))/2.0d0/sin(pi*nu_y) * k1l_at_source
    write(40,'(i6,2es14.5,"   ",a)') i, amp_factor_x, amp_factor_y, ele_at_source%name
  endif
enddo
close(40)

open(40,FILE='amp_factors_roll.dat')
write(40,'(a,f12.4)') "# roll amplification factors from source to s= ", eval_s
write(40,'(a)') "# source index, amp_factor_x, amp_factor_y, source name"
do i=1,lat%n_ele_track
  !ele_at_source => lat%ele(i)
  call twiss_and_track_at_s(lat,lat%ele(i)%s-lat%ele(i)%value(l$)/2.0d0,ele_at_source)
  if(has_attribute(ele_at_source, 'ANGLE')) then
    dphi_x = ele_at_eval%a%phi - ele_at_source%a%phi
    dphi_y = ele_at_eval%b%phi - ele_at_source%b%phi
    angle_at_source = value_of_attribute(ele_at_source, 'ANGLE', err_print_flag = .true.)
    amp_factor_x = sqrt(ele_at_source%a%beta*ele_at_eval%a%beta)*cos(pi*nu_x-abs(dphi_x))/2.0d0/sin(pi*nu_x) * angle_at_source
    amp_factor_y = sqrt(ele_at_source%b%beta*ele_at_eval%b%beta)*cos(pi*nu_y-abs(dphi_y))/2.0d0/sin(pi*nu_y) * angle_at_source
    write(40,'(i6,2es14.5,"   ",a)') i, amp_factor_x, amp_factor_y, ele_at_source%name
  endif
enddo
close(40)

end program












