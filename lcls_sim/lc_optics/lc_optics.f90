program lc_optics

use bmad
use write_lat_file_mod
use bmad_parser_mod

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s

integer i,j
integer ix
integer status

real(rp) Ha, Hb, e_tot
real(rp) s, delta_s
real(rp) emitx, emity, sigpop
real(rp) dummy

bp_com%always_parse = .true.
bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

allocate(orb(0:lat%n_ele_track))
orb(0)%vec(6) = 0.0d0
call twiss_and_track(lat,orb,status)

emitx = 50.0e-12
emity = 50.0e-12
sigpop = 1.0e-4

! Twiss by ele end

open(45,file='twiss.out')

write(45,*) "!   a-mode phase advance:    ", lat%ele(lat%n_ele_track)%a%phi /2./pi
write(45,*) "!   b-mode phase advance:    ", lat%ele(lat%n_ele_track)%b%phi /2./pi
write(45,*) "!   total length:            ", lat%param%total_length

write(45,'(a11,a15,a11,19a16,"   ",a14)') '!        ix', 's (m)', 'l (m)', 'beta_a', 'alpha_a', 'gamma_a', 'eta_a', 'etap_a', 'H_a', 'eta_x', 'etap_x', 'phi_x', &
                      'beta_b', 'alpha_b', 'gamma_b', 'eta_b', 'etap_b', 'H_b', 'eta_y', 'etap_y', 'phi_y', 'energy', 'name'
do i=1,lat%n_ele_track
  Ha = ( (lat%ele(i)%a%eta**2) + ( lat%ele(i)%a%beta*lat%ele(i)%a%etap + lat%ele(i)%a%alpha*lat%ele(i)%a%eta )**2)/lat%ele(i)%a%beta
  Hb = ( (lat%ele(i)%b%eta**2) + ( lat%ele(i)%b%beta*lat%ele(i)%b%etap + lat%ele(i)%b%alpha*lat%ele(i)%b%eta )**2)/lat%ele(i)%b%beta
  write(45,'(i11,f15.8,f11.4,19es16.7,"   ",a14)') i, &
        lat%ele(i)%s, lat%ele(i)%value(l$), &
        lat%ele(i)%a%beta, lat%ele(i)%a%alpha, lat%ele(i)%a%gamma, lat%ele(i)%a%eta, lat%ele(i)%a%etap, Ha, &
        lat%ele(i)%x%eta, lat%ele(i)%x%etap, lat%ele(i)%a%phi, &
        lat%ele(i)%b%beta, lat%ele(i)%b%alpha, lat%ele(i)%b%gamma, lat%ele(i)%b%eta, lat%ele(i)%b%etap, Hb, &
        lat%ele(i)%y%eta, lat%ele(i)%y%etap, lat%ele(i)%b%phi, value_of_attribute(lat%ele(i),'e_tot'), &
        lat%ele(i)%name
enddo
close(45)

! Twiss by s

open(55,file='twiss_by_s.out')

write(55,'(a11,2a11,19a14,"   ",a20,a14)') '!        ix', 's (m)', 'l (m)', 'beta_a', 'alpha_a', 'gamma_a', 'eta_a', 'etap_a', 'H_a', 'eta_x', 'etap_x', 'phi_x', &
                      'beta_b', 'alpha_b', 'gamma_b', 'eta_b', 'etap_b', 'H_b', 'eta_y', 'etap_y', 'phi_y', 'energy', 'name', 'beam_size'
s=0.0
delta_s = 0.01
do while (s .lt. lat%param%total_length)
  call twiss_and_track_at_s(lat, s, ele_at_s)
  Ha = ( (ele_at_s%a%eta**2) + ( ele_at_s%a%beta*ele_at_s%a%etap + ele_at_s%a%alpha*ele_at_s%a%eta )**2)/ele_at_s%a%beta
  Hb = ( (ele_at_s%b%eta**2) + ( ele_at_s%b%beta*ele_at_s%b%etap + ele_at_s%b%alpha*ele_at_s%b%eta )**2)/ele_at_s%b%beta
  write(55,'(i11,2f11.3,19es14.4,"   ",a20,2es14.4)') 0, &
        ele_at_s%s, delta_s, &
        ele_at_s%a%beta, ele_at_s%a%alpha, ele_at_s%a%gamma, ele_at_s%a%eta, ele_at_s%a%etap, Ha, &
        ele_at_s%x%eta, ele_at_s%x%etap, ele_at_s%a%phi, &
        ele_at_s%b%beta, ele_at_s%b%alpha, ele_at_s%b%gamma, ele_at_s%b%eta, ele_at_s%b%etap, Hb, &
        ele_at_s%y%eta, ele_at_s%y%etap, ele_at_s%b%phi, value_of_attribute(ele_at_s,'e_tot'), &
        ele_at_s%name, sqrt(emitx*ele_at_s%a%beta + sigpop**2 * ele_at_s%a%eta**2), sqrt(emity*ele_at_s%b%beta)
  s = s + delta_s
enddo
close(55)

!  ! Twiss midpoints
!  
!  open(55,file='twiss_at_midpoints.out')
!  
!  write(55,'(a11,a13,19a14,"   ",a20)') '!        ix', 's (m)', 'beta_a', 'alpha_a', 'gamma_a', 'eta_a', 'etap_a', 'H_a', 'eta_x', 'etap_x', 'phi_x', &
!                        'beta_b', 'alpha_b', 'gamma_b', 'eta_b', 'etap_b', 'H_b', 'eta_y', 'etap_y', 'phi_y', 'energy', 'name'
!  do i=1,lat%n_ele_track
!    if(lat%ele(i)%value(l$) .ge. 0) then
!      if(lat%ele(i)%value(l$) .gt. 0) then
!        s = lat%ele(i)%s - lat%ele(i)%value(l$)/2.0
!        if( s .lt. 0) cycle
!        call twiss_and_track_at_s(lat, s, ele_at_s)
!      else
!        ele_at_s = lat%ele(i)
!        s = ele_at_s%s
!        dummy = value_of_attribute(ele_at_s,'e_tot')
!      endif
!      Ha = ( (ele_at_s%a%eta**2) + ( ele_at_s%a%beta*ele_at_s%a%etap + ele_at_s%a%alpha*ele_at_s%a%eta )**2)/ele_at_s%a%beta
!      Hb = ( (ele_at_s%b%eta**2) + ( ele_at_s%b%beta*ele_at_s%b%etap + ele_at_s%b%alpha*ele_at_s%b%eta )**2)/ele_at_s%b%beta
!      ix = element_at_s(lat,s,.false.)
!      e_tot = value_of_attribute(ele_at_s,'e_tot')
!      write(55,'(i11,f13.5,19es14.4,"   ",a20)') ix, ele_at_s%s, &
!            ele_at_s%a%beta, ele_at_s%a%alpha, ele_at_s%a%gamma, ele_at_s%a%eta, ele_at_s%a%etap, Ha, &
!            ele_at_s%x%eta, ele_at_s%x%etap, ele_at_s%a%phi, &
!            ele_at_s%b%beta, ele_at_s%b%alpha, ele_at_s%b%gamma, ele_at_s%b%eta, ele_at_s%b%etap, Hb, &
!            ele_at_s%y%eta, ele_at_s%y%etap, ele_at_s%b%phi, e_tot, ele_at_s%name
!    endif
!  enddo
!  close(55)



end program












