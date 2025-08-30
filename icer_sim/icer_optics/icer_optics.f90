program icer_optics

use bmad
use bmad_parser_mod
use mode3_mod

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(coord_struct)  orb_at_s
type(ele_struct) ele_at_s
real(rp), allocatable :: sbend_k2_state(:)
real(rp), allocatable :: sextupole_state(:)
real(rp), allocatable :: multipole_state(:)
real(rp), allocatable :: partial_a(:)

integer i,j
integer status
integer nslices
logical err_flag

real(rp) Ha, Hb
real(rp) nat_chrom_x, nat_chrom_y
real(rp) cor_chrom_x, cor_chrom_y
real(rp) s, delta_s
real(rp) oor
real(rp) total_bend_angle
real(rp) emitx, emity, sigpop
real(rp) rho

bp_com%always_parse = .true.
bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

call set_on_off(rfcavity$,lat,off$)

!call write_bmad_lattice_file('test_lat.bmad',lat)

!call set_on_off(sextupole$, lat, off_and_save$, saved_values=sextupole_state)
!call set_on_off(multipole$, lat, off_and_save$, saved_values=multipole_state)
!call set_on_off(sbend$, lat, off_and_save$, saved_values=sbend_k2_state, attribute='K2')
!call twiss_and_track(lat,orb,status)
!call chrom_calc(lat,1.0D-6,nat_chrom_x,nat_chrom_y)
!
!call set_on_off(sextupole$, lat, restore_state$, saved_values=sextupole_state)
!call set_on_off(multipole$, lat, restore_state$, saved_values=multipole_state)
!call set_on_off(sbend$, lat, restore_state$, saved_values=sbend_k2_state, attribute='K2')
!call twiss_and_track(lat,orb,status)
!call chrom_calc(lat,1.0D-6,cor_chrom_x,cor_chrom_y)

allocate(orb(0:lat%n_ele_track))
orb(0)%vec(6) = 0.0d0
call twiss_and_track(lat,orb,status)

! call chrom_calc(lat,0.0D-6,cor_chrom_x,cor_chrom_y)
 
 call calc_z_tune(lat%branch(0))

total_bend_angle = 0.0d0
do i=1,lat%n_ele_track
  if(lat%ele(i)%key == sbend$ .or. lat%ele(i)%key == rbend$) then
    total_bend_angle = total_bend_angle + abs(lat%ele(i)%value(angle$))
  endif
enddo

open(45,file='closed_orbit.out')
do i=0,lat%n_ele_track
  write(45,'(i11,f11.5,6es14.4)') i, lat%ele(i)%s, orb(i)%vec(1:6)
enddo
close(45)

!
!
! Twiss by ele
!

emitx = 70.0e-12
emity = 70.0e-12
sigpop = 9.4e-4

open(45,file='twiss.out')

write(45,*) "!   a-mode tune:          ", lat%ele(lat%n_ele_track)%a%phi /2./pi
write(45,*) "!   b-mode tune:          ", lat%ele(lat%n_ele_track)%b%phi /2./pi
write(45,*) "!   Synchrotron tune:     ", lat%z%tune /2./pi
write(45,*) "!   nat chrom x:          ", nat_chrom_x
write(45,*) "!   nat chrom y:          ", nat_chrom_y
write(45,*) "!   cor chrom x:          ", cor_chrom_x
write(45,*) "!   cor chrom y:          ", cor_chrom_y
write(45,*) "!   total bend angle:     ", total_bend_angle / twopi * 360.0
write(45,*) "!   total length:         ", lat%param%total_length
write(45,*) "!   10^9 per mA:             ", 0.001 * lat%param%total_length / e_charge / c_light / 1.0e9

write(45,'(a11,2a11,19a16,"   ",a14)') '!        ix', 's (m)', 'l (m)', 'beta_a', 'alpha_a', 'gamma_a', 'eta_a', 'etap_a', 'H_a', 'eta_x', 'etap_x', 'phi_x', &
                      'beta_b', 'alpha_b', 'gamma_b', 'eta_b', 'etap_b', 'H_b', 'eta_y', 'etap_y', 'phi_y', 'energy', 'name'
do i=1,lat%n_ele_track
  Ha = ( (lat%ele(i)%a%eta**2) + ( lat%ele(i)%a%beta*lat%ele(i)%a%etap + lat%ele(i)%a%alpha*lat%ele(i)%a%eta )**2)/lat%ele(i)%a%beta
  Hb = ( (lat%ele(i)%b%eta**2) + ( lat%ele(i)%b%beta*lat%ele(i)%b%etap + lat%ele(i)%b%alpha*lat%ele(i)%b%eta )**2)/lat%ele(i)%b%beta
  write(45,'(i11,f11.5,f11.4,19es16.7,"   ",a14)') i, &
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
open(65,file='closed_orbit_by_a.out')

write(55,'(a11,2a11,19a14,"   ",a20,a14)') '!        ix', 's (m)', 'l (m)', 'beta_a', 'alpha_a', 'gamma_a', 'eta_a', 'etap_a', 'H_a', 'eta_x', 'etap_x', 'phi_x', &
                      'beta_b', 'alpha_b', 'gamma_b', 'eta_b', 'etap_b', 'H_b', 'eta_y', 'etap_y', 'phi_y', 'energy', 'name', 'beam_size'
s=0.0
delta_s = 0.01
i=1
do while (s .lt. lat%param%total_length)
  call twiss_and_track_at_s(lat, s, ele_at_s, orb, orb_at_s)
  write(65,'(i11,f11.5,6es14.4)') i, s, orb_at_s%vec(1:6)
  Ha = ( (ele_at_s%a%eta**2) + ( ele_at_s%a%beta*ele_at_s%a%etap + ele_at_s%a%alpha*ele_at_s%a%eta )**2)/ele_at_s%a%beta
  Hb = ( (ele_at_s%b%eta**2) + ( ele_at_s%b%beta*ele_at_s%b%etap + ele_at_s%b%alpha*ele_at_s%b%eta )**2)/ele_at_s%b%beta
  write(55,'(i11,2f11.3,19es14.4,"   ",a20,2es14.4)') element_at_s(lat,s,.true.), &
        ele_at_s%s, delta_s, &
        ele_at_s%a%beta, ele_at_s%a%alpha, ele_at_s%a%gamma, ele_at_s%a%eta, ele_at_s%a%etap, Ha, &
        ele_at_s%x%eta, ele_at_s%x%etap, ele_at_s%a%phi, &
        ele_at_s%b%beta, ele_at_s%b%alpha, ele_at_s%b%gamma, ele_at_s%b%eta, ele_at_s%b%etap, Hb, &
        ele_at_s%y%eta, ele_at_s%y%etap, ele_at_s%b%phi, value_of_attribute(ele_at_s,'e_tot'), &
        ele_at_s%name, sqrt(emitx*ele_at_s%a%beta + sigpop**2 * ele_at_s%a%eta**2), sqrt(emity*ele_at_s%b%beta)
  s = s + delta_s
  i = i + 1
enddo
close(65)
close(55)

s=0.0
delta_s = 0.01
nslices = floor(lat%param%total_length / delta_s)
i=1
allocate(partial_a(nslices))
do while (s .lt. lat%param%total_length)
  call twiss_and_track_at_s(lat, s, ele_at_s, orb, orb_at_s)
  rho = value_of_attribute(ele_at_s,"RHO",err_flag,.false.)
  if (.not. err_flag) then
    partial_a(i) = ele_at_s%a%eta / rho * delta_s
  endif
  s = s + delta_s
enddo

do i=1,nslices
  
enddo

end program












