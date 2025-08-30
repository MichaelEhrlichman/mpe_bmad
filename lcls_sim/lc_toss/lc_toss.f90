program lc_toss
  use bmad
  use bmad_parser_mod, only: bp_com
  use slice_mod

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:)
  type(coord_struct), allocatable ::  slorb(:)

  integer i,j,k
  integer n_found, bpm_ix, lc_ix, track_state
  integer nsteps

  character(20) set_string

  logical err_flag

  real(rp), allocatable :: slices(:)
  real(rp) vec_offset(6)
  real(rp) initial_x
  real(rp) s, delta_s, sf, si
  real(rp) phi0, phi
  real(rp) phi_min, phi_max, delta_phi
  real(rp) E0, V, KE

  character*100 lat_file

  bp_com%always_parse = .true.

  call getarg(1,lat_file)
  call bmad_parser(lat_file, lat)

  lc_ix = 1
  phi0 = value_of_attribute(lat%ele(lc_ix),'PHI0')
  write(*,*) "Initial cavity phase: ", phi0, " (rad./2pi)"

!  E0 = value_of_attribute(lat%ele(0),'E_TOT')
!  KE = sqrt(E0*E0-m_electron*m_electron)
!  V = value_of_attribute(lat%ele(lc_ix),'VOLTAGE')
!  write(*,*) "Initial Energy: ", E0
!  write(*,*) "Initial Kinetic Energy: ", KE
!  write(*,*) "Cavity Voltage: ", V
!  delta_phi = acos(-KE/V)/twopi
!  write(*,*) "delta_phi: ", delta_phi
  phi_min = -0.25 !-delta_phi
  phi_max =  0.25 !delta_phi
  write(*,*) "phi_min: ", phi_min
  write(*,*) "phi_max: ", phi_max

  nsteps=51

  bpm_ix = 1
  write(*,'(a,i6,a,f14.4,a)') "BPM is located at index ", bpm_ix, ".  (s = ", lat%ele(bpm_ix)%s, " (m) )"

  allocate(orb(0:lat%n_ele_track))

!
! vs phi
!

!initial_x = 50e-6
!vec_offset = 0.0d0
!vec_offset(1) = initial_x
vec_offset = (/ initial_x, 0.0d0, 0.0000d0, 0.0d0, 0.07143d0, -0.00612d0 /)
open(20,file='lc_toss_vs_phi.dat')
write(20,'(a,es14.5,a)') "#Nominal cavity phase is ", phi0*twopi, " radians."
write(20,'(a8,a19,a15,6a15,a15)') "# step", "phase (rad)", "xi", "x", "xp", "y", "yp", "z", "zp", "p0c"
do i=2,nsteps-1
  phi = phi_min + (phi_max-phi_min)/(nsteps-1.0)*(i-1.0)
  write(set_string,'(a,es14.4)') 'phi0=', phi
  call set_ele_attribute(lat%ele(lc_ix),set_string,err_flag)
  if(err_flag) then
    write(*,*) "Error in set_ele_attribute."
  endif
  call lattice_bookkeeper(lat,err_flag)
  if(err_flag) then
    write(*,*) "Error in lattice bookkeeper."
    stop
  endif
  call init_coord(orb(0),vec_offset,lat%ele(0),element_end=upstream_end$)
  call track_all(lat,orb)
  write(20,'(i8,f19.4,es15.6,6es15.6,es15.6)') i, phi*twopi, initial_x, orb(bpm_ix)%vec, orb(bpm_ix)%p0c
enddo
close(20)

  !
  ! vs xi
  !

  phi = phi0
  write(set_string,'(a,es14.4)') 'phi0=', phi
  write(*,*) "set_string: ", set_string
  call set_ele_attribute(lat%ele(lc_ix),set_string,err_flag)
  if(err_flag) then
    write(*,*) "Error in set_ele_attribute."
  endif
  call lattice_bookkeeper(lat,err_flag)
  if(err_flag) then
    write(*,*) "Error in lattice bookkeeper."
    stop
  endif

  open(20,file='lc_toss_vs_x.dat')
  write(20,'(a,es14.5,a)') "# Cavity phase is ", phi*twopi, " radians."
  write(20,'(a8,a19,a15,6a15,a15)') "# step", "phase (rad)", "xi", "x", "xp", "y", "yp", "z", "zp", "p0c"
  do i=1,nsteps
    initial_x = 500.0e-6 * (i-1.0) / (nsteps-1.0)
    vec_offset = (/ initial_x, 0.0d0, 0.0000d0, 0.0d0, 0.0d0, 0.0d0 /)
    call init_coord(orb(0),vec_offset,lat%ele(0),element_end=upstream_end$)
    call track_all(lat,orb, track_state=track_state)
    if(track_state .ne. moving_forward$) then
      write(*,*) "particle lost at initial_x = ", initial_x
    else
      write(20,'(i8,f19.4,es15.6,6es15.6,es15.6)') i, phi*twopi, initial_x, orb(bpm_ix)%vec, orb(bpm_ix)%p0c
    endif
  enddo
  close(20)

!
! trajectory through cavity
!

nsteps=20

allocate(slices(nsteps+1))
allocate(slorb(nsteps+1))
initial_x = 50.0d-6
vec_offset = (/ initial_x, 0.0d0, 0.0000d0, 0.0d0, 0.0d0, 0.0d0 /)
call init_coord(slorb(1),vec_offset,lat%ele(lc_ix-1),element_end=downstream_end$)
si = lat%ele(lc_ix-1)%s
sf = lat%ele(lc_ix)%s
do i=1,nsteps+1
  slices(i) = si + (sf-si)/(nsteps)*(i-1)
enddo
call track_s_to_s(lat,slices,1,nsteps+1,slorb)
open(20, file='lc_toss_s.dat')
write(20,'(a,es14.5,a)') "# Cavity phase is ", phi0*twopi, " radians."
write(20,'(a8,a11,a15,6a19,a19)') "# step", "s", "x", "xp", "y", "yp", "z", "zp", "p0c"
do i=1,nsteps+1
  write(20,'(i8,f11.4,6es19.10,es19.10)') i, slices(i), slorb(i)%vec, slorb(i)%p0c
enddo
close(20)

end program





