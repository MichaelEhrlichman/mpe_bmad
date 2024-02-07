program lc_toss
  use bmad
  use bmad_parser_mod, only: bp_com

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:)
  type(coord_struct) orb1
  type(ele_pointer_struct), allocatable :: lc_ele(:)
  type(ele_pointer_struct), allocatable :: bpm_ele(:)

  integer i,j,k
  integer n_found, bpm_ix, lc_ix, track_state
  integer, parameter :: nsteps=50

  character(20) set_string

  logical err_flag

  real(rp) vec_offset(6)
  real initial_x
  real(rp) s, delta_s
  real(rp) phi0, phi
  real(rp) phi_min, phi_max, delta_phi
  real(rp) E0, V, KE

  character*100 lat_file

  bp_com%always_parse = .true.

  call getarg(1,lat_file)
  call bmad_parser(lat_file, lat)

  call lat_ele_locator('CAVL011', lat, lc_ele, n_found)
  if(n_found .ne. 1) then
    write(*,*) "Cavity element not found"
    stop
  endif
  phi0 = value_of_attribute(lc_ele(1)%ele,'PHI0')
  lc_ix = lc_ele(1)%loc%ix_ele
  write(*,*) "Initial cavity phase: ", phi0

  E0 = value_of_attribute(lat%ele(lc_ix-1),'E_TOT')
  KE = sqrt(E0*E0-m_electron*m_electron)
  V = value_of_attribute(lc_ele(1)%ele,'VOLTAGE')
  write(*,*) "Initial Energy: ", E0
  write(*,*) "Initial Kinetic Energy: ", KE
  write(*,*) "Cavity Voltage: ", V
  delta_phi = acos(-KE/V)/twopi
  write(*,*) "delta_phi: ", delta_phi
  phi_min = -delta_phi
  phi_max =  delta_phi
  write(*,*) "phi_min: ", phi_min
  write(*,*) "phi_max: ", phi_max

  call lat_ele_locator('BPM', lat, bpm_ele, n_found)
  if(n_found .ne. 1) then
    write(*,*) "BPM not found"
    stop
  endif
  bpm_ix = bpm_ele(1)%loc%ix_ele
  write(*,'(a,i6,a,f14.4,a)') "BPM is located at index ", bpm_ix, ".  (s = ", bpm_ele(1)%ele%s, " (m) )"

  allocate(orb(0:lat%n_ele_track))

  initial_x = 50e-6
  vec_offset = (/ initial_x, 0.0, 0.0000, 0.0, 0.0, 0.0 /)

  open(20,file='lc_toss_vs_phi.dat')
  write(20,'(a8,a19,a15,6a15,a15)') "# step", "phase (.rad/2pi)", "xi", "x", "xp", "y", "yp", "z", "zp", "p0c"
  do i=2,nsteps-1
    phi = phi_min + (phi_max-phi_min)/(nsteps-1.0)*(i-1.0)
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
    call init_coord(orb(lc_ix-1),vec_offset,lat%ele(lc_ix-1),element_end=downstream_end$)
    call track_many(lat,orb,lc_ix-1,bpm_ix,1,track_state=track_state)
    write(20,'(i8,f19.4,es15.6,6es15.6,es15.6)') i, phi, initial_x, orb(bpm_ix)%vec, orb(bpm_ix)%p0c
  enddo
  close(20)


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
  write(20,'(a8,a19,a15,6a15,a15)') "# step", "phase (.rad/2pi)", "xi", "x", "xp", "y", "yp", "z", "zp", "p0c"
  do i=1,nsteps
    initial_x = 500.0e-6 * (i-1.0) / (nsteps-1.0)
    vec_offset = (/ initial_x, 0.0, 0.0000, 0.0, 0.0, 0.0 /)
    call init_coord(orb(lc_ix-1),vec_offset,lat%ele(lc_ix-1),element_end=downstream_end$)
    call track_many(lat,orb,lc_ix-1,bpm_ix,1,track_state=track_state)
    if(track_state .ne. moving_forward$) then
      write(*,*) "particle lost at initial_x = ", initial_x
    else
      write(20,'(i8,f19.4,es15.6,6es15.6,es15.6)') i, phi, initial_x, orb(bpm_ix)%vec, orb(bpm_ix)%p0c
    endif
  enddo
  close(20)

  !open(20, file='lc_toss_s.dat')
  !initial_x = 500.0e-6 * (i-1.0) / (nsteps-1.0)
  !vec_offset = (/ initial_x, 0.0, 0.0000, 0.0, 0.0, 0.0 /)
  !call init_coord(orb(lc_ix-1),vec_offset,lat%ele(lc_ix-1),element_end=downstream_end$)
  !x=0
  !do i=1,nsteps
  !  call track_s_to_s(lat,orb,lc_ix-1,bpm_ix,1,track_state=track_state)
  !enddo
  !close(20)

end program





