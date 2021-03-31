program footprint_simple

use bmad

implicit none

character(100) lat_file
real(rp) pz_min, pz_max
integer n_pz

character*60 in_file

type(lat_struct) lat
type(coord_struct), allocatable :: co(:), co_prev(:), co0(:)
real(rp) detector, detector_y
real(rp) orbit_change_rms, orbit_change_rms_prev
real(rp) orbit_change_y_rms, orbit_change_y_rms_prev
real(rp) :: detector_limit = 1.0d-6
real(rp), allocatable :: eta0(:), eta_y0(:)

integer i_prev
integer i, j, k
integer status

real(rp) dpz, pz
real(rp) nu_a, nu_b

namelist /fp/ lat_file, pz_min, pz_max, n_pz

call getarg(1,in_file)

open(10, file = in_file, action='read')
read(10, nml = fp)
close(10)

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)
allocate(co(0:lat%n_ele_track))
allocate(co0(0:lat%n_ele_track))
allocate(co_prev(0:lat%n_ele_track))
allocate(eta0(0:lat%n_ele_track))
allocate(eta_y0(0:lat%n_ele_track))

call set_on_off(rfcavity$,lat,off$)
bmad_com%radiation_damping_on=.false.

call twiss_and_track(lat,co,status)
eta0=lat%ele(:)%x%eta
eta_y0=lat%ele(:)%y%eta
co0 = co

open(21,file='footprint.dat')
write(21,'(a1,a10,2a14)') "#", "pz", "nu_x", "nu_y"
open(45,file='co.dat')
open(55,file='detector.dat')
do k=1,2
  if(k .eq. 1) then
    !assume pz_max is positive
    dpz = pz_max/n_pz
  else
    !assume pz_min is negative
    dpz = pz_min/n_pz
  endif
  i_prev = 0
  co_prev=co0
  orbit_change_rms_prev = 0.0d0
  orbit_change_y_rms_prev = 0.0d0
  do i=1,n_pz
    pz = dpz*i
    co(0)%vec = 0.0d0
    co(0)%vec(6) = pz
    call twiss_and_track(lat,co,status)
    if(status .eq. ok$) then
      do j=1,lat%n_ele_track
        write(45,'(i6,6es14.4)') j, co(j)%vec(:)
      enddo
      write(45,*)
      write(45,*)
      orbit_change_rms   = sqrt( sum( (co(:)%vec(1) - co_prev(:)%vec(1) + eta0(:)*(i-i_prev)*dpz)**2)/size(co(:)) )
      orbit_change_y_rms = sqrt( sum( (co(:)%vec(3) - co_prev(:)%vec(3) + eta_y0(:)*(i-i_prev)*dpz)**2)/size(co(:)) )
      if(i .gt. 1) then
        detector = abs(orbit_change_rms_prev - orbit_change_rms)/(i-i_prev)
        detector_y = abs(orbit_change_y_rms_prev - orbit_change_y_rms)/(i-i_prev)
        write(55,*) pz, detector, detector_y
      else
        detector = 0.0d0
        detector_y = 0.0d0
      endif
      if( (detector.lt.detector_limit) .and. (detector_y.lt.detector_limit)) then
        co_prev = co
        orbit_change_rms_prev = orbit_change_rms 
        orbit_change_y_rms_prev = orbit_change_y_rms 
        i_prev = i
      else
        status = no_closed_orbit$
      endif
    endif
    if(status .eq. ok$) then
      nu_a = lat%ele(lat%n_ele_track)%a%phi/twopi
      nu_b = lat%ele(lat%n_ele_track)%b%phi/twopi
      write(21,'(f18.11,2f14.6,2es14.4)') pz, nu_a, nu_b
    endif
  enddo
  call clear_lat_1turn_mats(lat)
enddo
close(21)
close(45)
close(55)

end program












