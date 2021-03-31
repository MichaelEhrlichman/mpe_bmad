program phase_space_program

use bmad
use sls_lib

implicit none

character lat_file*200
character de_str*5

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(coord_struct), allocatable :: co(:)
type(coord_struct), allocatable :: tbt(:)

integer h, i, j, k
integer Nic, Nturns
integer track_state
integer plane
integer status

real(rp) de
real(rp) xi, xf, dx

Nic = 100
xi = 0.0001
xf = 0.005
Nturns = 500

dx = (xf-xi)/(Nic-1)

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
if( iargc() == 2 ) then
  call getarg(2, de_str)
  read(de_str,*) de
else
  de_str = '0'
  de = 0.0d0
endif

call bmad_parser(lat_file, lat)
call set_on_off(rfcavity$, lat, off$)

allocate(co(0:lat%n_ele_track))
co(0)%vec = 0.0d0
co(0)%vec(6) = de
call twiss_and_track(lat,co,status)

allocate(orb(0:lat%n_ele_track))
allocate(tbt(Nturns))

!FOO do h=1,2
do h=2,2
  if(h == 1) then
    open(45,file='phase_space_'//trim(de_str)//'_x.out')
    plane = 1
  elseif(h == 2) then
    open(45,file='phase_space_'//trim(de_str)//'_y.out')
    plane = 3
  endif
  write(45,'(a6,a6,a14,6a14)') "#", "turn", "delta", "x", "px", "y", "py", "z", "pz"

  do i=1,Nic
    orb(0) = co(0)
    orb(0)%vec(plane) = orb(0)%vec(plane) + xi + dx*(i-1)
    write(*,'(a,6f12.7)') "Tracking ", orb(0)%vec
    do j=1,Nturns
      call track_all(lat,orb,track_state=track_state)
      if( track_state .eq. moving_forward$ ) then
        tbt(j) = orb(0)
        orb(0) = orb(lat%n_ele_track)
      else
        write(*,'(a,2i8)') "Particle lost at (turn, ele): ", j, track_state
        exit
      endif
    enddo
!    if( track_state .eq. moving_forward$ ) then
      do k=1,Nturns
        write(45,'(i6,i6,f14.5,6es14.4)') i, k, xi+dx*(i-1), tbt(k)%vec
      enddo
      write(45,*)
      write(45,*)
!    endif
  enddo
  close(45)
enddo


end program






