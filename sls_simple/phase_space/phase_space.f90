program phase_space_program

use bmad
use sls_lib
use all_phase_fft

implicit none

character lat_file*200
character de_str*5

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(coord_struct), allocatable :: co(:)
type(coord_struct), allocatable :: tbt(:)
type(coord_struct), allocatable :: norm_coords(:)
type(coord_struct) orb_in

integer h, i, j, k
integer Nic, Nturns
integer track_state
integer main_plane, secondary_plane
integer status
integer fft_a, fft_b

real(rp) de
real(rp) xi, xf, dx

logical ok

!apfft
real(rp) x_bounds(2), y_bounds(2)
real(rp) x_phase, x_amp, x_freq
real(rp) y_phase, y_amp, y_freq

real(rp) small_offset !small offset to sample tune

Nic = 100
xi = 0.0001
xf = 0.003
Nturns = 1000
fft_a = 500
fft_b = 1000

x_bounds(1) = 0.5
x_bounds(2) = 1.0
y_bounds(1) = 0.5
y_bounds(2) = 1.0

small_offset = 1e-6

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
allocate(norm_coords(Nturns))

do h=1,2
  if(h == 1) then
    open(45,file='phase_space_'//trim(de_str)//'_x.out')
    open(46,file='tunes_'//trim(de_str)//'_x.out')
    open(55,file='phase_space_'//trim(de_str)//'_x_norm.out')
    main_plane = 1
    secondary_plane = 3
  elseif(h == 2) then
    open(45,file='phase_space_'//trim(de_str)//'_y.out')
    open(46,file='tunes_'//trim(de_str)//'_y.out')
    open(55,file='phase_space_'//trim(de_str)//'_y_norm.out')
    main_plane = 3
    secondary_plane = 1
  endif
  write(45,'(a6,a6,a14,6a14)') "#", "turn", "delta", "x", "px", "y", "py", "z", "pz"
  write(46,'(a12,6a14)') "#particle", "x_phase", "x_amp", "x_freq", "y_phase", "y_amp", "y_freq"
  write(55,'(a6,a6,a14,4a14)') "#", "turn", "delta", "a", "pa", "b", "pb"

  do i=1,Nic
    orb_in = co(0)
    orb_in%vec(main_plane) = orb_in%vec(main_plane) + xi + dx*(i-1)
    orb_in%vec(secondary_plane) = orb_in%vec(secondary_plane) + small_offset
    call init_coord(orb(0), orb_in)
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
    if( track_state .eq. moving_forward$ ) then
      do k=1,Nturns
        write(45,'(i6,i6,f14.5,6es14.4)') i, k, xi+dx*(i-1), tbt(k)%vec
        call xy_to_action(lat, 0, tbt(k)%vec, norm_coords(k)%vec, ok) 
        write(55,'(i6,i6,f14.5,4es14.4)') i, k, xi+dx*(i-1), norm_coords(k)%vec(1:4)
      enddo
      write(45,*)
      write(45,*)
      write(55,*)
      write(55,*)
      call apfft_corr(norm_coords(fft_a:fft_b)%vec(1), x_bounds, 'han', x_phase, x_amp, x_freq)
      call apfft_corr(norm_coords(fft_a:fft_b)%vec(3), y_bounds, 'han', y_phase, y_amp, y_freq)
      write(46,'(i12,6f14.5)') i, x_phase, x_amp, x_freq, y_phase, y_amp, y_freq
    endif
  enddo
  close(45)
  close(46)
  close(55)
enddo


end program






