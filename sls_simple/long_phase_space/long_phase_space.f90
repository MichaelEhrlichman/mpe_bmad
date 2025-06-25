program long_phase_space_program

use bmad
use omp_lib

implicit none

character lat_file*200

type(lat_struct) lat
type(lat_struct), allocatable :: omp_lat(:)
type(coord_struct), allocatable :: orb(:)
type(coord_struct), allocatable :: co(:)
type(coord_struct), allocatable :: tbt(:)

real(rp), allocatable :: data(:,:,:,:)

integer h, i, j, k, l, m
integer Np_a, Np_b, Nturns
integer track_state
integer plane
integer status
integer alive
integer nok

real(rp) zi, zf, dz
real(rp) di, df, dd
real(rp) d_pt, z_pt

logical printed

integer omp_n, omp_i
logical first
integer n_ele_track

Np_a = 400
Np_b = 400
zi = -0.50
zf =  0.50
di = -0.08
df =  0.04
Nturns = 750

dz = (zf-zi)/(Np_a-1)
dd = (df-di)/(Np_b-1)

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)

call bmad_parser(lat_file, lat)
call set_on_off(rfcavity$, lat, on$)

n_ele_track = lat%n_ele_track

allocate(co(0:n_ele_track))
co(0)%vec = 0.0d0
call twiss_and_track(lat,co,status)


open(45,file='phase_space_z.out')
write(45,'(3a6,3a14)') "#", "i", "j", "z0", "p0", 'alive'

open(47,file='dead_or_alive.out')
write(47,'(3a6,3a14)') "#", "i", "j", "z0", "p0", 'alive'

allocate(data(Np_a,Np_b,Nturns,7))

omp_n = omp_get_max_threads()
allocate(omp_lat(omp_n))
do i=1,omp_n
  omp_lat(i) = lat
enddo
write(*,*) 'omp_get_max_threads(): ', omp_n

first = .true.

!$OMP PARALLEL DO &
!$OMP SCHEDULE(dynamic,1), &
!$OMP DEFAULT(SHARED), &
!$OMP FIRSTPRIVATE(first), &
!$OMP PRIVATE(i,j,k,l,z_pt,d_pt,nok,track_state,orb,tbt,printed,alive,omp_i)
do m=0, Np_a*Np_b-1
  omp_i = omp_get_thread_num()+1
  if (first) then
    allocate(orb(0:n_ele_track))
    allocate(tbt(Nturns))
    first = .false.
  endif
  i = m/Np_b + 1
  j = mod(m,Np_b) + 1
  ! effectively nexted i=1,Np_a, j=1,Np_b
    z_pt = zi + dz * (i-1)
    d_pt = di + dd * (j-1)
    orb(0) = co(0)
    orb(0)%vec(5) = orb(0)%vec(5) + z_pt
    orb(0)%vec(6) = orb(0)%vec(6) + d_pt
    write(*,'(a,i3,a,6f12.7)') "thread ", omp_i, " is tracking: ", orb(0)%vec
    nok = 0
    do k=1,Nturns
      call track_all(omp_lat(omp_i),orb,track_state=track_state)
      if( track_state .eq. moving_forward$ ) then
        nok = nok + 1
        tbt(k) = orb(0)
        orb(0) = orb(n_ele_track)
      else
        write(*,'(a,2i8)') "Particle lost at (turn, ele): ", k, track_state
        exit
      endif
    enddo
    if( track_state .eq. moving_forward$ ) then
      write(47,'(a6,2i6,2f14.5,i6)') '', i, j, z_pt, d_pt, 1
    else
      write(47,'(a6,2i6,2f14.5,i6)') '', i, j, z_pt, d_pt, 0
    endif
    printed=.false.
    do l=1,nok-1
      data(i,j,l,1:6) = tbt(l)%vec
      if( track_state .eq. moving_forward$ ) then
        alive = 1
      else
        alive = 0
      endif
      data(i,j,l,7) = alive
      write(45,'(a6,3i6,2f14.5,6es15.4,i6)') '', i, j, l, z_pt, d_pt, tbt(l)%vec, alive
      printed=.true.
    enddo
    if(printed) then
      write(45,*)
      write(45,*)
    endif
  !effective do loop
enddo
!$OMP END PARALLEL DO
close(45)
close(47)

!open(46,file='phase_space_z_transpose.out')
!write(46,'(4a6,9a14)') "#", "i", "j", "turn", "z0", "p0", "x", "px", "y", "py", "z", "pz", 'alive'
!do k=1,Nturns
!  do i=1,Np_a
!    do j=1,Np_b
!      z_pt = zi + dz*(j-1)
!      d_pt = di + dd*(i-1)
!      write(46,'(a6,3i6,2f14.5,6es15.4,i6)') '', i, j, k, z_pt, d_pt, data(i,j,k,1:6), int(data(i,j,k,7))
!    enddo
!  enddo
!  write(46,*)
!  write(46,*)
!enddo
!close(46)

deallocate(data)

end program






