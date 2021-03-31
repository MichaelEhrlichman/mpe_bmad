program grid
  use bmad
  use bmad_parser_mod, only: bp_com

  implicit none

  type(lat_struct) ring
  type(coord_struct), allocatable :: orb(:) 
  character*100 lat_file
  integer i,j,k
  integer nx, ny
  real(rp) x0, y0
  real(rp) xmin, xmax
  real(rp) ymin, ymax

  call getarg(1,lat_file)

  bp_com%always_parse = .true.
  call bmad_parser(lat_file, ring)
  bmad_com%radiation_damping_on = .false.
  bmad_com%radiation_fluctuations_on = .false.

  allocate(orb(0:ring%n_ele_track))

  nx = 11
  ny = 11 
  xmin = -0.01
  xmax =  0.01
  ymin = 0.00
  ymax = 0.01

  open(20,file='ought.dat')
  open(21,file='fin.dat')
  do j=1,nx
    x0 = xmin + (xmax-xmin)/(nx-1)*(j-1)
    do k=1,ny
      y0 = ymin + (ymax-ymin)/(ny-1)*(k-1)
      orb(0)%vec = 0.0d0
      orb(0)%vec(1) = x0
      orb(0)%vec(3) = y0
      call track_all(ring,orb)
      write(20,'(6es17.8)') orb(0)%vec
      write(21,'(6es17.8)') orb(ring%n_ele_track)%vec
    enddo
  enddo
  close(20)
  close(21)

!for making bucket
!   n_turns = 1000
!   n_pts = 60
!   pz_min = -0.07
!   pz_max = 0.07
!   allocate(tbt_coords(1:n_turns))
! 
!   call calc_ring(ring,dims,co,err)
!   open(21,file='z_pz_dead.dat')
!   open(22,file='z_pz_alive.dat')
!   do j=1,n_pts
!     pz0 = pz_min + (pz_max-pz_min)/(n_pts-1)*(j-1)
! 
!     orb(0) = co(0)
!     orb(0)%vec(6) = orb(0)%vec(6) + pz0
!     tbt_coords(1) = orb(0)
!     dead = .false.
!     do k=1,n_turns-1
!       call track_all(ring,orb,track_state=track_state)
!       tbt_coords(k+1) = orb(ring%n_ele_track)
!       if (track_state /= moving_forward$) then
!         dead = .true.
!         exit
!       endif
!       orb(0) = orb(ring%n_ele_track)
!     enddo
! 
!     do k=1,n_turns
!       if(dead) then
!         write(21,'(6es17.8)') tbt_coords(k)%vec
!       else
!         write(22,'(6es17.8)') tbt_coords(k)%vec
!       endif
!     enddo
! 
!     if(dead) then
!       write(21,*)
!       write(21,*)
!     else
!       write(21,*)
!       write(22,*)
!     endif
!   enddo
!   close(21)
!   close(22)

end program





