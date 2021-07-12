program inj_dist_tao

use bmad

implicit none

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)

integer i,j,k,ix
integer status
integer Nx, Ny, Nz
integer, parameter :: grid_nx = 200
integer, parameter :: grid_ny = 200
integer, parameter :: nsigx = 3
integer, parameter :: nPx = 400

#define ORDINARY

!integer, parameter :: Nmc = 3000
#ifdef SIGM_RING
integer, parameter :: Nmc = nsigx * nPx + 1
character(9), parameter :: dist_type = 'sigm_ring'
#endif
#ifdef GRID_GRID
integer, parameter :: Nmc = grid_nx*grid_nx
character(9), parameter :: dist_type = 'grid_grid'
#endif
#ifdef ORDINARY
integer, parameter :: Nmc = 9000
character(9), parameter :: dist_type = 'flat_wted'
!dist_type = 'true_gaus'
!dist_type = 'flat_wted'
!dist_type = 'flat_flat'
#endif

real(rp) x(Nmc), y(Nmc), z(Nmc)
real(rp) xp(Nmc), yp(Nmc), zp(Nmc)
real(rp) Npart(Nmc)
real(rp) grid_x_min, grid_x_max, grid_y_min, grid_y_max

real(rp) initial_offset(6)
real(rp) r
real(rp) phix, phiy, phiz
real(rp) Jx, Jy, Jz
real(rp) emit_x, emit_y, emit_z
real(rp) long_beta, long_gamma
real(rp) Ntot

real(rp) betax, alphax, betay, alphay
real(rp) etax,  etaxp, etay,  etayp  

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

ix = 0

initial_offset = 0.0d0

grid_x_min = -0.010   !-0.020
grid_x_max = 0.010   !0.00015625
grid_y_min = -0.002  !-0.0040
grid_y_max = 0.002  !0.0040

!old
!betax  = 14.66225186
!alphax = -0.09219613
!betay  =  5.23212891
!alphay = -0.27565469
!etax   =  0.11628667
!etaxp  =  0.0 
!etay   =  0.0
!etayp  =  0.0

!  !Booster
!  ! betax  = 14.81922113
!  betax  = 8.0
!  alphax = -0.13891292
!  betay  = 5.70144704
!  alphay = -0.41533194
!  etax   = 0.11628667
!  etaxp  =  0.0 
!  etay   =  0.0
!  etayp  =  0.0
!  emit_x = 300.0d-9  !ALS-U CDR Booster emittance
!  !emit_x = 1200.0d-9  !4x ALS-U CDR Booster emittance
!  emit_y = 30.0d-9
!  emit_z = 2.500e-5 !1.0e-3 * 25 mm
!  long_beta = 12.5
!  long_gamma = 0.08

!Unit
!betax  = 1.0
!alphax = 0.0
!betay  = 1.0
!alphay = 0.0
!etax   = 0.0
!etaxp  =  0.0 
!etay   =  0.0
!etayp  =  0.0
!emit_x = 300.0d-9  !ALS-U CDR Booster emittance
!emit_y = 30.0d-9
!emit_z = 2.500e-5 !1.0e-3 * 25 mm
!long_beta = 12.5
!long_gamma = 0.08

!!AR Equilibrium at element zero
betax  = 14.53867132
alphax = 0.0
betay  = 4.86263882
alphay = 0.0
etax   =  0.11628667
etaxp  =  0.0 
etay   =  0.0
etayp  =  0.0
emit_x = 1.806E-09
emit_y = 1.806E-11
emit_z = 4.367e-6
long_beta = 6.0
long_gamma = 0.167

Nx = 16 
Ny = 16
Nz = 16

!Produce monte carlo distributions
call random_seed()

Ntot = 0.0d0

if(dist_type == 'grid_grid') then
  do j=1,grid_nx
    do k=1,grid_ny
      ix = grid_ny*(j-1) + k
      x(ix) = grid_x_min + (grid_x_max - grid_x_min)/(grid_nx-1.0d0)*(j-1.0d0)
      xp(ix) = grid_y_min + (grid_y_max - grid_y_min)/(grid_ny-1.0d0)*(k-1.0d0)
      y(ix) = 0.0d0
      yp(ix) = 0.0d0
      z(ix) = 0.0d0
      zp(ix) = 0.0d0
      Npart(ix) = 1.0
    enddo  
  enddo
elseif(dist_type == 'sigm_ring') then
  x(1) = 0; xp(1) = 0
  y(1) = 0; yp(1) = 0
  z(1) = 0; zp(1) = 0
  Npart(1) = 1.0
  do i=1,nsigx
    do k=1,nPx
      ix = nPx*(i-1) + k + 1
      phix = twopi*(k-1.0)/(nPx-1.0)
       x(ix) = (2.0d0*i)/2.0d0*sqrt(emit_x*betax)*sin(phix)
      xp(ix) = (2.0d0*i)/2.0d0*sqrt(emit_x/betax)*( -alphax*sin(phix) + cos(phix) )
       y(ix) = 0.0d0  !sqrt(2*Jy*betay)*sin(phiy)
      yp(ix) = 0.0d0  !sqrt(2*Jy/betay)*( -alphay*sin(phiy) + cos(phiy) )
       z(ix) = 0.0d0  !sqrt(2*Jz*long_beta)*sin(phiz)
      zp(ix) = 0.0d0  !sqrt(2*Jz/long_beta)*cos(phiz) !longitudinal alpha is assumed very small
      Npart(ix) = 1.0
    enddo
  enddo
else
  do i=1, Nmc
    if(dist_type == 'true_gaus') then
      call random_number(r); Jx = -emit_x*log(1.0d0-r) !inverse transform sampling to get action J distrubition.
      call random_number(r); Jy = -emit_y*log(1.0d0-r)
      call random_number(r); Jz = -emit_z*log(1.0d0-r)

      call random_number(r); phix = r*twopi
      call random_number(r); phiy = r*twopi
      call random_number(r); phiz = r*twopi

      Npart(i) = 1.0
    elseif(dist_type == 'flat_wted') then
      call random_number(r); Jx = emit_x*r*Nx !flat distribution, weighted by exponential in J
      call random_number(r); Jy = emit_y*r*Ny
      call random_number(r); Jz = emit_z*r*Nz

      call random_number(r); phix = r*twopi
      call random_number(r); phiy = r*twopi
      call random_number(r); phiz = r*twopi

      Npart(i) = exp(-Jx/emit_x)*exp(-Jy/emit_y)*exp(-Jz/emit_z) / emit_x / emit_y / emit_z
    elseif(dist_type == 'flat_flat') then
      call random_number(r); Jx = emit_x*r*Nx !flat distribution, weighted by exponential in J
      call random_number(r); Jy = emit_y*r*Ny
      call random_number(r); Jz = emit_z*r*Nz

      call random_number(r); phix = r*twopi
      call random_number(r); phiy = r*twopi
      call random_number(r); phiz = r*twopi

      Npart(i) = 1.0
    elseif(dist_type == 'ring_wted') then
      call random_number(r); Jy = 0.0 !0.01*emit_y  !0.1*emit_y*r*Ny
      call random_number(r); Jz = 0.0 !0.01*emit_z  !0.1*emit_z*r*Nz
      call random_number(r); phiy = r*twopi
      call random_number(r); phiz = r*twopi
      k = (i-1) / nsigx + 1
      j = i - nsigx*(k-1)
      Jx = 4.0*emit_x*(1.0d0*j)/(1.0d0*nsigx) 
      phix = twopi*(k-1.0d0)/nPx
      Npart(i) = exp(-Jx/emit_x)*exp(-Jy/emit_y)*exp(-Jz/emit_z) / emit_x / emit_y / emit_z
    else
      write(*,*) "bomb"
      stop
    endif

     x(i) = sqrt(2*Jx*betax)*sin(phix)
    xp(i) = sqrt(2*Jx/betax)*( -alphax*sin(phix) + cos(phix) )
     y(i) = sqrt(2*Jy*betay)*sin(phiy)
    yp(i) = sqrt(2*Jy/betay)*( -alphay*sin(phiy) + cos(phiy) )
     z(i) = sqrt(2*Jz*long_beta)*sin(phiz)
    zp(i) = sqrt(2*Jz/long_beta)*cos(phiz) !longitudinal alpha is assumed very small
  enddo
endif
Ntot = sum(Npart)

open(45,file="dist_for_tao_mc.dat")
!tao header
write(45,'(i6,a)') ix, "   !ix_ele ignored"
write(45,'(i6,a)') 1, "   !number of bunches"
write(45,'(i6,a)') Nmc, "   !number of particles per bunch"
write(45,'(a)') "BEGIN_BUNCH"
write(45,'(a)') "ELECTRON"
write(45,'(es9.1,a)') 1.0,  "   !charge per bunch"
write(45,'(f3.1,a)') 0.0, "   !z position at center of bunch"
write(45,'(f3.1,a)') 0.0, "   !t position at center of bunch"
do i=1,Nmc
  write(45,'(6es12.3,es11.3,i3)') initial_offset(1) + zp(i)*etax  + x(i), &
                                initial_offset(2) + zp(i)*etaxp + xp(i), &
                                initial_offset(3) + zp(i)*etay  + y(i), &
                                initial_offset(4) + zp(i)*etayp + yp(i), &
                                initial_offset(5) + z(i), &
                                initial_offset(6) + zp(i), &
                                Npart(i)/Ntot, 1
enddo
write(45,'(a)') "END_BUNCH"
close(45)

end program




