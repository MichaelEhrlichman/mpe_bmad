
program profile_simulator

implicit none

integer, parameter :: dp = 8

real(dp), parameter :: two_pi = 2.0d0 * 3.14159265359d0
real(dp), parameter :: clight = 299792458.0d0 ! m/s
real(dp), parameter :: rc = 2.817d-15 ! m
real(dp), parameter :: hbar = 6.582119569d-16 ! eV s
real(dp), parameter :: me = 511000.0d0 ! eV / c^2

real z
real d
real C0
real g0
real a_kick

real gbend(3)
real L0(3)
real L0_

integer n
integer i, j, k
integer nslices

n=10e6
nslices=1
C0 = 119.587
g0 = 1957.0 ! relativistic gamma

       ! bend_sup, bend_r1, bend_main
gbend = (/ 1.0/1.88, 1.0/12.2, 1.0/2.5 /)
L0 = (/ 2.0, 14.4, 16.0 /)

!plot alpha(d) for diagnostic output
open(10,file='alpha.dat')
do i=1,500
  d = -0.05 + (i-1.0)*(0.05+0.05)/(500.0-1.0)
  write(10,'(es14.4,es14.4)') d, alpha(d)
enddo
close(10)

d = 0
z = 0
open(10,file='z.dat')
do i=1,n
  write(10,'(2es14.4)') d,z
  do j=1,3
    do k=1,nslices
      L0_ = L0(j) / nslices
      a_kick = kick(d, gbend(j), L0_)
      d = d + a_kick
    enddo
  enddo
  z = z + alpha(d) * C0
enddo
close(10)

contains

  function alpha(d)
    real d, alpha
    real a1, a2, a3
    !alpha = 1e-6 * d - 2e-2 * d**2 + 1e-0 * d**3
    a1 = 2.6112753e-07
    a2 = -3.5277239e-07
    a3 = 6.6528554e-05
    alpha = a1*d + a2*d**2 + a3*d**3
  end function

  function kick(d,gbend,L0)
    real kick, d
    real kick_d, kick_f
    real gbend, L0
    real kd, kf

    kf = g0**5 * 55.0*rc*hbar/24.0/sqrt(3.0)/(me/clight/clight)/clight
    kd = g0**3 * 2.0d0 * rc / 3.0d0

    kick_d = -1.0d0 * kd * (gbend**2) * L0 * ( (1.0d0+d)**2 - 1.0d0 ) * (1.569d0/2.0d0)
    kick_f = -sqrt(kf*(gbend**3)*L0) * xi() * (1.0+d)**2
    kick = kick_d + kick_f
  end function

  function xi()
    real(dp) u1, u2, z1, z2, xi
    call random_seed()
    call random_number(u1)
    call random_number(u2)
    z1 = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
    !z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2)
    xi = z1
  end function
end program
