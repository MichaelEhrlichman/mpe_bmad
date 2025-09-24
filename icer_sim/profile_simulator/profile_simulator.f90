
program profile_simulator

implicit none

real, parameter :: two_pi = 2.0 * 3.14159265359
real, parameter :: clight = 299792458.0 ! m/s
real, parameter :: rc = 2.828e-15 ! m
real, parameter :: hbar = 6.582119569e-16 ! eV s
real, parameter :: me = 511000.0 ! eV / c^2

real z
real d
real C0, L0
real g0
real gbend

integer n
integer i

n=2000000
C0 = 119.0
g0 = 2000.0 ! relativistic gamma

gbend = 1.0 / 2.5  ! gbend = 1/bend_radius
L0 = 40.  ! sum length over all bends

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
  d = d + kick(d)
  z = z + alpha(d) * C0
enddo
close(10)

contains

  function alpha(d)
    real d, alpha
    alpha = 1e-6 * d - 2e-2 * d**2 + 1e-0 * d**3
  end function

  function kick(d)
    real kick, d
    real kf

    kf = g0**5 * 55.0*rc*hbar/24.0/sqrt(3.0)/(me/clight/clight)/clight
    kick = -sqrt(kf*(gbend**3)*L0) * xi() * (1.0+d)**2
    
  end function

  function xi()
    real u1, u2, z1, z2, xi
    call random_seed()
    call random_number(u1)
    call random_number(u2)
    z1 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2)
    !z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2)
    xi = z1
  end function
end program
