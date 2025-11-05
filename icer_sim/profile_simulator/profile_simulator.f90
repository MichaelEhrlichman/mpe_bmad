program profile_simulator

implicit none

integer, parameter :: dp = 8

real(dp), parameter :: two_pi = 2.0d0 * 3.14159265359d0
real(dp), parameter :: clight = 299792458.0d0 ! m/s
real(dp), parameter :: rc = 2.817d-15 ! m
real(dp), parameter :: hbar = 6.582119569d-16 ! eV s
real(dp), parameter :: me = 511000.0d0 ! eV / c^2

real z, z0
real pz, pz0
real C0
real g0
real a_kick
real Jz
real kd, kf

real gbend(3)
real L0(3)
real L0_
real sigma_pz
real n_turns_in !namelist doesn't support exponential integers
real diffusion_parameter

integer n_report
integer(8) n_turns
integer i, j, k
integer particle_id
integer ierr

character(20) in_file, out_file
character(3) ix_str

logical diags

namelist /params/ z0, pz0, sigma_pz, n_turns_in, diags, n_report

call random_seed()

!beam description
n_turns_in=20.0e6
z0 = 0.0
pz0 = 0.0 !0.01 ! initial pz
sigma_pz = -1
diags = .false.
n_report = 1

call get_command_argument(1, in_file)
call get_command_argument(2, ix_str)
read(ix_str,'(i3)') particle_id

open(unit=10, file=in_file, status='old', action='read') !, iostat=ierr)
read(10, nml=params) !, iostat=ierr)
close(10)

n_turns = n_turns_in

if(pz0 .lt. 0.0d0) then
  if(sigma_pz > 0) then
    pz0 = xi() * sigma_pz 
  else
    write(*,*) "pz0 < 0, but sigma_pz not set in .in file"
    stop
  endif
endif

!lattice description
C0 = 119.587
g0 = 1957.0 ! relativistic gamma
Jz = 1.569 ! damping partition number
       ! bend_sup, bend_r1, bend_main
gbend = (/ 1.0/1.88, 1.0/12.2, 1.0/2.5 /)
L0 = (/ 2.0, 14.4, 16.0 /)

kd = g0**3 * 2.0d0 * rc / 3.0d0 * (Jz/2.0d0)
kf = g0**5 * 55.0*rc*hbar/24.0/sqrt(3.0)/(me/clight/clight)/clight

if(diags) then
  diffusion_parameter = 0.0d0
  do i=1,3
    diffusion_parameter = diffusion_parameter + sqrt(kf*(gbend(i)**3)*L0(i))
  enddo
  write(*,*) "diffusion parameter: ", diffusion_parameter
  !plot alpha(pz) for diagnostic output
  open(10,file='sample_zpz.dat')
  do i=1,500
    pz = -0.05 + (i-1.0)*(0.05+0.05)/(500.0-1.0)
    write(10,'(es14.4,es14.4)') pz, pz*alpha(pz)*C0
  enddo
  close(10)
endif

pz = pz0
z = z0
write(out_file,'(a,i3.3,a)') 'z_',particle_id,'.dat'
open(10,file=out_file)
write(10,'(a14,a14,a14)') "# turn", "z", "pz"
do i=1,n_turns
  if (mod(i,n_report) == 0) then
    write(10,'(i14,2es14.4)') i, z, pz
  endif
  do j=1,3
    L0_ = L0(j)
    a_kick = kick(pz, gbend(j), L0_)
    pz = pz + a_kick
  enddo
  z = z + alpha(pz)*pz * C0 +3.30237e-07
enddo
close(10)

contains

  function alpha(pz)
    real pz, alpha
    real a1, a2, a3, a4
    !alpha = 1e-6 * pz - 2e-2 * pz**2 + 1e-0 * pz**3
    !tru zero
    ! a1 = 2.6112753e-07
    ! a2 = -3.5277239e-07
    ! a3 = 6.6528554e-05
    ! gp fit to v6
    a1 = -8.39e-5
    a2 = -9.605e-6
    a3 = 1.37e-4
    a4 = -2.07e-1
    alpha = a1 + a2*pz + a3*pz**2 + a4*pz**3
  end function

  function kick(pz,gbend,L0)
    real kick, pz
    real kick_d, kick_f
    real gbend, L0

    kick_d = -1.0d0 * kd * (gbend**2) * L0 * ( (1.0d0+pz)**2 - 1.0d0 )
    kick_f = -sqrt(kf*(gbend**3)*L0) * xi() * (1.0+pz)**2
    kick = (kick_d + kick_f)
  end function

  function xi()
    real(dp) u1, u2, z1, z2, xi
    call random_number(u1)
    call random_number(u2)
    z1 = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
    !z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2)
    xi = z1
  end function
end program
