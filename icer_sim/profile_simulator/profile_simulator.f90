program profile_simulator
implicit none

integer, parameter :: dp = 8

integer, parameter :: N_BEND_MAX = 5

type :: lat_struct
  real(dp) C0  !storage ring circumference
  real(dp) E0  !beam energy
  real(dp) Jz  !longitudinal damping partition
  real(dp) :: B0(N_BEND_MAX) = -1.0
  real(dp) :: L0(N_BEND_MAX) = -1.0
  real(dp) :: a(5) = (/ 0.0, 0.0, 0.0, 0.0, 0.0 /)
end type 

type(lat_struct) lat

real(dp), parameter :: two_pi = 2.0d0 * 3.14159265359d0
real(dp), parameter :: clight = 299792458.0d0 ! m/s
real(dp), parameter :: rc = 2.817d-15 ! m, classical e- radius
real(dp), parameter :: hbar = 6.582119569d-16 ! eV s
real(dp), parameter :: me = 511000.0d0 ! eV / c^2
real(dp), parameter :: Cq = 3.832e-13 ! m for e-

real(dp) gbend
real(dp) rrr
real(dp) z
real(dp) pz, pz0
real(dp) g0
real(dp) a_kick
real(dp) kd, kf
real(dp) cell_kick
real(dp) rf_amp, rf_lambda
real(dp) I2, I3
real(dp) sigma_p
real(dp) barrier_z, barrier_kick

real(dp) sigma_pz
real(dp) n_turns_in !namelist doesn't support exponential integers

integer n_report
integer(8) i, n_turns
integer j, k
integer particle_id
integer ierr
integer progress, last_progress

character(20) in_file, out_file
character(5) ix_str

namelist /params/ pz0, rf_amp, rf_lambda, sigma_pz, n_turns_in, n_report, barrier_z, barrier_kick, lat

call random_seed()

last_progress = 0

!beam description
n_turns_in=5.0e6
pz0 = 0.0d0 !0.01 ! initial pz
sigma_pz = -1
n_report = 1
rf_amp = 1.0 ! eV
rf_lambda = 10.0 ! m
barrier_z = 40.0 !C0 * 0.90
barrier_kick = 100.0

call get_command_argument(1, in_file)
call get_command_argument(2, ix_str)
read(ix_str,'(i5)') particle_id

open(unit=10, file=in_file, status='old', action='read') !, iostat=ierr)
read(10, nml=params) !, iostat=ierr)
close(10)

n_turns = n_turns_in

write(*,'(a,i12,a)') "Tracking for ", n_turns, " turns."

if(pz0 .lt. 0.0d0) then
  if(sigma_pz > 0) then
    pz0 = xi() * sigma_pz 
  else
    write(*,*) "pz0 < 0, but sigma_pz not set in .in file"
    stop
  endif
endif

!lattice description
g0 = lat%E0 / me  !391.39 ! 200 MeV

kd = g0**3 * 2.0d0 * rc / 3.0d0 * (lat%Jz/2.0d0)
kf = g0**5 * 55.0*rc*hbar/24.0/sqrt(3.0)/(me/clight/clight)/clight

!calculate nominal restoring kick
cell_kick = 0.0d0
I2 = 0.0d0
I3 = 0.0d0
do j=1,N_BEND_MAX
  if(lat%B0(j) > 0) then 
    gbend = 1.0d0/lat%B0(j)
    I2 = I2 + gbend**2 * lat%L0(j)
    I3 = I3 + gbend**3 * lat%L0(j)
    cell_kick = cell_kick + kd * gbend**2 * lat%L0(j)
  endif
enddo
write(*,*) "Cell kick: ", cell_kick
sigma_p = sqrt(Cq * g0**2 * I3 / lat%Jz / I2)
write(*,'(a,es10.3,a,f10.3,a)') "Rad Int Energy Spread: ", sigma_p*100, "% (", sigma_p*lat%E0/1000, " keV)"

!modulation amplitude
write(*,'(a,es10.4,a)') "Modulation amplitude is ", rf_amp, " eV"
rf_amp = rf_amp / lat%E0
write(*,'(a,es10.4,a,es10.4,a)') "Modulation period is ", rf_lambda, " meters (", clight / rf_lambda / 1e6, " MHz)"

pz = 0 !FOO pz0
call random_number(rrr)
z = barrier_z - 2*barrier_z * rrr ! uniform distribution in the barrier
!z=0
write(out_file,'(a,i5.5,a)') 'z_',particle_id,'.dat'
open(10,file=out_file)
write(10,'(a14,a14,a14)') "# turn", "z", "pz"
do i=1,n_turns
  if (mod(i,n_report) == 0) then
    write(10,'(i14,2es16.7)') i, z, pz
  endif
  progress = int(100.0 * i / n_turns)
  if (progress >= last_progress + 10) then
    write(*,'(i12,a)') progress, "% complete"
    last_progress = progress
  endif
  do j=1,N_BEND_MAX
    if(lat%B0(j) > 0) then 
      a_kick = kick(pz, lat%B0(j), lat%L0(j), z)
      pz = pz + a_kick
    endif
  enddo
  pz = pz + barrier(z)
  z = z + alpha(pz)*pz*lat%C0
  z = modulo(z+lat%C0/2.0d0, lat%C0) - lat%C0/2.0
enddo
close(10)
write(*,*) "Complete!"

contains

  function alpha(pz)
    real(dp) alpha, pz
    
    alpha = lat%a(1) + lat%a(2)*pz + lat%a(3)*pz**2 + lat%a(4)*pz**3
  end function

  function barrier(z)
    real(dp) barrier, z
    real(dp) barrier_width
    real(dp) reset_kick

    barrier_width = barrier_z * 0.05 
    reset_kick = -10.0 * barrier_kick

    if (abs(z) .lt. barrier_z) then
      barrier = 0.0d0
    elseif (abs(z) .lt. barrier_z+barrier_width) then
      write(*,*) "BARRIER KICK"
      barrier = barrier_kick * sign(1.0d0,z)
    else
      write(*,*) "RESET KICK"
      barrier = reset_kick
    endif
  end function

  function kick(pz,B0,L0,z)
    real(dp) kick, pz, z
    real(dp) kick_d, kick_f, kick_rest_cell, kick_rest_rf
    real(dp) B0, gbend, L0
    !kick_d = -1.0d0 * kd * (gbend**2) * L0 * ( (1.0d0+pz)**2 - 1.0d0 )  !E restored
    !kick_d = -1.0d0 * kd * (gbend**2) * L0 * ( 1.0d0 + 2.0d0*pz + pz*pz )

    gbend = 1.0/B0
    kick_d = -1.0d0 * kd * gbend**2 * L0 * (1.0d0+pz)**2   ! Damping
    kick_f = -sqrt(kf*(gbend**3)*L0) * xi() * (1.0d0+pz)**2  ! Quantum excitation
    kick_rest_cell = kd * gbend**2 * L0  ! Restoring kick from the cell (no modulation)
    kick_rest_rf = rf_amp * sin(two_pi * z / rf_lambda)  ! Kick from RF modulation

    kick = kick_d + kick_f + kick_rest_cell !+ kick_rest_rf
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











