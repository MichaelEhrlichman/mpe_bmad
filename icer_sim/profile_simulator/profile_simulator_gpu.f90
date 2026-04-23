program profile_simulator

use openacc
use cudafor
use curand
use openacc

implicit none

integer, parameter :: dp = 8

real(dp), parameter :: two_pi = 2.0d0 * 3.14159265359d0
real(dp), parameter :: clight = 299792458.0d0 ! m/s
real(dp), parameter :: rc = 2.817d-15 ! m, classical e- radius
real(dp), parameter :: hbar = 6.582119569d-16 ! eV s
real(dp), parameter :: me = 511000.0d0 ! eV / c^2
real(dp), parameter :: Cq = 3.832e-13 ! m for e-

integer, parameter :: n_particles = 1000
integer, parameter :: BATCH_SIZE = 100 ! Tune based on GPU memory
integer, parameter :: n_batches = (n_particles + BATCH_SIZE - 1) / BATCH_SIZE
integer batch_id, local_id, global_id


real(dp) :: z_batch(BATCH_SIZE), pz_batch(BATCH_SIZE)
real(dp) :: z_history(BATCH_SIZE, n_turns/n_report)
real(dp) :: pz_history(BATCH_SIZE, n_turns/n_report)
integer :: save_counter
type(curandGenerator) :: gen

real(dp) C0
real(dp) g0, E0
real(dp) a_kick
real(dp) Jz
real(dp) kd, kf
real(dp) cell_kick
real(dp) mod_amp
real(dp) I2, I3
real(dp) sigma_p
real(dp) barrier_z
real(dp) mod_z

real(dp) gbend(3)
real(dp) L0(3)
real(dp) sigma_pz
real(dp) n_turns_in !namelist doesn't support exponential integers
real(dp) diffusion_parameter

integer n_report
integer(8) i, n_turns
integer j, k
integer particle_id
integer ierr
integer progress, last_progress

character(20) in_file, out_file
character(4) ix_str

namelist /params/ z0, pz0, sigma_pz, n_turns_in, n_report

call setup_random_gpu(gen)


last_progress = 0

!beam description
n_turns_in=20.0e6
z0 = 0.0d0
pz0 = 0.0d0 !0.01 ! initial pz
sigma_pz = -1
n_report = 1


call get_command_argument(1, in_file)

open(unit=10, file=in_file, status='old', action='read') !, iostat=ierr)
read(10, nml=params) !, iostat=ierr)
close(10)

n_turns = n_turns_in

!lattice description
C0 = 119.587d0
g0 = 1957.0d0 ! relativistic gamma
E0 = g0 * me
Jz = 1.569d0 ! damping partition number
       ! bend_sup, bend_r1, bend_main
gbend = (/ 1.0d0/1.88d0, 1.0d0/12.2d0, 1.0d0/2.5d0 /)
L0 = (/ 2.0d0, 14.4d0, 16.0d0 /)

kd = g0**3 * 2.0d0 * rc / 3.0d0 * (Jz/2.0d0)
kf = g0**5 * 55.0*rc*hbar/24.0/sqrt(3.0)/(me/clight/clight)/clight

!calculate nominal restoring kick
barrier_z = 40.0 !C0 * 0.90
mod_z = 10.0
cell_kick = 0.0d0
I2 = 0.0d0
I3 = 0.0d0
do j=1,3
  I2 = I2 + gbend(j)**2*L0(j)
  I3 = I3 + gbend(j)**3*L0(j)
  cell_kick = cell_kick + kd * L0(j) * gbend(j)**2
enddo
!write(*,*) "Cell kick: ", cell_kick
sigma_p = sqrt(Cq * g0**2 * I3 / Jz / I2)
write(*,'(a,es10.3,a,f10.3,a)') "Rad Int Energy Spread: ", sigma_p*100, "% (", sigma_p*E0/1000, " keV)"


!modulation amplitude
!mod_amp = 0.0000050d0
mod_amp = 0.0000200d0
write(*,'(a,f9.6,a,f8.3,a)') "Modulation amplitude is ", mod_amp*100, "% (", E0*mod_amp/1000, " keV)"
write(*,'(a,f8.5,a,f9.5,a)') "Modulation period is ", mod_z, " meters (", clight / mod_z / 1e6, " MHz)"

!$acc data copyin(gbend, L0, kd, kf, barrier_z, mod_z, mod_amp, C0)

do batch_id = 1, n_batches
  !$acc parallelize loop
  do local_id = 1, BATCH_SIZE
    global_id = (batch_id - 1) * BATCH_SIZE + local_id
    if (global_id > n_particles) exit

    write(out_file,'(a,i6.6,a)') 'z_', global_id, '.dat'
    file_units(local_id) = 1000 + local_id
    open(file_units(local_id), file=out_file)
    write(file_units(local_id),'(a14,a14,a14)') "# turn", "z", "pz"
  enddo

  !initialize batch
  !$acc parallel loop
  do local_id = 1, BATCH_SIZE
    global_id = (batch_id - 1) * BATCH_SIZE + local_id
    z_batch(local_id) = 0.0d0
    pz_batch(local_id) = 0.0d0
  enddo

  !$acc data copy(z_batch, pz_batch)

  do i=1, n_turns
    !$acc parallel loop private(z, pz, a_kick, j)
    do local_id = 1, BATCH_SIZE
      z = z_batch(local_id)
      pz = pz_batch(local_id)

      do j=1,3
        a_kick = kick_acc(pz, gbend(j), L0(j), z, &
                          (batch_id-1)*BATCH_SIZE + local_id, i, j, gen)
        pz = pz + a_kick
      enddo
      pz = pz + barrier_acc(z)
      z = z + alpha_acc(pz) * pz * C0
      z = modulo(z + C0/2.0d0, C0) - C0/2.0d0

      z_batch(local_id) = z
      pz_batch(local_id) = z
    enddo
    !$acc end parallel loop

    !Periodically write to disk
    if (mod(i, n_report) == 0) then
      !$acc update host(z_batch, pz_batch)

      do local_id = 1, BATCH_SIZE
        global_id = (batch_id - 1) * BATCH_SIZE + local_id
        if (global_id > n_particles) exit
        write(file_units(local_id),'(i14,2es16.7)') i, &
                z_batch(local_id), pz_batch(local_id)
      enddo
    endif
  enddo
  !$acc end data

contains

  function alpha_gpu(pz) result(alpha)
    real(dp), intent(in) :: pz
    real(dp) alpha
    real(dp), parameter: a1 = -8.39e-5
    real(dp), parameter: a2 = -9.605e-6
    real(dp), parameter: a3 = 1.37e-4
    real(dp), parameter: a4 = -2.07e-1
    !alpha = 1e-6 * pz - 2e-2 * pz**2 + 1e-0 * pz**3
    !tru zero
    ! a1 = 2.6112753e-07
    ! a2 = -3.5277239e-07
    ! a3 = 6.6528554e-05
    ! gp fit to v6
    alpha = a1 + a2*pz + a3*pz**2 + a4*pz**3
  end function

  function barrier_gpu(z) result(barrier)
    real(dp), intent(in) :: z
    real(dp) barrier
    real(dp) barrier_width, barrier_kick, reset_kick
    real(dp) U0, pz0

    barrier_width = 5.0
    U0 = 45.3e3 
    pz0 = U0 / E0
    barrier_kick = 0.01 * pz0
    reset_kick = -10.0 * pz0

    if (abs(z) .lt. barrier_z) then
      barrier = 0.0d0
    elseif (abs(z) .lt. barrier_z+barrier_width) then
      barrier = barrier_kick * sign(1.0d0,z)
    else
      barrier = reset_kick
    endif
  end function

  function kick_gpu(pz,gbend,L0,z, particle_id, turn, j, gen) result(kick)
    real(dp), intent(in) :: pz, gbend, L0, z
    integer, intent(in) :: particle_id, turn, j
    real(dp) :: kick, kick_d, kick_f, kick_restore
    real(dp) gbend, L0
    type(curandGenerator) :: gen
    !kick_d = -1.0d0 * kd * (gbend**2) * L0 * ( (1.0d0+pz)**2 - 1.0d0 )  !E restored

    !kick_d = -1.0d0 * kd * (gbend**2) * L0 * ( 1.0d0 + 2.0d0*pz + pz*pz )
    kick_d = -1.0d0 * kd * (gbend**2) * L0 * (1.0d0+pz)**2
    kick_f = -sqrt(kf*(gbend**3)*L0) * xi(gen) * (1.0d0+pz)**2

    kick_restore = kd * (gbend**2) * L0 * (1.0d0 + mod_amp*sin(two_pi * z / mod_z))

    kick = kick_d + kick_f + kick_restore
  end function

  subroutine setup_random_gpu(gen)
    type(curandGenerator), intent(out) :: gen
    integer :: istat
    integer :: time_values(8)
    integer(8) :: seed

    ! Create the cuRAND generator
    istat = curandCreateGenerator(gen, CURAND_RNG_PSEUDO_DEFAULT)

    call date_and_time(values=time_values)
    seed = int(time_values(8) + time_values(7)*1000 + time_values(6)*60000 + &
               time_values(5)*3600000, 8)

    istat = curandSetPseudoRandomGeneratorSeed(gen, seed)
  end subroutine

  subroutine generate_uniform_gpu(gen, n, d_array)
    type(curandGenerator) :: gen
    integer :: n
    real(dp), device :: d_array(*)
    integer :: istat
    
    istat = curandGenerateUniformDouble(gen, d_array, n)
  end subroutine

  function xi(gen)
    type(curandGenerator), intent(inout) :: gen
    real(dp) xi
    real(dp) u1, u2, z1, z2, istat
    real(dp), device :: d_random(2)

    call generate_uniform_gpu(gen, 2, d_random)
    u1 = d_random(1)
    u2 = d_random(2)
    z1 = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
    !z2 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2)
    xi = z1
  end function
end program
