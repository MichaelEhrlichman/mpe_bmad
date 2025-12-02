program diffusion_simulator

use mpi

implicit none

integer, parameter :: dp = 8

real(dp), parameter :: two_pi = 2.0d0 * 3.14159265359d0
real(dp), parameter :: clight = 299792458.0d0 ! m/s
real(dp), parameter :: rc = 2.817d-15 ! m
real(dp), parameter :: hbar = 6.582119569d-16 ! eV s
real(dp), parameter :: me = 511000.0d0 ! eV / c^2

real(dp), allocatable :: beam(:,:,:)
integer, allocatable :: block_tracker(:)

real z
real d, d0
real C0
real g0
real a_kick
real Jz

real gbend(3)
real L0(3)
real L0_

integer n_turns, n_turns_rpt, n_rpts
integer n_part, part_ix
integer i, j, k

integer :: seed_size, un
integer, allocatable :: seed(:)

!mpi housekeeping
integer num_workers, my_worker_num
integer task_at_hand, n_tasks
integer myrank, from_id
integer cluster_size
integer mpierr
integer mpistatus(MPI_STATUS_SIZE)
integer id, worker_status
integer n_sent, n_recv
integer mpi_i, incoming_ix
integer tag, source

logical tasker

! mpi stuff
call mpi_init(mpierr)                             ! Introduce yourself to the MPI daemon
call mpi_comm_rank(MPI_COMM_WORLD,myrank,mpierr)  ! Get your rank number, store in myrank.  Master is rank 0.
if(myrank .eq. 0) then
  tasker=.true.
  call mpi_comm_size(MPI_COMM_WORLD,cluster_size,mpierr)
  num_workers=cluster_size-1
  write(*,*) "Number of slaves: ", num_workers
else
  tasker=.false.
endif

! seed random number generator
call random_seed(size=seed_size)
allocate(seed(seed_size))
open(newunit=un, file="/dev/urandom", access="stream", &
     form="unformatted", action="read", status="old")
read(un) seed
close(un)
call random_seed(put=seed)

n_part=25
n_turns=5e6
n_turns_rpt=1000  !workers report back every . turns
n_rpts = n_turns / n_turns_rpt
C0 = 119.587
g0 = 1957.0 ! relativistic gamma
Jz = 1.569 ! damping partition number
d0 = 0.0 !0.01 ! initial pz

allocate(beam(n_part,n_rpts,2))  ! dim3 is z, pz
allocate(block_tracker(n_part))  ! dim3 is z, pz
block_tracker(:) = 0

       ! bend_sup, bend_r1, bend_main
gbend = (/ 1.0/1.88, 1.0/12.2, 1.0/2.5 /)
L0 = (/ 2.0, 14.4, 16.0 /)

!plot alpha(d) for diagnostic output
open(10,file='sample_zpz.dat')
do i=1,500
  d = -0.05 + (i-1.0)*(0.05+0.05)/(500.0-1.0)
  write(10,'(es14.4,es14.4)') d, d*alpha(d)*C0
enddo
close(10)

if(tasker) then
  !initialize beam
  do i=1,n_part
    beam(i,1,1) = d0
    beam(i,1,2) = 0.0d0
  enddo

  !mpi tags
  ! 1: sending particle index from tasker to worker
  ! 2: sending z,pz from tasker to worker
  ! 3: receiving particle index from worker to tasker
  ! 4: receiving z,pz from worker to tasker
  ! 5: termination signal from taster to worker

  do mpi_i=1,min(num_workers,n_part)
    n_sent = n_sent + 1
                !data,        , ndata, type,  to_who, tag
    call mpi_send(mpi_i, 1, MPI_INTEGER, mpi_i, 1, MPI_COMM_WORLD, mpierr)
    call mpi_send(beam(i,1,1:2), 2, MPI_DOUBLE_PRECISION, mpi_i, 2, MPI_COMM_WORLD, mpierr)
  enddo

  do while(n_sent < n_part)
    call mpi_recv(part_ix, 1, MPI_INTEGER, MPI_ANY_SOURCE, 3, MPI_COMM_WORLD, mpistatus, mpierr)  !blocking
    source = mpistatus(MPI_SOURCE)
    block_tracker(part_ix) = block_tracker(part_ix) + 1
    call mpi_recv(beam(part_ix,block_tracker(part_ix),1:2), 2, MPI_DOUBLE_PRECISION, source, 4, MPI_COMM_WORLD, mpistatus, mpierr)
    n_sent = n_sent + 1
    call mpi_send(n_sent, 1, MPI_INTEGER, source, 1, MPI_COMM_WORLD, mpierr)
    call mpi_send(beam(n_sent,1,1:2), 2, MPI_DOUBLE_PRECISION, source, 2, MPI_COMM_WORLD, mpierr)
  enddo

  ! send termination signal to workers
  do mpi_i=1,num_workers                                                                                                                                                                                                                                       
    call mpi_send(-1, 1, MPI_INTEGER, mpi_i, 5, MPI_COMM_WORLD, mpierr)                                                                                                                                                                                        
  enddo                                                                                                                                                                                                                                                        
  call mpi_finalize(mpierr)       
else
  !get z0,d0 from mpi
  write(filename,'(a,i3.3,a)') 'z_', myrank, '.dat'
  open(10,file=filename)
  do i=1,n_turns_rpt
    write(10,'(2es14.4)') d,z
    do j=1,3
      L0_ = L0(j)
      a_kick = kick(d, gbend(j), L0_)
      d = d + a_kick
    enddo
    z = z + alpha(d)*d * C0 +3.30237e-07
  enddo
close(10)
endif

contains

  function alpha(d)
    real d, alpha
    real a1, a2, a3, a4
    !alpha = 1e-6 * d - 2e-2 * d**2 + 1e-0 * d**3
    !tru zero
    ! a1 = 2.6112753e-07
    ! a2 = -3.5277239e-07
    ! a3 = 6.6528554e-05
    ! gp fit to v6
    a1 = -8.39e-5
    a2 = -9.605e-6
    a3 = 1.37e-4
    a4 = -2.07e-1
    alpha = a1 + a2*d + a3*d**2 + a4*d**3
  end function

  function kick(d,gbend,L0)
    real kick, d
    real kick_d, kick_f
    real gbend, L0
    real kd, kf

    kf = g0**5 * 55.0*rc*hbar/24.0/sqrt(3.0)/(me/clight/clight)/clight
    kd = g0**3 * 2.0d0 * rc / 3.0d0 * (Jz/2.0d0)

    kick_d = -1.0d0 * kd * (gbend**2) * L0 * ( (1.0d0+d)**2 - 1.0d0 )
    kick_f = -sqrt(kf*(gbend**3)*L0) * xi() * (1.0+d)**2
    kick = (kick_d + kick_f)
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
