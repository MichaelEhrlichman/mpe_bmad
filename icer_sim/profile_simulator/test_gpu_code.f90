program test_random_numbers

  use cudafor
  use curand
  implicit none

  integer, parameter :: dp = kind(1.0d0)
  type(curandGenerator) :: gen
  integer i, istat
  real(dp), parameter :: two_pi = 2 * 3.141592653

  call setup_random_gpu(gen)

  do i=1,10
    write(*,*) xi(gen)
  enddo

  ! Clean up
  istat = curandDestroyGenerator(gen)

contains

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
