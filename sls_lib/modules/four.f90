subroutine xy_to_action(ring, ix, X, J, ok)

  use bmad
  use sim_utils

  type(lat_struct) ring
  integer ix, i
  real(rp) X(6)
  real(rp) J(6)
  logical ok

  real(rp) Afour(6,6)
  real(rp) Afour_inv(6,6)
  type(ele_struct), pointer :: ele

  real(rp) sqrtabeta
  real(rp) sqrtbbeta

  ele => ring%ele(ix)

  sqrtabeta = sqrt(ele%a%beta)
  sqrtbbeta = sqrt(ele%b%beta)
  !write(*,*) "FOO: ", sqrtabeta, sqrtbbeta

  Afour = 0.0d0
  Afour(1,1) = sqrtabeta
  Afour(2,1) = -ele%a%alpha/sqrtabeta
  Afour(2,2) = 1.0d0/sqrtabeta
  Afour(1,6) = ele%x%eta
  Afour(2,6) = ele%x%etap
  Afour(3,3) = sqrtbbeta
  Afour(4,3) = -ele%b%alpha/sqrtbbeta
  Afour(4,4) = 1.0d0/sqrtbbeta
  Afour(3,6) = ele%y%eta
  Afour(4,6) = ele%y%etap
  Afour(5,5) = 1.0d0
  Afour(6,6) = 1.0d0

  ! do i=1,6
  !   write(*,'(6ES14.4)') Afour(i,:)
  ! enddo

  !call mat_inverse(Afour,Afour_inv, ok)
  Afour_inv = my_mat_symp_conj(Afour)
  ok = .true.

  J = matmul(Afour_inv,X)

contains

function my_mat_symp_conj(mat) result (mat_conj)

use output_mod, only: rp, out_io, s_fatal$, global_com

implicit none

integer i, j, nn, dims

real(rp) mat(:,:)
real(rp) :: mat_conj(size(mat, 1), size(mat, 1))
real(rp), allocatable :: S(:,:)

character(*), parameter :: r_name = 'mat_symp_conj'

! Check bounds

nn = size(mat, 1)

if (mod(nn, 2) /= 0 .or. nn /= size(mat, 2)) then
  call out_io (s_fatal$, r_name, 'ARRAY SIZE IS NOT EVEN!')
  if (global_com%exit_on_error) call err_exit
endif

dims = nn/2
allocate(S(nn,nn))
S=0.0d0
do i=1,nn,2
  S(i,i+1) =  1.0d0
  S(i+1,i) = -1.0d0
enddo

mat_conj = -1.0d0*matmul(S,matmul(transpose(mat),S))

deallocate(S)

end function


end subroutine

