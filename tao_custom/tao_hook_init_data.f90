!+
! Subroutine tao_hook_init_data (do_standard_setup) 
!
! Hook routine to initialize data.
!-

subroutine tao_hook_init_data ()

use tao_init_data_mod, dummy => tao_hook_init_data

implicit none

integer i,j

!
do i=1,size(s%u)
  do j=1,size(s%u(i)%data)
    if(trim(s%u(i)%data(j)%data_type) == 'tune_fp') then
      s%u(i)%data(j)%exists = .true.
    endif
  enddo
enddo

end subroutine




