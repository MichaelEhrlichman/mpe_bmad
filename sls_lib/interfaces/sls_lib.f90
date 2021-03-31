module sls_lib

interface

  subroutine binary_search(lost,delta_m,accuracy,reset)
    use bmad
    implicit none
    logical lost
    real(rp) delta_m
    real(rp) accuracy
    logical reset
  end subroutine

  subroutine locate_elements(lat,mask,n_elements,elements)
    use bmad
    implicit none
    type(lat_struct) lat
    character(10) mask
    integer n_elements
    integer, allocatable :: elements(:)
  end subroutine

end interface

end module
