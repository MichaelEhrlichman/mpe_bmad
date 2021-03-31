module namelist_fp

! provides:
! pz_min, pz_max, n_pz

use bmad

implicit none

real(rp) pz_min, pz_max
integer n_pz

namelist / fp /         pz_min, &
                        pz_max, &
                        n_pz

end module
