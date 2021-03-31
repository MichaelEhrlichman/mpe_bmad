module namelist_adts

! provides:
! adts: x_min, x_max, y_min, y_max, use_bounds_file, n_steps

use bmad

implicit none

real(rp) x_min, x_max, y_min, y_max
logical use_linear_bounds
integer n_steps

namelist / adts / x_min, &
                  x_max, &
                  y_min, &
                  y_max, &
                  use_linear_bounds, &
                  n_steps

end module
