module namelist_raster

! provides:
! linear_bouds_file, x_min, x_max, y_min, y_max, nx, ny, n_naff

use bmad

implicit none

!character*10 calc_tunes
logical use_linear_bounds
real(rp) x_min, x_max, y_min, y_max
real(rp) pz_min, pz_max
integer nx, ny, n_naff, n_turns
integer npz
integer, parameter :: max_dE = 3
real(rp) dE(max_dE)
logical radiation

namelist / da_raster /  use_linear_bounds, &
                        radiation, &
                        dE, &
                        x_min, &
                        x_max, &
                        y_min, &
                        y_max, &
                        nx, &
                        ny, &
                        n_turns, &
                        n_naff, &
                        pz_min, &
                        pz_max, &
                        npz
end module
