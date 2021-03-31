module namelist_fm

! provides:
! dE, xmax, ymax, n_turns, fft_turns, nx, ny

use precision_def

implicit none

real(rp) dE
real(rp) xmax
real(rp) ymax
integer pre_turns
integer naff_turns
integer nx, ny

namelist / fm / dE, &
                xmax, &
                ymax, &
                nx, &
                ny, &
                pre_turns, &
                naff_turns

end module
