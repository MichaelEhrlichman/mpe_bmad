module namelist_fm

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
