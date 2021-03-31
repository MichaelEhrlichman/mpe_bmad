program tunes

use bmad
use bmad_parser_mod, only: bp_com

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: co(:)
integer status

call getarg(1, lat_file)
bp_com%always_parse = .true.
call bmad_parser(lat_file, lat)
call set_on_off(rfcavity$, lat, off$)

call twiss_and_track(lat,co,status)

write(*,'(2f14.5,es14.5)') lat%ele(lat%n_ele_track)%a%phi/2./pi, lat%ele(lat%n_ele_track)%b%phi/2./pi, lat%ele(0)%a%eta

end program












