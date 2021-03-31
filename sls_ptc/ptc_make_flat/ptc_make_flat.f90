program ptc_make_flat

use bmad
use ptc_layout_mod
use madx_ptc_module

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: co(:)

integer status

call getarg(1, lat_file)

call bmad_parser(lat_file, lat)
call twiss_and_track(lat,co,status)

call lat_to_ptc_layout(lat)

call write_ptc_flat_file_lattice("ptc.lat",lat%branch(0))

end program ptc_make_flat












