program for_dcs

use bmad

implicit none

character lat_file*200

type(lat_struct) lat
character(40) loc
type (ele_pointer_struct), allocatable :: eles(:)
integer n_loc

loc = 'quad::sec1:sec2'

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call lat_ele_locator(loc, lat, eles, n_loc)

end program











