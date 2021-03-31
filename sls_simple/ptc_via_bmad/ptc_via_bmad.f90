program ptc_via_bmad

use bmad
use ptc_layout_mod
use madx_ptc_module

implicit none

character lat_file*200

type(lat_struct), target :: lat
type(coord_struct), allocatable :: co(:)

type (normal_form_struct), pointer :: normal_form

integer i
integer j(6)
integer status

complex(rp) c_val

call getarg(1, lat_file)

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)
call set_on_off(rfcavity$,lat,off$)
call twiss_and_track(lat, co, status)

normal_form => lat%branch(0)%normal_form_no_rf
normal_form%ele_origin => lat%branch(0)%ele(0)

call lat_to_ptc_layout(lat)

call ptc_one_turn_map_at_ele (normal_form%ele_origin, normal_form%m, pz = 0.0_rp )
call alloc(normal_form%vb)

call normal_form_rd_terms(normal_form, .false.)

do i=1,size(normal_form%rd_term)
  if(normal_form%rd_term(i)%F_index == 0) then
    write(*,'(i1," h",6i1,2es14.5,es14.5)') normal_form%rd_term(i), abs(normal_form%rd_term(i)%c_val)
  endif
enddo

end program












