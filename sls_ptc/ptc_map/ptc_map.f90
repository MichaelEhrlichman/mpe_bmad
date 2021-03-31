PROGRAM ptc_map

USE bmad
USE ptc_layout_mod

IMPLICIT none

CHARACTER lat_file*200
CHARACTER(3) order_str
CHARACTER(5) ele_str

TYPE(lat_struct) lat
TYPE(taylor_struct) one_turn_map(6)
TYPE(taylor_struct) dhdj(6)
TYPE(taylor_struct) A(6), A_inverse(6)
TYPE(complex_taylor_struct) A1(6), A1_inverse(6)
TYPE(taylor_struct) :: h, h6(6)
TYPE(complex_taylor_struct) F(6), L(6)

INTEGER i
INTEGER ele_ix

!parameters
REAL(rp) pz

CALL GETARG(1, lat_file)
CALL GETARG(2, ele_str)

READ(ele_str,'(I5)') ele_ix

WRITE(*,*) "Preparing lattice..."

CALL bmad_parser(lat_file, lat)
bmad_com%radiation_damping_on=.false.
CALL lat_to_ptc_layout(lat)

IF( (ele_ix .lt. 0) .or. (ele_ix .gt. lat%n_ele_track) ) THEN
  ele_ix = 0
ENDIF
WRITE(*,*) "Calculating one-turn-map at ele: ", ele_ix

pz = 0.0
call set_ptc(no_cavity=.true.)
CALL ptc_one_turn_map_at_ele(lat%ele(ele_ix), one_turn_map, pz=pz, order=2)
CALL type_taylors(one_turn_map)

!CALL normal_form_taylors(one_turn_map, .false., dhdj, A, A_inverse)
!CALL type_taylors(dhdj)

!CALL normal_form_complex_taylors(one_turn_map,.false., F, L, A1, A1_inverse, order)

!CALL write_ptc_flat_file_lattice("ptc.lat",lat%branch(0))

END PROGRAM ptc_map












