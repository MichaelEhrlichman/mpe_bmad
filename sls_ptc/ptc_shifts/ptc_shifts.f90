PROGRAM ptc_shifts

USE bmad
USE ptc_layout_mod
USE madx_ptc_module

IMPLICIT none

CHARACTER lat_file*200

TYPE(lat_struct) lat
TYPE(taylor_struct) one_turn_map(6)
TYPE(taylor_struct) dhdj(6)
TYPE(complex_taylor_struct) :: F(6), L(6)
!TYPE(complex_taylor_struct) :: H_lie_poly_bas(6)
TYPE(coord_struct), ALLOCATABLE :: co(:)

!parameters
REAL(rp) pz
LOGICAL rf_on

CALL GETARG(1, lat_file)

WRITE(*,*) "Preparing lattice..."

CALL bmad_parser(lat_file, lat)
CALL set_on_off(rfcavity$,lat,on$)

bmad_com%radiation_damping_on=.false.
CALL lat_to_ptc_layout(lat)

rf_on = .false.

pz = 0.0
WRITE(*,'(A,F10.4)') "Calculating terms for pz = ", pz
CALL ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, rf_on, pz)
CALL normal_form_taylors(one_turn_map, rf_on, dhdj)

!CALL normal_form_complex_taylors(one_turn_map, rf_on, H_lie_poly=H_lie_poly_bas)
open(21,file="tune_shifts.dat")
write(21,'(A,11ES14.4)') "Qx_pz   ", &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 0]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 1]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 2]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 3]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 4]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 5]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 6]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 7]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 8]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 9]), &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 10])
write(21,'(A,11ES14.4)') "Qy_pz   ", &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 0]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 1]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 2]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 3]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 4]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 5]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 6]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 7]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 8]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 9]), &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 10])
write(21,'(A,6ES14.4)') "Qx_Jx   ", &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 0]), &
                      taylor_coef(dhdj(1), [2, 0, 0, 0, 0, 0])*2.0, &
                      taylor_coef(dhdj(1), [4, 0, 0, 0, 0, 0])*4.0, &
                      taylor_coef(dhdj(1), [6, 0, 0, 0, 0, 0])*6.0, &
                      taylor_coef(dhdj(1), [8, 0, 0, 0, 0, 0])*24.0, &
                      taylor_coef(dhdj(1), [10, 0, 0, 0, 0, 0])*120.0
WRITE(21,'(A,6ES14.4)') "Qx_Jy   ", &
                      taylor_coef(dhdj(1), [0, 0, 0, 0, 0, 0]), &
                      taylor_coef(dhdj(1), [0, 0, 2, 0, 0, 0])*2.0, &
                      taylor_coef(dhdj(1), [0, 0, 4, 0, 0, 0])*4.0, &
                      taylor_coef(dhdj(1), [0, 0, 6, 0, 0, 0])*6.0, &
                      taylor_coef(dhdj(1), [0, 0, 8, 0, 0, 0])*24.0, &
                      taylor_coef(dhdj(1), [0, 0, 10, 0, 0, 0])*120.0
write(21,'(A,6ES14.4)') "Qy_Jy   ", &
                      taylor_coef(dhdj(2), [0, 0, 0, 0, 0, 0]), &
                      taylor_coef(dhdj(2), [0, 0, 2, 0, 0, 0])*2.0, &
                      taylor_coef(dhdj(2), [0, 0, 4, 0, 0, 0])*4.0, &
                      taylor_coef(dhdj(2), [0, 0, 6, 0, 0, 0])*6.0, &
                      taylor_coef(dhdj(2), [0, 0, 8, 0, 0, 0])*24.0, &
                      taylor_coef(dhdj(2), [0, 0, 10, 0, 0, 0])*120.0
write(21,'(A,7ES14.4)') "z_pz   ", &
                      taylor_coef(dhdj(6), [0, 0, 0, 0, 0, 1]), &
                      taylor_coef(dhdj(6), [0, 0, 0, 0, 0, 2]), &
                      taylor_coef(dhdj(6), [0, 0, 0, 0, 0, 3]), &
                      taylor_coef(dhdj(6), [0, 0, 0, 0, 0, 4]), &
                      taylor_coef(dhdj(6), [0, 0, 0, 0, 0, 5]), &
                      taylor_coef(dhdj(6), [0, 0, 0, 0, 0, 6]), &
                      taylor_coef(dhdj(6), [0, 0, 0, 0, 0, 7])
close(21)

END PROGRAM ptc_shifts












