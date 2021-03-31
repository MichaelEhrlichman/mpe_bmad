PROGRAM ptc_h

USE bmad
USE ptc_layout_mod
USE madx_ptc_module

IMPLICIT none

CHARACTER lat_file*200

TYPE(lat_struct) lat
TYPE(taylor_struct) one_turn_map(6)
TYPE(taylor_struct) dhdj(6)
TYPE(complex_taylor_struct) :: F(6), L(6)
TYPE(complex_taylor_struct) :: H_lie_poly_bas(6)
TYPE(complex_taylor_struct) :: H_lie_poly_neg(6)
TYPE(complex_taylor_struct) :: H_lie_poly_pos(6)
TYPE(coord_struct), ALLOCATABLE :: co(:)

INTEGER i, j
INTEGER order
INTEGER hix(6)

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
pz = -0.03
WRITE(*,'(A,F10.4)') "Calculating terms for pz = ", pz
CALL ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, rf_on, pz, order)
CALL normal_form_taylors(one_turn_map, rf_on, dhdj)
CALL normal_form_complex_taylors(one_turn_map, rf_on, H_lie_poly=H_lie_poly_neg)
pz = 0.0
WRITE(*,'(A,F10.4)') "Calculating terms for pz = ", pz
CALL ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, rf_on, pz, order)
CALL normal_form_taylors(one_turn_map, rf_on, dhdj)
CALL normal_form_complex_taylors(one_turn_map, rf_on, H_lie_poly=H_lie_poly_bas)
pz = 0.03
WRITE(*,'(A,F10.4)') "Calculating terms for pz = ", pz
CALL ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, rf_on, pz, order)
CALL normal_form_taylors(one_turn_map, rf_on, dhdj)
CALL normal_form_complex_taylors(one_turn_map, rf_on, H_lie_poly=H_lie_poly_pos)

open(20,FILE="ptc_h_bas.out")
!WRITE(20,'(A10,4A20)') "# Term", "Resonance", "Real Part", "Imaginary Part", "Absolute Value"
hix = [1, 1, 0, 0, 0, 1]
WRITE(20,'(A10,A20,3F20.10)') "h11001", "--",      real(complex_taylor_coef(H_lie_poly_bas(1), hix))/2, aimag(complex_taylor_coef(H_lie_poly_bas(1), hix))/2, abs(complex_taylor_coef(H_lie_poly_bas(1), hix))/2
hix = [0, 0, 1, 1, 0, 1]
WRITE(20,'(A10,A20,3F20.10)') "h00111", "--",      real(complex_taylor_coef(H_lie_poly_bas(1), hix))/2, aimag(complex_taylor_coef(H_lie_poly_bas(1), hix))/2, abs(complex_taylor_coef(H_lie_poly_bas(1), hix))/2
!3rd order terms
hix = [2, 1, 0, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h21000", "Qx",        real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [3, 0, 0, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h30000", "3Qx",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [1, 0, 1, 1, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h10110", "Qx",        real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [1, 0, 0, 2, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h10020", "Qx-2Qy",    real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [1, 0, 2, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h10200", "Qx+2Qy",    real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [2, 0, 0, 0, 0, 1]
WRITE(20,'(A10,A20,3F20.10)') "h20001", "2Qx",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [0, 0, 2, 0, 0, 1]
WRITE(20,'(A10,A20,3F20.10)') "h00201", "2Qy",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
!4th order terms
write(20,*)
write(20,*)
hix = [2, 2, 0, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h22000", "--     ",   real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [0, 0, 2, 2, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h00220", "--     ",   real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [1, 1, 1, 1, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h11110", "--     ",   real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [3, 1, 0, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h31000", "2Qx",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [4, 0, 0, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h40000", "4Qx",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [2, 0, 1, 1, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h20110", "2Qx",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [1, 1, 2, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h11200", "2Qy",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [2, 0, 0, 2, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h20020", "2Qx-2Qy",   real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [2, 0, 2, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h20200", "2Qx+2Qy",   real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [0, 0, 3, 1, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h00310", "2Qy",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
hix = [0, 0, 4, 0, 0, 0]
WRITE(20,'(A10,A20,3F20.10)') "h10002", "4Qy",       real(complex_taylor_coef(H_lie_poly_bas(1), hix)), aimag(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_bas(1), hix))
! write(20,*) "Special"
! hix = [3, 0, 1, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h30100", "3Qx-Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 0, 1, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h10100", "Qx+Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [6, 0, 2, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3ES20.10)') "h60200", "6Qx+2Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
close(20)

! open(20,FILE="ptc_h_pdep.out")
! WRITE(20,'(A10,4A20)') "# Term", "Resonance", "On Momentum", "Off Momentum (-)", "Off Momentum (+)"
! hix = [1, 1, 0, 0, 0, 1]
! WRITE(20,'(A10,A20,3F20.10)') "h11001", "--",        abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [0, 0, 1, 1, 0, 1]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h00111", "--",        abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! !3rd order terms                                                                                                                                                                                         
! hix = [2, 1, 0, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h21000", "Qx",        abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [3, 0, 0, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h30000", "3Qx",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 0, 1, 1, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h10110", "Qx",        abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 0, 0, 2, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h10020", "Qx-2Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 0, 2, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h10200", "Qx+2Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [2, 0, 0, 0, 0, 1]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h20001", "2Qx",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [0, 0, 2, 0, 0, 1]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h00201", "2Qy",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 0, 0, 0, 0, 2]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h10002", "Qx",        abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! write(20,*) "Special"
! hix = [3, 0, 1, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h30100", "3Qx-Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 0, 1, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h10100", "Qx+Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! !4th order terms                                                                                                                                                                                         
! write(20,*)                                                                                                                                                                                              
! write(20,*)                                                                                                                                                                                              
! hix = [2, 2, 0, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h22000", "--     ",   abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [0, 0, 2, 2, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h00220", "--     ",   abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 1, 1, 1, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h11110", "--     ",   abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [3, 1, 0, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h31000", "2Qx",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [4, 0, 0, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h40000", "4Qx",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [2, 0, 1, 1, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h20110", "2Qx",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 1, 2, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h11200", "2Qy",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [2, 0, 0, 2, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h20020", "2Qx-2Qy",   abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [2, 0, 2, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h20200", "2Qx+2Qy",   abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [0, 0, 3, 1, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h00310", "2Qy",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [0, 0, 4, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h10002", "4Qy",       abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! write(20,*) "Special"
! hix = [3, 0, 1, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h30100", "3Qx-Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! hix = [1, 0, 1, 0, 0, 0]                                                                                                                                                                                 
! WRITE(20,'(A10,A20,3F20.10)') "h10100", "Qx+Qy",    abs(complex_taylor_coef(H_lie_poly_bas(1), hix)), abs(complex_taylor_coef(H_lie_poly_neg(1), hix)), abs(complex_taylor_coef(H_lie_poly_pos(1), hix))
! close(20)

END PROGRAM ptc_h












