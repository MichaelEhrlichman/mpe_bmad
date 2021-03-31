PROGRAM test_h

USE bmad
USE ptc_layout_mod
USE madx_ptc_module

IMPLICIT none

CHARACTER lat_file*200
CHARACTER(3) order_str

TYPE(lat_struct) lat
TYPE(taylor_struct) one_turn_map(6)
TYPE(taylor_struct) A(6), A_inverse(6), dhdj(6)
TYPE(taylor_struct) :: h, h6(6)
TYPE(complex_taylor_struct) :: F(6), L(6)
TYPE(complex_taylor_struct) :: H_lie_poly(6)
TYPE(coord_struct), ALLOCATABLE :: co(:)
TYPE(normal_modes_struct) mode

INTEGER i,j
INTEGER order
INTEGER hix(6)

!parameters
REAL(rp) pz

REAL(rp) chromx
REAL(rp) chromy
COMPLEX(rp) h20001
COMPLEX(rp) h00201
COMPLEX(rp) h21000
COMPLEX(rp) h30000
COMPLEX(rp) h10110
COMPLEX(rp) h10020
COMPLEX(rp) h10200
COMPLEX(rp) h22000
COMPLEX(rp) h31000
COMPLEX(rp) h21001
COMPLEX(rp) h40000

LOGICAL rf_on

CALL GETARG(1, lat_file)
CALL GETARG(2, order_str)

READ(order_str,'(I3)') order

IF( (order .lt. 2) .or. (order .gt. 8) ) THEN
  order = 4
ENDIF
WRITE(*,*) "Calculating map to order: ", order

WRITE(*,*) "Preparing lattice..."

CALL bmad_parser(lat_file, lat)
CALL set_on_off(rfcavity$,lat,off$)
CALL closed_orbit_calc(lat,co,4)
CALL lat_make_mat6(lat,-1,co)
CALL twiss_at_start(lat)
CALL twiss_propagate_all(lat)
call chrom_calc(lat,1.0d-4,chromx,chromy)
write(*,*) "chrom_calc ", chromx, chromy

chromx = 0
chromy = 0
h20001 = 0
h00201 = 0
h21000 = 0
h30000 = 0
h10110 = 0
h10020 = 0
h10200 = 0
do i=1,lat%n_ele_track
  if( any(lat%ele(i)%key .eq. (/sbend$,quadrupole$,sextupole$,octupole$,rbend$,multipole$/))) then
    chromx   = chromx   + lat%ele(i)%value(k1$)*lat%ele(i)%value(l$)*lat%ele(i)%a%beta
    chromy   = chromy   + lat%ele(i)%value(k1$)*lat%ele(i)%value(l$)*lat%ele(i)%b%beta
    h20001   = h20001   + lat%ele(i)%value(k1$)*lat%ele(i)%value(l$)*lat%ele(i)%a%beta * exp( cmplx(0.0d0,2.0d0*lat%ele(i)%a%phi) )
    h00201   = h00201   + lat%ele(i)%value(k1$)*lat%ele(i)%value(l$)*lat%ele(i)%b%beta * exp( cmplx(0.0d0,2.0d0*lat%ele(i)%b%phi) )
    h21000   = h21000   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(i)%a%beta**(3./2.) * exp( cmplx(0.0d0,1.0d0*lat%ele(i)%a%phi) )
    h30000   = h30000   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(i)%a%beta**(3./2.) * exp( cmplx(0.0d0,3.0d0*lat%ele(i)%a%phi) )
    h10110   = h10110   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(i)%a%beta**(1./2.)*lat%ele(i)%b%beta * exp( cmplx(0.0d0,1.0d0*lat%ele(i)%a%phi) )
    h10020   = h10020   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(i)%a%beta**(1./2.)*lat%ele(i)%b%beta * exp( cmplx(0.0d0,lat%ele(i)%a%phi-2.0*lat%ele(i)%b%phi) )
    h10200   = h10200   + lat%ele(i)%value(k2$)*lat%ele(i)%value(l$)*lat%ele(i)%a%beta**(1./2.)*lat%ele(i)%b%beta * exp( cmplx(0.0d0,lat%ele(i)%a%phi+2.0*lat%ele(i)%b%phi) )
  endif
enddo

chromx = -1.0 * chromx / 4.0 / 3.1415926
chromy = chromy / 4.0 / 3.1415926
h20001 = h20001 / 8.0
h00201 = -1.0 * h00201 / 8.0
h21000 = -1.0 * h21000 / 8.0
h30000 = -1.0 * h30000 / 24.0
h10110 =  h10110 / 4.0
h10020 =  h10020 / 8.0
h10200 =  h10200 / 8.0
WRITE(*,*) "----------------------------"
WRITE(*,*) "Calculations from Analytic Formulas"
WRITE(*,'(A10,4A20)') "Term", "Resonance", "Re", "Im", "Abs"
WRITE(*,*) "----------------------------"
WRITE(*,'(A10,3A20,F20.10)') "chromx", "", "", "", chromx
WRITE(*,'(A10,3A20,F20.10)') "chromy", "", "", "", chromy
WRITE(*,'(A10,A20,3F20.10)') "h21000", "Qx",     h21000, ABS(h21000)
WRITE(*,'(A10,A20,3F20.10)') "h30000", "3Qx",    h30000, ABS(h30000)
WRITE(*,'(A10,A20,3F20.10)') "h10110", "Qx",     h10110, ABS(h10110)
WRITE(*,'(A10,A20,3F20.10)') "h10020", "Qx-2Qy", h10020, ABS(h10020)
WRITE(*,'(A10,A20,3F20.10)') "h10200", "Qx+2Qy", h10200, ABS(h10200)
WRITE(*,'(A10,A20,3F20.10)') "h20001", "2Qx",    h20001, ABS(h20001)
WRITE(*,'(A10,A20,3F20.10)') "h00201", "2Qy",    h00201, ABS(h00201)

bmad_com%radiation_damping_on=.false.
CALL lat_to_ptc_layout(lat)

pz = 0.0
rf_on = .false.
CALL ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, rf_on, 0.0d0, order)
CALL normal_form_taylors(one_turn_map, rf_on, dhdj)
CALL normal_form_complex_taylors(one_turn_map, rf_on, H_lie_poly=H_lie_poly)

! CALL ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, rf_on, 0.03d0, order)
! CALL normal_form_taylors(one_turn_map, rf_on, dhdj)
! CALL normal_form_complex_taylors(one_turn_map, rf_on, H_lie_poly=H_lie_poly)
! 
! CALL ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, rf_on, -0.03d0, order)
! CALL normal_form_taylors(one_turn_map, rf_on, dhdj)
! CALL normal_form_complex_taylors(one_turn_map, rf_on, H_lie_poly=H_lie_poly)

WRITE(*,*) "----------------------------"
WRITE(*,*) "Calculations from PTC"
WRITE(*,'(A10,4A20)') "Term", "Resonance", "Real Part", "Imaginary Part", "Absolute Value"
WRITE(*,*) "----------------------------"

hix = [0, 0, 0, 0, 0, 1]
WRITE(*,'(A10,3A20,F20.10)') "chromx:", "",        "", "", taylor_coef(dhdj(1), hix)
hix = [0, 0, 0, 0, 0, 1]
WRITE(*,'(A10,3A20,F20.10)') "chromy:", "",        "", "", taylor_coef(dhdj(2), hix)
hix = [2, 1, 0, 0, 0, 0]
WRITE(*,'(A10,A20,3F20.10)') "h21000:", "Qx",      real(complex_taylor_coef(H_lie_poly(1), hix)), aimag(complex_taylor_coef(H_lie_poly(1), hix)), abs(complex_taylor_coef(H_lie_poly(1), hix))
hix = [3, 0, 0, 0, 0, 0]
WRITE(*,'(A10,A20,3F20.10)') "h30000:", "3Qx",     real(complex_taylor_coef(H_lie_poly(1), hix)), aimag(complex_taylor_coef(H_lie_poly(1), hix)), abs(complex_taylor_coef(H_lie_poly(1), hix))
hix = [1, 0, 1, 1, 0, 0]
WRITE(*,'(A10,A20,3F20.10)') "h10110:", "Qx",      real(complex_taylor_coef(H_lie_poly(1), hix)), aimag(complex_taylor_coef(H_lie_poly(1), hix)), abs(complex_taylor_coef(H_lie_poly(1), hix))
hix = [1, 0, 0, 2, 0, 0]
WRITE(*,'(A10,A20,3F20.10)') "h10020:", "Qx-2Qy",  real(complex_taylor_coef(H_lie_poly(1), hix)), aimag(complex_taylor_coef(H_lie_poly(1), hix)), abs(complex_taylor_coef(H_lie_poly(1), hix))
hix = [1, 0, 2, 0, 0, 0]
WRITE(*,'(A10,A20,3F20.10)') "h10200:", "Qx+2Qy",  real(complex_taylor_coef(H_lie_poly(1), hix)), aimag(complex_taylor_coef(H_lie_poly(1), hix)), abs(complex_taylor_coef(H_lie_poly(1), hix))
hix = [2, 0, 0, 0, 0, 1]
WRITE(*,'(A10,A20,3F20.10)') "h20001:", "2Qx",     real(complex_taylor_coef(H_lie_poly(1), hix)), aimag(complex_taylor_coef(H_lie_poly(1), hix)), abs(complex_taylor_coef(H_lie_poly(1), hix))
hix = [0, 0, 2, 0, 0, 1]
WRITE(*,'(A10,A20,3F20.10)') "h00201:", "2Qy",     real(complex_taylor_coef(H_lie_poly(1), hix)), aimag(complex_taylor_coef(H_lie_poly(1), hix)), abs(complex_taylor_coef(H_lie_poly(1), hix))
WRITE(*,*)
hix = [2, 2, 0, 0, 0, 0]
write(*,'(a10,a20,3f20.10)') "h22000:", "   ",     real(complex_taylor_coef(h_lie_poly(1), hix)), aimag(complex_taylor_coef(h_lie_poly(1), hix)), abs(complex_taylor_coef(h_lie_poly(1), hix))
hix = [3, 1, 0, 0, 0, 0]
write(*,'(a10,a20,3f20.10)') "h31000:", "2Qx",     real(complex_taylor_coef(h_lie_poly(1), hix)), aimag(complex_taylor_coef(h_lie_poly(1), hix)), abs(complex_taylor_coef(h_lie_poly(1), hix))
hix = [4, 0, 0, 0, 0, 0]
write(*,'(a10,a20,3f20.10)') "h40000:", "   ",     real(complex_taylor_coef(h_lie_poly(1), hix)), aimag(complex_taylor_coef(h_lie_poly(1), hix)), abs(complex_taylor_coef(h_lie_poly(1), hix))
hix = [2, 1, 0, 0, 0, 1]
write(*,'(a10,a20,3f20.10)') "h21001:", "   ",     real(complex_taylor_coef(h_lie_poly(1), hix)), aimag(complex_taylor_coef(h_lie_poly(1), hix)), abs(complex_taylor_coef(h_lie_poly(1), hix))

WRITE(*,*)


END PROGRAM test_h












