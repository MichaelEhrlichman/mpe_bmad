PROGRAM ptc_ndts

use bmad
use ptc_layout_mod
use madx_ptc_module
!use dynap_pisa_mod

implicit none

character*60 in_file
character lat_file*200

type(lat_struct) lat
type(taylor_struct) one_turn_map(6)
type(complex_taylor_struct) :: H_lie_poly(6)

integer i
real(rp) pz_array(3)

!parameters
real(rp) pz
logical rf_on

!ndts
real(rp) h21000, h30000, h10110, h10020, h10200, h20001, h00201, h10002
real(rp) h31000, h40000, h20110, h11200, h20020, h20200, h00310, h00400
real(rp) ndt3sum, ndt4sum

namelist / general /    lat_file

call getarg(1, in_file)
open(10, file = in_file)
read(10, nml = general)
close(10)

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)

call set_on_off(rfcavity$,lat,off$)
rf_on = .false.
bmad_com%radiation_damping_on=.false.

call lat_to_ptc_layout(lat)

open(22,file='ptc_ndts.dat')
pz_array = (/0.0d0, -0.03d0, 0.03d0/)
do i=1,3
  pz = pz_array(i)
  call ptc_one_turn_map_at_ele(lat%ele(0), one_turn_map, rf_on,pz, order=5)
  call normal_form_complex_taylors(one_turn_map, rf_on, H_lie_poly=H_lie_poly)
  h21000 = abs(complex_taylor_coef(H_lie_poly(1), [2,1,0,0,0,0]))
  h30000 = abs(complex_taylor_coef(H_lie_poly(1), [3,0,0,0,0,0]))
  h10110 = abs(complex_taylor_coef(H_lie_poly(1), [1,0,1,1,0,0]))
  h10020 = abs(complex_taylor_coef(H_lie_poly(1), [1,0,0,2,0,0]))
  h10200 = abs(complex_taylor_coef(H_lie_poly(1), [1,0,2,0,0,0]))
  h20001 = abs(complex_taylor_coef(H_lie_poly(1), [2,0,0,0,0,1]))
  h00201 = abs(complex_taylor_coef(H_lie_poly(1), [0,0,2,0,0,1]))
  h10002 = abs(complex_taylor_coef(H_lie_poly(1), [1,0,0,0,0,2]))
  ndt3sum = h21000+h30000+h10110+h10020+h10200+h20001+h00201+h10002
  h31000 = abs(complex_taylor_coef(H_lie_poly(1), [3,1,0,0,0,0]))
  h40000 = abs(complex_taylor_coef(H_lie_poly(1), [4,0,0,0,0,0]))
  h20110 = abs(complex_taylor_coef(H_lie_poly(1), [2,0,1,1,0,0]))
  h11200 = abs(complex_taylor_coef(H_lie_poly(1), [1,1,2,0,0,0]))
  h20020 = abs(complex_taylor_coef(H_lie_poly(1), [2,0,0,2,0,0]))
  h20200 = abs(complex_taylor_coef(H_lie_poly(1), [2,0,2,0,0,0]))
  h00310 = abs(complex_taylor_coef(H_lie_poly(1), [0,0,3,1,0,0]))
  h00400 = abs(complex_taylor_coef(H_lie_poly(1), [0,0,4,0,0,0]))
  ndt4sum = h31000+h40000+h20110+h11200+h20020+h20200+h00310+h00400

  write(22,'(f12.4,a,es14.6)') pz, "  h21000:  ", h21000
  write(22,'(f12.4,a,es14.6)') pz, "  h30000:  ", h30000
  write(22,'(f12.4,a,es14.6)') pz, "  h10110:  ", h10110
  write(22,'(f12.4,a,es14.6)') pz, "  h10020:  ", h10020
  write(22,'(f12.4,a,es14.6)') pz, "  h10200:  ", h10200
  write(22,'(f12.4,a,es14.6)') pz, "  h20001:  ", h20001
  write(22,'(f12.4,a,es14.6)') pz, "  h00201:  ", h00201
  write(22,'(f12.4,a,es14.6)') pz, "  h10002:  ", h10002
  write(22,'(f12.4,a,es14.6)') pz, "  ndt3sum: ", ndt3sum
  write(22,'(f12.4,a,es14.6)') pz, "  h31000:  ", h31000
  write(22,'(f12.4,a,es14.6)') pz, "  h40000:  ", h40000
  write(22,'(f12.4,a,es14.6)') pz, "  h20110:  ", h20110
  write(22,'(f12.4,a,es14.6)') pz, "  h11200:  ", h11200
  write(22,'(f12.4,a,es14.6)') pz, "  h20020:  ", h20020
  write(22,'(f12.4,a,es14.6)') pz, "  h20200:  ", h20200
  write(22,'(f12.4,a,es14.6)') pz, "  h00310:  ", h00310
  write(22,'(f12.4,a,es14.6)') pz, "  h00400:  ", h00400
  write(22,'(f12.4,a,es14.6)') pz, "  ndt4sum: ", ndt4sum
  write(22,*)
  write(22,*)
enddo
close(22)

end program












