program analytic_h

use bmad
use srdt_mod

implicit none

character lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: co(:)
type(summation_rdt_struct) srdt
type(summation_rdt_struct), allocatable :: per_ele_rdt(:)

integer i
integer status
integer n_gen_slice, n_sxt_slice

complex(rp), allocatable :: cache(:,:,:)

call getarg(1, lat_file)

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)
call set_on_off(rfcavity$,lat,off$)
call twiss_and_track(lat, co, status)


n_gen_slice = 60
n_sxt_slice = 120

call srdt_calc(lat, srdt, 2, n_gen_slice, n_sxt_slice, per_ele_rdt)

!call make_srdt_cache(lat, n_gen_slice, n_sxt_slice, cache)
!call srdt_calc_with_cache(lat, srdt, 2, n_gen_slice, n_sxt_slice, cache)

open(20, file='srdt_by_ele.out')
write(20,'(a14,a,a40,16a14)') 's', '   ', 'name', 'h20001.r', 'h20001.i', 'h00201.r', 'h00201.i', 'h10002.r', 'h10002.i', &
                                   'h21000.r', 'h21000.i', 'h30000.r', 'h30000.i', 'h10110.r', 'h10110.i', &
                                   'h10020.r', 'h10020.i', 'h10200.r', 'h10200.i'
do i=1,size(per_ele_rdt)
  write(20,'(i14,a,a,16es14.4)') i, "   ", lat%ele(i)%name, per_ele_rdt(i)%h20001, per_ele_rdt(i)%h00201, per_ele_rdt(i)%h10002, per_ele_rdt(i)%h21000, &
                                                            per_ele_rdt(i)%h30000, per_ele_rdt(i)%h10110, per_ele_rdt(i)%h10020, per_ele_rdt(i)%h10200
enddo
close(20)

open(20, file='srdt_sum_by_ele.out')
write(20,'(a6,a12,a,a20,16a14)') 'ix', 's', '   ', 'name', 'h20001.r', 'h20001.i', 'h00201.r', 'h00201.i', 'h10002.r', 'h10002.i', &
                                   'h21000.r', 'h21000.i', 'h30000.r', 'h30000.i', 'h10110.r', 'h10110.i', &
                                   'h10020.r', 'h10020.i', 'h10200.r', 'h10200.i'
do i=1,size(per_ele_rdt)
  write(20,'(i6,f12.4,a,a20,16es14.4)') i, lat%ele(i)%s, "   ", adjustl(lat%ele(i)%name), sum(per_ele_rdt(1:i)%h20001), sum(per_ele_rdt(1:i)%h00201), sum(per_ele_rdt(1:i)%h10002), &
                                                            sum(per_ele_rdt(1:i)%h21000), sum(per_ele_rdt(1:i)%h30000), sum(per_ele_rdt(1:i)%h10110), &
                                                            sum(per_ele_rdt(1:i)%h10020), sum(per_ele_rdt(1:i)%h10200)
enddo
close(20)

write(*,*) "----------------------------"
write(*,*) "Calculations from Analytic Formulas"
write(*,'(A10,4A20)') "Term", "Resonance", "Re", "Im", "Abs"
write(*,*) "----------------------------"
write(*,'(A10,A20,3F20.10)') "h11001", "chi_x",  srdt%h11001, abs(srdt%h11001)
write(*,'(A10,A20,3F20.10)') "h00111", "chi_y",  srdt%h00111, abs(srdt%h00111)
write(*,*)
write(*,'(A10,A20,3F20.10)') "h21000", "Qx",     srdt%h21000, abs(srdt%h21000)
write(*,'(A10,A20,3F20.10)') "h30000", "3Qx",    srdt%h30000, abs(srdt%h30000)
write(*,'(A10,A20,3F20.10)') "h10110", "Qx",     srdt%h10110, abs(srdt%h10110)
write(*,'(A10,A20,3F20.10)') "h10020", "Qx-2Qy", srdt%h10020, abs(srdt%h10020)
write(*,'(A10,A20,3F20.10)') "h10200", "Qx+2Qy", srdt%h10200, abs(srdt%h10200)
write(*,'(A10,A20,3F20.10)') "h20001", "2Qx",    srdt%h20001, abs(srdt%h20001)
write(*,'(A10,A20,3F20.10)') "h00201", "2Qy",    srdt%h00201, abs(srdt%h00201)
write(*,'(A10,A20,3F20.10)') "h10002", "eta_a'", srdt%h10002, abs(srdt%h10002)
write(*,*)
write(*,'(A10,A20,3F20.10)') "h31000", "2Qx",     srdt%h31000, abs(srdt%h31000)
write(*,'(A10,A20,3F20.10)') "h00310", "2Qy",     srdt%h00310, abs(srdt%h00310)
write(*,'(A10,A20,3F20.10)') "h11200", "2Qy",     srdt%h11200, abs(srdt%h11200)
write(*,'(A10,A20,3F20.10)') "h20110", "2Qx",     srdt%h20110, abs(srdt%h20110)
write(*,'(A10,A20,3F20.10)') "h20200", "2Qy-2Qy", srdt%h20200, abs(srdt%h20200)
write(*,'(A10,A20,3F20.10)') "h20020", "2Qy+2Qy", srdt%h20020, abs(srdt%h20020)
write(*,'(A10,A20,3F20.10)') "h40000", "4Qx",     srdt%h40000, abs(srdt%h40000)
write(*,'(A10,A20,3F20.10)') "h00400", "4Qy",     srdt%h00400, abs(srdt%h00400)
write(*,*)
write(*,'(A20,A10,3F20.10)') "h22000", "   ",    srdt%h22000
write(*,'(A20,A10,3F20.10)') "h00220", "   ",    srdt%h00220
write(*,'(A20,A10,3F20.10)') "h11110", "   ",    srdt%h11110
!write(*,'(A20,A10,3F20.10)') "nux_Jx from h22000", "   ",    srdt%h22000 * -2.0/pi
!write(*,'(A20,A10,3F20.10)') "nux_Jy from h11110", "   ",    srdt%h11110 * -1.0/pi
!write(*,'(A20,A10,3F20.10)') "nuy_Jy from h00220", "   ",    srdt%h00220 * -2.0/pi

end program












