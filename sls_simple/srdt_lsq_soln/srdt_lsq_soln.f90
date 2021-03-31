program srdt_lsq_soln

use bmad
use srdt_mod

implicit none

character in_file*200, lat_file*200

type(lat_struct) lat
type(coord_struct), allocatable :: co(:)
type(summation_rdt_struct) srdt
type(summation_rdt_struct), allocatable :: per_ele_rdt(:)

real(rp), allocatable :: ls_soln(:)
real(rp) chrom_x, chrom_y

integer i, j, k, nVar
integer status
integer gen_slices, sxt_slices
integer, allocatable :: var_indexes(:)

character(40) var_names(200)

namelist /srdt_lsq/ lat_file, var_names, gen_slices, sxt_slices, chrom_x, chrom_y

gen_slices = 60
sxt_slices = 120
chrom_x = 1.0
chrom_y = 1.0
var_names = ''

call getarg(1, in_file)
open (unit = 10, file = in_file, action='read')
read (10, nml = srdt_lsq)
close (10)

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)
call set_on_off(rfcavity$,lat,off$)
call twiss_and_track(lat, co, status)

do i=1,size(var_names)
  var_names(i) = upcase(var_names(i))
enddo
call name_to_list(lat,var_names)
nVar = count(lat%ele(:)%select)
allocate(var_indexes(nVar))
allocate(ls_soln(nVar))
var_indexes = pack([(i,i=0,lat%n_ele_track)],lat%ele(:)%select)

call srdt_lsq_solution(lat, var_indexes, ls_soln, gen_slices, sxt_slices, chrom_x, chrom_y)

open(45,file='srdt_lsq_soln.dat')
do i=1,nVar
  write(45,'(a,a,es16.8)') trim(lat%ele(var_indexes(i))%name), '[k2]=', ls_soln(i)
enddo
close(45)

end program


