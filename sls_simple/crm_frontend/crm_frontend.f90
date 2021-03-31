program crm_frontend

use bmad
use crm_mod

implicit none

character lat_file*200
character mag_file*200
type(lat_struct) lat
type(coord_struct), allocatable :: co(:)
type(crm_struct) crm
type(mag_struct), pointer :: c_mags(:)
integer i, status
integer j, n_chrom
logical err_flag
real(rp) crmat(2)

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,co,status)

call getarg(2, mag_file)
open(20,file='chrom_mags.list', action='read', status='old')
n_chrom = 0
do
  read(20,*,end=10)
  n_chrom = n_chrom + 1
enddo
10 close(20)
write(*,*) "Number of chrom mags in file: ", n_chrom
call crm_alloc(crm,n_chrom)
allocate(c_mags(n_chrom))
crm%c_mags => c_mags
crm%set_chrom_x = 1.0d0
crm%set_chrom_y = 1.0d0

open(20,file='chrom_mags.list')
do i=1,n_chrom
  read(20,*) crm%c_mags(i)%name, crm%c_mags(i)%property
  crm%c_mags(i)%type = 'c'
enddo
close(20)

call crm_build(lat, crm, err_flag)
if(err_flag) then
  write(*,*) "could not build chromaticity matrices at program start.  aborting."
  stop
endif

open(45,file='crm.dat')
do i=1,n_chrom
  write(45,'(a,f16.7)') trim(crm%c_mags(i)%name)//'[k2]=',crm%ApC(i)
enddo
close(45)

end program












