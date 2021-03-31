program footprint

use bmad
use sls_lib

use namelist_general !general: lat_file, use_hybrid
use namelist_fp !fp: pz_min, pz_max, n_pz

implicit none

character*60 in_file
character*100 lat_file_override

type(lat_struct) lat
type(coord_struct), allocatable :: co(:)

integer iargs
integer i, i0
integer j
integer status
integer pz_min_int, pz_max_int

real(rp) dpz
real(rp), allocatable :: pz(:), nu_a(:), nu_b(:), tr_a(:), tr_b(:)
real(rp) vec0(6)

logical err

iargs = iargc()
lat_file_override = ''
if( iargs == 1 ) then
  call getarg(1,in_file)
elseif (iargs == 2 ) then
  call getarg(1,in_file)
  call getarg(2,lat_file_override)
endif

open(10, file = in_file, action='read')
read(10, nml = general)
read(10, nml = fp)
close(10)

if( lat_file_override .ne. '' ) then
  lat_file = lat_file_override
endif

!pz_min_int = floor(100*pz_min)
!pz_max_int = floor(100*pz_max)
pz_min_int = 100*pz_min
pz_max_int = 100*pz_max

write(*,*) "Preparing lattice..."

call bmad_parser(lat_file, lat)
allocate(co(0:lat%n_ele_track))

call set_on_off(rfcavity$,lat,off$)
bmad_com%radiation_damping_on=.false.

allocate(nu_a(n_pz))
allocate(nu_b(n_pz))
allocate(tr_a(n_pz))
allocate(tr_b(n_pz))
allocate(pz(n_pz))

dpz = (pz_max-pz_min)/(n_pz-1)
i0 = (n_pz+1)/2
do i=i0,1,-1
  pz(i) = pz_min + dpz*(i-1)
  vec0 = 0.0d0
  vec0(6) = pz(i)
  call init_coord(co(0), vec0, lat%ele(0), upstream_end$)
  !call calc_ring(lat,4,co,err)
  call twiss_and_track(lat,co,status)
  if(status .ne. ok$) then
    nu_a(i) = -1.0
    nu_b(i) = -1.0
  else
    nu_a(i) = lat%ele(lat%n_ele_track)%a%phi/twopi
    nu_b(i) = lat%ele(lat%n_ele_track)%b%phi/twopi
  endif
  tr_a(i) = lat%param%t1_no_RF(1,1)+lat%param%t1_no_RF(2,2)
  tr_b(i) = lat%param%t1_no_RF(3,3)+lat%param%t1_no_RF(4,4)
enddo
call clear_lat_1turn_mats(lat)
do i=i0,n_pz
  pz(i) = pz_min + dpz*(i-1)
  vec0 = 0.0d0
  vec0(6) = pz(i)
  call init_coord(co(0), vec0, lat%ele(0), upstream_end$)
  !call calc_ring(lat,4,co,err)
  call twiss_and_track(lat,co,status)
  if(status .ne. ok$) then
    nu_a(i) = -1.0
    nu_b(i) = -1.0
  else
    nu_a(i) = lat%ele(lat%n_ele_track)%a%phi/twopi
    nu_b(i) = lat%ele(lat%n_ele_track)%b%phi/twopi
  endif
  tr_a(i) = lat%param%t1_no_RF(1,1)+lat%param%t1_no_RF(2,2)
  tr_b(i) = lat%param%t1_no_RF(3,3)+lat%param%t1_no_RF(4,4)
enddo

open(21,file='footprint.dat')
write(21,'(a1,a10,2a14)') "#", "pz", "nu_x", "nu_y"
do i=1,n_pz
  if( nu_a(i).lt.0 .or. nu_b(i).lt.0 ) then
    write(21,'(f18.11,2a14,2es14.4)') pz(i), "*", "*", tr_a(i), tr_b(i)
  else
    write(21,'(f18.11,2f14.6,2es14.4)') pz(i), nu_a(i), nu_b(i), tr_a(i), tr_b(i)
  endif
enddo
close(21)

deallocate(nu_a)
deallocate(nu_b)
deallocate(tr_a)
deallocate(tr_b)
deallocate(pz)
n_pz = pz_max_int-pz_min_int+1
i0 = (n_pz+1)/2
allocate(nu_a(n_pz))
allocate(nu_b(n_pz))
allocate(tr_a(n_pz))
allocate(tr_b(n_pz))
allocate(pz(n_pz))

do i=i0,1,-1
  pz(i) = 0.01*(pz_min_int + (i-1))
  vec0 = 0.0d0
  vec0(6) = pz(i)
  call init_coord(co(0), vec0, lat%ele(0), upstream_end$)
  !call calc_ring(lat,4,co,err)
  call twiss_and_track(lat,co,status)
  if(status .ne. ok$) then
    nu_a(i) = -1.0
    nu_b(i) = -1.0
  else
    nu_a(i) = lat%ele(lat%n_ele_track)%a%phi/twopi
    nu_b(i) = lat%ele(lat%n_ele_track)%b%phi/twopi
  endif
  tr_a(i) = lat%param%t1_no_RF(1,1)+lat%param%t1_no_RF(2,2)
  tr_b(i) = lat%param%t1_no_RF(3,3)+lat%param%t1_no_RF(4,4)
enddo
call clear_lat_1turn_mats(lat)
do i=i0,n_pz
  pz(i) = 0.01*(pz_min_int + (i-1))
  vec0 = 0.0d0
  vec0(6) = pz(i)
  call init_coord(co(0), vec0, lat%ele(0), upstream_end$)
  !call calc_ring(lat,4,co,err)
  call twiss_and_track(lat,co,status)
  if(status .ne. ok$) then
    nu_a(i) = -1.0
    nu_b(i) = -1.0
  else
    nu_a(i) = lat%ele(lat%n_ele_track)%a%phi/twopi
    nu_b(i) = lat%ele(lat%n_ele_track)%b%phi/twopi
  endif
  tr_a(i) = lat%param%t1_no_RF(1,1)+lat%param%t1_no_RF(2,2)
  tr_b(i) = lat%param%t1_no_RF(3,3)+lat%param%t1_no_RF(4,4)
enddo

open(21,file='footprint_1pct.dat')
write(21,'(a1,a10,2a14)') "#", "pz", "nu_x", "nu_y"
do i=1,n_pz
  if( nu_a(i).lt.0 .or. nu_b(i).lt.0 ) then
    write(21,'(f18.11,2a14,2f14.6)') pz(i), "*", "*", tr_a(i), tr_b(i)
  else
    write(21,'(f18.11,4f14.6)') pz(i), nu_a(i), nu_b(i), tr_a(i), tr_b(i)
  endif
enddo
close(21)

deallocate(nu_a)
deallocate(nu_b)
deallocate(tr_a)
deallocate(tr_b)

end program












