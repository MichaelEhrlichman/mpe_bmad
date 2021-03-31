program amp_fact2

use bmad
!use sls_lib

implicit none

character(200) lat_file
character(200) data_file
character(200) line

type(lat_struct) lat
type(coord_struct), allocatable :: orb(:)
type(ele_struct) ele_at_s
type(ele_struct), pointer :: ele_at_eval_p
type(ele_struct), target :: ele_at_eval

integer i,j
integer ikick
integer stat
integer ndata
integer eval_ix

real(rp) eval_s
real(rp) s, delta_s
real(rp) k1l_at_kick
real(rp) nu, dphi
real(rp) amp_factor, filter, PSD
real(rp), allocatable :: freq(:), beam_cas2(:)
real(rp), allocatable :: cas_a(:), cas_b(:), filtered_beam_psd(:), integrated_filtered_beam_cas(:)

bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

call getarg(1, lat_file)
call getarg(2, data_file)
call bmad_parser(lat_file, lat)
call twiss_and_track(lat,orb,stat)

!count data
ndata = -1 !account for header line
open(2020,file=data_file)
do while(.true.)
  read(2020,'(A)',iostat=stat) line
  if(stat < 0) exit
  ndata = ndata + 1
enddo
close(2020)
write(*,*) "Number of data: ", ndata
allocate(freq(ndata))
allocate(cas_a(ndata))
allocate(cas_b(ndata))
allocate(filtered_beam_psd(ndata-1))
allocate(integrated_filtered_beam_cas(ndata-1))
allocate(beam_cas2(ndata))
!read data
open(2020,file=data_file)
read(2020,*) line !skip header
do i=1,ndata
  read(2020,*) freq(i), cas_a(i), cas_b(i)
enddo
close(2020)

eval_s = 13.06  !0.4 m past end of SHF
call twiss_and_track_at_s(lat, eval_s, ele_at_eval)
ele_at_eval_p => ele_at_eval

! eval_ix = 6
! ele_at_eval_p => lat%ele(eval_ix)

nu = lat%ele(lat%n_ele_track)%a%phi / twopi
!FOOy nu = lat%ele(lat%n_ele_track)%b%phi / twopi

beam_cas2 = 0
do i=1,lat%n_ele_track
  if(has_attribute(lat%ele(i), 'K1')) then
    dphi = ele_at_eval_p%a%phi - lat%ele(i)%a%phi
    !y dphi = ele_at_eval_p%b%phi - lat%ele(i)%b%phi
    k1l_at_kick = lat%ele(i)%value(l$) * value_of_attribute(lat%ele(i), 'K1', err_print_flag = .true.)
    amp_factor = sqrt(lat%ele(i)%a%beta*ele_at_eval_p%a%beta)*cos(pi*nu-abs(dphi))/2.0d0/sin(pi*nu) * k1l_at_kick
    !y amp_factor = sqrt(lat%ele(i)%b%beta*ele_at_eval_p%b%beta)*cos(pi*nu-abs(dphi))/2.0d0/sin(pi*nu) * k1l_at_kick
    do j=1,ndata
      if((i .gt. 157) .and. (i .lt. 470)) then
        beam_cas2(j) = beam_cas2(j) + (amp_factor*cas_b(j))**2
      else
        beam_cas2(j) = beam_cas2(j) + (amp_factor*cas_a(j))**2
      endif
    enddo
  endif
enddo

open(45,file='unfiltered_CAS_at_ele.out')
write(45,'(a14,a14)') "# Freq (Hz)", "Summed CAS(nm)"
do i=1,ndata
  write(45,'(f14.6,2es14.5)') freq(i), sqrt(beam_cas2(i))
enddo
close(45)

open(45,file='filtered_CAS_at_ele.out')
write(45,'(a14,a14)') "# Freq (Hz)", "Filtered Summed CAS(nm)"
! Apply filter to PSD obtained by taking derivative of beam CAS
do i=1,ndata-1
  filter = PSD_attenuation_ratio((freq(i)+freq(i+1))/2.0d0)
  PSD = ( sqrt(beam_cas2(i)) - sqrt(beam_cas2(i+1)) ) / ( freq(i) - freq(i+1) )
  filtered_beam_psd(i) = -1.0d0*filter*PSD
enddo
! Backwards integrate PSD to get CAS
integrated_filtered_beam_cas(ndata-1) = filtered_beam_psd(ndata-1)*0.046875
do i=ndata-2,1,-1
  integrated_filtered_beam_cas(i) = integrated_filtered_beam_cas(i+1) + filtered_beam_psd(i)*0.046875
enddo
do i=1,ndata-1
  write(45,'(f14.6,es14.5)') freq(i), integrated_filtered_beam_cas(i)
enddo
close(45)

contains

  function PSD_attenuation_ratio(f) result(atten)
    real(rp) f
    real(rp) atten
    real(rp) :: fbw = 40.0d0 !200.0d0

    atten = (f/fbw)**2 / (1.0d0+(f/fbw)**2)
  end function

end program












