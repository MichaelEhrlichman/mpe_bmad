program dynap_slim

  use bmad
  use sls_lib
  use measurement_mod

  implicit none

  type (lat_struct) ring_ideal
  type (lat_struct) ring_ma
  type (coord_struct), allocatable :: co_ideal(:)
  type (coord_struct), allocatable :: co_ma(:)
  type (coord_struct) co_ma_at_s
  type(ele_struct) ele_ideal
  type(ele_struct) ele_ma
  type(bpm_phase_coupling_struct) k_and_cbar
  type(bpm_phase_coupling_struct) k_and_cbar_rms

! type bpm_phase_coupling_struct
!   real(rp) K_22a  ! In-phase y/x for a-mode oscillations.
!   real(rp) K_12a  ! Out-of-phase y/x for a-mode oscillations.
!   real(rp) K_11b  ! In-phase x/y for b-mode oscillations.
!   real(rp) K_12b  ! Out-of-phase x/y for b-mode oscillations.
!   real(rp) Cbar22_a ! Cbar22 as calculated from K_22a.
!   real(rp) Cbar12_a ! Cbar12 as calculated from K_12a.
!   real(rp) Cbar11_b ! Cbar11 as calculated from K_11b.
!   real(rp) Cbar12_b ! Cbar12 as calculated from K_12b.
!   real(rp) phi_a    ! a-mode betatron phase.
!   real(rp) phi_b    ! b-mode betatron phase.
! end type


  integer i,j,n
  integer status

  real(rp) s, mean_percent_x, mean_percent_y
  real(rp), parameter :: delta_s = 0.01
  real(rp) cbar_ma(2,2)
  real(rp) cbar_rms(2,2)
  real(rp) gamma_c_rms
  real(rp) x_rms, y_rms
  real(rp) x_meas, y_meas

  logical err

  character(11) mode

  character*200 lat_1, lat_2

  call getarg(1,lat_1)
  call getarg(2,lat_2)
  call getarg(3,mode)
    
  call bmad_parser(lat_1, ring_ideal)
  call bmad_parser(lat_2, ring_ma)
  
  call twiss_and_track(ring_ideal,co_ideal,status)
  call twiss_and_track(ring_ma,co_ma,status)
  ! call calc_ring(ring_ideal,4,co_ideal,err)
  ! call calc_ring(ring_ma,4,co_ma,err)

  s = 0.0
  n = 0
  mean_percent_x = 0.0d0
  mean_percent_y = 0.0d0
  cbar_rms = 0.0d0
  k_and_cbar_rms%K_22a = 0.0d0
  k_and_cbar_rms%K_12a = 0.0d0
  k_and_cbar_rms%K_11b = 0.0d0
  k_and_cbar_rms%K_12b = 0.0d0
  k_and_cbar_rms%Cbar22_a = 0.0d0
  k_and_cbar_rms%Cbar12_a = 0.0d0
  k_and_cbar_rms%Cbar11_b = 0.0d0
  k_and_cbar_rms%Cbar12_b = 0.0d0
  x_rms = 0.0d0
  y_rms = 0.0d0
  gamma_c_rms = 0.0d0
  if(mode == 'instruments') then
    write(*,*) "Instruments-only mode"
    do i=1, ring_ideal%n_ele_track
      if( ring_ideal%ele(i)%key == instrument$ ) then
        ele_ideal = ring_ideal%ele(i)
        ele_ma = ring_ma%ele(i)
        co_ma_at_s = co_ma(i)
        !call twiss_and_track_at_s(ring_ideal,s,ele_ideal,co_ideal)
        !call twiss_and_track_at_s(ring_ma,s,ele_ma,co_ma,co_ma_at_s)
        call c_to_cbar(ele_ma, cbar_ma)
        cbar_rms = cbar_rms + cbar_ma**2
        call to_orbit_reading(co_ma_at_s, ele_ma, x_plane$, .false., x_meas, err)
        call to_orbit_reading(co_ma_at_s, ele_ma, y_plane$, .false., y_meas, err)
        call to_phase_and_coupling_reading (ele_ma, .false., k_and_cbar, err)
        k_and_cbar_rms%Cbar22_a = k_and_cbar_rms%Cbar22_a + k_and_cbar%Cbar22_a**2
        k_and_cbar_rms%Cbar12_a = k_and_cbar_rms%Cbar12_a + k_and_cbar%Cbar12_a**2
        k_and_cbar_rms%Cbar11_b = k_and_cbar_rms%Cbar11_b + k_and_cbar%Cbar11_b**2
        k_and_cbar_rms%Cbar12_b = k_and_cbar_rms%Cbar12_b + k_and_cbar%Cbar12_b**2
        x_rms = x_rms + x_meas**2
        y_rms = y_rms + y_meas**2
        gamma_c_rms = gamma_c_rms + ele_ma%gamma_c**2
        mean_percent_x = mean_percent_x + abs(ele_ideal%a%beta-ele_ma%a%beta)/ele_ideal%a%beta
        mean_percent_y = mean_percent_y + abs(ele_ideal%b%beta-ele_ma%b%beta)/ele_ideal%b%beta
        s = s + delta_s
        n = n + 1
      endif
    enddo
  else
    do while (s .lt. ring_ideal%param%total_length)
      call twiss_and_track_at_s(ring_ideal,s,ele_ideal,co_ideal)
      call twiss_and_track_at_s(ring_ma,s,ele_ma,co_ma,co_ma_at_s)
      call c_to_cbar(ele_ma, cbar_ma)
      cbar_rms = cbar_rms + cbar_ma**2
      x_rms = x_rms + co_ma_at_s%vec(1)**2
      y_rms = y_rms + co_ma_at_s%vec(3)**2
      gamma_c_rms = gamma_c_rms + ele_ma%gamma_c**2
      mean_percent_x = mean_percent_x + abs(ele_ideal%a%beta-ele_ma%a%beta)/ele_ideal%a%beta
      mean_percent_y = mean_percent_y + abs(ele_ideal%b%beta-ele_ma%b%beta)/ele_ideal%b%beta
      s = s + delta_s
      n = n + 1
    enddo
  endif
  mean_percent_x = mean_percent_x / n
  mean_percent_y = mean_percent_y / n
  cbar_rms = sqrt(cbar_rms/n)
  x_rms = sqrt(x_rms/n)
  y_rms = sqrt(y_rms/n)
  gamma_c_rms = sqrt(gamma_c_rms/n)

  write(*,'(a,f11.3)') "Mean x beta beat %: ", mean_percent_x * 100.
  write(*,'(a,f11.3)') "Mean y beta beat %: ", mean_percent_y * 100.
  write(*,'(a,f9.3,a)') "Second lattice RMS of the cbar matrix (rms gamma_c = ", gamma_c_rms, "):"
  write(*,'(2f11.5)') cbar_rms(1,1), cbar_rms(1,2)
  write(*,'(2f11.5)') cbar_rms(2,1), cbar_rms(2,2)
  write(*,'(a,f9.3,a)') "Second lattice RMS of the BPM cbar matrix:"
  write(*,'(2f11.5)') sqrt(k_and_cbar_rms%Cbar11_b/n), sqrt(k_and_cbar_rms%Cbar12_b/n)
  write(*,'(2f11.5)') sqrt(k_and_cbar_rms%Cbar12_a/n), sqrt(k_and_cbar_rms%Cbar22_a/n)
  write(*,'(a,es14.4)') "x_rms: ", x_rms
  write(*,'(a,es14.4)') "y_rms: ", y_rms
  open(45,file='avg_beta_beats.dat')
  write(45,*) mean_percent_x*100.
  write(45,*) mean_percent_y*100.
  close(45)

end program
                                                            










