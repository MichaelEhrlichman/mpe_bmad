program rad_int

use bmad
!use sls_lib

implicit none

procedure(track1_custom_def) :: track1_custom

character lat_file*200

type(lat_struct) lat
type(normal_modes_struct) mode
type(coord_struct), allocatable :: co(:)
type(rad_int_all_ele_struct) rad_int_by_ele

integer ix_cache, i
integer status

track1_custom_ptr => track1_custom

!bmad_com%convert_to_kinetic_momentum = .true. !needed for LGB tracking
bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

WRITE(*,*) "Preparing lattice..."
call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

call twiss_and_track(lat,co,status)
ix_cache = -1
!call radiation_integrals(lat, co, mode, ix_cache)
call radiation_integrals (lat, co, mode, ix_cache, rad_int_by_ele=rad_int_by_ele)
call calc_z_tune(lat%branch(0))

call write_output(6)
open(60,file='rad_int.out')
call write_output(60)
close(60)

open(61,file='rad_int_s.out')
do i=1,lat%n_ele_track
  write(61,'(i6,f14.6,es14.4,a,a)') i, lat%ele(i)%s, rad_int_by_ele%branch(0)%ele(i)%i4a, "     ", lat%ele(i)%name
enddo
close(61)

contains

  subroutine write_output(lun)
    integer lun
    write(lun,*) "!   a-mode tune:          ", lat%ele(lat%n_ele_track)%a%phi /2./pi
    write(lun,*) "!   b-mode tune:          ", lat%ele(lat%n_ele_track)%b%phi /2./pi
    write(lun,*) "!   Synchrotron tune:     ", lat%z%tune /2./pi
    write(lun,*) 
    write(lun,*) "!   a-mode Emittance:     ", mode%a%emittance
    write(lun,*) "!   b-mode Emittance:     ", mode%b%emittance
    write(lun,*) "!   Energy Spread:        ", mode%sigE_E
    write(lun,*) "!   Bunch Length:         ", mode%sig_z
    write(lun,*) 
    write(lun,*) "!   Energy Loss/Turn:     ", mode%e_loss
    write(lun,*) "!   Momentum Compaction   ", mode%synch_int(1)/lat%param%total_length
    write(lun,*) "!   Slip Factor           ", mode%synch_int(1)/lat%param%total_length - (mass_of(lat%param%particle)/lat%ele(0)%value(E_tot$))**2
    write(lun,*) "!   I0 =                  ", mode%synch_int(0)
    write(lun,*) "!   I1 =                  ", mode%synch_int(1)
    write(lun,*) "!   I2 =                  ", mode%synch_int(2)
    write(lun,*) "!   I3 =                  ", mode%synch_int(3)
    write(lun,*) "!   I4a =                 ", mode%a%synch_int(4)
    write(lun,*) "!   I5a =                 ", mode%a%synch_int(5)
    write(lun,*) "!   I6a =                 ", mode%a%synch_int(6)
    write(lun,*) "!   I4b =                 ", mode%b%synch_int(4)
    write(lun,*) "!   I5b =                 ", mode%b%synch_int(5)
    write(lun,*) "!   I6b =                 ", mode%b%synch_int(6)
    write(lun,*) "!   I4z =                 ", mode%z%synch_int(4)
    write(lun,*) "!   I5z =                 ", mode%z%synch_int(5)
    write(lun,*) "!   I6z =                 ", mode%z%synch_int(6)
    write(lun,*) "!    Ja =                 ", mode%a%j_damp
    write(lun,*) "!    Jb =                 ", mode%b%j_damp
    write(lun,*) "!    Jz =                 ", mode%z%j_damp
    write(lun,*) "! damping tau_a =         ", lat%param%total_length / c_light / mode%a%alpha_damp
    write(lun,*) "! damping tau_b =         ", lat%param%total_length / c_light / mode%b%alpha_damp
    write(lun,*) "! damping tau_z =         ", lat%param%total_length / c_light / mode%z%alpha_damp
  end subroutine
end program





