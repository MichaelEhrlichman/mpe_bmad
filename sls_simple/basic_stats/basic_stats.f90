program basic_stats

use bmad
use bmad_parser_mod, only: bp_com

implicit none

character lat_file*200

type(lat_struct) lat
type(normal_modes_struct) mode
type(coord_struct), allocatable :: co(:)
real(rp) cor_chrom_x, cor_chrom_y
real(rp) nat_chrom_x, nat_chrom_y
integer ix_cache
integer status

bp_com%always_parse = .true.
bmad_com%radiation_damping_on = .false.
bmad_com%radiation_fluctuations_on = .false.

WRITE(*,*) "Preparing lattice for basic_stats ..."
call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

call twiss_and_track(lat,co,status)
call chrom_calc(lat,1.0D-6,cor_chrom_x,cor_chrom_y)
ix_cache = -1
CALL radiation_integrals(lat, co, mode, ix_cache)

call set_on_off(sextupole$,lat,off$)
call chrom_calc(lat,1.0D-6,nat_chrom_x,nat_chrom_y)

write(*,*) "!   a-mode tune:          ", lat%ele(lat%n_ele_track)%a%phi /2./pi
write(*,*) "!   b-mode tune:          ", lat%ele(lat%n_ele_track)%b%phi /2./pi
write(*,*) "!   a-mode Emittance:     ", mode%a%emittance
write(*,*) "!   straight beta x:      ", lat%ele(0)%a%beta
write(*,*) "!   straight beta y:      ", lat%ele(0)%b%beta
write(*,*) "!   straight eta x:       ", lat%ele(0)%a%eta
write(*,*) "!   chromaticity x:       ", cor_chrom_x
write(*,*) "!   chromaticity y:       ", cor_chrom_y
write(*,*) "!     nat chrom  x:       ", nat_chrom_x
write(*,*) "!     nat chrom  y:       ", nat_chrom_y
write(*,*) "!               Jx:       ", mode%a%j_damp
write(*,*) "!               Jy:       ", mode%b%j_damp
write(*,*) "!               Jz:       ", mode%z%j_damp

open(60,file='basic_stats.out')
write(60,*) "!   a-mode tune:          ", lat%ele(lat%n_ele_track)%a%phi /2./pi
write(60,*) "!   b-mode tune:          ", lat%ele(lat%n_ele_track)%b%phi /2./pi
write(60,*) "!   a-mode Emittance:     ", mode%a%emittance
write(60,*) "!   straight beta x:      ", lat%ele(0)%a%beta
write(60,*) "!   straight beta y:      ", lat%ele(0)%b%beta
write(60,*) "!   straight eta x:       ", lat%ele(0)%a%eta
write(60,*) "!   chromaticity x:       ", cor_chrom_x
write(60,*) "!   chromaticity y:       ", cor_chrom_y
write(60,*) "!     nat chrom  x:       ", nat_chrom_x
write(60,*) "!     nat chrom  y:       ", nat_chrom_y
write(60,*) "!               Jx:       ", mode%a%j_damp
write(60,*) "!               Jy:       ", mode%b%j_damp
write(60,*) "!               Jz:       ", mode%z%j_damp
close(60)

end program





