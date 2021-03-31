program special_chrom_set
  use bmad
  use bmad_parser_mod, only: bp_com
  use dynap_mod, only: set_magnet_strengths
  use moga_struct_mod

  implicit none

  type (lat_struct) ring
  type (coord_struct), allocatable :: co(:)
  type (mag_struct) chrom_mags(2)

  integer i
  integer status

  real(rp) init_chrom_x, init_chrom_y 
  real(rp) chrom_rep_x, chrom_rep_y
  real(rp) set_chrom_x, set_chrom_y
  real(rp) A(2,2), Ainv(2,2), chrom_vec(2), chrom_vars(2)

  logical err_flag, ok

  character(30) lat_file
  character(30) in_file
  character(10) chrom_mag_1, chrom_mag_2

  namelist /special_chrom_set_nl/ lat_file, set_chrom_x, set_chrom_y, chrom_mag_1, chrom_mag_2

  chrom_mags(1:2)%property = 'k2'

  call getarg(1,in_file)

  open (unit = 10, file = in_file, action='read')
  read (10, nml = special_chrom_set_nl)
  close (10)

  chrom_mags(1)%name = chrom_mag_1
  chrom_mags(2)%name = chrom_mag_2

  bp_com%always_parse = .true.

  call bmad_parser(lat_file, ring)

  call twiss_and_track(ring,co,status)

  call set_magnet_strengths(chrom_mags,ring,(/0.0d0,0.0d0/))
  call chrom_calc(ring, 1.0d-6, init_chrom_x, init_chrom_y, err_flag)

  call set_magnet_strengths(chrom_mags,ring,(/1.0d0,0.0d0/))
  call twiss_and_track(ring,co,status)
  call chrom_calc(ring, 1.0d-6, chrom_rep_x, chrom_rep_y, err_flag)
  A(1,1) = (chrom_rep_x-init_chrom_x)/1.0d0
  A(2,1) = (chrom_rep_y-init_chrom_y)/1.0d0

  call set_magnet_strengths(chrom_mags,ring,(/0.0d0,1.0d0/))
  call twiss_and_track(ring,co,status)
  call chrom_calc(ring, 1.0d-6, chrom_rep_x, chrom_rep_y, err_flag)
  A(1,2) = (chrom_rep_x-init_chrom_x)/1.0d0
  A(2,2) = (chrom_rep_y-init_chrom_y)/1.0d0

  call mat_inverse(A,Ainv,ok,.true.)
  if(.not. ok) then
    write(*,*) "Chromaticity response matrix inversion failed."
  endif

  chrom_vec(1) = set_chrom_x - init_chrom_x
  chrom_vec(2) = set_chrom_y - init_chrom_y
  chrom_vars = matmul(Ainv,chrom_vec)

  open(10,file='special_chrom_set.bmad')
  write(10,'(a)') 'call, file = '//lat_file
  write(10,'(A,es14.5)') trim(adjustl(chrom_mag_1))//'[k2]=', chrom_vars(1)
  write(10,'(A,es14.5)') trim(adjustl(chrom_mag_2))//'[k2]=', chrom_vars(2)
  close(10)

end program
                                                            










