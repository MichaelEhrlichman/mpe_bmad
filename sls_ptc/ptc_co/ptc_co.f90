PROGRAM ptc_closed_orbit

USE bmad
USE sim_utils
!USE ptc_layout_mod

IMPLICIT none

CHARACTER lat_file*200

INTEGER i
integer radcache
real(rp) t6(6,6)

TYPE(lat_struct) lat
TYPE(coord_struct), ALLOCATABLE :: co(:)
TYPE(coord_struct), ALLOCATABLE :: ptc_co(:)
TYPE(coord_struct), ALLOCATABLE :: orb(:)
TYPE(coord_struct), ALLOCATABLE :: ptc_orb(:)
type(coord_struct) init_guess
type(normal_modes_struct) mode

CALL GETARG(1, lat_file)

WRITE(*,*) "Preparing lattice..."
bmad_com%radiation_damping_on = .true.
bmad_com%radiation_fluctuations_on = .false.
CALL bmad_parser(lat_file, lat)

CALL set_on_off(rfcavity$, lat, on$)
CALL closed_orbit_calc(lat,co,6)
CALL lat_make_mat6(lat, -1, co)
WRITE(*,'(A,6ES14.4)') "bmad closed orbit calc: ", co(0)%vec(1:6)

do i=1,lat%n_ele_track
  write(98,'(I6,6ES14.4)') i, co(i)%vec(1:6)
enddo

! CALL lat_make_mat6(lat, -1, co)
! CALL transfer_matrix_calc (lat, .true., t6, ix1=0, one_turn=.TRUE.)
! write(*,*) "Determinant: ", determinant(t6)

CALL twiss_at_start(lat)
CALL twiss_propagate_all(lat)
call radiation_integrals(lat, co, mode)
write(*,*) "bunch length: ", mode%sig_z

! CALL lat_to_ptc_layout(lat)
! call reallocate_coord(ptc_co, lat)
! call ptc_closed_orbit_calc (lat%branch(0), ptc_co)
! WRITE(*,'(A,6ES14.4)') "ptc closed orbit:  ", ptc_co(0)%vec(1:6)

! call reallocate_coord(ptc_orb, lat)
! ptc_orb(0)%vec(:) = 0.0d0
! ptc_orb(0)%vec(6) = 5.0e-2
! call ptc_track_all (lat%branch(0), ptc_orb)
! 
! call reallocate_coord(orb, lat)
! orb(0)%vec(:) = 0.0d0
! orb(0)%vec(6) = 5.0e-2
! call track_all (lat, orb)
! 
! OPEN(45,FILE='orb.out')
! do i=1,lat%n_ele_track
!   WRITE(45,'(I6,6ES14.4)') i, orb(i)%vec(1:6)
! enddo
! CLOSE(45)
! 
! OPEN(45,FILE='ptc_orb.out')
! do i=1,lat%n_ele_track
!   WRITE(45,'(I6,6ES14.4)') i, ptc_orb(i)%vec(1:6)
! enddo
! CLOSE(45)

END PROGRAM ptc_closed_orbit












