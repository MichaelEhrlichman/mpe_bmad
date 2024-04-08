program rmat

  use bmad
  use bmad_parser_mod, only: bp_com

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: orb(:)

  real(rp) mat(6,6)

  integer i,j
  logical exists

  real(rp) E_tot_start, E_tot_end
  real(rp) gamma0, gammaf, check_dg

  character*100 lat_file

  bp_com%always_parse = .true.

  call getarg(1,lat_file)
  call bmad_parser(lat_file, lat)

  !do i=1,lat%n_ele_track
  !  write(*,'(a,a)') "element: ", lat%ele(i)%name
  !  do j=1,6
  !    write(*,'(6es14.5)') lat%ele(i)%mat6(j,:)
  !  enddo
  !  write(*,*)
  !enddo

  mat = lat%ele(1)%mat6

  E_tot_start = lat%ele(1)%value(e_tot_start$)
  call convert_total_energy_to(E_tot_start,particle=electron$,gamma=gamma0)
  E_tot_end = lat%ele(1)%value(e_tot$)
  call convert_total_energy_to(E_tot_end,particle=electron$,gamma=gammaf)
  check_dg = (gammaf-gamma0)/gamma0

  inquire(file='rmat.out',exist=exists)
  if(.not. exists) then
    open(4,file='rmat.out',status='new',action='write')
    write(4,'(10a14)') "E0", "mat11", "mat12", "mat21", "mat22", "mat55", "mat56", "mat65", "mat66", 'dg/g0'
  else
    open(4,file='rmat.out',access='append',status='old',action='write')
  endif
  write(4,'(10es14.5)') lat%ele(0)%value(E_tot$), mat(1,1), mat(1,2), mat(2,1), mat(2,2), mat(5,5), mat(5,6), mat(6,5), mat(6,6), check_dg

  !do i=1,1 !lat%n_ele_track
  !  write(*,'(a,a)') "element: ", lat%ele(i)%name
  !  write(*,'(a)') "x,x' block:"
  !  do j=1,2
  !    write(*,'(2es14.5)') lat%ele(i)%mat6(j,1:2)
  !  enddo
  !  write(*,*)
  !enddo

  !do i=1,1 !lat%n_ele_track
  !  write(*,'(a)') "z,z' block:"
  !  do j=5,6
  !    write(*,'(2es14.5)') lat%ele(i)%mat6(j,5:6)
  !  enddo
  !  write(*,*)
  !enddo
end program rmat
