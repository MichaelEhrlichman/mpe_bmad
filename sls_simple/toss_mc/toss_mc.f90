program toss_mc
  use bmad
  use bmad_parser_mod, only: bp_com

  implicit none

  type (lat_struct) lat
  type (coord_struct), allocatable :: traj(:)
  type(coord_struct) orb1

  integer i,j,k
  integer ios, npart
  integer ele_in, ele_out

  real(rp) vec_offset(6)
  real(rp) s, delta_s
  real(rp) tot_particles, tot_good, kDz_mu0
  real(rp), allocatable :: part_in(:,:), part_out(:,:)

  character*500 line
  character*100 param_file
  character*100 lat_file
  character*100 part_file
  character*100 out_file
  
  namelist /toss_mc_params/ lat_file, part_file, ele_in, ele_out

  bp_com%always_parse = .true.

  call getarg(1,param_file)
  open(10,file=param_file)
  read(10,nml=toss_mc_params)
  close(10)

  call file_suffixer(part_file, out_file, 'tossed', .true.)

  call bmad_parser(lat_file, lat)

  allocate(traj(0:lat%n_ele_track))

  open(55,file=part_file)
  npart = 0
  do while(.true.)
    read(55,*,iostat=ios) line
    if(ios /= 0) exit
    npart = npart + 1
  enddo
  rewind(55)
  allocate(part_in(npart,7))  ! 6-vec and npart
  do i=1,npart
    read(55,*) part_in(i,:)
  enddo
  close(55)

  tot_particles = sum(part_in(:,7))

  write(*,*) "Total particles: ", tot_particles
  write(*,*) "First ele: ", ele_in, lat%ele(ele_in)%name
  write(*,*) "Last ele:  ", ele_out, lat%ele(ele_out)%name

  open(65,file=out_file)
  tot_good = 0
  do i=1,npart
    !call init_coord(traj(ele_in),part_in(i,:),lat%ele(ele_in),element_end=upstream_end$)
    traj(ele_in)%vec = part_in(i,1:6)
    call track_many(lat, traj, ele_in, ele_out,1)
    !kDz_mu0 = traj(ele_out)%vec(5) * 2.31199e6  ! z * k / mu0, k=1um, mu0=2.4
    kDz_mu0 = traj(ele_out)%vec(5) * 3.26568e6  ! z * k / mu0
    if(abs(kDz_mu0) .le. 1.0) tot_good = tot_good + part_in(i,7)
    write(65,'(6es12.3,es12.3)') traj(ele_out)%vec(:), part_in(i,7)/tot_particles
  enddo
  close(65)

  write(*,'(a,f11.4)') "fraction in OSC dampable envelope (kDz < mu0): ", tot_good / tot_particles

end program





