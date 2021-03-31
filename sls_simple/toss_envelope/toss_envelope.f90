program toss_envelope
  use bmad
  use bmad_parser_mod, only: bp_com
  use mode3_mod

  implicit none

  type(lat_struct) lat
  type(coord_struct), allocatable :: traj(:)
  type(coord_struct), allocatable :: orb(:)

  integer i,j,k
  integer, parameter :: nsteps=200
  integer npart
  integer ele_in, ele_out
  integer status

  logical error

  real(rp) vec_offset(6)
  real(rp) s, delta_s
  real(rp), allocatable :: part_in(:,:), part_out(:,:)
  real(rp) x(6), emit(3)
  real(rp) klight
  real(rp) r, nsig

  character*100 param_file
  character*100 lat_file

  namelist /toss_env_params/ lat_file, ele_in, ele_out, emit, klight, npart, nsig

  bp_com%always_parse = .true.

  call getarg(1,param_file)
  open(10,file=param_file)
  read(10,nml=toss_env_params)
  close(10)

  !call file_suffixer(part_file, out_file, 'tossed', .true.)

  call bmad_parser(lat_file, lat)
  call twiss_and_track(lat,orb,status)
  call twiss3_at_start(lat,error)
  call twiss3_propagate_all(lat)

  allocate(traj(0:lat%n_ele_track))

  write(*,*) "First ele: ", ele_in, lat%ele(ele_in)%name
  write(*,*) "Last ele:  ", ele_out, lat%ele(ele_out)%name

  call random_seed()

  do i=1,npart
    call random_number(r); r = 2*r-1
    x(5) = 0.0d0
    x(6) = ssqrt(2*nsig*emit(3)*r*lat%ele(ele_in)%z%beta) 

    call random_number(r); r = 2*r-1
    x(1) = orb(ele_in)%vec(1) + lat%ele(ele_in)%a%eta*x(6) + ssqrt(2*nsig*emit(1)*r*lat%ele(ele_in)%a%beta) 
    x(2) = orb(ele_in)%vec(2) + lat%ele(ele_in)%a%etap*x(6) + lat%ele(ele_in)%a%alpha*ssqrt(2*nsig*emit(1)*r/lat%ele(ix)%a%beta)

    call random_number(r); r = 2*r-1
    x(3) = orb(ele_in)%vec(3) + lat%ele(ele_in)%b%eta*x(6) + ssqrt(2*nsig*emit(2)*r*lat%ele(ele_in)%b%beta) 
    x(4) = orb(ele_in)%vec(4) + lat%ele(ele_in)%b%etap*x(6) + lat%ele(ele_in)%b%alpha*ssqrt(2*nsig*emit(2)*r/lat%ele(ix)%b%beta)

    !call init_coord(traj(ele_in),part_in(i,:),lat%ele(ele_in),element_end=upstream_end$)
    traj(ele_in)%vec = part_in(i,:)
    call track_many(lat, traj, ele_in, ele_out,1)
    write(65,'(6es12.3)') traj(ele_out)%vec(:)
  enddo

  contains

  function rho(x,ele,emit)
    real(rp) x(6) !x, x', y, y', z, z'
    type(ele_struct) ele
    real(rp) emit(3) !emit_x, emit_y, emit_z
    real(rp) rho

    real(rp) Jx, Jy, Jz

    Jx = ele%a%gamma*x(1)**2 + 2*ele%a%alpha*x(1)*x(2) + ele%a%beta*x(2)**2
    Jy = ele%b%gamma*x(2)**2 + 2*ele%b%alpha*x(3)*x(4) + ele%b%beta*x(4)**2
    Jz = ele%z%gamma*x(5)**2 + 2*ele%z%alpha*x(5)*x(6) + ele%z%beta*x(6)**2

    rho = exp(-Jx/emit(1) -Jy/emit(2) - Jz/emit(3)) / twopi**3 / emit(1) / emit(2) / emit(3)

  end function

  function ssqrt(x)
    real(rp) x, ssqrt
    ssqrt = sign(1.0d0,x)*sqrt(abs(x))
  end function


end program





