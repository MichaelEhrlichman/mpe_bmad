program ptc_jac

USE pointer_lattice
USE orm_mod

implicit none

integer, parameter :: max_fam = 500
integer, parameter :: max_terms = 100

type :: terms_struct
  integer :: stage = 0
  integer F_index  !used to index complex vector field.  Only for dhdj
  character(5) name
  character(10) str
  real(dp) target
  real(dp) weight
end type
type(terms_struct) termdefs(max_terms)

character(100) lat_file
character(nlp) magname(max_fam)
character(50) jac_file_name
character(30) format_spec
character(2) nfam_str
character(1) F_index_str
character(7) current_terms_str(max_terms)

type(layout), pointer :: r
TYPE(FIBRE),POINTER :: p

!PTC datatypes - Not allocatable
type(pol_block) polyb
type(probe) ray
type(internal_state), target :: state

!PTC datatypes - Must be allocated
type(c_damap) cmap, L
type(c_vector_field) F
type(c_normal_form) cnorm
type(damapspin) maps
type(damapspin) ids   ! Initialized to identity map
type(c_taylor) t
type(c_taylor) h
type(probe_8) ray8
type(normalform) normal
type(damap) rdamap

integer n_terms
integer nfam
integer magtype(max_fam)
integer i, j, k, ks
integer nvars

logical is_new

real(dp) x(6)
real(dp) prec
real(dp) kmag(max_fam)
real(dp) mag
real(dp) der
real(dp), allocatable :: jac_buf(:)
real(dp), allocatable :: jac(:,:)

! Set global parameters
use_info=.true.
my_estate=>state

! Set local parameters
prec=1.d-8

!----------------------------------------------Initialize PTC and read in universe
call ptc_ini_no_append
call getarg(1,lat_file)
call READ_AND_APPEND_VIRGIN_general(M_U,lat_file)

c_verbose=.false.

r=>m_u%start

! state=delta0
state = default + nocavity0

use_complex_in_ptc=my_true

call init(state,5,1)

call alloc(ids)
call alloc(maps)
call alloc(ray8)
call alloc(cmap)
call alloc(cnorm)
call alloc(L)
call alloc(F)
call alloc(t)
call alloc(h)
call alloc(rdamap)
call alloc(normal)

              !terms_struct(stage, name, target, weight)
termdefs(1)  = terms_struct(1, 0, '21000', 'h21000', 0.0, 0.1)         ! Qx
termdefs(2)  = terms_struct(1, 0, '30000', 'h30000', 0.0, 0.1)         ! 3Qx
termdefs(3)  = terms_struct(1, 0, '10110', 'h10110', 0.0, 0.1)         ! Qx
termdefs(4)  = terms_struct(1, 0, '10020', 'h10020', 0.0, 0.1)         ! Qx-2Qy
termdefs(5)  = terms_struct(1, 0, '10200', 'h10200', 0.0, 0.1)         ! Qx+2Qy
termdefs(6)  = terms_struct(1, 0, '20001', 'h20001', 0.0, 0.1)         ! Synchro-betatron resonances ???
termdefs(7)  = terms_struct(1, 0, '00201', 'h00201', 0.0, 0.1)         ! Momentum-dependence of beta functions ???
termdefs(8)  = terms_struct(1, 0, '10002', 'h10002', 0.0, 0.1)         ! 2nd order dispersion

! termdefs(9)  = terms_struct(2, 0, '31000', 0.0, 0.0001)         ! 2Qx
! termdefs(10) = terms_struct(2, 0, '00310', 0.0, 0.0001)         ! 2Qy
! termdefs(11) = terms_struct(2, 0, '11200', 0.0, 0.0001)         ! 2Qy
! termdefs(12) = terms_struct(2, 0, '20110', 0.0, 0.0001)         ! 2Qx
! 
! termdefs(13) = terms_struct(3, 0, '22000', 0.0, 0.01)         ! dQx/dJx
! termdefs(14) = terms_struct(3, 0, '00220', 0.0, 0.01)         ! dQy/dJy
! termdefs(15) = terms_struct(3, 0, '11110', 0.0, 0.01)         ! dQx/dJy

! termdefs(16) = terms_struct(3, 0, '20200', 0.0, 0.0001)         ! 2Qx+2Qy
! termdefs(17) = terms_struct(3, 0, '20020', 0.0, 0.0001)         ! 2Qx-2Qy
! termdefs(18) = terms_struct(3, 0, '40000', 0.0, 0.0001)         ! 4Qx
! termdefs(19) = terms_struct(3, 0, '00400', 0.0, 0.0001)         ! 4Qy

termdefs(9)   = terms_struct(1, 1, '00001', 'chrom_x',  1.0, 1.)     ! chrom x
termdefs(10)  = terms_struct(2, 1, '00002', 'chrom_x2', 0.0, 1.0)     ! 2nd order chrom x
termdefs(11)  = terms_struct(3, 1, '00003', 'chrom_x3', 0.0, 1.0)     ! 3rd order chrom x
termdefs(12)  = terms_struct(1, 2, '00001', 'chrom_y',  1.0, 1.)     ! chrom y
termdefs(13)  = terms_struct(2, 2, '00002', 'chrom_y2', 0.0, 1.0)     ! 2nd order chrom y
termdefs(14)  = terms_struct(3, 2, '00003', 'chrom_y3', 0.0, 1.0)     ! 3rd order chrom y
termdefs(15)  = terms_struct(6, 6, '00001', 'dz/dp',   -0.1, 1.0)     ! dz/dp: first order momentum compaction (related by negative sign)
termdefs(16)  = terms_struct(6, 6, '00002', 'd2z/dp2',  0.0, 1.0)     ! d2z/dp2: second order momentum compaction (related by negative sign)
termdefs(17)  = terms_struct(6, 6, '00003', 'd3z/dp3', -1.0, 1.0)     ! d3z/dp3: third order momentum compaction (related by negative sign)
termdefs(18)  = terms_struct(6, 6, '00004', 'd4z/dp4', -5.0, 1.0)     ! d4z/dp4: third order momentum compaction (related by negative sign)

n_terms = 18

!--------------------------------------------- Locate multipole families
k=0
p=>r%start
do i=1,r%n
  if(p%mag%p%nmul > 2) then
    if(p%mag%bn(3) .ne. 0.d0) then
      is_new = .true.
      do j=1,k
        if(p%mag%name == magname(j)) then  !already found this family
          is_new = .false.
          exit
        endif
      enddo
      if(is_new) then 
        write(6,*) i,p%mag%name,p%mag%bn(3)
        k=k+1
        magname(k)=p%mag%name
        magtype(k)=3
      endif
      !call add(p,3,0,0.d0)   ! sets multipole field to zero
    elseif(p%mag%bn(4) .ne. 0.d0) then
      is_new = .true.
      do j=1,k
        if(p%mag%name == magname(j)) then  !already found this family
          is_new = .false.
          exit
        endif
      enddo
      if(is_new) then 
        write(6,*) i,p%mag%name,p%mag%bn(4)
        k=k+1
        magname(k)=p%mag%name
        magtype(k)=4
      endif
      !call add(p,3,0,0.d0)   ! sets multipole field to zero
    endif
  endif
p=>p%next
enddo
nfam=k ! Total number of multipoles
write(nfam_str,'(I2)') nfam
write(6,*) "Number of multitupole families found: ", nfam

allocate(jac_buf(n_terms))
allocate(jac(n_terms,nfam))

jac=0.0_dp
jac_buf = 0.0_dp

!------------------------------------- 
! For each multipole family i, get dh_xxxxx/dSi, for each driving term h_xxxxx.
!------------------------------------- 
do i=1,nfam
  ks = i
  polyb=0                     ! pol_block:  polymorphic_block.  sets a knob.
  polyb%ibn(magtype(ks))=1    ! set normal multipole strength as first parameter
  polyb%name=magname(ks)      ! set name of family
  r=polyb                     ! Overloaded to scan_for_polymorphs(layout,pol_block).  
                              ! For all magnets in r with name == polyb%name, bn(magtype(ks)) becomes
                              ! a knob.

  ids = 1
  x=0.d0 
  ray=0  
  ray=x  
  ray8=ray+ids  ! ray8 is a probe_8

  call track_probe(r,ray8,+state,fibre1=1)  ! get the 1-turn map by mapping the identity map.
                                            ! + symbol activates the knob
  call kill_para(r)  ! removes all polymorphic blocks (pol_block) from a layout

  ! Process 1-turn map into dhdj, which contains the amplitude-dependent tune shifts.
  rdamap = ray8    ! damap <- probe_8
  normal = rdamap  ! normalform <- damap

  ! Process 1-turn map into h, which contains the h_xxxxx driving terms.
  maps=ray8   ! damapspin <- probe_8
  cmap=maps   ! c_damap <- damapspin
  call c_normal(cmap,cnorm,dospin=my_false,no_used=1)  ! c_normal_form <- c_damap
  cmap=cnorm%n
  cmap=from_phasor(-1)*cmap*from_phasor(1)
  call c_factor_map(cmap,L,F,1)
  call c_clean_vector_field(F,F,prec)
  h=cgetpb(F,-1)   ! cgetpb includes 1/(-2*i_)   ! TODO: figure out crash if cgetpb replaced by getpb

  nvars = 0
  do j=1, n_terms
    call process_term(h, normal%dhdj%v, termdefs(j)%F_index, termdefs(j)%name, mag, der)

    nvars = nvars + 1
    jac_buf(nvars) = der
    write(F_index_str,'(I1)') termdefs(j)%F_index
    current_terms_str(nvars) = F_index_str//'_'//termdefs(j)%name
  enddo

  write(*,*) "Completed Family Number ",ks, " of ", nfam
  jac(1:nvars,ks) = jac_buf(1:nvars)
enddo

open(111,file='ptc_jac.out')
format_spec = '(A10,'//nfam_str//'A14)'
write(111,format_spec) '#         ', magname(1:nfam)
format_spec = '(A10,'//nfam_str//'ES14.4)'
do j=1, nvars
  write(111,format_spec) termdefs(j)%str, jac(j,1:nfam)
enddo
close(111)

deallocate(jac_buf)
deallocate(jac)

call ptc_end

contains

  subroutine process_term(h_local, v_local, F_index, lie_index, mag, der)
  
    implicit none

    type(c_taylor) h_local
    type(taylor) v_local(ndim2)
    integer F_index
    character(5) lie_index
    real(dp) mag
    real(dp), optional :: der

    !working space
    type(taylor) t_local
    type(c_taylor) c_t_local
    real(dp) hr, hi, dhrdS, dhidS


    if( F_index .gt. 0 ) then
      call alloc(t_local)
      t_local = v_local(F_index)
      t_local = t_local.par.lie_index
      t_local=t_local<=6  ! shift indexes down by 6, which moves the knob index 7 into
                          ! the first position.

      mag = t_local.sub.'0'
      if( present(der) ) then
        der = t_local.sub.'1'
      endif
      call kill(t_local)
    else
      call alloc(c_t_local)
      c_t_local=h_local.par.lie_index  ! pull desired polynomial coefficient from h.
      c_t_local=c_t_local<=6  ! shift indexes down by 5, which moves the knob index 6 into
                          ! the first position.

      dhrdS = real(c_t_local).sub.'1'
      dhidS = aimag(c_t_local).sub.'1'
      hr = real(c_t_local).sub.'0'
      hi = aimag(c_t_local).sub.'0'
      mag = sqrt(hr**2 + hi**2)
      if( present(der) ) then
        der = (hr*dhrdS+hi*dhidS) / mag
      endif
      call kill(c_t_local)
    endif

  end subroutine process_term

end program ptc_jac








