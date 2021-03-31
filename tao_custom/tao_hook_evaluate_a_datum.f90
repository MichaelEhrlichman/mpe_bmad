!+
! Subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value)
!
! See the Programmer's manual for how to add custom data types here.
!
! Input:
!   datum        -- tao_data_struct: The current datum to evaluate
!   u            -- tao_universe_struct: Universe this datum is in
!   tao_lat      -- tao_lattice_struct: Lattice to use.
!
! Output:
!   found        -- Logical: True if  this datum is evaluated in this subroutine.
!   datum_value  -- real(rp): Which datum value to compute (model_value, design_value, etc...)
!   valid_value  -- Logical: Set False when there is a problem. Set True otherwise.
!   why_invalid  -- Character(*), optional: Tells why datum value is invalid.
!-

subroutine tao_hook_evaluate_a_datum (found, datum, u, tao_lat, datum_value, valid_value, why_invalid)

use tao_data_and_eval_mod, dummy => tao_hook_evaluate_a_datum
use, intrinsic :: iso_c_binding
!$use omp_lib
!test change

implicit none

type (tao_universe_struct), target :: u
type (tao_data_struct) datum
type (tao_lattice_struct), target :: tao_lat

real(rp) datum_value
logical found, valid_value

character(*), optional :: why_invalid
character(*), parameter :: r_name = 'tao_hook_evaluate_a_datum'

!my stuff
integer k, counter
type(lat_struct), pointer :: lat
type(branch_struct), pointer :: branch

type(bmad_normal_form_struct), pointer :: normal_form
integer i, order_for_normal_form
integer order
logical :: rf_on, c_verbose_save

!for 'survival' data_type
integer j, saved_ele
integer nturns, N, ix1, track_state
real(rp) nsurvive, ndoa
type (tao_model_branch_struct), pointer :: model_branch
type(coord_struct), allocatable :: orbit(:)
integer ix_branch

logical dump_dist
type vec_struct
  logical alive
  real(rp) vec(6)
end type
type(vec_struct), allocatable :: saved_coords(:,:)
character(1) ix_uni_str
logical loud
character(1), allocatable :: status_array(:)

ix_branch = datum%ix_branch
lat => tao_lat%lat
branch => tao_lat%lat%branch(ix_branch)

found = .false.
datum_value = 0
dump_dist = .true.
saved_ele = 157  !ap_start
loud = .false.

if(trim(datum%data_type) == 'survival') then
  if(datum%good_user) then
    model_branch => u%model_branch(ix_branch)
    ix1 = datum%ix_ele
    N = size(model_branch%ele(ix1)%beam%bunch(1)%particle(:))
    allocate(status_array(1:N))
    status_array(:) = "!"
    nturns = 300 !250
    allocate(orbit(0:branch%n_ele_track))
    if(dump_dist) then
      write(ix_uni_str,'(i1)') u%ix_uni
      open(1010,FILE='beam_coords_'//ix_uni_str//'.dat')
      allocate(saved_coords(nturns,N))
      saved_coords(:,:)%alive = .false.
    endif

    nsurvive = model_branch%ele(ix1)%beam%bunch(1)%charge_tot
    ndoa = 0.0d0
    write(*,*) "Custom survival tracking starting ..."
    !$OMP PARALLEL DO &
    !$OMP DEFAULT(FIRSTPRIVATE), &
    !$OMP SHARED(model_branch,N,nturns,lat,ix1,ix_branch,dump_dist,saved_coords,branch,saved_ele,loud,status_array), &
    !$OMP PRIVATE(orbit,track_state,j), &
    !$OMP REDUCTION(-:nsurvive), REDUCTION(+:ndoa)
    do i=1,N
      track_state = moving_forward$
      orbit(ix1) = model_branch%ele(ix1)%beam%bunch(1)%particle(i)
      if(loud) write(*,'(2(a,i6))') "Particle ", i, " of ", N
      if(orbit(ix1)%state == alive$) then
        do j=1,nturns
          call track_many (lat, orbit, ix1, ix1, 1, track_state=track_state, ix_branch=ix_branch)
          if(dump_dist .and. (track_state == moving_forward$)) then
            saved_coords(j,i)%alive = .true.
            saved_coords(j,i)%vec(1:6) = orbit(saved_ele)%vec(1:6)
          endif
          if(track_state .ne. moving_forward$) then
            if(loud) then
              write(*,'(a,i6,a,a,a,i6)') "tao custom survival tracking.  particle ", i, " lost at ",  branch%ele(track_state)%name, " on turn ", j
            else
              status_array(i) = 'X'
            endif
            nsurvive = nsurvive - model_branch%ele(ix1)%beam%bunch(1)%particle(i)%charge
            exit
          endif
        enddo
      else
        if(loud) then
          write(*,'(a,i6,a)') "Particle ", i, " is dead on arrival (probably lost in upstream branch)."
        else
          status_array(i) = 'O'
        endif
        ndoa = ndoa + model_branch%ele(ix1)%beam%bunch(1)%particle(i)%charge
        nsurvive = nsurvive - model_branch%ele(ix1)%beam%bunch(1)%particle(i)%charge
      endif    
      if(j==nturns+1) then
        status_array(i) = '.'
      endif
    enddo
    !$OMP END PARALLEL DO
    if(dump_dist) then
      do j=1,nturns
        do i=1,N
          if(saved_coords(j,i)%alive) then
            write(1010,'(2i6,6es14.4)') j, i, saved_coords(j,i)%vec(1:6)
          endif
        enddo
        write(1010,*)
        write(1010,*)
      enddo
    endif
    deallocate(orbit)
    if(dump_dist) then
      deallocate(saved_coords)
      close(1010)
    endif
    datum_value = nsurvive / model_branch%ele(ix1)%beam%bunch(1)%charge_tot
    valid_value = .true.
    found=.true.
    write(*,'(10000a1)') status_array(1:N)
    write(*,'(a,f14.5)') "Percentage of charge dead on arrival: ", ndoa / model_branch%ele(ix1)%beam%bunch(1)%charge_tot * 100.0
    write(*,'(a,f14.5)') "Percentage of charge surviving (DOA + lost during tracking): ", datum_value * 100.0
    deallocate(status_array)
  else
    datum_value = -1
    valid_value = .false.
    why_invalid = 'disabled by user (good_user)'
    found=.true.
  endif
endif

end subroutine tao_hook_evaluate_a_datum













