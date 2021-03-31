module set_on_off_value_mod
  use bmad

  implicit none

  type value_struct
    character(40) :: ele_name = ''
    real(rp) :: value = 0.0d0
  end type

contains

!+
! Subroutine set_on_off_value(lat,ele_key,value_ix,switch,store)
!
! If called with switch==off$, sets value(value_ix) to zero for each element
! in the lattice of type key_name(ele_key), and stores the original value in
! store.
!
! If called with switch==on$, sets value(value_ix) to the stored value in store
! for each element in the lattice of type key_name(ele_key).
!
! Input:
!   lat                  -- type(lat_struct): lattice
!   ele_key              -- integer:  matched to ele_struct%key.  Identified elements to switch on or off.
!   value_ix             -- integer:  index of ele_struct%value to switch on or off.
!   switch               -- integer:  off$ or on$
!   store                -- type(value_struct):  Preallocated array.  Stores previous values.  Keyed by element name.
! Output:
!   lat%ele(:)%value     -- real(rp): values set to either zero or restored
!-
subroutine set_on_off_value(lat,ele_key,value_ix,switch,store)
  implicit none

  type(lat_struct) lat
  integer ele_key
  integer value_ix
  integer switch
  type(value_struct) store(:)

  integer i, j
  integer store_counter
  integer n_loc
  logical error
  logical stored, found
  type (ele_pointer_struct), allocatable :: eles(:)

  store_counter = 0

  call lat_ele_locator(key_name(ele_key)//'::*', lat, eles, n_loc, error) 
  
  select case(switch)
    case(off$)
      store(:)%ele_name = ''
      store(:)%value = 0.0d0
      do i=1,n_loc
        stored = .false.
        do j=1,store_counter
          if( eles(i)%ele%name .eq. store(j)%ele_name ) then
            stored = .true.
            exit
          endif
        enddo
        if( .not. stored ) then
          store_counter = store_counter + 1
          store(store_counter)%ele_name = eles(i)%ele%name
          store(store_counter)%value = eles(i)%ele%value(value_ix)
        endif
        eles(i)%ele%value(value_ix) = 0.0d0
        call set_flags_for_changed_attribute(eles(i)%ele,eles(i)%ele%value(value_ix))
      enddo
    case(on$)
      do i=1,n_loc
        found = .false.
        j=1
        do while(trim(store(j)%ele_name) .ne. '')
          if( eles(i)%ele%name .eq. store(j)%ele_name ) then
            eles(i)%ele%value(value_ix) = store(j)%value
            found = .true.
            exit
          endif
          j=j+1
        enddo
        if( .not. found) then
          write(*,'(A)') "set_on_off_value unable to find element in store.  This is bad and your code is broken."
          if(global_com%exit_on_error) call err_exit 
          return
        endif
        call set_flags_for_changed_attribute(eles(i)%ele,eles(i)%ele%value(value_ix))
      enddo
    end select
  deallocate(eles)
end subroutine

end module
