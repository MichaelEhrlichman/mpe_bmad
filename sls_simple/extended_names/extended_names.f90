program extended_names

use bmad

implicit none

integer, allocatable :: order_index(:)
integer, allocatable :: nt2(:)
logical, allocatable :: more_than_one(:)
integer counter
character(100) curstr
character(6) str

character lat_file*200
type(lat_struct) lat
integer i, j
integer ixa

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

curstr = lat%nametable%name(lat%nametable%indexx(lat%nametable%n_min))
ixa = 0
allocate(nt2(lat%nametable%n_min:lat%nametable%n_max))
allocate(more_than_one(lat%nametable%n_min:lat%nametable%n_max))
more_than_one = .false.
do i=lat%nametable%n_min+1,lat%nametable%n_max
  if(trim(lat%nametable%name(lat%nametable%indexx(i))) .ne. trim(curstr)) then
    curstr = lat%nametable%name(lat%nametable%indexx(i))
    nt2(ixa:i-1) = lat%nametable%indexx(ixa:i-1)
    call shell_sort_i(nt2(ixa:i-1))
    if ( ixa .lt. i-1 ) then
      more_than_one(ixa:i-1) = .true. 
    endif
    ixa = i
  endif
enddo
nt2(ixa:) = lat%nametable%indexx(ixa:)
call shell_sort_i(nt2(ixa:))
if ( ixa .lt. lat%nametable%n_max ) then
  more_than_one(ixa:lat%nametable%n_max) = .true. 
endif

curstr = ''
counter = 0
allocate(order_index(lat%nametable%n_min:lat%nametable%n_max))
do i=lat%nametable%n_min,lat%nametable%n_max
  if(trim(lat%nametable%name(nt2(i))) .ne. trim(curstr)) then
    counter = 0
    curstr = lat%nametable%name(nt2(i))
  endif
  counter = counter + 1
  order_index(nt2(i)) = counter
enddo

open(45,file='extended_names.dat')
do i=0,lat%n_ele_track
  do j=lat%nametable%n_min,lat%nametable%n_max
    if(i == nt2(j)) exit
  enddo
  write(str,'(i6)') order_index(nt2(j))
  if( more_than_one(j) ) then
    write(45,'(i6,"   ",a20)') i, trim(lat%ele(i)%name)//"##"//adjustl(trim(str))
  else
    write(45,'(i6,"   ",a20)') i, trim(lat%ele(i)%name)
  endif
enddo

contains

  subroutine shell_sort_i(arr)
    implicit none
    integer, dimension(:), intent(inout) :: arr
    integer i,j,inc,n
    integer v
    n=size(arr)
    inc=1
    do
      inc=3*inc+1
      if (inc > n) exit
    end do
    do 
      inc=inc/3
      do i=inc+1,n
        v=arr(i)
        j=i
        do
          if (arr(j-inc) <= v) exit
          arr(j)=arr(j-inc)
          j=j-inc
          if (j <= inc) exit
        end do
        arr(j)=v
      end do
      if (inc <= 1) exit
    end do
  end subroutine
end program












