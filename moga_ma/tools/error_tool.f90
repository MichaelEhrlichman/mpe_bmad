program error_tool

use bmad
use bmad_parser_mod
use extended_names_mod

implicit none

  type error_struct
    integer keys(10) 
    character(20) mask
    character(20) property
    real(rp) cutoff
    real(rp) rms
  end type

  ! namelist general
  character(100) lat_file
  character(100) lat_file_override
  logical use_hybrid
  integer iargs

  type (lat_struct) ring
  integer i,j,k
  integer attribute_ix
  integer, parameter :: MAX_ERR_SETS=40
  real(rp) rnum
  real(rp) magnet_error, magnet_error2
  real(rp) center_offset, angle
  character(60) in_file
  character(10) seed_str, num_str
  character(20), allocatable :: extended_names(:)
  logical attribute_is_free

  ! namelist moga_errors
  integer magnet_error_seed
  type(error_struct) error(MAX_ERR_SETS)

  namelist / general /    lat_file, use_hybrid
  namelist / moga_errors / magnet_error_seed, error

  call getarg(1,in_file)
  lat_file_override = ''
  iargs = iargc()
  if (iargs == 2 ) then
    call getarg(2,lat_file_override)
  endif

  error(:)%mask = 'null'
  do i=1,MAX_ERR_SETS
    error(i)%keys(:) = -1
  enddo

  open (unit = 10, file = in_file, action='read')
  read (10, nml = general)
  read (10, nml = moga_errors)

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  bp_com%always_parse = .true.
  call bmad_parser(lat_file,ring)
  allocate(extended_names(0:ring%n_ele_max))
  call make_extended_names(ring,extended_names)

  write(seed_str,'(i6.6)') magnet_error_seed
  !open(45, file='moga_errors_'//trim(seed_str)//'.lat') 
  open(45, file='errors.lat') 
  !write(45,'(a)') 'call, file = '//lat_file
  write(45,'(a)') 'expand_lattice'
  call ran_seed_put(magnet_error_seed)
  do i=1,MAX_ERR_SETS
    if(error(i)%mask .ne. 'null') then
      do j=1,ring%n_ele_max
        if( ( str_match_wild(ring%ele(j)%name,error(i)%mask) .and. any(ring%ele(j)%key==error(i)%keys(:)) ) .or. &
            ( str_match_wild(ring%ele(j)%name,error(i)%mask) .and. all(error(i)%keys(:) == -1) ) ) then
          attribute_ix = attribute_index(ring%ele(j),error(i)%property)
          attribute_is_free = attribute_free(ring%ele(j), error(i)%property, err_print_flag=.false.)
          if( (attribute_ix .gt. 0 .and. attribute_is_free) .or. str_match_wild(error(i)%property,'*_ENDS*') ) then
            do while(.true.)
              call ran_gauss(rnum)
              if(abs(rnum) .lt. error(i)%cutoff) exit
            enddo
            magnet_error = rnum * error(i)%rms
            !write(num_str,'(i6.6)') j
            !if(ring%ele(j)%key == 37) then  !girder
            if( str_match_wild(error(i)%property,'*_ENDS*') ) then  !girder
              do while(.true.)
                call ran_gauss(rnum)
                if(abs(rnum) .lt. error(i)%cutoff) exit
              enddo
              magnet_error2 = rnum * error(i)%rms
              center_offset = (magnet_error+magnet_error2)/2.0d0
              angle = (magnet_error2-magnet_error)/value_of_attribute(ring%ele(j),'L')  ! sin theta = (x2-x1)/l
              if(trim(error(i)%property) == 'X_OFFSET_ENDS') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[X_OFFSET]=',center_offset,'   ! '//trim(ring%ele(j)%name)
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[X_PITCH]=',angle,'   ! '//trim(ring%ele(j)%name)
              elseif(trim(error(i)%property) == 'Y_OFFSET_ENDS') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[Y_OFFSET]=',center_offset,'   ! '//trim(ring%ele(j)%name)
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[Y_PITCH]=',angle,'   ! '//trim(ring%ele(j)%name)
              endif
            else
              if(trim(error(i)%property) == 'K1') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[K1]=',(1.0+magnet_error)*value_of_attribute(ring%ele(j),'K1'),'   ! '//trim(ring%ele(j)%name)
              elseif(trim(error(i)%property) == 'A1') then  !skew gradient errors
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[A1]=',magnet_error*value_of_attribute(ring%ele(j),'K1')*value_of_attribute(ring%ele(j),'L'),'   ! '//trim(ring%ele(j)%name)
              elseif(trim(error(i)%property) == 'X_OFFSET') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[X_OFFSET]=',magnet_error,'   ! '//trim(ring%ele(j)%name)
              elseif(trim(error(i)%property) == 'Y_OFFSET') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[Y_OFFSET]=',magnet_error,'   ! '//trim(ring%ele(j)%name)
              elseif(trim(error(i)%property) == 'TILT') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[TILT]=',magnet_error,'   ! '//trim(ring%ele(j)%name)
              elseif(trim(error(i)%property) == 'ROLL') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[ROLL]=',magnet_error,'   ! '//trim(ring%ele(j)%name)
              elseif(trim(error(i)%property) == 'X_GAIN_ERR') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[X_GAIN_ERR]=',magnet_error,'   ! '//trim(ring%ele(j)%name)
              elseif(trim(error(i)%property) == 'Y_GAIN_ERR') then
                write(45,'(a,es15.6,a)') trim(extended_names(j))//'[Y_GAIN_ERR]=',magnet_error,'   ! '//trim(ring%ele(j)%name)
              else
                write(*,*) "Error property not yet implemented: ", trim(error(i)%property)
              endif
            endif
          elseif( attribute_ix .le. 0) then
            write(*,*) "When applying errors, element does not have requested property."
            write(*,*) "Name and property: ", ring%ele(j)%name, error(i)%property
          endif
        endif
      enddo
    endif
  enddo
  close(45)
end program
