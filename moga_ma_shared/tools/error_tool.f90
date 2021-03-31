program error_tool

use bmad

implicit none

  ! namelist general
  character(100) lat_file
  character(100) lat_file_override
  logical use_hybrid
  integer iargs

  type (lat_struct) ring
  integer i,j,k
  real(rp) rnum
  real(rp) magnet_error
  integer attribute_ix
  integer, parameter :: MAX_ERR_SETS=10
  character(60) in_file
  character(10) seed_str, num_str

  ! namelist moga_errors
  integer magnet_error_seed
  character(20) error_mask(MAX_ERR_SETS)
  character(10) error_property(MAX_ERR_SETS)
  real(rp) error_cutoff(MAX_ERR_SETS)
  real(rp) error_rms(MAX_ERR_SETS)

  namelist / general /    lat_file, use_hybrid
  namelist / moga_errors / magnet_error_seed, &
                           error_mask, error_property, error_cutoff, error_rms

  call getarg(1,in_file)
  lat_file_override = ''
  iargs = iargc()
  if (iargs == 2 ) then
    call getarg(2,lat_file_override)
  endif

  error_mask = 'null'

  open (unit = 10, file = in_file, action='read')
  read (10, nml = general)
  read (10, nml = moga_errors)

  if( lat_file_override .ne. '' ) then
    lat_file = lat_file_override
  endif

  bp_com%always_parse = .true.
  call bmad_parser(lat_file,ring)

  write(seed_str,'(i6.6)') magnet_error_seed
  !open(45, file='moga_errors_'//trim(seed_str)//'.lat') 
  open(45, file='errors.lat') 
  write(45,'(a)') 'call, file = '//lat_file
  write(45,'(a)') 'expand_lattice'
  call ran_seed_put(magnet_error_seed)
  do i=1,MAX_ERR_SETS
    if(error_mask(i) .ne. 'null') then
      do j=1,ring%n_ele_track
        if(str_match_wild(ring%ele(j)%name,error_mask(i))) then
          attribute_ix = attribute_index(ring%ele(j),error_property(i))
          if( attribute_ix .gt. 0 ) then
            do while(.true.)
              call ran_gauss(rnum)
              if(abs(rnum) .lt. error_cutoff(i)) exit
            enddo
            magnet_error = rnum * error_rms(i)
            write(num_str,'(i6.6)') j
            if(trim(error_property(i)) == 'K1') then
              write(45,'(aes15.6,a)') trim(num_str)//'[K1]=',(1.0+magnet_error)*value_of_attribute(ring%ele(j),'K1'),'   ! '//trim(ring%ele(j)%name)
            elseif(trim(error_property(i)) == 'A1') then  !skew gradient errors
              write(45,'(aes15.6,a)') trim(num_str)//'[A1]=',magnet_error*value_of_attribute(ring%ele(j),'K1')*value_of_attribute(ring%ele(j),'L'),'   ! '//trim(ring%ele(j)%name)
            elseif(trim(error_property(i)) == 'X_OFFSET') then
              write(45,'(aes15.6,a)') trim(num_str)//'[X_OFFSET]=',magnet_error,'   ! '//trim(ring%ele(j)%name)
            elseif(trim(error_property(i)) == 'Y_OFFSET') then
              write(45,'(aes15.6,a)') trim(num_str)//'[Y_OFFSET]=',magnet_error,'   ! '//trim(ring%ele(j)%name)
            else
              write(*,*) "Error property not yet implemented: ", trim(error_property(i))
            endif
          else
            write(*,*) "When applying errors, element does not have requested property."
            write(*,*) "Name and property: ", ring%ele(j)%name, error_mask(i)
          endif
        endif
      enddo
    endif
  enddo
  close(45)
end program
