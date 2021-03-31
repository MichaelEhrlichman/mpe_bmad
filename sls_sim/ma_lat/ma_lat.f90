program ma_lat

use bmad
use bmad_parser_mod, only: bp_com
use sls_lib
use sim_utils_interface
use random_mod
!use ifport

implicit none

type(lat_struct) lat
type(coord_struct), allocatable :: co(:)
type(normal_modes_struct) mode

real(rp), allocatable :: results(:,:)
real(rp) avg_emit_x, avg_emit_y, avg_sigmaE_E, avg_sigma_z
real(rp) rms_etay, avg_rms_etay

character(100) in_file
character(200) lat_file
character(4) seed_str
character(6) num_str
character(30) dir

integer radcache
integer attribute_ix
integer i,j,k
integer n_seeds
integer dirstatus

logical err
logical ma_cor

integer, parameter :: MAX_SETS = 40
character(10) mask(MAX_SETS)
character(10) property(MAX_SETS)
character(3) abs_or_rel(MAX_SETS)
real(rp) rms(MAX_SETS)
real(rp) cutoff(MAX_SETS)

integer ran_seed
real(rp) rnum
real(rp) ma

namelist /parameters/  lat_file,&
                       ran_seed,&
                       ma_cor, &  !misalign correctors?
                       mask,&
                       property,&
                       rms,&
                       abs_or_rel,&
                       cutoff,&
                       n_seeds

mask = 'null'
ran_seed = 0
abs_or_rel = 'abs'

call getarg(1,in_file)
open(10,file=in_file,status='old')
read(10,nml=parameters)
close(10)

bp_com%always_parse = .true.
call bmad_parser(lat_file, lat)

call ran_seed_put(ran_seed)

do k=1, n_seeds
  write(seed_str,'(i4.4)') k
  dir = 'seed_'//seed_str
  !dirstatus = makedirqq(trim(dir))
  call system('mkdir -p ' // adjustl(trim( dir ) ) )
  open(45,file=trim(dir)//'/seed_results_'//seed_str//'.out')
  write(45,'(a)') 'call, file = ../'//lat_file
  write(45,'(a)') 'expand_lattice'
  do i=1,MAX_SETS
    if(mask(i) .ne. 'null') then
      j = 0
      do while(j.lt.lat%n_ele_max)
        j=j+1
        if(str_match_wild(lat%ele(j)%name,mask(i))) then
          attribute_ix = attribute_index(lat%ele(j),property(i))
          if( attribute_ix .gt. 0 ) then
            do while(.true.)
              call ran_gauss(rnum) 
              if(abs(rnum) .lt. cutoff(i)) exit
            enddo
            ma = rnum * rms(i)
            write(num_str,'(i6.6)') j
            if(trim(property(i)) == 'K1') then
              if(abs_or_rel(i) .eq. 'abs') then
                write(45,'(a,es14.5,a)') num_str//'['//trim(property(i))//']=',ma+lat%ele(j)%value(k1$),'   ! '//trim(lat%ele(j)%name)
              elseif(abs_or_rel(i) .eq. 'rel') then
                write(45,'(a,es14.5,a)') num_str//'['//trim(property(i))//']=',(1.0+ma)*lat%ele(j)%value(k1$),'   ! '//trim(lat%ele(j)%name)
              else
                write(*,*) "Inconsistent definition in .in file.  Terminating."
                stop
              endif
            elseif(trim(property(i)) == 'A1') then  !skew gradient errors
              write(45,'(a,es14.5,a)') num_str//'['//trim(property(i))//']=',ma*lat%ele(j)%value(k1$)*lat%ele(j)%value(l$),'   ! '//trim(lat%ele(j)%name)
            else
              write(45,'(a,es14.5,a)') num_str//'['//trim(property(i))//']=',ma,'   ! '//trim(lat%ele(j)%name)
            endif
            if( (lat%ele(j+1)%key==hkicker$ .and. lat%ele(j+2)%key==vkicker$) .or. &
                (lat%ele(j+1)%key==vkicker$ .and. lat%ele(j+2)%key==hkicker$) ) then

              if(trim(lat%ele(j+3)%name)==trim(lat%ele(j)%name)) then
                !we have two sextupole halves with a steering corrector in the middle
                if(ma_cor) then
                  write(num_str,'(i6.6)') j+1
                  write(45,'(a,es14.5,a)') num_str//'['//trim(property(i))//']=',ma,'   ! '//trim(lat%ele(j+1)%name)
                  write(num_str,'(i6.6)') j+2
                  write(45,'(a,es14.5,a)') num_str//'['//trim(property(i))//']=',ma,'   ! '//trim(lat%ele(j+2)%name)
                endif
                write(num_str,'(i6.6)') j+3
                write(45,'(a,es14.5,a)') num_str//'['//trim(property(i))//']=',ma,'   ! '//trim(lat%ele(j+3)%name)
                j=j+3
              endif
            endif
          else
            write(*,*) "bad: attribute ", property(i), " not found for element ", lat%ele(j)%name
          endif
        endif
      enddo
    endif
  enddo
  close(45)
enddo

close(45)

end program ma_lat



