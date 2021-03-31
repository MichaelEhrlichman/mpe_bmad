program element_list

use bmad

implicit none

character lat_file*200
type(lat_struct) lat
integer i, counter
real(rp) newang, prevang
character*8 elename
character*8 eledone(25)

call getarg(1, lat_file)
call bmad_parser(lat_file, lat)

!  !write out combined bends
!  eledone = ''
!  counter = 0
!  do i=1,lat%n_ele_track
!    elename = lat%ele(i)%name
!    call str_downcase(elename,elename)
!    if( .not. any(eledone(:) == trim(elename)) ) then
!      !if( elename(1:2) == 'vb' .and. lat%ele(i)%key == sbend$) then
!      if( lat%ele(i)%key == sbend$) then
!        write(*,'(a8,a,f10.6,a,f10.6,a,f10.6,a)') adjustl(elename)," : combined, l = ", lat%ele(i)%value(l$), ", t = ",  &
!                                                  r2d(lat%ele(i)%value(angle$)), ", k = ", lat%ele(i)%value(k1$), ","
!        write(*,'(a,f10.6,a,f10.6,a,f10.6,a)') "       t1 = ", r2d(lat%ele(i)%value(e1$)), ", t2 = ", r2d(lat%ele(i)%value(e2$)), &
!                                               ", m = ", lat%ele(i)%value(k2$)/2.0, ", ax = 10.0, ay = 10.0;"
!        counter = counter+1
!        eledone(counter) = trim(elename)
!      endif
!    endif
!  enddo
!  write(*,*) 

!write out drifts
eledone = ''
counter = 0
do i=1,lat%n_ele_track
  elename = lat%ele(i)%name
  call str_downcase(elename,elename)
  if( .not. any(eledone(:) == trim(elename)) ) then
    if( elename(1:1) == 'd' .and. lat%ele(i)%key == drift$) then
      write(*,'(a8,a,f10.6,a)') adjustl(elename), " : drift, l = ", lat%ele(i)%value(l$), ", ax = 10.0, ay = 10.0;"
      counter = counter+1
      eledone(counter) = trim(elename)
    endif
  endif
enddo
write(*,*) 

!write out quads
eledone = ''
counter = 0
do i=1,lat%n_ele_track
  elename = lat%ele(i)%name
  call str_downcase(elename,elename)
  if( .not. any(eledone(:) == trim(elename)) ) then
    if( elename(1:1) == 'q' .and. lat%ele(i)%key == quadrupole$) then
      write(*,'(a8,a,f10.6,a,f10.6,a)') adjustl(elename), " : quadrupole, l = ", lat%ele(i)%value(l$), ", k = ", lat%ele(i)%value(k1$), ", ax = 10.0, ay = 10.0;"
      counter = counter+1
      eledone(counter) = trim(elename)
    endif
  endif
enddo
write(*,*) 

!write out sextupoles
eledone = ''
counter = 0
do i=1,lat%n_ele_track
  elename = lat%ele(i)%name
  call str_downcase(elename,elename)
  if( .not. any(eledone(:) == trim(elename)) ) then
    if( elename(1:1) == 's' .and. lat%ele(i)%key == sextupole$) then
      write(*,'(a8,a,f10.6,a,f13.6,a)') adjustl(elename), " : sextupole, l = ", lat%ele(i)%value(l$), ", k = ", lat%ele(i)%value(k2$)/2.0, ","
      write(*,'(a)') "                  n =8, ax = 10.0, ay = 10.0;"
      counter = counter+1
      eledone(counter) = trim(elename)
    endif
  endif
enddo
write(*,*) 

!write out octupoles
eledone = ''
counter = 0
do i=1,lat%n_ele_track
  elename = lat%ele(i)%name
  call str_downcase(elename,elename)
  if( .not. any(eledone(:) == trim(elename)) ) then
    if( elename(1:1) == 'o' .and. lat%ele(i)%key == multipole$) then
      write(*,'(a8,a,f13.6,a)') adjustl(elename), " : multipole, n = 4, k = ", lat%ele(i)%a_pole(3)/6.0, ","
      write(*,'(a)') "                  ax = 10.0, ay = 10.0;"
      counter = counter+1
      eledone(counter) = trim(elename)
    endif
  endif
enddo
write(*,*) 

! write out LGBs
eledone = ''
prevang = 0.0
newang = 0.0
counter = 0
do i=1,lat%n_ele_track
  elename = lat%ele(i)%name
  call str_downcase(elename,elename)
  if( .not. any(eledone(:) == trim(elename)) ) then
    if( elename(1:2) == 'bn' .and. lat%ele(i)%key == sbend$) then
      newang = prevang + lat%ele(i)%value(angle$)/twopi*360.0
      write(*,'(a8,a,f10.6,a,f10.6,a,f10.6,a)') adjustl(elename)," : bending, l = ", lat%ele(i)%value(l$), ", t = ",  &
                                                r2d(lat%ele(i)%value(angle$)), ", k = ", lat%ele(i)%value(k1$), ","
      write(*,'(a,f10.6,a,f10.6,a)') "       t1 = ", r2d(lat%ele(i)%value(e1$)), ", t2 = ", r2d(lat%ele(i)%value(e2$)), ", ax = 10.0, ay = 10.0;"
      prevang = newang
      counter = counter+1
      eledone(counter) = trim(elename)
    endif
  endif
enddo
write(*,*) 

!write out combined anti-bends
eledone = ''
counter = 0
do i=1,lat%n_ele_track
  elename = lat%ele(i)%name
  call str_downcase(elename,elename)
  if( .not. any(eledone(:) == trim(elename)) ) then
    if( elename(1:2) == 'an' .and. lat%ele(i)%key == sbend$) then
      write(*,'(a8,a,f10.6,a,f10.6,a,f10.6,a)') adjustl(elename)," : combined, l = ", lat%ele(i)%value(l$), ", t = ",  &
                                                r2d(lat%ele(i)%value(angle$)), ", k = ", lat%ele(i)%value(k1$), ","
      write(*,'(a,f10.6,a,f10.6,a)') "       t1 = ", r2d(lat%ele(i)%value(e1$)), ", t2 = ", r2d(lat%ele(i)%value(e2$)), ", ax = 10.0, ay = 10.0;"
      counter = counter+1
      eledone(counter) = trim(elename)
    endif
  endif
enddo

!write out combined anti-bends
eledone = ''
counter = 0
do i=1,lat%n_ele_track
  elename = lat%ele(i)%name
  call str_downcase(elename,elename)
  if( .not. any(eledone(:) == trim(elename)) ) then
    if( elename(1:2) == 'vb' .and. lat%ele(i)%key == sbend$) then
      write(*,'(a8,a,f10.6,a,f10.6,a,f10.6,a)') adjustl(elename)," : combined, l = ", lat%ele(i)%value(l$), ", t = ",  &
                                                r2d(lat%ele(i)%value(angle$)), ", k = ", lat%ele(i)%value(k1$), ","
      write(*,'(a,f10.6,a,f10.6,a)') "       t1 = ", r2d(lat%ele(i)%value(e1$)), ", t2 = ", r2d(lat%ele(i)%value(e2$)), ", ax = 10.0, ay = 10.0;"
      counter = counter+1
      eledone(counter) = trim(elename)
    endif
  endif
enddo

contains
  function r2d(r)
    real(rp) r, r2d
    r2d = r/twopi*360.0
  end function

end program












