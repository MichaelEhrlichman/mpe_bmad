!+
!-
subroutine binary_search(lost,delta_m,accuracy,reset) !updates delta_m according to one iteration of binary search
  use bmad

  logical lost
  real(rp) delta_m
  real(rp) accuracy
  logical reset

  logical, save :: high_found
  logical, save :: low_found
  logical, save :: bracket_found
  real(rp), save :: last_high
  real(rp), save :: last_low

  if(reset) then
    high_found = .false.
    low_found = .false.
    bracket_found = .false.
    accuracy = 10.0
  else
    ! if an upper and lower bound for the aperture have have not yet been found, double or half delta_m
    if(.not. bracket_found) then
      if(lost) then
        high_found = .true.
        last_high = delta_m
        if(.not. low_found) then
          delta_m = delta_m / 2.0_rp
        endif
      else
        low_found = .true.
        last_low = delta_m
        if(.not. high_found) then
          delta_m = delta_m * 2.0_rp
        endif
      endif
      bracket_found = (high_found .and. low_found)
      accuracy = 10.0
    endif

    ! if and upper and lower bound have been found, move half the distance
    if(bracket_found) then
      if(lost) then
        ! decrease delta_m
        last_high = delta_m
        delta_m = delta_m - abs(delta_m-last_low)/2.0_rp
      else
        ! increase delta_m
        last_low = delta_m
        delta_m = delta_m + abs(last_high-delta_m)/2.0_rp
      endif
      accuracy = abs(last_high-last_low)/min(abs(last_high),abs(last_low))
    endif

    ! check for extreme values.  if the aperture is larger then 100%, record 100% and move on.
    if(abs(delta_m) .gt. 2.0) then
      accuracy = 0.0
      !the next two lines result in 1.0 being recorded as the aperture for this location.
      last_high = 1.1_rp
      last_low = 0.9_rp
    endif
  endif
end subroutine binary_search

