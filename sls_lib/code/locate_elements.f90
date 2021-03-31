SUBROUTINE locate_elements(lat,mask,n_elements,elements)
  USE bmad

  IMPLICIT NONE

  TYPE(lat_struct) lat
  CHARACTER(10) mask
  INTEGER n_elements
  INTEGER, ALLOCATABLE :: elements(:)

  INTEGER i

  n_elements = 0
  DO i=1,lat%n_ele_track
    IF(str_match_wild(lat%ele(i)%name, mask)) THEN
      n_elements = n_elements + 1
    ENDIF
  ENDDO
  ALLOCATE(elements(n_elements))
  n_elements = 0
  DO i=1,lat%n_ele_track
    IF(str_match_wild(lat%ele(i)%name, mask)) THEN
      n_elements = n_elements + 1
      elements(n_elements) = i
    ENDIF
  ENDDO
END SUBROUTINE locate_elements

