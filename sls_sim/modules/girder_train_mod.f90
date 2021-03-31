MODULE girder_train_mod

implicit none

CONTAINS
!+
! Subroutine girder_train(girders,survey,ma)
!
! Input:
!   girders(:)      - INTEGER, consecutive element indexes of girder train.
!   survey(:)       - REAL(rp), offset and angle for
!                     beginning and end of train.  relative offset and angle between
!                     girders.
! Output:
!   ma(:)           - REAL(rp), offset and angle for each 
!                     girder and train relative to lab frame.
!-
SUBROUTINE girder_train_convert(survey,ma)

USE bmad
USE f95_lapack

IMPLICIT NONE

REAL(rp) survey(:)
REAL(rp) ma(:)
REAL(rp), allocatable :: lapack_vec(:,:)
INTEGER n_meas
INTEGER i

INTEGER M, N, lwork, info
INTEGER, allocatable :: lapack_working_space(:)

REAL(rp), allocatable :: mat(:,:)

n_meas = SIZE(survey)

 allocate(lapack_vec(n_meas,1))
 allocate(lapack_working_space(2*n_meas*SIZE(ma)))
 allocate(mat(n_meas,SIZE(ma)))
 
 IF( SIZE(ma) .NE. SIZE(survey)-2 ) THEN
   WRITE(*,*) "FAILING ERROR: size of ma must be exactly two less than the size of survey."
   WRITE(*,*) "STOPPING PROGRAM"
   STOP
 ENDIF
 
 ! make linear system
 mat = 0
 mat(1,1) = 1
 mat(2,1:2) = (/-1,1/)
 DO i=3,n_meas-3,2
   !mat( i   , i-1:i ) = (/-1,1/)
   !mat( i+1 , i-2:i+1 ) = (/1,-1,-1,1/)
   mat( i   , i-1 ) = -1
   mat( i   , i ) = 1
   mat( i+1 , i-2 ) = 1
   mat( i+1 , i-1 ) = -1
   mat( i+1 , i) = -1
   mat( i+1 , i+1 ) = 1
 ENDDO
 mat( n_meas-1 , n_meas-2 ) = 1
 mat( n_meas, n_meas-3:n_meas-2 ) = (/-1,1/)

 DO i=1,n_meas
   WRITE(*,'(8F4.0)') mat(i,:)
 ENDDO

 M = n_meas
 N = n_meas-2
 lwork = SIZE(lapack_working_space)
 lapack_vec = 0
 lapack_vec(:,1) = survey
 CALL dgels('N', M, N, 1, mat, M, lapack_vec, M, lapack_working_space, lwork, info)
 ma = lapack_vec(1:N,1)

END SUBROUTINE girder_train_convert

END MODULE girder_train_mod









