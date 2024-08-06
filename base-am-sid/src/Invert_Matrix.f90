! -------------------------------------------------------------------------
!
SUBROUTINE INVERT_MATRIX(MA,INV,N)
!
! -------------------------------------------------------------------------
!     
!     MATRIX INVERSION ALGORITHM OBTAINED FROM GOTOP FORTRAN 90 TEXT,
!     ISBN:957-566-172-9
!
!     COPIED FROM THE COMPANION CD TO THAT TEXT BY MIKE, 1 APRIL 2005

IMPLICIT NONE

INTEGER,INTENT(IN) :: N
REAL(8),INTENT(IN) :: MA(N,N)
REAL(8),INTENT(OUT) :: INV(N,N)
REAL(8),dimension(:,:),allocatable :: temp
INTEGER I,J
 
allocate(temp(n,n))

DO I=1,N
  DO J=1,N
    TEMP(I,J)=MA(I,J)
    INV(I,J)=0.0d0
  END DO
  INV(I,I)=1.0d0
END DO

CALL UPPER(TEMP,INV,N)
CALL LOWER(TEMP,INV,N)

DO I=1,N
  DO J=1,N
    INV(I,J)=INV(I,J)/TEMP(I,I)
  END DO
END DO

deallocate(temp)

CONTAINS

SUBROUTINE UPPER(M,S,N)
INTEGER,INTENT(IN):: N
INTEGER I,J,K
REAL(8) E
REAL(8),INTENT(INOUT):: M(N,N)
REAL(8),INTENT(INOUT):: S(N,N)

DO I=1,N-1
  DO J=I+1,N            
    E=M(J,I)/M(I,I)
    DO K=1,N
      M(J,K)=M(J,K)-M(I,K)*E
      S(J,K)=S(J,K)-S(I,K)*E
    END DO
  END DO
END DO

END SUBROUTINE UPPER

SUBROUTINE LOWER(M,S,N)
INTEGER,INTENT(IN):: N
REAL(8),INTENT(INOUT):: M(N,N)
REAL(8),INTENT(INOUT):: S(N,N)
INTEGER I,J,K
REAL(8) E

DO I=N,2,-1
  DO J=I-1,1,-1         
    E=M(J,I)/M(I,I)
    DO K=1,N
      M(J,K)=M(J,K)-M(I,K)*E
      S(J,K)=S(J,K)-S(I,K)*E 
    END DO
  END DO
END DO

END SUBROUTINE LOWER

END SUBROUTINE INVERT_MATRIX
