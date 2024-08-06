! -------------------------------------------------------------------------
!
SUBROUTINE EPSW(F,T,EW)
!
! -------------------------------------------------------------------------
!
! REVISED CODE FROM MEMLS3&A, WRITTEN BY JINMEI
! INPUT:
! F - Frequency [GHz]
! T - Temperature [K]

IMPLICIT NONE

REAL(8),INTENT(IN) :: F,T
COMPLEX(8),INTENT(OUT) :: EW
REAL(8) TETA,TK,e0,e1,e2,f1,f2,fGHz
COMPLEX(8) :: i0


TK=T
fGHz=F;

TETA=1.0d0-300d0/TK;
e0=77.66d0-103.3d0*TETA;
e1=0.0671d0*e0;
f1=20.2d0+146.4d0*TETA+316.0d0*TETA*TETA;
e2=3.52d0+7.52d0*TETA;

f2=39.8d0*f1;
i0=dcmplx(0.0d0,1.0d0)
EW=e2+(e1-e2)/(1.0-i0*fGHz/f2)+(e0-e1)/(1.0-i0*fGHz/f1);



END SUBROUTINE EPSW
