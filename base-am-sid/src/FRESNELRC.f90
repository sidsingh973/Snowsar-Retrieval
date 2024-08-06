! -------------------------------------------------------------------------
!
SUBROUTINE FRESNELRC(TEI,EPSI,FH,FV,NUM)
!
! -------------------------------------------------------------------------
!
!     CODE ORIGINALLY OBTAINED IN MATLAB FROM MATZLER
!
!
!   FRESNEL REFLECTION COEFFICIENTS (ASSUMING EPS'' = 0)
!     (LAYER N+1 IS THE AIR ABOVE THE SNOWPACK)
!
!   [FH,FV] = FRESNELRC(TEI,EPSR)
!       FH:   FRESNEL REFLECTION COEFFICIENT AT H POL
!       FV:   FRESNEL REFLECTION COEFFICIENT AT V POL
!       TEI:  LOCAL INCIDENCE ANGLE
!       EPSR: (REAL PART) DIELECTRIC PERMITTIVITY
!
!   VERSION HISTORY:
!      1.0    WI 15.7.95
!      2.0    MD 1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!
!   USES:
!
!
!
!   COPYRIGHT (C) 1997 BY THE INSTITUTE OF APPLIED PHYSICS,
!   UNIVERSITY OF BERN, SWITZERLAND

IMPLICIT NONE

INTEGER,INTENT(IN) :: NUM
REAL(8),INTENT(IN) :: TEI(NUM+1),EPSI(NUM+1)
REAL(8),INTENT(OUT) :: FH(NUM),FV(NUM)
INTEGER N
REAL(8) EPSN,TEIN,SINQ,QEPS,WURZ,WSUB,ND

DO N=1,NUM
    EPSN=EPSI(N)/EPSI(N+1)
    TEIN=TEI(N+1)
    SINQ=SIN(TEIN)**2d0
    QEPS=SINQ/EPSN
    WURZ=(1-QEPS)**0.5d0
    WSUB=EPSN-SINQ
    ND=EPSN**0.5d0

    FH(N)=((ND*WURZ-COS(TEIN))/(ND*WURZ+COS(TEIN)))
    FV(N)=((WURZ-ND*COS(TEIN))/(WURZ+ND*COS(TEIN)))
END DO

END SUBROUTINE FRESNELRC
