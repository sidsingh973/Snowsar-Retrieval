! -------------------------------------------------------------------------
!
SUBROUTINE MIXMOD(FREQ,TI,WIFR,EPSI,EPSII,NUM)
!
! -------------------------------------------------------------------------
!
!     CODE ORIGINALLY OBTAINED IN MATLAB FROM MATZLER
!
!   CALCULATES THE PERMITTIVITY FOR WETNESS > 0
!      PHYSICAL MIXING MODEL WEISE 97 AFTER MATZLER 1987 (CORRECTED)
!      WATER TEMPERATURE IS ASSUMED CONSTANT AT 273.15 K
!
!   [EPSI,EPSII] = MIXMOD(F,TI,WI,EPSI,EPSII)
!       EPSI:  REAL PART OF THE PERMITTIVITY
!       EPSII: IMAGINARY PART OF THE PERMITTIVITY
!       F:     FREQUENCY [GHZ]
!       TI:    PHYSICAL SNOW TEMPERATURE [K]
!       WI:    WETNESS [FRAC]
!       EPSI:  REAL PART OF DRY SNOW PERM.
!       EPSII: IMAGINARY PART OF DRY SNOW PERM.
!
!   VERSION HISTORY:
!      1.0    WI 15.7.95
!      2.0    MD  1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!      2.1    MD 21 NOV 05 MADE ALL LOCALS DYNAMICALLY ALLOCATABLE
!
!   USES: - NONE
!
!   COPYRIGHT (C) 1997 BY THE INSTITUTE OF APPLIED PHYSICS,
!   UNIVERSITY OF BERN, SWITZERLAND

IMPLICIT NONE

INTEGER, INTENT(IN) :: NUM
REAL(8), INTENT(IN) :: FREQ, TI(NUM), WIFR(NUM)
REAL(8), INTENT(INOUT) :: EPSI(NUM), EPSII(NUM)

REAL(8) :: AA,AB,Wi
COMPLEX(8) :: EW,epsd,Ka,Kb,K,epsz,epsn,eps
COMPLEX(8) :: i0
INTEGER :: i

AA=0.005d0
AB=0.4975d0

i0=dcmplx(0.0d0,1.0d0)

DO i=1,NUM
    Call EPSW(FREQ,TI(i),ew)
    epsd=EPSI(i)+EPSII(i)*i0
    Wi=WIFR(i)

    Ka=epsd/(epsd+AA*(ew-epsd))
    Kb=epsd/(epsd+AB*(ew-epsd))
    K =(Ka+2.0d0*Kb)/3.0d0
    epsz=(1.0d0-Wi)*epsd+Wi*ew*K
    epsn=1.0d0-Wi*(1.0d0-K)
    eps=epsz/epsn   !Maxwell-Garnett Mixing of water in dry snow

    EPSI(i)=real(eps)
    EPSII(i)=imag(eps)
ENDDO

END SUBROUTINE MIXMOD
