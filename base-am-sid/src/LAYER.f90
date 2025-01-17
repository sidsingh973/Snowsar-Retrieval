! -------------------------------------------------------------------------
!
SUBROUTINE LAYER(RI,SI,TRI,TI,TGND,TSKY,D,NUM)
!
! -------------------------------------------------------------------------
!
!     CODE ORIGINALLY OBTAINED IN MATLAB FROM MATZLER
!
!   CALCULATES THE UPWELLING BRIGHTNESS TEMPERATURES D (SEE NOTE 6)
!
!   D = LAYER(RI,SI,TI,TI,TGND,TSKY)
!       D:    UPWELLING BRIGHTNESS TEMPERATURE
!       RI:   LAYER REFLECTIVITY
!       SI:   INTERFACE REFLECTIVITY
!       TI:   LAYER TRANSMISSIVITY
!       TI:   PHYSICAL TEMPERATURE [K]
!       TGND: BRIGHTNESS TEMPERATURE OF THE SOIL BELOW THE SNOWPACK
!       TSKY: BRIGHTNESS TEMPERATURE OF THE SKY
!
!   VERSION HISTORY:
!      1.0    WI 15.7.95
!      1.1    WI 26.9.97  HANDLES ALSO THE SPECIAL CASE OF A SINGLE LAYER NOW
!      1.2    WI 02.03.99 FIXED ERROR IN 1 LAYER HANDLING
!      2.0    MD 1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!
!   USES: - INVERT_MATRIX SUBROUTINE BY MIKE,ZEROSI INTERNAL FUNCTION AND
!              EYEI INTERNAL FUNCTION BY MIKE
!
!   COPYRIGHT (C) 1997 BY THE INSTITUTE OF APPLIED PHYSICS,
!   UNIVERSITY OF BERN, SWITZERLAND

IMPLICIT NONE

INTEGER,INTENT(IN) :: NUM
REAL(8),INTENT(IN) :: RI(NUM),SI(NUM+1),TRI(NUM),TI(NUM),TGND,TSKY
REAL(8),INTENT(OUT) :: D(NUM)
REAL(8) K1
REAL(8),DIMENSION(:),ALLOCATABLE ::  EI,F,E
REAL(8),DIMENSION(:,:),ALLOCATABLE :: M1,H,EYE,M2,M3,M4,M5,INVIM1,INVIM5,EYEM1,&
EYEM5
INTEGER ROW,COL,I

!integer,intent(in) :: pixel

!PRINT *, ' '
!PRINT *, 'RI',RI,'SI',SI,'TRI',TRI,'TI',TI,'TGND',TGND,'TSKY',TSKY,'NUM',NUM


ALLOCATE( EI(NUM),M1(NUM,NUM),H(NUM-1,NUM-1),EYE(NUM,NUM),M2(NUM,NUM),&
M3(NUM,NUM),M4(NUM,NUM),E(NUM),F(NUM),M5(NUM,NUM),&
INVIM1(NUM,NUM),INVIM5(NUM,NUM),EYEM1(NUM,NUM),EYEM5(NUM,NUM) )

EI=1-RI-TRI


IF (NUM==1) THEN
    !FOR ONE LAYER CASE, THIS IS THE ENTIRE CALCULATION; THE MATRIX INVERSION
    ! ROUTINE IS NOT NEEDED.  THEREFORE, A RETURN STATEMENT IS USED HERE...
    K1=(1d0-RI(1)*SI(1))*(1d0-RI(1)*SI(2))-TRI(1)*SI(1)*TRI(1)*SI(2)
    D=TRI(1)*SI(1)*((1d0-SI(1))*RI(1)*TGND+(1d0-SI(2))*TSKY*TRI(1)+&
    EI(1)*TI(1))/K1 + (1d0-RI(1)*SI(1))*((1d0-SI(1))*TGND*TRI(1)+&
    (1d0-SI(2))*TSKY*RI(1)+EI(1)*TI(1))/K1

    !PRINT *, 'THIS IS D', D
    !PRINT *, 'TSKY', TSKY
    !PRINT *,TRI(1)*SI(1)*((1-SI(1))*RI(1)*TGND+(1-SI(2))*TSKY*TRI(1)+EI(1)*TI(1))/K1
    !PRINT *, (1-RI(1)*SI(1))*((1-SI(1))*TGND*TRI(1)+(1-SI(2))*TSKY*RI(1)+EI(1)*TI(1))/K1
    !PRINT *,'RI(1)',RI(1),'SI(1)',SI(1),'TGND',TGND,'TRI(1)',TRI(1),'SI(2)',SI(2),'EI(1)',EI(1)
    !PRINT *,'TI(1)',TI(1),'K1',K1
    !PRINT *,(1-RI(1)*SI(1)),(1-SI(1))*TGND*TRI(1),(1-SI(2))*TSKY*RI(1),EI(1)*TI(1)
    !PRINT *,(1-SI(1))*TGND*TRI(1),(1-SI(1)),TGND,TRI(1)

    RETURN
ELSE
    ! INITIALIZE ARRAYS
    M1=ZEROSI(NUM,NUM)
    H=ZEROSI(NUM-1,NUM-1)
    M2=ZEROSI(NUM,NUM)
    M3=ZEROSI(NUM,NUM)
    M4=ZEROSI(NUM,NUM)

    DO I=1,NUM
        M1(I,I)=RI(I)*SI(I)
    END DO
    DO I=1,NUM-1
        H(I,I)=TRI(I)*(1-SI(I+1))
    END DO
    M1(1:NUM-1:1,2:NUM:1)=M1(1:NUM-1:1,2:NUM:1)+H

    EYE=EYEI(NUM)
    H=ZEROSI(NUM-1,NUM-1)
    DO I=1,NUM
        M2(I,I)=TRI(I)*SI(I+1)
    END DO
    DO I=2,NUM
        H(I-1,I-1)=RI(I)*(1-SI(I))
    END DO
    M2(2:NUM:1,1:NUM-1:1)=M2(2:NUM:1,1:NUM-1:1)+H

    H=ZEROSI(NUM-1,NUM-1)
    DO I=1,NUM
        M3(I,I)=TRI(I)*SI(I)
    END DO
    DO I=1,NUM-1
        H(I,I)=RI(I)*(1-SI(I+1))
    END DO
    M3(1:NUM-1:1,2:NUM:1)=M3(1:NUM-1:1,2:NUM:1)+H

    H=ZEROSI(NUM-1,NUM-1)
    DO I=1,NUM
        M4(I,I)=RI(I)*SI(I+1)
    END DO
    DO I=2,NUM
        H(I-1,I-1)=TRI(I)*(1-SI(I))
    END DO
    M4(2:NUM:1,1:NUM-1:1)=M4(2:NUM:1,1:NUM-1:1)+H

    E=EI*TI
    E(1)=E(1)+RI(1)*(1-SI(1))*TGND
    E(NUM)=E(NUM)+TRI(NUM)*(1-SI(NUM+1))*TSKY

    F=EI*TI
    F(1)=F(1)+TRI(1)*(1-SI(1))*TGND
    F(NUM)=F(NUM)+RI(NUM)*(1-SI(NUM+1))*TSKY

END IF

EYEM1=EYE-M1
CALL INVERT_MATRIX(EYEM1,INVIM1,NUM)

M5=MATMUL(M3,MATMUL(INVIM1,M2))+M4
EYEM5=EYE-M5
CALL INVERT_MATRIX(EYEM5,INVIM5,NUM)

D=MATMUL(INVIM5,(MATMUL(MATMUL(M3,INVIM1),E)+F))

!PRINT *, 'THIS IS D',D


DEALLOCATE( EI,M1,H,EYE,M2,M3,M4,E,F,M5,INVIM1,INVIM5,EYEM1,EYEM5 )

CONTAINS

FUNCTION ZEROSI(ROW,COL)

! BY MIKE, 1 APRIL 2005
!
! THIS FUNCTION IS INTERNAL TO LAYER SUBROUTINE
! MANY VALUES IN THE ABOVE ARRAYS ARE NEVER A VALUE OTHER THAN ZERO, SO INSTEAD
! OF USING WHATEVER VALUES ARE RANDOMLY ALLOCATED TO THE ARRAYS WHEN THEY ARE
! ALLOCATED, I SPECIFICALLY SET EACH POSITION TO ZERO.  THIS FUNCTION SHOULD
! BE ENTIRELY unneceSSARY, SINCE THE DEFAULT VALUES IN THE ARRAY SHOULD BE ZERO. NONETHELESS, I DECIDED TO EXPLICITLY ZERO OUT 2-D ARRAYS BEFORE USING THEM

IMPLICIT NONE

INTEGER,INTENT(IN):: ROW,COL
REAL(8) :: ZEROSI(ROW,COL)
INTEGER I,J

DO I=1,ROW
    DO J=1,COL
        ZEROSI(I,J)=0.0d0
    END DO
END DO

END FUNCTION ZEROSI

FUNCTION EYEI(NUM)


! BY MIKE, 1 APRIL 2005
!
! THIS FUNCTION IS INTERNAL TO LAYER SUBROUTINE
! THIS FUNCTION DEFINES AN IDENTITY MATRIX OF SIZE NUM

INTEGER, INTENT(IN) :: NUM
REAL(8) :: EYEI(NUM,NUM)
INTEGER I,J
    DO I=1,NUM
        DO J=1,NUM
            IF (I==J) THEN
                EYEI(I,J)=1d0
            ELSE
                EYEI(I,J)=0d0
            END IF
        END DO
    END DO
END FUNCTION EYEI

END SUBROUTINE LAYER
