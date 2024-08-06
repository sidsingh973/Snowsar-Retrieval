! -------------------------------------------------------------------------
!
SUBROUTINE BORN2(K0,VFI,PCIMM,EPSI,EICE,GB6,GC6,GF6,GS6,NUM)
!
! -------------------------------------------------------------------------
!
!   CODE ORIGINALLY OBTAINED IN MATLAB FROM MATZLER, MODIFIED BY MIKE... SEE
!     VERSION HISTORY
!
!   CALCULATES THE SCATTERING COEFFICIENT USING BORN APPROXIMATION
!
!   [GB6,GC6,GF6,GS6] = BORNA(K0,VFI,PCI,EPSI,EICE,EPSEFF,KP)
!       GB6: 6-FLUX BACK SCATTERING COEFFICIENT
!       GC6: 6-FLUX CROSS SCATTERING COEFFICIENT
!       GF6: 6-FLUX FORWARD SCATTERING COEFFICIENT
!       GS6: 6-FLUX SCATTERING COEFFICIENT
!       K0:   WAVE NUMBER. K0=FREQ*2.*PI/0.299793.
!       VFI: VOLUME FRACTION OF ICE [FRAC]
!       PCIMM: CORRELATION LENGTH [mm]
!       EPSI: DIELECTRIC CONSTANT OF SNOW
!       EICE: DIELECTRIC CONSTANT OF SNOW
!       EPSEFF: EFFECTIVE PERMITTIVITY
!
!   VERSION HISTORY:
!      1.0     WI 27.05.98
!      2.0     WI 31.03.05.  TRANSLATED TO FORTRAN BY MIKE.  BORNA AND
!                             BORNSNK COMBINED INTO ONE.  CALL TO POLDER
!                             IS USED TO COMPUTE EFFECTIVE PERMITTIVITY
!      3.0    MD  1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!      3.1    MD 22 DEC 11 Added Unity,NegUnity variables to prevent compiler warnings
!
!   USES:
!
!   COPYRIGHT (C) 1998 BY THE INSTITUTE OF APPLIED PHYSICS,
!   UNIVERSITY OF BERN, SWITZERLAND

IMPLICIT NONE

INTEGER, INTENT(IN) :: NUM
REAL(8), INTENT(IN) :: K0,VFI(NUM),PCIMM(NUM)
REAL(8), INTENT(IN) :: EPSI(NUM),EICE(NUM)
REAL(8), INTENT(OUT) :: GB6(NUM),GC6(NUM),GF6(NUM),GS6(NUM)
REAL(8),DIMENSION(:),ALLOCATABLE :: PCI,A,A3,EA,EA3,K1,K3,KP
Real(8),Dimension(:),Allocatable :: MUC,AA,XX,EPSEFF,X2
Real(8),Dimension(:),Allocatable :: BB,BT,BF,BTOT
INTEGER STEPS
!Real,Parameter,Dimension(1) :: Unity=1.,NegUnity=-1.
!Integer :: ARG_LENGTH(2)

!New Variables Jinmei added
Integer I,J,K,STEPS3
Real(8) :: maxi, mini, maxo, mino, dmu, delta
Real(8) :: sii
Real(8) :: PI
Real(8),Dimension(:),Allocatable ::mui, muo_b, muo_t, muo_f, phi
Real(8),Dimension(:,:),Allocatable :: Grid_muo, Grid_phi
Real(8),Dimension(:,:),Allocatable :: sio, cofi, si2fi, coste, si2chi, func_phi
Real(8),Dimension(:),Allocatable :: inte_phi,inte_muo_b,inte_muo_t,inte_muo_f


ALLOCATE( PCI(NUM),A(NUM),A3(NUM),EA(NUM),EA3(NUM),K1(NUM),K3(NUM),KP(NUM),&
MUC(NUM),AA(NUM),XX(NUM),BB(NUM),BT(NUM),BF(NUM),BTOT(NUM),EPSEFF(NUM),&
X2(NUM))




! 0)  CONSTANTS AND CONVERSION
PI=3.14159265358979d0
STEPS=11
PCI=PCIMM*0.001d0

! 1)  COMPUTE FIELD FACTOR AND DEPOLARIZATION RATIO, FROM BORNSNK
!     'A' AFTER NOTE 10, MATZLER 1997, COMMENT AND CALL FROM BORNSNK
CALL SNOWAO(VFI,A,NUM)

! 2)  COMPUTE EFFECTIVE PERMITTIVITY USING NEW FUNCTION POLDER.M
!CALL POLDER(VFI,A,EICE,EPSI,EPSEFF,NUM)  !Jinmei commented on Dec22,2015, because MEMLS3 code does not contain this
EPSEFF=EPSI

! 3)  COMPUTE KP, THE SQUARED RATIO BETWEEN INTERNAL / EXTERNAL FIELDS
!     THIS CODE ORIGINALLY IN BORNSNK, STARTING AT LINE 41
A3=1-2d0*A
EA=EPSEFF*(1d0-A)+A
EA3=EPSEFF*(1d0-A3)+A3
K1=(EA/(EA+A*(EICE-1d0)))**2d0
K3=(EA3/(EA3+A3*(EICE-1)))**2d0
KP=(2d0*K1+K3)/3d0

! 4)  COMPUTE SCATTERING COEFFICIENTS, FROM BORNA, STARTING AT LINE 34

!NOTE: THE NEXT LINE CONSISTENT WITH MATZLER AND WIESMANN 99 PAPER EQUATION (7).
!  BASED ON CORRESPONDENCE WITH MATZLER, EPSEFF=N^2

MUC=((EPSEFF-1)/EPSEFF)**0.5d0

AA=2d0*(PCI*K0)**3d0*K0*VFI*(1-VFI)*(EICE-1)**2d0*KP
XX=PCI*K0*EPSI**0.5d0
X2 = 2.0d0*XX**2.0d0


! Jinmei Revised for Borna2.f90
! 5) Calculate integration
Do i=1,NUM

    ! Calculate integration steps of thetai
    STEPS3=STEPS*3
    Allocate(mui(STEPS),muo_b(STEPS), muo_t(STEPS), muo_f(STEPS),phi(STEPS))
    Allocate(Grid_phi(STEPS3,STEPS),Grid_muo(STEPS3,STEPS))
    Allocate(inte_muo_b(STEPS),inte_muo_t(STEPS),inte_muo_f(STEPS))

    mini=MUC(i)
    maxi=1.0d0
    dmu = maxi - mini
    delta = dmu/real(STEPS)
    Do j=1,STEPS
        mui(j)=mini+delta*0.5d0+(j-1)*delta
    End Do

    mino=-1.0d0
    maxo=-1.0d0*MUC(i)
    dmu = maxo - mino
    delta = dmu/real(STEPS)
    Do j=1,STEPS
        muo_b(j)=mino+delta*0.5d0+(j-1)*delta
    End Do


    mino=-1.0d0*MUC(i)
    maxo=MUC(i)
    dmu = maxo - mino
    delta = dmu/real(STEPS)
    Do j=1,STEPS
        muo_t(j)=mino+delta*0.5d0+(j-1)*delta
    End Do


    mino=MUC(i)
    maxo=1.0d0
    dmu = maxo - mino
    delta = dmu/real(STEPS)
    Do j=1,STEPS
        muo_f(j)=mino+delta*0.5d0+(j-1)*delta
    End Do

    delta=PI/real(STEPS)
    Do j=1,STEPS
        phi(j)=delta*0.5d0+(j-1)*delta
    End Do


    Do j=1,STEPS
        Grid_phi(1:STEPS3,j)=phi(j)
    End Do

    Do j=1,STEPS
        Grid_muo(j,1:STEPS)=muo_b(j)
        Grid_muo(STEPS+j,1:STEPS)=muo_t(j)
        Grid_muo(STEPS*2+j,1:STEPS)=muo_f(j)
    End Do

    ! 5) Calculate integration
    Do j=1,STEPS

        Allocate(sio(STEPS3,STEPS),cofi(STEPS3,STEPS),si2fi(STEPS3,STEPS),&
            coste(STEPS3,STEPS),si2chi(STEPS3,STEPS),func_phi(STEPS3,STEPS))

        sii=sqrt(1.0d0-mui(j)**2.0d0)
        sio=sqrt(1.0d0-Grid_muo**2.0d0)
        cofi=cos(Grid_phi)
        si2fi=1d0-cofi**2.0d0
        coste=mui(j)*Grid_muo + sii * sio * cofi
        si2chi=0.5d0*(1.0d0+coste**2.0d0)
        func_phi=si2chi / (1.0d0+(1.0d0-coste)*X2(i))**2.0d0

        !integration over phi
        Allocate(inte_phi(STEPS3))
        Do k=1,STEPS3
            inte_phi(k)=sum(func_phi(k,1:STEPS))/Real(STEPS)
        End Do

        !integration over muo
        inte_muo_b(j)=sum(inte_phi(1:STEPS)) *(1.0d0-MUC(i))/STEPS/2.0d0
        inte_muo_t(j)=sum(inte_phi(STEPS+1:STEPS*2)) * (2.0d0*MUC(i))/STEPS/2.0d0
        inte_muo_f(j)=sum(inte_phi(STEPS*2+1:STEPS*3)) * (1.0d0-MUC(i))/STEPS/2.0d0

        Deallocate(sio,cofi,si2fi,coste,si2chi,func_phi)
        Deallocate(inte_phi)
    End Do

    !integration over mui
    BB(i)=sum(inte_muo_b)/STEPS
    BT(i)=sum(inte_muo_t)/STEPS
    BF(i)=sum(inte_muo_f)/STEPS
    Deallocate(mui,muo_b, muo_t, muo_f,phi)
    Deallocate(Grid_phi,Grid_muo)
    Deallocate(inte_muo_b,inte_muo_t,inte_muo_f)
End Do


! 6)  CALCULATION OF SCATTERING COEFFICIENTS
BTOT=BB+BT+BF
GB6=AA*BB
GC6=0.25d0*AA*BT
GF6=AA*BF
GS6=AA*BTOT


DEALLOCATE(PCI,A,A3,EA,EA3,K1,K3,KP,MUC,AA,XX,BB,BT,BF,BTOT,EPSEFF,X2)

END SUBROUTINE BORN2
