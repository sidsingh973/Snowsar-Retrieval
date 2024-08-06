! -------------------------------------------------------------------------
!
SUBROUTINE EPSSOIL_DOBSON(MV0,TD,F,RHOS_KG,SAND,SILT,CLAY,EPS_SOIL)
!
! -------------------------------------------------------------------------
!
!     CODE ORIGINALLY OBTAINED IN MATLAB FROM PULLIAINEN
!
!     FUNCTION FOR CALCULATING EPSILON FOR SOIL USING FREQUENCY,
!     TEMPERATURE AND VOLUMETRIC SOIL MOSITURE. USES EPSW.M FOR
!     DIELECTRICITY OF WATER.
!
!     BY J. PULLIAINEN (MOD. BY K. TIGERSTEDT)
!
!     MV - SOIL VOLUMETRIC MOISTURE[FRAC]
!     TD - SOIL TEMPERATURE [degC]
!     F - FREQUENCY [GHZ]
!     RHOS_KG - SOIL BULK DENSITY [KG/M3]
!     SAND - SAND CONTENT[%]
!     SILT - SILT CONTENT[%]
!     CLAY - CLAY CONTENT[%]
!
!     NOTE: SOME FINNISH COMMENTS WERE NOT COPIED IN ENTIRETY -MD
!   VERSION HISTORY:
!      1.0    JP ?.?.?
!      1.1    KT ?.?.?
!      2.0    MD 1 APR 05 TRANSLATED TO FORTRAN FROM MATLAB
!      3.0    JINMEI, REVISED ACCORDING TO THE DOBSON MODEL

IMPLICIT NONE

REAL(8), INTENT(IN) :: MV0,TD,F,RHOS_KG,SAND,SILT,CLAY
COMPLEX(8),INTENT(OUT) :: EPS_SOIL

INTEGER OPT_CONSIDER_ICE
REAL(8) RHOS,MV,TK
REAL(8) SIGMAE
REAL(8) BETA_R, BETA_I
REAL(8) ER_W,EI_W,ER_ICE,EI_ICE,ER_SOIL, EI_SOIL
Real(8) S,A,B
Real(8) er_dsoil, rho_ice, alpha, mvu, mvi, Temp


!1. parameters
er_dsoil=4.7d0
alpha=0.65d0
rho_ice=0.917d0


!2. process data
RHOS=RHOS_KG/1000.0d0
TK=TD+273.15d0


IF(MV0==0.0d0)Then
    !July,2015, revised to prevent Inf values
    MV=0.00001d0
Else
    MV=MV0
End if


!3. treat frozen soil
! OPT_CONSIDER_ICE=0 means, MV is the equivalent unfrozen water content (used by passive Tb MCMC before)
! OPT_CONSIDER_ICE=1 means, the water is totaly frozen
! OPT_CONSIDER_ICE=2 means, use the revised Dobson model by Zhang et al.
OPT_CONSIDER_ICE=2

Select Case(OPT_CONSIDER_ICE)
    Case(0)
        mvu=mv
        mvi=0.0d0

    Case(1)
        If(TD.lt.0.0d0)Then
            mvi=mv/rho_ice
            mv=0.0d0
        Else
            mvu=mv
            mvi=0.0d0
        End if

    Case(2)
        If(TD.lt.0.0d0)Then
            S=0.042d0 +4.23d0*clay/100.0d0 +1.12d0*silt/100.0d0 -1.16d0*sand/100.0d0
            A=0.2618d0+0.5519d0*log(S)
            B=-0.264d0*log(S)+0.3711d0
            mvu=exp(A)*(abs(TD)**(-exp(B))) *RHOS/100.0d0

            If(S.lt.0.0d0)Then
                mvu=0.0d0
            End if

            If (mvu.lt.0.0d0)Then
                mvu=0.0d0
            End if

            If(mvu.lt.mv)Then
                mvi=(mv-mvu)/rho_ice
                mv=mvu
            Else
                mvu=mv
                mvi=0.0d0
            End if

        Else
            mvu=mv
            mvi=0.0d0
        EndIf
End Select



!4. Begin calculation
CALL EPSW(F,TK,ER_W,EI_W)

IF (F>=1.4) THEN
    SIGMAE=-1.645d0 +1.9390d0*RHOS -0.0225622d0*SAND +0.01594d0*CLAY
ELSE
    SIGMAE=0.0467d0 +0.2204d0*RHOS -0.004111d0*SAND +0.006614d0*CLAY
ENDIF
EI_W=EI_W+SIGMAE*(1-RHOS/2.65d0)/(0.05563132d0*F*MV)



!MATZLER AND WEGMULLER 1987
ER_ICE=3.1884d0+9.1d-4*TD
CALL EPSIICE(TK,F,EI_ICE,1)


BETA_R=1.2748 - 0.00519 * SAND - 0.00152 * CLAY
BETA_I=1.33797 - 0.00603 * SAND - 0.00166 * CLAY


Temp=1.0 + RHOS/2.65*(er_dsoil**alpha-1.0) + MV**BETA_R*(ER_W**alpha) -MV &
        + mvi*(er_ice**alpha)- mvi
ER_SOIL=Temp**(1.0/alpha)

IF(F<1.4)THEN
    ER_SOIL=ER_SOIL*1.15-0.68
ENDIF

Temp=MV**(BETA_I)*EI_W + MVi**(BETA_I)*EI_ICE
EI_SOIL=Temp**(1.0/alpha)


EPS_SOIL=complex(ER_SOIL,-EI_SOIL)


END SUBROUTINE EPSSOIL_DOBSON
