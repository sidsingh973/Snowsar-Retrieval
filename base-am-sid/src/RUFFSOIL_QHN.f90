! -------------------------------------------------------------------------
!
SUBROUTINE RUFFSOIL_QHN(F,MV,TK,KSIGMA,GND_SIG,THETA,EPS_TOP,&
RHO, SAND, SILT, CLAY,EPS_SOIL,R_H_MOD,R_V_MOD,RS_H_MOD,RS_V_MOD)
!
! -------------------------------------------------------------------------

!     F = FREQUENCY [GHZ]
!     MV = VOLUMETRIC MOISTURE [FRAC]
!     TK = TEMP [K]
!     KSIGMA = NORMALIZED SURFACE SDEV [m]
!     GND_SIG=Gound roughness [m]
!     THETA = NADIR ANGLE [DEG]
!     RHO = SOIL BULK DENSITY [KG/M3]
!     EPS_TOP = EPSILON OF OVERLYING MEDIUM (Complex)
!     EPS_SOIL = DEFINE EPSILON EXPLICITLY (Complex)
!
!     THIS IS THE QHN MODEL REFITTED BY BENOIT
!     BASED ON FREQUENCY-INDEPEDENT PARAMETERS

IMPLICIT NONE
REAL(8),INTENT(IN) :: F,MV,TK,KSIGMA,THETA
REAL(8),INTENT(IN) :: RHO,SAND,SILT,CLAY
REAL(8),INTENT(IN) ::  GND_SIG
COMPLEX(8),INTENT(IN) :: EPS_TOP
REAL(8),INTENT(OUT) :: R_H_MOD,R_V_MOD
REAL(8),INTENT(OUT) :: RS_H_MOD,RS_V_MOD


REAL(8),PARAMETER :: PI=3.14159d0, MJU0=1.2566d-006, EPS0=8.8542d-012
REAL(8) :: TD,A1,A2,A3,QR,NRV,NRH,HR
REAL(8) :: RMS_G_MM,THETA_R
REAL(8) :: FRESNEL_H,FRESNEL_V
Real(8) :: r_v_mod0, r_h_mod0
COMPLEX(8) :: EPS_SOIL,EPS_EFF


TD=TK-273.15d0 !CONVERT DEGK TO DEGC


If(.False.)Then
    IF(F.LT.15d0)THEN
        CALL EPSSOIL_MBSDM_MIRONOV(MV,TD,F,RHO,SAND,SILT,CLAY,EPS_SOIL)
    ELSE
        CALL EPSSOIL_WS80_TUNED_CALVET(MV,TD,F,RHO,SAND,SILT,CLAY,EPS_SOIL)
    END IF
else
    !CALL EPSSOIL_DOBSON(MV,TD,F,RHO,SAND,SILT,CLAY,EPS_SOIL)
    Call Epssoil_Pan(MV,TD,F,RHO,SAND,SILT,CLAY,EPS_SOIL)
end if



!EPS_SOIL=complex(real(EPS_SOIL),0)
EPS_EFF=EPS_SOIL/EPS_TOP
!Print *,'EPS_SOIL',EPS_SOIL
!print *,'EPS_TOP',EPS_TOP
!print *,'EPS_EFF',EPS_EFF


!CALCULATE THE SPECULAR PART OF REFLECTIVITY
CALL GAMMAH(EPS_EFF,THETA,FRESNEL_H)
CALL GAMMAV(EPS_EFF,THETA,FRESNEL_V)

THETA_R=THETA/180.0d0*PI
!Original, rough soil surface models
!RS_H_MOD=FRESNEL_H*EXP(-(2.0*KSIGMA*COS(THETA_R))**2.0)
!RS_V_MOD=FRESNEL_V*EXP(-(2.0*KSIGMA*COS(THETA_R))**2.0)

!Print *,'Print Ruffsoil_QHN!'
!Print *,'FRESNEL_H',FRESNEL_H,'KSIGMA',KSIGMA,'THETA_R',THETA_R
!print *,'RS_H_MOD',RS_H_MOD
!Print *,'FRESNEL_V',FRESNEL_V,'KSIGMA',KSIGMA,'THETA_R',THETA_R
!print *,'RS_V_MOD',RS_V_MOD


!Calculate the total reflectivity using QHN MODEL
a1=0.887d0
a2=0.796d0
a3=3.517d0
QR=0.075d0
NRv=1.503d0
NRh=0.131d0

RMS_G_MM = KSIGMA/REAL(2d0*PI*F*1.0d9*(MJU0*EPS0*EPS_TOP)**0.5d0)*1000.0d0

HR= ((a1* RMS_G_MM)/(a2 * RMS_G_MM + a3))**6.0d0
r_v_mod0 = (1-QR) * fresnel_v + QR * fresnel_h;
r_v_mod = r_v_mod0 * exp(-HR * cos(THETA_R)**NRv )

r_h_mod0 = (1-QR) * fresnel_h + QR * fresnel_v;
r_h_mod = r_h_mod0 * exp(-HR * cos(THETA_R)**NRh )

rs_v_mod = r_v_mod0 * EXP(-(2.0d0*KSIGMA*COS(THETA_R))**2.0d0)
rs_h_mod = r_h_mod0 * EXP(-(2.0d0*KSIGMA*COS(THETA_R))**2.0d0)

If(.False.)Then
    Print *,'QR&HR&NRv',QR,HR,NRv
    Print *, 'fresnel_v',fresnel_v
    Print *, 'RMS_G_MM',RMS_G_MM
    Print *, 'THETA_R',THETA_R
    Print *,'Eps_Soil',EPS_SOIL
    Print *,'F',F
    Print *,'rv',r_v_mod
    Print *,'rh',r_h_mod
EndIf


END SUBROUTINE RUFFSOIL_QHN
