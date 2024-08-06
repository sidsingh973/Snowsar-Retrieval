SUBROUTINE EPSSOIL_MBSDM_MIRONOV(MV,TD,FGHZ,RHOKG,SAND,SILT,CLAY,EPS_SOIL)

!     PARAMETERS:
!     MV - SOIL VOLUMETRIC MOISTURE[FRAC]
!     TD - SOIL TEMPERATURE [degC]
!     FGHZ - FREQUENCY [GHZ]
!     RHOKG - SOIL BULK DENSITY [KG/M3]
!     SAND - SAND CONTENT[%]
!     SILT - SILT CONTENT[%]
!     CLAY - CLAY CONTENT[%]

!     THIS IS A MODEL USED FOR 20-22 DEGC. THE FROZEN SOIL IS NOT CONSIDERED.
!     THERFORE, THE MV HERE IS AN EQUIVALENT VALUE.
!     TD, RHO, SAND, SILT ARE NOT USED!


IMPLICIT NONE

REAL(8), INTENT(IN) :: MV,TD,FGHZ,RHOKG,SAND,SILT,CLAY
COMPLEX(8),INTENT(OUT) :: EPS_SOIL


REAL(8), PARAMETER :: PI=3.1415926536
complex(8), parameter :: i0=complex(0,1)
REAL(8) :: C, F, OMEGA,nd,kd,mvt,e0b,tb,sb,su,e0u,tu
real(8) :: einf,eps0
real(8) :: nb,kb,nu,ku
real(8) :: nm,km, emr,emi
complex(8) :: eb,eu,nbc,nuc



C=CLAY/100.0d0
F=FGHZ*1.0d9
OMEGA=2d0*PI*F



!fit of MBSDM spectroscopic parameters, correlated to C
nd = 1.634d0 - 0.539d0*C + 0.2748d0*C**2.0d0
kd = 0.03952d0 - 0.04038d0*C
mvt = 0.02863d0 + 0.30673d0*C
e0b = 79.8d0 - 85.4d0*C + 32.7d0*C**2.0d0
tb = 1.062d-11 + 3.450d-12*C
sb = 0.3112d0 + 0.467d0*C
su = 0.3631d0 + 1.217d0*C
e0u = 100d0;
tu = 8.5d-12

!Debye relaxation equations for bound and free water components
einf = 4.9d0;
eps0 = 8.854d-12;
eb = einf + (e0b - einf)/(1 - i0*omega*tb) + i0*sb/omega/eps0;
eu = einf + (e0u - einf)/(1 - i0*omega*tu) + i0*su/omega/eps0;

nbc = sqrt(eb)
nb = real(nbc)
kb = imag(nbc)
nuc = sqrt(eu)
nu = real(nuc)
ku = imag(nuc)

!refractive mixing dielectric model
if(mv .lt. mvt)then
    nm = nd + (nb - 1)*mv;
    km = kd + kb*mv;
else
    nm = nd + (nb - 1)*mvt + (nu - 1)*(mv - mvt);
    km = kd + kb*mvt + ku*(mv - mvt);
end if

emr = nm**2.0d0 - km**2.0d0;
emi = 2d0*nm*km;

EPS_SOIL = complex(emr,-emi);




END SUBROUTINE EPSSOIL_MBSDM_MIRONOV
