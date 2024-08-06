SUBROUTINE EPSSOIL_PAN(MV,TD,FGHZ,RHOKG,SAND,SILT,CLAY,EPS_SOIL)

!     PARAMETERS:
!     MV - SOIL VOLUMETRIC MOISTURE[FRAC]
!     TD - SOIL TEMPERATURE [degC]
!     FGHZ - FREQUENCY [GHZ]
!     RHOKG - SOIL BULK DENSITY [KG/M3]
!     SAND - SAND CONTENT[%]
!     SILT - SILT CONTENT[%]
!     CLAY - CLAY CONTENT[%]

!     T_limit is set as default value, 0degC.


IMPLICIT NONE

REAL(8), INTENT(IN) :: MV,TD,FGHZ,RHOKG,SAND,SILT,CLAY
COMPLEX(8),INTENT(OUT) :: EPS_SOIL


REAL(8), PARAMETER :: PI=3.1415926536d0
complex(8), parameter :: i0=complex(0d0,1d0)
REAL(8) :: C, F, OMEGA, eps0
Real(8) :: mvt, mvu
Real(8) :: A_C, A_mv, B_C, B_mv, A, B
Real(8) :: nd, kd
Real(8) :: einfu,e0u,tu,su, nu, ku
Real(8) :: einfb,e0b,tb,sb, nb, kb
Real(8) :: einft,e0t,tt,st, nt, kt
Real(8) :: epsr_ice, epsi_ice, ni, ki
Real(8) :: nm,km, emr,emi
Complex(8) :: eu, nuc, eb, nbc, et, ntc, nic



C=CLAY/100.0d0
F=FGHZ*1.0d9
OMEGA=2d0*PI*F
eps0 = 8.854d-12
A_C=0.00306d0
A_mv=0.394d0
B_C=0.00582d0
B_mv=-1.073d0


!fit of MBSDM spectroscopic parameters, correlated to C
nd = 1.634d0 - 0.539d0*C + 0.2748d0*C**2.0d0
kd = 0.03952d0 - 0.04038d0*C
mvt = 0.02863d0 + 0.30673d0*C


!Debye relaxation equations for free water components
einfu = 4.9d0
e0u=88.045d0 - 0.4147d0*TD + 6.295d-4*TD**2.0d0 + 1.075d-5*TD**3.0d0
tu = (17.68d0 - 0.6086d0*TD + 0.01104d0*TD**2.0d0 - 8.111d-5*TD**3.0d0) * 1d-12
su = 0.3631d0 + 1.217d0*C
eu = einfu + (e0u - einfu)/(1.0d0 - i0*omega*tu) + i0*su/omega/eps0
nuc = sqrt(eu)
nu = real(nuc)
ku = imag(nuc)


!Debye relaxation equations for bound water components
einfb=24.9d0+0.0685d0*TD
e0b=80.8d0+0.715d0*TD
tb=26.04d-12
sb = 0.3112d0 + 0.467d0*C
eb = einfb + (e0b - einfb)/(1.0d0 - i0*omega*tb) + i0*sb/omega/eps0
nbc = sqrt(eb)
nb = real(nbc)
kb = imag(nbc)


!Debye relaxation equations for transient water components
!T_limit is 0.0 degC
If(TD.le.0.0d0)Then
    !Related to unfrozen water content
    A = A_C*CLAY + A_mv*mv
    B = B_C*CLAY + B_mv*mv

    mvu=A*abs(TD)**B

    If(mvu.gt.mv)Then
        mvu=mv
    EndIf

    If(mvu.lt.mvt)Then
        mvu=mvt
    EndIf

    einft = 4.78d0
    e0t=77.90d0
    tt = 24.0*1d-12
    st = 0.3913d0
    et = einft + (e0t - einft)/(1.0d0 - i0*omega*tt) + i0*st/omega/eps0
    ntc = sqrt(et)
    nt = real(ntc)
    kt = imag(ntc)

    !Debye relaxation equations for ice components
    epsr_ice = 3.1884d0 + 9.1d-4*(TD+273.15d0-273d0)
    epsi_ice=0.0d0
    nic=sqrt(complex(epsr_ice,epsi_ice))
    ni = real(nic)
    ki = imag(nic)
EndIf




!refractive mixing dielectric model
if(TD.gt.0)then
    if(mv.lt.mvt)then
        nm = nd + (nb - 1.0d0)*mv
        km = kd + kb*mv
    else
        nm = nd + (nb - 1.0d0)*mvt + (nu - 1.0d0)*(mv - mvt)
        km = kd + kb*mvt + ku*(mv - mvt)
    endif
else
    if(mv.lt.mvt)then
        nm = nd + (nb - 1.0d0)*mv
        km = kd + kb*mv
    else
        if(mv.lt.mvu)then
            nm = nd + (nb - 1.0d0)*mvt + (nt - 1.0d0)*(mv - mvt)
            km = kd + kb*mvt + kt*(mv - mvt)
        else
            nm = nd + (nb - 1.0d0)*mvt + (nt - 1.0d0)*(mvu - mvt) + (ni - 1.0d0)*(mv - mvu)
            km = kd + kb*mvt + kt*(mvu - mvt) + ki*(mv-mvu)
        endif
    endif
endif


emr = nm**2.0d0 - km**2.0d0;
emi = 2d0*nm*km;

EPS_SOIL = complex(emr,-emi);


If(.False.)Then
    Print *,'(MV,TD,FGHZ,RHOKG,SAND,SILT,CLAY,EPS_SOIL)'
    Print *, MV,TD,FGHZ,RHOKG,SAND,SILT,CLAY,EPS_SOIL
EndIf

END SUBROUTINE EPSSOIL_PAN
