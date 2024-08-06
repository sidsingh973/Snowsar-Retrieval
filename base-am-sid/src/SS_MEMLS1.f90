!Note that CommonVars module is not used in this function, therefore, the same variable
!names are allowed.
!Variables Description:
!Freq: Frequencies
!Angle: Observation Angle
!Y: snow parameters, there is one row for each snow layer, where the Y(1) is the bottom layer
!       the variables in order are: LAYER THICKNESS [M], LAYER DENSITY [KG/M3],
!       LAYER GRAIN DIAMETER [mm], LAYER LIQUID WATER CONTENT [FRAC],
!       LAYER TEMPERATURE [K]
!Aux_Soil: Soil parameters:
!       Aux_Soil(1): Soil Temperature [K]
!       Aux_Soil(2): Soil Volumbetric water content [frac], Theta(Ntheta-1)
!       Aux_Soil(3): Ground roughness [m], Theta(Ntheta)
!       Aux_Soil(4): Soil density [kg/m^3]
!       Aux_Soil(5): Soil sand content [%]
!       Aux_Soil(6): Soil silt content [%]
!       Aux_Soil(7): Soil clay content [%]
! Aux_Model vector speficies 3 auxliliary observation model parameter
!       Aux_Model(1): M
!       Aux_Model(2): Q
!       Aux_Model(3): SR
!Tsky: the downwelling Tb upon the snow surface. The dimension is Tsky(Nfreq,2).
!       Tsky(:,1) is for vertical polarization, Tsky(:,2) is for horizontal polarization
!Pol_passive: Polarizations for Passive Measurement(Pol_Passive(Np_passive)): 0 (horizontal), 1 (vertica)
!Pol_active: Polarizations for Active Measurement(Pol_Active(Np_active)): 1 (vv), 2 (hh), 3 (vh), 4 (hv)
!Nfreq: Dimension of Freq
!Nangle: Dimension of Angle
!Np_passive: Dimension of Pol_passive
!Np_active: Dimension of Pol_active
!Np: Np_passive+Np_active+Np_other
!Num: number of snow layers, i.e. NLYR in the outside function
!ScatOpt: Scattering Coefficient Option (ScatOpt): 1 (Empirical MEMLS, Hallikainen-HUT), 2 (MEMLS-IBA, combined HUT), or 3 (Roy-HUT)


SUBROUTINE SS_MEMLS1(FREQ,ANGLE,Y,Aux_Soil,Aux_Model,TSKY,&
            POL_PASSIVE,POL_ACTIVE, &
            NFREQ,NANGLE,NP_PASSIVE,NP_ACTIVE,NP,NUM,SCATOPT,EstimateP_Q,Obs,ObsOut)

!
! ----------------------------------------------------------------------------
!
!     1.DESCRIPTION.  THE CODE IS A TRANSLATION OF THE MEMLS2 CODE RECEIVED
!     THROUGH CORRESPONDENCE WITH MATZLER AND A SOIL RADIATIVE TRANSFER
!     MODEL RECEIVED THROUGH CORRESPONDENCE WITH PULLIAINEN SOMETIME IN 2003.
!     THE MEMLS CODE USES THE GROUND REFLECTIVITY AND THE SNOWPACK DATA TO
!     COMPUTE THE BRIGHTNESS TEMPERATURE OF THE SNOWPACK. THE CODE WAS
!     TRANSLATED TO FORTRAN IN ORDER  TO SAVE COMPUTATION TIME IN A DATA
!     ASSIMILATION SCHEME IN WHICH MILLIONS OF CALLS MUST BE MADE TO THE
!     RADIATIVE TRANSFER MODEL .  IN TESTS, ONE HUNDRED BRIGHTNESS TEMPERATURE
!     CALCULATIONS TOOK TEN TIMES LESS PROCESSOR TIME IN THIS FORTRAN VERSION
!     THAN IN THE ORIGINAL MATLAB.  EXTENSIVE VERIFICATION WAS PERFORMED USING
!     SNOWPIT DATA FOR NINE DAYS IN ONE WINTER DURING THE CLPX PROGRAM AS WELL
!     AS 3-LAYER SNOWMODEL RESULTS FROM THE MODIFIED SAST (SUN ET AL, 99 IN
!     IN JGR) + SSIB MODEL RUN AT MAMMOTH MOUNTAIN, TO ENSURE THAT THE FORTRAN
!     CODE PRODUCES IDENTICAL RESULTS AS THE MATLAB CODE.
!
!     2. MODIFICATION.  THE MAIN CHANGE MADE TO THE MEMLS CODE IN THE
!     TRANSLATION PROCESS WAS THAT THE VAN POLDER APPROXIMATION (EFFECTIVE
!     MIXING THEORY) WAS  NO LONGER BEING USED IN THE MATLAB CODES WHEN I
!     RECEIVED THEM.  IN FACT, THE REAL PART OF THE SNOW DIELECTRIC WAS USED
!     INSTEAD IN SEVERAL SUBROUTINES. THIS WAS SLIGHTLY INCONSISTENT WITH THE
!     PAPERS, SO (FOLLOWING MATZLER'S ADVICE VIA EMAIL) I WROTE A NEW 'POLDER'
!     SUBROUTINE WHICH USES THE NEWTON-RAPHSON APPROXIMATION TO SOLVE FOR THE
!     SNOW EFFECTIVE PERMITTIVITY - SEE MATZLER'S 1996 PAPER IN IEEE FOR THE
!     EQUATIONS AND THE SUBROUTINE COMMENTS BELOW.  ACCORDING TO COMPARISON
!     WITH THE SNOWPIT DATA RESULTS, THE EFFECT OF USING THE SNOW DIELECTRIC
!     INSTEAD OF EFFECTIVE PERMITTIVITY ON THE BRIGHTNESS TEMPREATURE WAS
!     MINIMAL (~0.1 DEGREE).  HOWEVER, BECAUSE I COULD NOT TEST WITH A MORE
!     EXTENSIVE DATASET, I LEFT THE SLIGHTLY MORE COMPUTATIONALLY EXPENSIVE
!     EFFECTIVE PERMITTIVITY CALCULATION IN PLACE FOR THE SAKE OF INTERNAL
!     CONSISTENCY.
!
!     ANOTHER CHANGE THAT I MADE TO MEMLS WAS TO SET UP THE PROGRAM TO PASS IN
!     SNOW GRAIN DIAMETER SINCE THIS IS USUALLY WHAT IS MEASUSRED IN SNOWPITS
!     AND MODELED BY PROGNOSTIC EQUATIONS IN SNOW MODELS.  THE CONSTANT OF
!     PROPORTIONALITY BETWEEN THE GRAIN DIAMETER AND THE CORRELATION LENGTH (SEE
!     MATLZER, 2002 IN JOURNAL OF GLACIOLOGY) IS ALSO PASSED IN.
!
!     A THIRD CHANGE I MADE TO MEMLS WAS TO USE A SWITCH TO DETERMINE WHETHER
!     THE BORN APPROXIMATION OR THE EMPIRICAL SCATTERING COEFFICIENT WOULD BE
!     USED.  IF THE MAXIMUM CORRELATION LENGTH IN THE SNOWPACK IS GREATER THAN
!     0.33 MM, THE BORN APPROXIMATION IS USED; OTHERWISE, THE EMPIRICAL FORMULA
!     IS USED.
!
!     JINMEI'S REVISION (MAY18,2015): NOW THE SNOW DIELECTRIC CONSTANTS FOLLOWS
!     MEMLS3 - THE NEWEST VERSION OF MEMLS
!
!     HISTORY
!       1.0  - ORIGINAL TRANSLATION - MD 11/2005
!       2.0  - ADDED SUPPORT FOR SUCCESSIVE COHERENT LAYERS - MD & JG 4/2011
!       3.0  - ADDED DOBSON MODEL TO SIMULATE SOIL DIELECTRIC CONSTANTS



IMPLICIT NONE


!  A. DECLARATIONS
!   A0.INPUTS AND OUTPUTS
INTEGER,INTENT(IN) :: NFREQ,NANGLE,NP_PASSIVE,NP_ACTIVE,NP,NUM,ScatOpt
Integer,intent(In) :: Pol_Passive(Np_Passive), Pol_Active(Np_Active)
REAL(8),INTENT(IN) :: Aux_soil(7),Tsky(NFreq,Nangle),Y(NUM,5),Aux_Model(3)
Real(8),Intent(In) :: Freq(Nfreq),Angle(Nangle)
Real(8),Intent(In) :: Obs(Nfreq,Nangle,Np)
Real(8),Intent(Out) :: ObsOut(Nfreq,Nangle,Np)

!   A1.Constants and Set-up
INTEGER I,J,K
Integer Is_data
!constants
REAL(8) :: MJU0,EPS0,PI
!Snow parameters in Y
REAL(8),DIMENSION(:), ALLOCATABLE :: DI,ROIKG,PCI,WIFR,TI
!Soil Parameters in Aux_Soil
REAL(8) GND_TEMP,GND_MV,GND_SIG,SOIL_ROIKG,SAND,SILT,CLAY,GND_EPSI,GND_EPSII

!   A2. snow radiative transfer property calculation
Integer :: RNUM
REAL(8) TETA
!frequency-depedent only
REAL(8), Dimension(:),Allocatable :: EPSI,EPSII,GAI,NS
!angle&frequency-depedent
REAL(8), Dimension(:),Allocatable :: TEI,DEI,SIH,SIV
!variables after treating the coherent effect for thin layers
REAL(8),DIMENSION(:), ALLOCATABLE :: XROI,XEPSI,XEPSII,XTEI,&
    XSIH,XSIV,XDI,XDEI,XTI,XPCI,XWIFR,XGAI
REAL(8),DIMENSION(:), ALLOCATABLE :: RROI,REPSI,REPSII,RTEI,&
    RSIH,RSIV,RDI,RDEI,RTI,RPCI,RWIFR,RGAI

! A3. variables to compute soil properties
COMPLEX(8) :: EPS_UPPER,KSIG,EPS_SOIL
REAL(8) :: SOH,SOV,TETAD_UPPER,S0H,S0V,SS0H,SS0V


! A4. VARIALBES TO compute extinction coefficients
REAL(8),DIMENSION(:),ALLOCATABLE :: GBIH,GBIV,GS6,GA2I,&
    TSCAT,RSIHLONG,RSIVLONG
REAL(8),DIMENSION(:),ALLOCATABLE :: RSIHLONG_OLD,RSIVLONG_OLD

! A5. Variable to calculate brightness temperature and backscattering coefficient
!     Note that TRI is T matrix in MEMLS
REAL(8),DIMENSION(:),ALLOCATABLE :: RI,TRI,DH,DV0
REAL(8) TSKYH,TSKYV,TBV,TBH
REAL(8) ESG_H,ESG_V,SIGMAVV,SIGMAHH,SIGMAVH,SIGMAHV

! A6. Parameter used for active MEMLS
! Oct-12, add pex and q scaler
Integer EstimateP_Q
REAL(8) :: Param_M, Param_Q, Param_SR
REAL(8) :: Pex_scaler,Q_scaler

!  B. CONSTANTS AND CONTROL STATEMENTS
MJU0=1.2566d-006
EPS0=8.8542d-012
PI=3.14159d0
ObsOut=0.0d0


!  C. ALLOCATE AND EXTRACT STATEMENTS
ALLOCATE(DI(1:NUM),ROIKG(1:NUM),PCI(1:NUM),WIFR(1:NUM),TI(1:NUM))

DI(1:NUM)=   Y(1:NUM,1)
ROIKG(1:NUM)=Y(1:NUM,2)
PCI(1:NUM)=  Y(1:NUM,3)
WIFR(1:NUM)= Y(1:NUM,4)
TI(1:NUM)=   Y(1:NUM,5)

GND_TEMP=AUX_SOIL(1)
GND_MV=AUX_SOIL(2)
GND_SIG=AUX_SOIL(3)
SOIL_ROIKG=AUX_SOIL(4)
SAND=AUX_SOIL(5)
SILT=AUX_SOIL(6)
CLAY=AUX_SOIL(7)

Param_M=AUX_Model(1)
Param_Q=AUX_Model(2)
Param_SR=AUX_Model(3)


! SET UP
If(Np_passive.Eq.0 .and. Np_active.eq.0) Then
    Print *,'No radar or radiometer measurements. Exit!'
    Return
EndIf


! D. COMPUTE BRIGHTNESS TEMPERATURE OF SNOWPACK AND SOIL
DO i=1,NFREQ

    !D1.COMPUTE RADIATIVE TRANSFER PROPERTIES OF SNOW
    Allocate(EPSI(1:NUM),EPSII(1:NUM),GAI(1:NUM),NS(1:NUM))

    CALL RO2EPSD(ROIKG,TI,FREQ(i),EPSI,EPSII,NUM)
    CALL MIXMOD(FREQ(i),TI,WIFR,EPSI,EPSII,NUM)
    CALL ABSCOEFF(EPSI,EPSII,TI,FREQ(i),WIFR,GAI,NUM)
    NS=EPSI**0.5d0

    !SET UP
    Pex_scaler=1.0d0
    If(Np_passive.eq.0 .and. Np_active.gt.0) Then
        !Pex_scaler=1.23
        Pex_scaler=1.0d0
    EndIf

    Do j=1,Nangle

        !Check is there is data for this frequency&angle
        Is_data=0
        Do k=1,Np_active+Np_passive
           If(Obs(i,j,k)<900) Then
             Is_data=1
           End If
        EndDo

        If(Is_data.Eq.0) Cycle

        !D2.Compute angle-depedent RT properties of snow
        Allocate(TEI(1:NUM+1),DEI(1:NUM),SIH(1:NUM),SIV(1:NUM))
        Allocate(XROI(1:NUM),XEPSI(1:NUM),XEPSII(1:NUM),XTEI(1:NUM+1),&
            XSIH(1:NUM),XSIV(1:NUM),XDI(1:NUM),XDEI(1:NUM),&
            XTI(1:NUM),XPCI(1:NUM),XWIFR(1:NUM),XGAI(1:NUM))

        TETA=ANGLE(J)*PI/180.0
        !Calculate the angle inside the i-th snow layer
        TEI(1:NUM)=ASIN(SIN(TETA)/NS)
        TEI(NUM+1)=TETA
        !Calculate actual path length, Dei=di/cos(tei)
        CALL PFADI(TEI,DI,DEI,NUM)
        !Calculate the reflectivity at smooth interfaces between snow layers
        CALL FRESNELC(TEI,EPSI,SIH,SIV,NUM)
        !Treat coherent effect due to thin snow layers
        !place Pex scaler here~!
        CALL SLRED(NUM,ROIKG,EPSI,EPSII,TEI,SIH,SIV,DI,DEI,TI,PCI*Pex_scaler,WIFR,GAI,&
            FREQ(i),RNUM,XROI,XEPSI,XEPSII,XTEI,XSIH,XSIV,XDI,XDEI,&
            XTI,XPCI,XWIFR,XGAI)

        ALLOCATE(RROI(RNUM),REPSI(RNUM),REPSII(RNUM),RTEI(RNUM+1),&
        RSIH(RNUM),RSIV(RNUM),RDI(RNUM),RDEI(RNUM),&
        RTI(RNUM),RPCI(RNUM),RWIFR(RNUM),RGAI(RNUM))

        DO k=1,RNUM
            RROI(k)=XROI(k)
            REPSI(k)=XEPSI(k)
            REPSII(k)=XEPSII(k)
            RSIH(k)=XSIH(k)
            RSIV(k)=XSIV(k)
            RDI(k)=XDI(k)
            RDEI(k)=XDEI(k)
            RTI(k)=XTI(k)
            RPCI(k)=XPCI(k)
            RWIFR(k)=XWIFR(k)
            RGAI(k)=XGAI(k)
        END DO
        DO k=1,RNUM+1
            RTEI(k)=XTEI(k)
        END DO

        !D3.COMPUTE GROUND REFLECTIVITIES
        TETAD_UPPER=RTEI(1)/pi*180.0d0
        EPS_UPPER=dCMPLX(REPSI(1),(-1*REPSII(1)))
        KSIG=REAL(2d0*PI*FREQ(i)*1.0d9*(MJU0*EPS0*EPS_UPPER)**0.5d0)*GND_SIG


        IF(.False.)Then
            CALL RUFFSOIL_WM(FREQ(i),GND_MV,GND_TEMP,KSIG,GND_SIG,TETAD_UPPER,EPS_UPPER,&
                SOIL_ROIKG,SAND,SILT,CLAY,EPS_SOIL,S0H,S0V,SS0H,SS0V)

        ELSE
            CALL RUFFSOIL_QHN(FREQ(i),GND_MV,GND_TEMP,KSIG,GND_SIG,TETAD_UPPER,EPS_UPPER,&
                SOIL_ROIKG,SAND,SILT,CLAY,EPS_SOIL,S0H,S0V,SS0H,SS0V)
            !Print *,'KSIG',KSIG

        END IF


        !D4.COMPUTE SCATTERING COEFFICIENTS, REFLECTIVITIES AND
        !     TRANSMISSIVITIES

        DEALLOCATE(RTEI,RDEI) ! BECAUSE THESE ARRAYS WILL CHANGE SIZES

        ALLOCATE(GBIH(RNUM),GBIV(RNUM),GS6(RNUM),GA2I(RNUM),&
            TSCAT(RNUM),RTEI(RNUM),RDEI(RNUM),RSIHLONG(RNUM+1),&
            RSIVLONG(RNUM+1))
        ALLOCATE(RSIHLONG_OLD(RNUM+1),RSIVLONG_OLD(RNUM+1))

        !Improved Born Approximation simplified
        CALL SCCOEFF(RROI,RTI,RPCI,FREQ(i),RWIFR,RGAI,ScatOpt,GBIH,GBIV,&
            GS6,GA2I,RNUM,REPSI,REPSII)


        CALL PFADC(TETA,RDI,REPSI,GS6,RDEI,RTEI,TSCAT,RNUM)


        RSIHLONG(1)=S0H
        RSIHLONG(2:RNUM+1)=RSIH
        RSIVLONG(1)=S0V
        RSIVLONG(2:RNUM+1)=RSIV

        RSIHLONG_OLD=RSIHLONG
        RSIVLONG_OLD=RSIVLONG

        !Calculate polarization mixing
        CALL POLMIX(TSCAT,RSIHLONG,RSIVLONG,RNUM)

        !D5.  COMPUTE BRIGHTNESS TEMPERATURES AND EMISSIVITIES OF SNOWPACK AND
        !       GROUND


        !!!!!!!!!
        !Print *, 'ALLOCATE(DV0(1:RNUM))'
        ALLOCATE(DV0(1:RNUM))
        !!!!!!!!!

        !Print *, 'ALLOCATE(DH(1:RNUM))'
        ALLOCATE(DH(1:RNUM))

        !Print *, 'ALLOCATE(RI(1:RNUM),TRI(1:RNUM))'
        ALLOCATE(RI(1:RNUM),TRI(1:RNUM))


        ! Jinmei 9/27/2016
        TskyV=Tsky(i,j)
        TskyH=Tsky(i,j)


        CALL RT(GA2I,GBIH,RDEI,RI,TRI,RNUM)

        CALL LAYER(RI,RSIHLONG,TRI,RTI,GND_TEMP,TSKYH,DH,RNUM)

        CALL RT(GA2I,GBIV,RDEI,RI,TRI,RNUM)

        CALL LAYER(RI,RSIVLONG,TRI,TI,GND_TEMP,TSKYV,DV0,RNUM)

        TBV=(1-RSIVLONG(RNUM+1))*DV0(RNUM)+RSIVLONG(RNUM+1)*TskyV
        TBH=(1-RSIHLONG(RNUM+1))*DH(RNUM)+RSIHLONG(RNUM+1)*TskyH


        Do k=1,Np_passive
            Select Case(Pol_passive(k))
            Case(1)
                ObsOut(i,j,k)=TBV
            Case(2)
                ObsOut(i,j,k)=TBH
            End Select
        End Do


        !Revised Sep 23, 2016
        If(Np_passive.gt.0 .and. Np_active .gt. 0 .and. Pex_scaler.eq.1.0d0)Then
            !Pex_scaler=1.23
            Pex_scaler=1.0d0

            !Recalculate for scaled pex
            !Redo from D4

            Deallocate(GBIH,GBIV,GS6,GA2I,TSCAT)
            Deallocate(RI,TRI,DH,DV0)

            ALLOCATE(GBIH(RNUM),GBIV(RNUM),GS6(RNUM),GA2I(RNUM),TSCAT(RNUM))
            ALLOCATE(RI(1:RNUM),TRI(1:RNUM),DH(1:RNUM),DV0(1:RNUM))


            !Improved Born Approximation simplified
            !Place Pex_scaler here ~!
            CALL SCCOEFF(RROI,RTI,RPCI*Pex_scaler,FREQ(i),RWIFR,RGAI,ScatOpt,GBIH,GBIV,&
                GS6,GA2I,RNUM,REPSI,REPSII)

            CALL PFADC(TETA,RDI,REPSI,GS6,RDEI,RTEI,TSCAT,RNUM)

            !Calculate polarization mixing
            CALL POLMIX(TSCAT,RSIHLONG,RSIVLONG,RNUM)

            !D5.  COMPUTE BRIGHTNESS TEMPERATURES AND EMISSIVITIES OF SNOWPACK AND
            !       GROUND
            CALL RT(GA2I,GBIH,RDEI,RI,TRI,RNUM)

            CALL LAYER(RI,RSIHLONG,TRI,RTI,GND_TEMP,TSKYH,DH,RNUM)

            CALL RT(GA2I,GBIV,RDEI,RI,TRI,RNUM)

            CALL LAYER(RI,RSIVLONG,TRI,TI,GND_TEMP,TSKYV,DV0,RNUM)

            TBV=(1-RSIVLONG(RNUM+1))*DV0(RNUM)+RSIVLONG(RNUM+1)*TskyV
            TBH=(1-RSIHLONG(RNUM+1))*DH(RNUM)+RSIHLONG(RNUM+1)*TskyH
        EndIf


        !D6. COMPUTE BACKSCATTERING COEFFIICENT
        CALL EMISSIVITY(GA2I,GBIH,GBIV,RDEI,RSIHLONG,RSIVLONG,GND_TEMP,RTI,&
            ESG_H,ESG_V,RNUM)

        !IF(FREQ(i).gt.10.0 .and. FREQ(i).lt.11.0)Then
        !    Q_scaler=1.0
        !EndIf

        !IF(FREQ(i).gt.13.0 .and. FREQ(i).lt.14.0)Then  !Q fitted for 13.3GHz
        !    Q_scaler=1.06
        !EndIf

        !IF(FREQ(i).gt.16.0 .and. FREQ(i).lt.17.0)Then  !Q fitted for 16.7GHz
        !    Q_scaler=1.28
        !EndIf
        Q_scaler=1.0d0

        CALL BACKSCATTERING(RGAI,GS6,RDEI,TETA,SS0H,SS0V,&
            RSIHLONG_OLD,RSIVLONG_OLD,ESG_H,ESG_V,&
            Param_M, Param_Q*Q_scaler, Param_SR,EstimateP_Q,&
            SIGMAVV,SIGMAHH,SIGMAVH,SIGMAHV,RNUM)

        Do k=1, Np_active
            Select Case(Pol_active(k))
            Case(1)
                ObsOut(i,j,Np_passive+k)=10*log10(SigmaVV)
            Case(2)
                ObsOut(i,j,Np_passive+k)=10*log10(SigmaHH)
            Case(3)
                ObsOut(i,j,Np_passive+k)=10*log10(SigmaVH)
            Case(4)
                ObsOut(i,j,Np_passive+k)=10*log10(SigmaHV)
            End Select
        End Do



        If(.False.)Then
            Print *,'Calculate Parameter!!'
            Print *,'DI',DI
            Print *,'ROIKG',ROIKG
            Print *,'PCI',PCI
            Print *,'RPCI',RPCI
            Print *,'TI',TI
            Print *,'Y',Y
            Print *,'Param_Q',Param_Q
            Print *,'GND_TEMP',GND_TEMP
            Print *,'GND_MV',GND_MV
            Print *,'GND_SIG',GND_SIG
            Print *,'EPS_SOIL',EPS_SOIL
            Print *,'S0V,S0H,SS0V,SS0H',S0V,S0H,SS0V,SS0H
            Print *,'RSIHLONG',RSIHLONG
            Print *,'RSIVLONG',RSIVLONG
            Print *,'EPSI,EPSII',EPSI,EPSII
            Print *,'GAI',GAI
            Print *,'GA2I',GA2I
            Print *,'GBIV',GBIV
            Print *,'GBIH',GBIH
            Print *, 'GS6',GS6
            !Print *, 'TBV,TBH',TBV,TBH
            !Print *, 'SigmaVV, SigmaHH', SigmaVV, SigmaHH
            !Print *, 'SigmaVH, SigmaHV', SigmaVH, SigmaHV
            Print *, 'ObsOut',ObsOut
        ENDIF


        ! F.  DEALLOCATE VARIABLE SPACE
        Deallocate(TEI,DEI,SIH,SIV)
        Deallocate(XROI,XEPSI,XEPSII,XTEI, &
            XSIH,XSIV,XDI,XDEI, &
            XTI,XPCI,XWIFR,XGAI)
        Deallocate(RROI,REPSI,REPSII,RTEI, &
            RSIH,RSIV,RDI,RDEI, &
            RTI,RPCI,RWIFR,RGAI)
        Deallocate(GBIH,GBIV,GS6,GA2I, &
            TSCAT,RSIHLONG,RSIVLONG)
        Deallocate(RSIHLONG_Old,RSIVLONG_Old)
        Deallocate(RI,TRI,DH,DV0)

    End Do !End Angle

    Deallocate(EPSI,EPSII,GAI,NS)

End Do !End Frequency


DEALLOCATE(DI,ROIKG,PCI,WIFR,TI)

!Print *, 'ObsOut',ObsOut

END SUBROUTINE SS_MEMLS1









! -------------------------------------------------------------------------
!
SUBROUTINE EMISSIVITY(GA2I,GBIH,GBIV,DEI,SIHLONG,SIVLONG,GND_TEMP,TI,&
EH,EV,NUM)
!
! -------------------------------------------------------------------------
!
!  CODE ORIGINALLY OBTAINED IN MATLAB FROM MATZLER, MODIFIED BY MIKE.
!    SEE VERSION HISTORY
!
!   CALCULATES THE SCATTERING COEFFICIENT USING BORN APPROXIMATION
!
!       GA2I: ABSORPTION COEFFICIENT
!       GBIH: 2-FLUX SCATTERING COEFFICIENT, H POLARIZATION
!       GBIV: 2-FLUX SCATTERING COEFFICIENT, V POLARIZATION
!       DEI:  EFFECTIVE PATH LENGTH [M]
!       SIHLONG: LAYER INTERFACE REFLECTIVITY, H POLARIZATION
!       SIVLONG: LAYER INTERFACE REFLECTIVITY, V POLARIZATION
!       EH:   EMISSIVITY, H POLARIZATION
!       EV:   EMISSIVITY, V POLARIZATION
!
!   VERSION HISTORY:
!      1.0     MD 1 APR 05 THIS CODE WAS PART OF LMAIN.  I TRANSLATED TO
!                 FORTRAN FROM MATLAB AND MOVED IT TO A SEPARATE SUBROUTINE.
!                 COMPARE WIESMANN AND MATZLER, 99 EQN (8)
!
!   USES: RT, LAYER
!
!   COPYRIGHT (C) 1998 BY THE INSTITUTE OF APPLIED PHYSICS,
!   UNIVERSITY OF BERN, SWITZERLAND

IMPLICIT NONE

INTEGER, INTENT(IN) :: NUM
REAL(8),INTENT(IN) :: GA2I(NUM),GBIH(NUM),GBIV(NUM),DEI(NUM),&
SIHLONG(NUM+1),SIVLONG(NUM+1),GND_TEMP,&
TI(NUM)
REAL(8),INTENT(OUT) :: EH,EV
REAL(8) RI(NUM),TRI(NUM),DH(NUM),DV(NUM),TBH0,TBH100,TBV0,TBV100,TSKY



! HORIZONTAL BRIGHTNESS TEMPERATURES UNDER DIFFERENT TSKY VALS

CALL RT(GA2I,GBIH,DEI,RI,TRI,NUM)

TSKY=0.0d0
CALL LAYER(RI,SIHLONG,TRI,TI,GND_TEMP,TSKY,DH,NUM)
TBH0=(1-SIHLONG(NUM+1))*DH(NUM)+SIHLONG(NUM+1)*TSKY

TSKY=100.0d0
CALL LAYER(RI,SIHLONG,TRI,TI,GND_TEMP,TSKY,DH,NUM)
TBH100=(1-SIHLONG(NUM+1))*DH(NUM)+SIHLONG(NUM+1)*TSKY

! VERTICAL BRIGHTNESS TEMPERATUERS UNDER DIFFERENT TSKY VALS

CALL RT(GA2I,GBIV,DEI,RI,TRI,NUM)

TSKY=0.0d0
CALL LAYER(RI,SIVLONG,TRI,TI,GND_TEMP,TSKY,DV,NUM)
TBV0=(1-SIVLONG(NUM+1))*DV(NUM)+SIVLONG(NUM+1)*TSKY

TSKY=100.0d0
CALL LAYER(RI,SIVLONG,TRI,TI,GND_TEMP,TSKY,DV,NUM)
TBV100=(1-SIVLONG(NUM+1))*DV(NUM)+SIVLONG(NUM+1)*TSKY

! COMPUTE EMISSIVITIES
EH=1-(TBH100-TBH0)/100.0d0
EV=1-(TBV100-TBV0)/100.0d0

END SUBROUTINE EMISSIVITY



! -------------------------------------------------------------------------
!
SUBROUTINE BACKSCATTERING(RGAI,GS6,RDEI,TETA,SS0H,SS0V,&
    RSIHLONG,RSIVLONG,ESG_H,ESG_V,&
    Param_M,Param_Q,Param_SR,EstimateP_Q, &
    SIGMAVV,SIGMAHH,SIGMAVH,SIGMAHV,RNUM)
!
! -------------------------------------------------------------------------
! CODED BY JINMEI, MAY 10,2016
! REFER TO MEMLS3&A - AMEMLSMAIN
! RGAI: ABSORPTION COEFFIICENT
! GS6:  KS, SCATTERING COEFFIICENT
! RDEI: DI/COS(THETA)
! TETA: INCIDENCE ANGLE AT THE SNOW SURFACE
! SS0H: SPECULAR PRART OF THE SOIL REFILECTIVITY, H POL.
! SS0V: SPECULAR PRART OF THE SOIL REFILECTIVITY, V POL.
! RSIHLONG: REFLECTIVITIES AT BOUNDARIES, H POL.
! RSIVLONG: REFLECTIVITIES AT BOUNDARIES, V POL.
! ESG_H: TOTAL EMISSIVITY OF THE SNOW
! EGS_V: TOTAL EMISSIVITY OF THE SNOW
! SIGMAVV, SIGMAHH, SIGMAVH, SIGMAHV: BACKSCATTERING COEFFIICENT
! RNUM: NUMBER OF SNOW LAYERS

IMPLICIT NONE
INTEGER,INTENT(IN) :: RNUM
REAL(8),INTENT(IN) :: RGAI(RNUM),GS6(RNUM),RDEI(RNUM),TETA,SS0H,SS0V
REAL(8),INTENT(IN) :: RSIHLONG(RNUM+1),RSIVLONG(RNUM+1),ESG_H,ESG_V
REAL(8),INTENT(IN) :: Param_M,Param_Q,Param_SR
Integer :: EstimateP_Q
REAL(8),INTENT(OUT) :: SIGMAVV,SIGMAHH,SIGMAVH,SIGMAHV

INTEGER I
REAL(8) M,Q
REAL(8) :: REXTC(RNUM),U(RNUM),U2(RNUM),RH(RNUM+1),RV(RNUM+1)
REAL(8) :: FORWA,R_SH,R_SV,R_S0,R_H,R_V,R_DH,R_DV
REAL(8) :: MU,MU2,M2,SIGMA_DV,SIGMA_DH,SIGMA_DVV,SIGMA_DHH,XPON,SIGMA_S

!CONSTANTS, TO BE REVISED
!M=0.1
!Q=0.1
M=Param_M
Q=Param_Q

!CALCULATE THE SPECULAR PART OF REFLECTIVITY
FORWA=4.0d0
REXTC=RGAI+FORWA*GS6
U=EXP(-REXTC*RDEI)
U2=U*U
!RH(1)=SS0H
!RV(1)=SS0V
RH(1)=SS0H*Param_SR
RV(1)=SS0V*Param_SR

DO I=2,RNUM+1
    RH(I)=RSIHLONG(I)+RH(I-1)*((1-RSIHLONG(I))*U(I-1))**2.0d0 &
            /(1-U2(I-1)*RSIHLONG(I)*RH(I-1))
    RV(I)=RSIVLONG(I)+RV(I-1)*((1-RSIVLONG(I))*U(I-1))**2.0d0 &
            /(1-U2(I-1)*RSIVLONG(I)*RV(I-1))
END DO
R_SH=RH(RNUM+1)
R_SV=RV(RNUM+1)
R_S0=0.5d0*(R_SH+R_SV)

!CALCULATE THE BACKSCATTERING COEFFICIENT
R_V=1.0d0-ESG_V
R_H=1.0d0-ESG_H

R_DV=R_V-R_SV
R_DH=R_H-R_SH

!IF(EstimateP_Q.eq.0)Then    !Revised Q, to be calcualted from DV/V
!    !Q=0.04605 * R_DV/R_V + 0.07935
!    Q=0.056171d0 * R_DV/R_V + 0.070302d0
!EndIf

MU=COS(TETA)
MU2=MU*MU
M2=M**2.0d0
SIGMA_DV=4.0d0*R_DV*MU2
SIGMA_DH=4.0d0*R_DH*MU2
SIGMA_DVV=(1-Q)*SIGMA_DV
SIGMA_DHH=(1-Q)*SIGMA_DH
SIGMAHV=Q*0.5d0*(SIGMA_DV+SIGMA_DH)
SIGMAVH=SIGMAHV

XPON= -(TAN(TETA))**2.0d0/(2.0d0*M2);
SIGMA_S=R_S0*EXP(XPON)/(2.0d0*M2*MU2*MU2)
SIGMAVV=SIGMA_DVV + SIGMA_S
SIGMAHH=SIGMA_DHH + SIGMA_S

If(.False.)Then
    Print *,'Check in BACKSCATTERING'
    Print *,'REXTC',REXTC
    Print *,'RDEI',RDEI
    Print *, 'U',U
    Print *, 'RH', RH
    Print *, 'RV', RV
    Print *,'R_V&H', R_V, R_H
    Print *,'R_SV&H', R_SV, R_SH
    Print *,'R_DV&H', R_DV,R_DH
    Print *,'MU', MU
    Print *,'M2', M2
    Print *,'SIGMA_DV',SIGMA_DV
    Print *,'SIGMA_DH',SIGMA_DH
    Print *,'Q',Q
    Print *,'SIGMA_DVV',SIGMA_DVV
    Print *,'SIGMA_S',SIGMA_S
End If

END SUBROUTINE BACKSCATTERING
