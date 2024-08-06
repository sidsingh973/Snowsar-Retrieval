!Revisd initial values for normal and uniforma distributions
!Revised, now ThetaU and ThetaV is in log(Theta), ThetaPost saved Theta
!revised, use normal distribution instead of lognormal !!

Subroutine Initialize

Use CommonVars
Implicit None
Integer :: i
Real(8) :: Depth

Do i=1,NVars
  Na(i)=0;
End Do

Do i=1,Nlyr
  iDz(i)=i
  iRho(i)=Nlyr+i
  iD(i)=2*Nlyr+i
  iTsnow(i)=3*Nlyr+i
End Do

iTsoil(1)=4*Nlyr+1
iMvS(1)=4*Nlyr+2
iGndSig(1)=4*Nlyr+3
!New parameters
iP_M(1)=4*Nlyr+NsoilVars+1
iP_Q(1)=4*Nlyr+NsoilVars+2
iP_SR(1)=4*Nlyr+NsoilVars+3


Do i=1,Nlyr

  !Initialize snow thickness
  Select Case(InitialDz)
  Case(-1)
    ThetaU(iDz(i))=DzMu(i)-dsqrt(DzCov(i))
  Case(0)
    ThetaU(iDz(i))=DzMu(i)
  Case(1)
    ThetaU(iDz(i))=DzMu(i)+dsqrt(DzCov(i))
  End Select

  If(EstimateDz.Eq.0)Then
     ThetaU(iDz(i))=DzMu(i)    !+0.5*DzCov(i)
  End If


  !Initialize snow density
  Select Case(InitialRho)
  Case(-1)
    ThetaU(iRho(i))=RhoMu(i)-dsqrt(RhoCov(i))
  Case(0)
    ThetaU(iRho(i))=RhoMu(i)
  Case(1)
    ThetaU(iRho(i))=RhoMu(i)+dsqrt(RhoCov(i))
  End Select

  If(EstimateRho.Eq.0)Then
     ThetaU(iRho(i))=RhoMu(i)  !+0.5*RhoCov(i)
  End If


  !Initialize grain size
  Select Case(InitialD)
   Case(-1)
     ThetaU(iD(i))=DMu(i)-dsqrt(DCov(i))
   Case(0)
     ThetaU(iD(i))=DMu(i)
   Case(1)
     ThetaU(iD(i))=DMu(i)+dsqrt(DCov(i))
   End Select

   If(EstimateD.Eq.0)Then
      ThetaU(iD(i))=DMu(i)      !+0.5*DCov(i)
   Endif


   !Initialize snow temperature
   Select Case(InitialTsnow)
   Case(-1)
     ThetaU(iTsnow(i))=TsnowMu(i)-dsqrt(TsnowCov(i))
   Case(0)
     ThetaU(iTsnow(i))=TsnowMu(i)
   Case(1)
     ThetaU(iTsnow(i))=TsnowMu(i)+dsqrt(TsnowCov(i))
   End Select

   If(EstimateTsnow.Eq.0)Then
      ThetaU(iTsnow(i))=TsnowMu(i)   !+0.5*TsnowCov(i)
   End If
End Do


!Initialize soil temperature
Select Case(InitialSoil)
  Case(-1)
    ThetaU(iTsoil(1))=TsoilMu(1)-dsqrt(TsoilCov(1))
    ThetaU(iMvS(1))=MvSMu(1)-dsqrt(MvSCov(1))
    ThetaU(iGndSig(1))=GndSigMu(1)-dsqrt(GndSigCov(1))
Case(0)
    ThetaU(iTsoil(1))=TsoilMu(1)
    ThetaU(iMvS(1))=MvSMu(1)
    ThetaU(iGndSig(1))=GndSigMu(1)
Case(1)
    ThetaU(iTsoil(1))=TsoilMu(1)+dsqrt(TsoilCov(1))
    ThetaU(iMvS(1))=MvSMu(1)+dsqrt(MvSCov(1))
    ThetaU(iGndSig(1))=GndSigMu(1)+dsqrt(GndSigCov(1))
End Select

If(EstimateSoil.Eq.0)Then
   ThetaU(iTsoil(1))=TsoilMu(1)   !+0.5*TsoilCov(1)
   ThetaU(iMvS(1))=MvSMu(1)       !+0.5*MvSCov(1)
   ThetaU(iGndSig(1))=GndSigMu(1) !+0.5*GndSigCov(1)
End If


!Initialize model parameter
ThetaU(iP_M(1))=P_MMu(1)+0.5*P_MCov(1)
ThetaU(iP_Q(1))=P_QMu(1)+0.5*P_QCov(1)
ThetaU(iP_SR(1))=P_SRMu(1)+0.5*P_SRCov(1)


!Print first theta
Print *,'First Initialzied theta'
Call Show_Theta(ThetaU)


! Calculate Jump function
Do i=1,Nlyr
    JmpDzStd(i) = 0.3d0*dsqrt(DzCov(i))

    JmpRhoStd(i) = 0.3d0*dsqrt(RhoCov(i))

    JmpDStd(i) = 0.3d0*dsqrt(DCov(i))

    JmpTsnowStd(i) = 0.3d0*dsqrt(TsnowCov(i))
End Do

JmpTsoilStd(1)= 0.3d0*dsqrt(TsoilCov(1))

JmpMvSStd(1)= 0.3d0*Sqrt(MvSCov(1))

JmpGndSigStd(1) = 0.3d0*Sqrt(GndSigCov(1))

JmpP_MStd(1) = 0.3d0*Sqrt(P_MCov(1))

JmpP_QStd(1) = 0.3d0*Sqrt(P_QCov(1))

JmpP_SRStd(1) = 0.3d0*Sqrt(P_SRCov(1))


Print *,'Initial Jumps on log(Theta):'
Print *,'Initialized jump dz (cm):',JmpDzStd
Print *,'Initialized jump rho (kg/m^3):',JmpRhoStd
Print *,'Initialized jump D (mm):',JmpDStd
Print *,'Initialized jump Tsnow (degC):',JmpTsnowStd
Print *,'Initialized jump Tsoil (degC):',JmpTsoilStd
Print *,'Initialized jump MvS (%):',JmpMvSStd
Print *,'Initialized jump GndSig (mm):',JmpGndSigStd
Print *,'Initialized jump M:',JmpP_MStd
Print *,'Initialized jump Q:',JmpP_QStd
Print *,'Initialized jump SR:',JmpP_SRStd


!Treat the Radiometric Observations
Call ObsModel(ThetaU,ObsU)
Print *,'Initialized ObsU:', ObsU


!Treat the Other Observations
If(Np_other.gt.0)Then
    Depth=0.0d0
    Do i=1,Nlyr
        Depth=Depth+dexp(ThetaU(iDz(i)))
    End Do

    ObsU(1,1,Np_passive+Np_active+1)=Depth
EndIf


Call NProbObs(ObsU,Fu)
Call NProb(ThetaU,iDz,DzMu,DzCov,PuDz,Nlyr)
Call NProb(ThetaU,iRho,RhoMu,RhoCov,PuRho,Nlyr)
Call NProb(ThetaU,iD,DMu,DCov,PuD,Nlyr)
Call NProb(ThetaU,iTsnow,TsnowMu,TsnowCov,PuTsnow,Nlyr)
Call NProb(ThetaU,iTsoil,TsoilMu,TsoilCov,PuTsoil,1)
Call NProb(ThetaU,iMvS,MvSMu,MvSCov,PuMvS,1)
Call NProb(ThetaU,iGndSig,GndSigMu,GndSigCov,PuGndSig,1)
Call NProb(ThetaU,iP_M,P_MMu,P_MCov,PuP_M,1)
Call NProb(ThetaU,iP_Q,P_QMu,P_QCov,PuP_Q,1)
Call NProb(ThetaU,iP_SR,P_SRMu,P_SRCov,PuP_SR,1)


End Subroutine Initialize
