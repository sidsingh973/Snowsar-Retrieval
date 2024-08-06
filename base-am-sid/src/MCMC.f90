Subroutine MCMC

Use CommonVars
Implicit None
Integer :: i,j,k

Real(8) :: Z1(Niter*Nlyr),Z2(Niter*Nlyr),Z3(Niter*Nlyr),Z4(Niter*Nlyr)
Real(8) :: Z5(Niter),Z6(Niter),Z7(Niter),Z8(Niter),Z9(Niter),Z10(Niter)
Real(8) :: U1(Niter),U2(Niter),U3(Niter),U4(Niter),U5(Niter), U6(Niter),U7(Niter)
Real(8) :: U8(Niter),U9(Niter),U10(Niter)
Real(8) :: Depth
integer :: ConstraintD, RndSeed

Allocate(ThetaU(Ntheta),ThetaV(Ntheta))
Allocate(JmpDzStd(Nlyr),JmpRhoStd(Nlyr),JmpDStd(Nlyr),JmpTsnowStd(Nlyr))
Allocate(JmpTsoilStd(1),JmpMvSStd(1),JmpGndSigStd(1))
Allocate(JmpP_MStd(1),JmpP_QStd(1),JmpP_SRStd(1))
Allocate(iDz(Nlyr),iRho(Nlyr),iD(Nlyr),iTsnow(Nlyr))
Allocate(iTsoil(1),iMvS(1),iGndSig(1))
Allocate(iP_M(1),iP_Q(1),iP_SR(1))


!Initialize, Set the first values of the state vector using the priors. Calculate the jump function before burn-in
Call Initialize


!Calcualte Intermediate Parameters
!Z*, to be multiplied by JmpStd
RndSeed=1
Call RandNormal(Niter*Nlyr,Z1,0+iPit+Nlyr+RndSeed)
Call RandNormal(Niter*Nlyr,Z2,100+iPit+Nlyr+RndSeed)
Call RandNormal(Niter*Nlyr,Z3,200+iPit+Nlyr+RndSeed)
Call RandNormal(Niter*Nlyr,Z4,300+iPit+Nlyr+RndSeed)
Call RandNormal(Niter,Z5,800+iPit+Nlyr+RndSeed)
Call RandNormal(Niter,Z6,900+iPit+Nlyr+RndSeed)
Call RandNormal(Niter,Z7,1000+iPit+Nlyr+RndSeed)
Call RandNormal(Niter,Z8,1100+iPit+Nlyr+RndSeed)
Call RandNormal(Niter,Z9,1200+iPit+Nlyr+RndSeed)
Call RandNormal(Niter,Z10,1300+iPit+Nlyr+RndSeed)

!U*, to be compared with likelihood ratio
Call RandUniform(Niter,U1,400+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U2,500+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U3,600+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U4,700+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U5,1100+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U6,1200+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U7,1300+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U8,1400+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U9,1500+iPit+Nlyr+RndSeed)
Call RandUniform(Niter,U10,1600+iPit+Nlyr+RndSeed)



!Start running MCMC iterations
Do i=1,Niter

  !Adjust the jump function during the burn-in period at every 100 iterations.
  If(i.le.Nburn)Then
    If(Mod(i,100).eq.0)Then
        Print *, '#Iter - ', i
        Call AdjustJump(i)
    EndIf
  Else
    If(Mod(i,1000).eq.0)Then
        Print *, '#Iter - ', i
        Call AdjustJump(i)
    EndIf
  EndIf


  !Update thickness
  If(EstimateDz.Eq.1)Then
    Call Iterate(Z1((i-1)*Nlyr+1:i*Nlyr),U1(i),iDz,&
       JmpDzStd,DzMu,DzCov,PuDz,1,&
       DzMinLim,DzMaxLim,Nlyr,0,i)
  End If


  !Update density
  If(EstimateRho.Eq.1)Then
    Call Iterate(Z2((i-1)*Nlyr+1:i*Nlyr),U2(i),iRho,&
       JmpRhoStd,RhoMu,RhoCov,PuRho,2,&
       RhoMinLim,RhoMaxLim,Nlyr,ConstrainRho,i)
  End If

  !Update grain size, Constrained D here
  ConstraintD=1
  If(EstimateD.Eq.1)Then
    Call Iterate(Z3((i-1)*Nlyr+1:i*Nlyr),U3(i),iD,&
       JmpDStd,DMu,DCov,PuD,3,&
       DMinLim,DmaxLim,Nlyr,ConstraintD,i)
  End If

  !Update temperature
  If(EstimateTsnow.Eq.1)Then
    Call Iterate(Z4((i-1)*Nlyr+1:i*Nlyr),U4(i),iTsnow,&
       JmpTsnowStd,TsnowMu,TsnowCov,PuTsnow,4,&
       TsnowMinLim,TsnowMaxLim,Nlyr,ConstrainTsnow,i)

    Call Iterate(Z5(i),U5(i),iTsoil,&
        JmpTsoilStd,TsoilMu,TsoilCov,PuTsoil,5,&
        TsoilMinLim,TsoilMaxLim,1,0,i)
  End If


  !Always estimate soil moisture, use EstimateSoil to decide estimating roughness or not~
  Call Iterate(Z6(i),U6(i),iMvS,&
    JmpMvSStd,MvSMu,MvSCov,PuMvS,6,&
    MvSMinLim,MvSMaxLim,1,0,i)


  !Update soil parameters
  If(EstimateSoil.Eq.1)Then
     Call Iterate(Z7(i),U7(i),iGndSig,&
        JmpGndSigStd,GndSigMu,GndSigCov,PuGndSig,7,&
        GndSigMinLim,GndSigMaxLim,1,0,i)
  End If

  !Update model parameter, M
  If(EstimateP_M.Eq.1)Then
     Call Iterate(Z8(i),U8(i),iP_M, &
        JmpP_MStd,P_MMu,P_MCov, PuP_M, 8, &
        P_MMinLim,P_MMaxLim,1,0,i)
  End If

  !Update model parameter, Q
  If(EstimateP_Q.Eq.1)Then
    Call Iterate(Z9(i),U9(i),iP_Q, &
        JmpP_QStd,P_QMu,P_QCov, PuP_Q, 9, &
        P_QMinLim,P_QMaxLim,1,0,i)
  End If

  !Update model parameter, SR (adjusting ratio of specular part of reflectivity)
  If(EstimateP_SR.Eq.1)Then
    Call Iterate(Z10(i),U10(i),iP_SR, &
        JmpP_SRStd,P_SRMu,P_SRCov, PuP_SR, 10, &
        P_SRMinLim,P_SRMaxLim,1,0,i)
  End If


  !Update ThetaPost
  !ThetaPost(1:Ntheta,i)=dexp(ThetaU)
  ThetaPost(1:Ntheta,i)=ThetaU

  !Jinmei, 2017/6/30 Use relative thickness for 2 to Nlyr
  !Update Thickness
  If(Nlyr.gt.1)Then
    Do j=2,Nlyr
        !ThetaPost(iDz(j),i)=dexp(ThetaU(iDz(j))) * dexp(ThetaU(iDz(1)))
        ThetaPost(iDz(j),i)=ThetaU(iDz(j)) * ThetaU(iDz(1))
    End Do
  EndIf


  !Update TbPost
  Do j=1,Nobs
     ObsPost(j,i)=ObsU(ObsI(j),ObsJ(j),ObsK(j))
  End Do


  !Print *, 'ThetaPost', ThetaPost(1:Ntheta,i)
  !Print *, 'ObsPost(j,i)',ObsPost(j,i)
End Do



!Deallocate parameters allocated in MCMC.f90
Deallocate(ThetaU,ThetaV)
Deallocate(JmpDzStd,JmpRhoStd,JmpDStd,JmpTsnowStd)
Deallocate(JmpTsoilStd,JmpMvSStd,JmpGndSigStd)
Deallocate(JmpP_MStd,JmpP_QStd,JmpP_SRStd)
Deallocate(iDz,iRho,iD,iTsnow,iTsoil,iMvS,iGndSig)
Deallocate(iP_M,iP_Q,iP_SR)


End Subroutine MCMC

