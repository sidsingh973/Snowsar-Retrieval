!Iterate function:
!Z: Random number to produce the actual jump
!U: Random number following uniform distribution
!iSel: index for this current variable in theta
!JmpStd: jump function
!Mu: mean value
!Cov: covariance
!Pu: prior function
!NaX: Acceptance counts
!ThetaMin: minimal value of Theta
!ThetaMax: maximal value of Theta
!NiSel: Number of layers
!ConstrainOpt: option to constrain the values
!iIter: current number of iterations in MCMC

Subroutine Iterate(Z,U,iSel,JmpStd,Mu,Cov,Pu,NaXi,ThetaMin,ThetaMax,NiSel,ConstrainOpt,iIter)

Use CommonVars
Implicit None

Integer,Intent(In) :: iSel(NiSel),NaXi,NiSel,ConstrainOpt,iIter
Integer :: i
Real(8),Intent(In) :: Z(NiSel),JmpStd(NiSel),Mu(NiSel),Cov(NiSel),U,ThetaMin,ThetaMax
Real(8),Intent(InOut) :: Pu
Real(8) :: R,Pv
Real(8) :: Temp1, Temp2, ThetaTemp(NiSel),Acept_Temp, Depth

ThetaV=ThetaU

!Update Theta
!Upate thickness here.
!Call Iterate(Z1((i-1)*Nlyr+1:i*Nlyr),U1(i),iDz,&
!JmpDzStd,DzMu,DzCov,PuDz,1,&
!DzMinLim,DzMaxLim,Nlyr,0,i)

Do i=1,NiSel
  ThetaTemp(i)=ThetaU(iSel(i))+Z(i)*JmpStd(i)
  !ThetaTemp(i)=Min(ThetaTemp(i),dlog(ThetaMax))
  !ThetaTemp(i)=Max(ThetaTemp(i),dlog(ThetaMin))
  ThetaTemp(i)=Min(ThetaTemp(i),ThetaMax)
  ThetaTemp(i)=Max(ThetaTemp(i),ThetaMin)
End Do


!Apply Layer Constrains. Note that 1 is bottom layer, NiSel is the surface layer
If(NiSel.Gt.1 .And. ConstrainOpt.Ne.0)Then

  Temp1=ThetaTemp(1)   !Temp1 is the reference

  Do i=2,NiSel

    Temp2=ThetaTemp(i)

    !If bottom>surface
    If(ConstrainOpt.Eq.1)Then
        If(Temp2.Gt.Temp1)Then
           ThetaTemp(i)=Temp1
        Else
           Temp1=Temp2
        End If
    End If

    !If bottom<surface
    If(ConstrainOpt.Eq.2)Then
        If(Temp2.Lt.Temp1)Then
          ThetaTemp(i)=Temp1
        Else
          Temp1=Temp2
        End If
    End If
  End Do
End If


Do i=1,NiSel
  ThetaV(iSel(i))=ThetaTemp(i)
End Do


Call ObsModel(ThetaV,ObsV)


!Treat the Other Observations when the function is updating the snow layer thickness
!Will it influence the simulation? Add Np_other.Gt.0 at Apr18,2017 (prevent rewriting ObsV wrongly)
If (iSel(1).Eq.1 .And. Np_other.Gt.0) Then
    Depth=ThetaV(iDz(1))
    Do i=2,NiSel
        !Depth=Depth+dexp(ThetaV(iDz(i)))

        Depth=Depth+ThetaV(iDz(i)) * ThetaV(iDz(1))

    End Do
    ObsV(1,1,Np_passive+Np_active+1)=Depth
End If

!If(iSel(1).Eq.1 .And. Np_other.Eq.0)Then
    !Print *,'Current ObsV after updating Dz only', ObsV
!EndIf


Call NProbObs(ObsV,Fv)


!Calculate Theta Probability Compared to Prior
!Use different distributions for different estimated parameters
!If(Naxi.Le.NsnowVars+1)Then  !Snow parameter + SoilT
!    Call LogNProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)
!EndIf

!If(NaXi.Eq.6 .Or. Naxi.Eq.8  .Or. Naxi.Eq.10)Then  !6=Mv,8=P_M,10=P_SR
!    Call NProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)
!End If

!If(Naxi.Eq.7 .Or. Naxi.Eq.9)Then   !7=Sig,9=P_Q
!    Call UProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)
!End If

Call NProb(ThetaV,iSel,Mu,Cov,Pv,NiSel)

!If(Naxi.Eq.1 .And. NiSel.eq.1)Then
!    Print *,'ThetaV(iSel)',ThetaV(iSel),'Mu',Mu,'Cov',Cov,'Pv',Pv
!EndIf


If (UsePrior.Eq.0) Pv=Pu

R=Pv/Pu*Exp(-0.5d0*(Fv-Fu))


If (R.Gt.U) Then
  Na(NaXi)=Na(NaXi)+1
  ThetaU=ThetaV
  ObsU=ObsV
  Fu=Fv
  Pu=Pv
End If


End Subroutine Iterate
