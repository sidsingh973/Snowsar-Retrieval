Subroutine ObsModel(Theta,ObsS)

Use CommonVars
Implicit None

Integer :: i
Real(8) :: Profile(Nlyr,5)
Real(8),Intent(In) :: Theta(Ntheta)
Real(8),Intent(Out) :: ObsS(Nfreq,Nangle,Np)


!Make soil property vector
!Aux_Soil(1)=279.0d0-dexp(Theta(iTsoil(1)))   !Soil temperature [K]
!Aux_Soil(2)=dexp(Theta(iMvS(1)))             !Soil water content [FRAC]
!Aux_Soil(3)=dexp(Theta(iGndSig(1)))          !Ground roughness [m]
Aux_Soil(1)=279.0d0-Theta(iTsoil(1))   !Soil temperature [K]
Aux_Soil(2)=Theta(iMvS(1))             !Soil water content [FRAC]
Aux_Soil(3)=Theta(iGndSig(1))          !Ground roughness [m]


!Make model parameter vector
!Aux_Model(1)=dexp(Theta(iP_M(1)))
!Aux_Model(2)=dexp(Theta(iP_Q(1)))
!Aux_Model(3)=dexp(Theta(iP_SR(1)))
Aux_Model(1)=Theta(iP_M(1))
Aux_Model(2)=Theta(iP_Q(1))
Aux_Model(3)=Theta(iP_SR(1))


Do i=1,Nlyr
  !Profile(i,1)=dexp(Theta(iDz(i)))  !Layer thickness [m]
  !Profile(i,2)=dexp(Theta(iRho(i))) !Density [kg/m3]
  !Profile(i,3)=dexp(Theta(iD(i)))   !Exponential correlation length [mm]
  Profile(i,1)=Theta(iDz(i))  !Layer thickness [m]
  Profile(i,2)=Theta(iRho(i)) !Density [kg/m3]
  Profile(i,3)=Theta(iD(i))   !Exponential correlation length [mm]
  Profile(i,4)=0.0d0                !Liquid fraction [m3/m3]
  !Profile(i,5)=274.0d0-dexp(Theta(iTsnow(i))) !Snow temperature [K]
  Profile(i,5)=274.0d0-Theta(iTsnow(i)) !Snow temperature [K]
End Do


!Jinmei, 2017/6/30, change dz to relative thickness for layer 2 to Nlyr
Do i=2,Nlyr
    Profile(i,1)=Profile(i,1)*Profile(1,1)
End Do


If(ModelOpt.Eq.1)Then

 If(.True.)Then
    !Revised Sep23, 2016, if EstimateP_Q==0, Use model to estimate Q (recommended)
    Call SS_MEMLS1(Freq,Angle,Profile,Aux_Soil,Aux_Model,Tsky,&
        Pol_Passive,Pol_Active,&
        Nfreq,Nangle,Np_passive,Np_active,Np,&
        Nlyr,ScatOpt,EstimateP_Q,Obs,ObsS)
 Else
    ObsS=Obs
 EndIf

EndIf



!print *,'worked inside ObsModel before sending out~'
End Subroutine ObsModel
