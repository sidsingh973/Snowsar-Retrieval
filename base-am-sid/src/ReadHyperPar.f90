Subroutine ReadHyperPar
! HyperPar revised to include two soil parameters, soil moisture [FRAC] and soil roughness rms-height [m]

Use CommonVars
Implicit None
Integer :: i
Real(8) :: RhoNMean,RhoNStd

Read(HyperFileUnit,*) !Skip for number of layers


! Read layer thickness prior, from snow bottom to snow surface
Read(HyperFileUnit,*)
Do i=1,Nlyr
  Read(HyperFileUnit,'(F12.8)') DzMu(i)
End Do

Do i=1,Nlyr
  Read(HyperFileUnit,'(F12.8)') DzCov(i)
End Do
Print *, 'Priors:'
Print *, 'DzMu=',DzMu,', DzCov=',DzCov


! Read snow grain size prior, from snow bottom to snow surface
Read(HyperFileUnit,*)
Do i=1,Nlyr
  Read(HyperFileUnit,'(F12.8)') DMu(i)
End Do
Do i=1,Nlyr
  Read(HyperFileUnit,'(F12.8)') DCov(i)
End Do
Print *, 'DMu=',DMu,', DCov=',DCov


! Read snow density prior, from snow bottom to snow surface
Read(HyperFileUnit,*)
Do i=1,Nlyr
  Read(HyperFileUnit,'(F12.8)') RhoMu(i)
End Do
Do i=1,Nlyr
  Read(HyperFileUnit,'(F12.8)') RhoCov(i)
End Do
Print *, 'RhoMu=',RhoMu,', RhoCov=',RhoCov


! Added Sep 23, 2016, Cut-off density at 2 std around mean(logX)
!Do i=1,Nlyr
!!  !RhoMaxLim(i)=Exp(RhoMu(i)+2.5*RhoCov(i))
!!  !RhoMinLim(i)=Exp(RhoMu(i)-1.5*RhoCov(i))
!  RhoNMean = exp(RhoMu(i) + RhoCov(i)/2)
!  RhoNStd = sqrt(exp(2*RhoMu(i) + RhoCov(i))*(exp(RhoCov(i))-1))
!  RhoMaxLim(i)=RhoNMean+3*RhoNStd
!  RhoMinLim(i)=RhoNMean-3*RhoNStd
!End Do
!Print *, 'RhoNMean',RhoNMean
!Print *, 'RhoNStd',RhoNStd

Do i=1,Nlyr
    RhoMaxLim(i)=RhoMaxLim0
    RhoMinLim(i)=RhoMinLim0
End Do
Print *, 'RhoMaxLim',RhoMaxLim
Print *, 'RhoMinLim',RhoMinLim



! Read snow temperature prior, from snow bottom to snow surface
Read(HyperFileUnit,*)
Do i=1,Nlyr
  Read(HyperFileUnit,'(F12.8)') TsnowMu(i)
End Do
Do i=1,Nlyr
  Read(HyperFileUnit,'(F12.8)') TsnowCov(i)
End Do
Print *, 'TsnowMu=',TsnowMu,', TsnowCov=',TsnowCov


! Read soil temprarature prior
Read(HyperFileUnit,*)
Do i=1,1
  Read(HyperFileUnit,'(F12.8)') TsoilMu(i)
  Read(HyperFileUnit,'(F12.8)') TsoilCov(i)
End Do
Print *, 'TsoilMu=',TsoilMu,', TsoilCov=',TsoilCov

! Read soil moisture prior
Read(HyperFileUnit,*)
Do i=1,1
    Read(HyperFileUnit,'(F12.8)') MvSMu(i)
    Read(HyperFileUnit,'(F12.8)') MvSCov(i)
End Do
Print *, 'MvSMu=',MvSMu,', MvSCov=',MvSCov

! Read soil roughness rms-height prior
Read(HyperFileUnit,*)
Do i=1,1
    Read(HyperFileUnit,'(F12.8)') GndSigMu(i)
    Read(HyperFileUnit,'(F12.8)') GndSigCov(i)
End Do
Print *, 'GndSigMu=',GndSigMu,', GndSigCov=',GndSigCov

! Added Nov 11, 2016, Cut-off soil roughness at 2 std around mean(X)
!GndSigMaxLim=Exp(GndSigMu(1)+2.5*GndSigCov(1))
!GndSigMinLim=Exp(GndSigMu(1)-1.5*GndSigCov(1))


! Read prior for m for active MEMLS
Read(HyperFileUnit,*)
Do i=1,1
    Read(HyperFileUnit,'(F12.8)') P_MMu(i)
    Read(HyperFileUnit,'(F12.8)') P_MCov(i)
End Do
Print *, 'P_MMu=',P_MMu,', P_MCov=',P_MCov


! Read prior for q for active MEMLS
Read(HyperFileUnit,*)
Do i=1,1
    Read(HyperFileUnit,'(F12.8)') P_QMu(i)
    Read(HyperFileUnit,'(F12.8)') P_QCov(i)
End Do
Print *, 'P_QMu=',P_QMu,', P_QCov=',P_QCov


! Read prior for adjusting ratio for specular part of reflectivity, active MEMLS
Read(HyperFileUnit,*)
Do i=1,1
    Read(HyperFileUnit,'(F12.8)') P_SRMu(i)
    Read(HyperFileUnit,'(F12.8)') P_SRCov(i)
End Do
Print *, 'P_SRMu=',P_SRMu,', P_SRCov=',P_SRCov


End Subroutine ReadHyperPar
