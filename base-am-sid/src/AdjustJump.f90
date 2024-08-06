!Adjust the jump function at each 100 step in the burn-in period

Subroutine AdjustJump(iIter)

Use CommonVars
Implicit None
Integer,Intent(In) :: iIter
Integer :: i, AddAdjust
Real(8) :: DzStd(Nlyr),RhoStd(Nlyr),DStd(Nlyr),TsnowStd(Nlyr)
Real(8) :: TsoilStd,MvSStd,GndSigStd
Real(8) :: P_MStd,P_QStd,P_SRStd
Real(8) :: JmpLimit
Print *, Na
!Cacluate Standard deviation
Do i=1,Nlyr
  DzStd(i)=1000000d0
  RhoStd(i)=1000000d0
  DStd(i)=1000000d0
  TsnowStd(i)=1000000d0
End Do
TsoilStd=1000000d0
MvSStd=1000000d0
GndSigStd=1000000d0
P_MStd=1000000d0
P_QStd=1000000d0
P_SRStd=1000000d0
!Calcuate current iteration

Do i=1,NVars
  Acept(i)=Real(Na(i),8)/Real(iIter,8)  !Convert to Real(8)
End Do


!Adjust jump function
!Every 100th iteration, I invoke the proportionality-based adjustment several of us have used:  
! (new jump standard deviation) / (current jump standard deviation) = (current acceptance rate) /
! (desired acceptance rate)
!At the fist 100th iteration, use the standard deviation first


Print *,'#Iter - ',iIter
Print *,'Number of Acceptance:'
Print *, Na
Print *, 'Acceptance Rate'
Print *, 'For Dz, Rho, D, T, Tsoil, Mv, GndSig, M, Q, SR'
Print *, Acept
Print *, '    Jump Steps Before constraining:'
Print *,'Jump dz (cm):',JmpDzStd
Print *,'Jump rho (kg/m^3):',JmpRhoStd
Print *,'Jump D (mm):',JmpDStd
Print *,'Jump Tsnow (degC):',JmpTsnowStd
Print *,'Jump Tsoil (degC):',JmpTsoilStd
Print *,'Jump MvS (%):',JmpMvSStd
Print *,'Jump GndSig (mm):',JmpGndSigStd
Print *,'Jump M:',JmpP_MStd
Print *,'Jump Q:',JmpP_QStd
Print *,'Jump SR:',JmpP_SRStd
Print *,'MvSStd:',MvSStd

If(iIter.eq.100)Then

  Do i=1,Nlyr
      
      JmpDzStd(i)=min(DzStd(i), 4.0d0*dsqrt(DzCov(i)))
      JmpRhoStd(i)=min(RhoStd(i), 4.0d0*dsqrt(RhoCov(i)))
      JmpDStd(i)=min(DStd(i), 4.0d0*dsqrt(DCov(i)))
      JmpTsnowStd(i)=min(TsnowStd(i), 4.0d0*dsqrt(TsnowCov(i)))
  End Do

  JmpTsoilStd(1)=min(TsoilStd, 4.0d0*dsqrt(TsoilCov(1)))
  JmpMvSStd(1)=min(MvSStd, 4.0d0*dsqrt(MvSCov(1)))
  JmpGndSigStd(1)=min(GndSigStd, 4.0d0*dsqrt(GndSigCov(1)))
  JmpP_MStd(1)=min(P_MStd, 4.0d0*dsqrt(P_MCov(1)))
  JmpP_QStd(1)=min(P_QStd, 4.0d0*dsqrt(P_QCov(1)))
  JmpP_SRStd(1)=min(P_SRStd, 4.0d0*dsqrt(P_SRCov(1)))

Else

  !Add code, the jump function can not be larger than half of the Max-Min range
  Do i=1,Nlyr
     !JmpLimit=dlog(DzMaxLim)-dlog(DzMinLim)
     JmpLimit=DzMaxLim-DzMinLim
     JmpDzStd(i)= min(Acept(1)/AceptGoal(1) * JmpDzStd(i), dabs(JmpLimit)/2.0d0)
     JmpDzStd(i)= min(JmpDzStd(i),4.0d0*dsqrt(DzCov(i)))

     !JmpLimit=dlog(RhoMaxLim0)-dlog(RhoMinLim0)
     Print *,'	   Rhostd:',dsqrt(RhoCov(i))
     Print *,'Jump Jmplim:',JmpLimit
     JmpLimit=(RhoMaxLim0)-(RhoMinLim0)
     JmpRhoStd(i)= min(Acept(2)/AceptGoal(2) * JmpRhoStd(i), dabs(JmpLimit)/2.0d0)
     Print *,'Jump Rhostd:',JmpRhoStd(i)
     JmpRhoStd(i)= min(JmpRhoStd(i),4.0d0*dsqrt(RhoCov(i)))
     Print *,'Jump Rhostd:',JmpRhoStd(i)
     
     !JmpLimit=dlog(DMaxLim)-dlog(DMinLim)
     JmpLimit=(DMaxLim)-(DMinLim)
     JmpDStd(i)= min(Acept(3)/AceptGoal(3) * JmpDStd(i), dabs(JmpLimit)/2.0d0)
     JmpDStd(i)= min(JmpDStd(i),4.0d0*dsqrt(DCov(i)))

     !JmpLimit=dlog(TsnowMaxLim)-dlog(TsnowMinLim)
     JmpLimit=(TsnowMaxLim)-(TsnowMinLim)
     JmpTsnowStd(i)= min(Acept(4)/AceptGoal(4) * JmpTsnowStd(i), dabs(JmpLimit)/2.0d0)
     JmpTsnowStd(i)= min(JmpTsnowStd(i),4.0d0*dsqrt(TsnowCov(i)))
  End Do

  !JmpLimit=dlog(TsoilMaxLim)-dlog(TsoilMinLim)
  JmpLimit=(TsoilMaxLim)-(TsoilMinLim)
  JmpTsoilStd(1)=min(Acept(5)/AceptGoal(5) * JmpTsoilStd(1), dabs(JmpLimit)/2.0d0)
  JmpTsoilStd(1)= min(JmpTsoilStd(1),4.0d0*dsqrt(TsoilCov(1)))

  !JmpLimit=dlog(MvSMaxLim)-dlog(MvSMinLim)
  JmpLimit=(MvSMaxLim)-(MvSMinLim)
  JmpMvSStd(1)=min(Acept(6)/AceptGoal(6) * JmpMvSStd(1), dabs(JmpLimit)/2.0d0)
  JmpMvSStd(1)= min(JmpMvSStd(1),4.0d0*dsqrt(MvSCov(1)))

  !JmpLimit=dlog(GndSigMaxLim)-dlog(GndSigMinLim)
  JmpLimit=(GndSigMaxLim)-(GndSigMinLim)
  JmpGndSigStd(1)=min(Acept(7)/AceptGoal(7) * JmpGndSigStd(1), dabs(JmpLimit)/2.0d0)
  JmpGndSigStd(1)= min(JmpGndSigStd(1),4.0d0*dsqrt(GndSigCov(1)))

  !JmpLimit=dlog(P_MMaxLim)-dlog(P_MMinLim)
  JmpLimit=(P_MMaxLim)-(P_MMinLim)
  JmpP_MStd(1)=min(Acept(8)/AceptGoal(8) * JmpP_MStd(1), dabs(JmpLimit)/2.0d0)
  JmpP_MStd(1)= min(JmpP_MStd(1),4.0d0*dsqrt(P_MCov(1)))

  !JmpLimit=dlog(P_QMaxLim)-dlog(P_QMinLim)
  JmpLimit=(P_QMaxLim)-(P_QMinLim)
  JmpP_QStd(1)=min(Acept(9)/AceptGoal(9) * JmpP_QStd(1), dabs(JmpLimit)/2.0d0)
  JmpP_QStd(1)= min(JmpP_QStd(1),4.0d0*dsqrt(P_QCov(1)))

  !JmpLimit=dlog(P_SRMaxLim)-dlog(P_SRMinLim)
  JmpLimit=(P_SRMaxLim)-(P_SRMinLim)
  JmpP_SRStd(1)=min(Acept(10)/AceptGoal(10) * JmpP_SRStd(1), dabs(JmpLimit)/2.0d0)
  JmpP_SRStd(1)= min(JmpP_SRStd(1),4.0d0*dsqrt(P_SRCov(1)))
EndIf



Print *, '    Jumpping Steps After Adjustment on log(Theta):'
Print *,'Adjusted jump dz (cm):',JmpDzStd
Print *,'Adjusted jump rho (kg/m^3):',JmpRhoStd
Print *,'Adjusted jump D (mm):',JmpDStd
Print *,'Adjusted jump Tsnow (degC):',JmpTsnowStd
Print *,'Adjusted jump Tsoil (degC):',JmpTsoilStd
Print *,'Adjusted jump MvS (%):',JmpMvSStd
Print *,'Adjusted jump GndSig (mm):',JmpGndSigStd
Print *,'Adjusted jump M:',JmpP_MStd
Print *,'Adjusted jump Q:',JmpP_QStd
Print *,'Adjusted jump SR:',JmpP_SRStd


End Subroutine AdjustJump

