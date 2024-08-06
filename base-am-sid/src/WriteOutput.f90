Subroutine WriteOutput

Use CommonVars
Implicit None
Integer :: i,j


Do i=1,Ntheta
  Write(ThetaOutFileUnit) (ThetaPost(i,j),j=1,Niter)
End Do

Do i=1,Nobs
   Write(TbOutFileUnit) (ObsPost(i,j),j=1,Niter)
End Do

Write(AcptFileUnit,'(I5,1X,I5)') iPit,Nlyr
Do i=1,NVars
  Write(AcptFileUnit,'(F5.3)') Real(Na(i))/Real(Niter)
EndDo

End Subroutine WriteOutput
