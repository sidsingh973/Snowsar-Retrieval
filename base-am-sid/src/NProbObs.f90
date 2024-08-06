!When Normal distribution was applied over the Observation vector,
!Exp is not calculated, because it could produce many zero prob due to bad simulations at
!a single band

Subroutine NProbObs(ObsS,F)

Use CommonVars
Implicit None
Integer :: i,j,k
Real(8),Intent(In) :: ObsS(Nfreq,Nangle,Np)
Real(8),Intent(Out) :: F
Real(8) :: Temp


!Probability of normal distribution. Here the constant multification factor 1/sqrt(2*pi) is neglected.
F=0.0d0

Do i=1,Nfreq
  Do j=1,Nangle
    Do k=1,Np
       Temp=Obs(i,j,k)
       If(Temp<900)Then
          F=F+(ObsS(i,j,k)-Obs(i,j,k))**2.0d0/(StdObs(k)**2.0d0)
       End If
    End Do
  End Do
End Do


!!!! Commented out !!!
!values for observations is large
!F=Exp(-0.5d0*F)

!Do k=1,Np
!   F=F/StdObs(k))
!End Do
!!!!!!!!!!!!!!!!!!!!!!




End Subroutine NProbObs
