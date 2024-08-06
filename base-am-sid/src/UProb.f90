Subroutine UProb(Theta,iSel,Mu,Cov,P,NiSel)

Use CommonVars
Implicit None
Integer,Intent(In) :: iSel(NiSel),NiSel
Integer :: i
Real(8),Intent(In) :: Theta(Ntheta),Mu(Nlyr),Cov(Nlyr)
Real(8),Intent(Out) :: P
Real(8) :: LowValue,HighValue

P=0.0

!Multiply sqrt(2*pi) to compensate for LogNProb & NProb
Do i=1,NiSel
    LowValue=Mu(i)-Sqrt(Cov(i))
    HighValue=Mu(i)+Sqrt(Cov(i))

    If(Theta(iSel(i)).lt.HighValue .And. Theta(iSel(i)).gt.LowValue)Then
        P=P+1.0/(HighValue-LowValue)
    End If
End Do


End Subroutine UProb
