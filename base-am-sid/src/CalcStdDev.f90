Subroutine CalcStdDev(X,N,StdDevX)

Implicit None
Integer,Intent(In) :: N
Integer :: i
Real(8),Intent(In) :: X(N)
Real(8),Intent(Out) :: StdDevX
Real(8) :: MuX

MuX=Sum(X)/Real(N,8)
StdDevX=0.0d0

Do i=1,N
   StdDevX=StdDevX+(dlog(X(i))-MuX)**2d0
End Do

StdDevX=Sqrt(StdDevX/Real(N-1,8))

End Subroutine CalcStdDev
