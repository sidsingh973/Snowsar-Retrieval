Subroutine RandUniform(Num,U,SeedNum)

Implicit None
Integer,Intent(In) :: Num,SeedNum
Integer,Dimension(:),Allocatable :: Seed
Integer :: Rsz,i
Real(8),Intent(Out) :: U(Num)

Call Random_Seed(Size=Rsz)

Allocate(Seed(Rsz))
Do i=1,Rsz
  Seed(i)=Seednum+i
End Do

Call Random_Seed(Put=Seed(1:Rsz))

Call Random_Number(U)

Deallocate(Seed)

End Subroutine RandUniform

!Random_Number: Returns a single pseudorandom number or an array of 
!pseudorandom numbers from the uniform distribution over the range 0 < x < 1.
