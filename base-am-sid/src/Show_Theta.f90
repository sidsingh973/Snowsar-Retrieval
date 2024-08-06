! -------------------------------------------------------------------------
!
SUBROUTINE Show_Theta(Theta)
!
! -------------------------------------------------------------------------

Use CommonVars
Implicit None

Real(8),Intent(In) :: Theta(Ntheta)

Print *,'Snow Params (Nlyr=',Nlyr,'): ', Theta(1:4*Nlyr)
Print *,'Soil Params:',Theta(4*Nlyr+1:4*Nlyr+3)
Print *,'Model Params:',Theta(4*Nlyr+NsoilVars+1:4*Nlyr+NsoilVars+3)

END SUBROUTINE Show_Theta
