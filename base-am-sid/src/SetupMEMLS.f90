Subroutine SetupMEMLS

Use CommonVars
Implicit None
Integer :: i, i_start, i_end

Allocate(Aux_Soil(8),Aux_Model(3))

!Calculate intermediate parameters
Np=Np_passive+Np_active+Np_other


!Fixed number of mode variables, yes
NmodelVars=3
NVars=NsnowVars+NsoilVars+NmodelVars


!Give Fixed soil parameters
!Aux_Soil(4)=1500d0   !Soil bulk density in g/cm3
!Aux_Soil(5)=51.5d0 !Sand content of soil texture [%]
!Aux_Soil(6)=35.1d0 !Silt content of soil texture [%]
!Aux_Soil(7)=13.4d0 !Clay content of soil texture [%]
Aux_Soil(4)=1300d0   !Soil bulk density in g/cm3
Aux_Soil(5)=70.0d0 !Sand content of soil texture [%]
Aux_Soil(6)=29.0d0 !Silt content of soil texture [%]
Aux_Soil(7)=1.0d0 !Clay content of soil texture [%]


MvSMaxLim=1.0d0-1.5d0/2.65d0 !Change the maximum volumetric water content from 1 to porosity

!Define Std for observations
Allocate(StdObs(Np))
Do i=1,Np_passive
   StdObs(i)=StdTb
End Do

i_start=Np_passive+1
i_end=Np_passive+Np_active
Do i=i_start,i_end
   StdObs(i)=StdSigma
End Do

i_start=Np_passive+Np_active+1
i_end=Np
Do i=i_start,i_end
   StdObs(i)=StdOther(i-Np_passive-Np_active)
End Do

Print *, 'StdObs', StdObs


!Define the acceptance to achieve
Allocate(AceptGoal(NVars))
AceptGoal=0.3d0


End Subroutine SetupMEMLS
