! Read Observations from tb_obs.txt file.
! The observations contains Nfreq*Nangle rows, Np lines.

Subroutine ReadObs

Use CommonVars
Implicit None
Integer :: i,j,k
Real(8) Temp


Read(ObsFileUnit,*)
Allocate(ObsTemp(Np))
Print *, 'Observations'

Obs=0.0d0
ObsI=0
ObsJ=0
ObsK=0

Nobs=0

Do i=1,Nfreq
  Do j=1,Nangle
    Read(ObsFileUnit,'(10F12.4)') ObsTemp
    Print *,ObsTemp
    Do k=1,Np
       Temp=ObsTemp(k)
       Obs(i,j,k)=Temp
       If(Temp<900)Then
          Nobs=Nobs+1
          ObsI(Nobs)=i
          ObsJ(Nobs)=j
          ObsK(Nobs)=k
       End If
    End Do
  End Do
End Do

Deallocate(ObsTemp)

Print *, 'Obs', Obs


End Subroutine ReadObs
