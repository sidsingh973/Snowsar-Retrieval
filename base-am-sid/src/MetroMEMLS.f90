Program MetroMEMLS

Use CommonVars

Implicit None
Integer :: CountI,CountF
Integer :: idxPath,iFile

Character*150 FullPath
Character*150 RootFolder
Character*150 TempFile



Call System_Clock(CountI)

Print *, "MetroMEMLS version 3.0, 23 Sep 2016"


!Read filenames
Call getarg(0,FullPath)
idxPath=index(FullPath,'/',.true.)
RootFolder=FullPath(1:idxPath)
Print *, RootFolder
Open(81,file=trim(RootFolder)//'FILENAME.txt',status='old')

Do iFile=1,6
    Read(81,'(a)') fnm(iFile)
    TempFile=trim(RootFolder)//trim(fnm(iFile))
    Print *,TempFile
    fnm(iFile)=TempFile
End Do
Close(81)


Open(Unit=ObsFileUnit,File=fnm(2),Status='Old')
Open(Unit=AcptFileUnit,File=fnm(4),Status='Unknown')
Open(Unit=TbOutFileUnit,File=fnm(5),Status='Unknown',Form='Unformatted')
Open(Unit=ThetaOutFileUnit,File=fnm(6),Status='Unknown',Form='Unformatted')


Call ReadRunParams

Call SetupMEMLS



Do iPit=1,Npits

  Print *, 'Running pit #',iPit,'/',Npits


  Allocate(Obs(Nfreq,Nangle,Np))
  Allocate(ObsI(Nfreq*Nangle*Np),ObsJ(Nfreq*Nangle*Np),ObsK(Nfreq*Nangle*Np))
  Call ReadObs

  Open(Unit=HyperFileUnit,File=fnm(3),Status='Old')

  !Run different layer plan
  Do Nlyr=1,Nlyr_total

    Print *, 'Layer plan #',Nlyr

    !Calculate number of variables in Theta vector
    !Revised 18/5/15, add two more params: soil moisture and roughness
    Ntheta=Nlyr*NsnowVars+NsoilVars+NmodelVars

    Allocate(DzMu(Nlyr),DzCov(Nlyr),RhoMu(Nlyr),RhoCov(Nlyr),&
      DMu(Nlyr),DCov(Nlyr),TsnowMu(Nlyr),TsnowCov(Nlyr))
    Allocate(RhoMinLim(Nlyr),RhoMaxLim(Nlyr))
    Allocate(TsoilMu(1),TsoilCov(1),MvSMu(1),MvSCov(1),GndSigMu(1),GndSigCov(1))
    Allocate(P_MMu(1),P_MCov(1),P_QMu(1),P_QCov(1),P_SRMu(1),P_SRCov(1))
    Allocate(ThetaPost(Ntheta,Niter),ObsPost(Nobs,Niter))
    Allocate(Na(NVars),Acept(NVars))
    Allocate(ObsU(Nfreq,Nangle,Np),ObsV(Nfreq,Nangle,Np))

    Call ReadHyperPar

    Call MCMC

    Call WriteOutput

    Deallocate(DzMu,DzCov,RhoMu,RhoCov,DMu,DCov,TsnowMu,TsnowCov)
    Deallocate(MvSMu,MvSCov,GndSigMu,GndSigCov,TsoilMu,TsoilCov)
    Deallocate(RhoMinLim,RhoMaxLim)
    Deallocate(P_MMu,P_MCov,P_QMu,P_QCov,P_SRMu,P_SRCov)
    Deallocate(ThetaPost,ObsPost)
    Deallocate(Na,Acept)
    Deallocate(ObsU,ObsV)

  End Do

  Deallocate(Obs,ObsI,ObsJ,ObsK)
  Close(HyperFileUnit)
End Do


Close(ObsFileUnit)
Close(ThetaOutFileUnit)
Close(AcptFileUnit)


!Deallocate variables in ReadRunParams
Deallocate(StdOther,Freq,Angle,Pol_Passive,Pol_Active,Tsky)
!Dealllocate variables in SetupMEMLS
Deallocate(Aux_Soil,Aux_Model,StdObs,AceptGoal)


Call System_Clock(CountF)
Print *, 'Job took',Real(CountF-CountI,8)/1000.0d0/60.0d0,' minutes.'


End Program MetroMEMLS


