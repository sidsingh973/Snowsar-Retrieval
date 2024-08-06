Subroutine ReadRunParams

Use CommonVars
Implicit None
Integer :: i,j
Real(8) :: TCMinLim,TCMaxLim
Real(8) :: TskyTemp

Open(Unit=ParFileUnit,File=fnm(1),Status='Unknown')


! Number of pits to run (Npits)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') Npits
Print *, 'Npits=', Npits

! Number of iterations in the Markov Chain (Niter)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I12)') Niter
Print *, 'Niter=', Niter

! Number of burn-in iterations in the Markov Chain (Nburn)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I12)') Nburn
Print *, 'Nburn=', Nburn

! Number of snow layers to predict (Nlyr_total)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') Nlyr_total
Print *, 'Nlyr_total=', Nlyr_total

! Number of observation frequencies (Nfreq)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') Nfreq
Print *, 'Nfreq=', Nfreq

! Number of observation angles (Nangle)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') Nangle
Print *, 'Nangle=', Nangle

! Number of polarizations of passive Tb (Np_passive)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') Np_passive
Print *, 'Np_passive=', Np_passive

! Number of polarizations of active backscattering coefficient (Np_active)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') Np_active
Print *, 'Np_active=', Np_active

! Number of other observations (Np_other)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') Np_other
Print *, 'Np_other=', Np_other

! Number of Theta variables (NsnowVars)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') NsnowVars
Print *, 'NsnowVars=', NsnowVars

! Number of soil parameters (NsoilVars)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') NsoilVars
Print *, 'NsoilVars=', NsoilVars

! Scattering Coefficient Option (ScatOpt): 1 (Empirical MEMLS, Hallikainen-HUT), 2 (MEMLS-IBA, combined HUT), or 3 (Roy-HUT)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') ScatOpt
Print *, 'ScatOpt=', ScatOpt

! Observation Model Option (ModelOpt): 1 (MEMLS), 2 (HUT), 3 (DMRT-ML), 4 (DMRT-QMS), 5 (Bi-Continuous)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') ModelOpt
Print *, 'ModelOpt=', ModelOpt


! Error standard deviation of Tb observations (StdTb)
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') StdTb
Print *, 'StdTb=', StdTb

! Error standard deviation of Backscattering coefficient observations (StdSigma)
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') StdSigma
Print *, 'StdSigma=', StdSigma

! Error standard deviation of Tb observations (StdOther)
Allocate(StdOther(Np_other))
Read(ParFileUnit,*)
Do i=1,Np_other
    Read(ParFileUnit,'(F10.4)') StdOther(i)
    Print *, 'StdOther=', StdOther
End Do

! Observation frequencies (Freq(Nfreq))
Allocate(Freq(Nfreq))
Read(ParFileUnit,*)
Do i=1,Nfreq
  Read(ParFileUnit,'(F10.4)') Freq(i)
End Do
Print *, 'Freq=', Freq

! Observation angles (Angle), Nf lines
Allocate(Angle(Nangle))
Read(ParFileUnit,*)
Do i=1,Nangle
  Read(ParFileUnit,'(F10.4)') Angle(i)
End Do
Print *, 'Angle=', Angle

!Polarizations for Passive Measurement(Pol_Passive(Np_passive)): 1 (vertical), 2 (horizontal)
Allocate(Pol_Passive(Np_passive))
Read(ParFileUnit,*)
Do i=1,Np_passive
   Read(ParFileUnit,'(I5)') Pol_Passive(i)
End Do
Print *, 'Pol_Passive=',Pol_Passive

! Polarizations for Active Measurement(Pol_Active(Np_active)): 1 (vv), 2 (hh), 3 (vh), 4 (hv)
Allocate(Pol_Active(Np_active))
Read(ParFileUnit,*)
Do i=1,Np_active
  Read(ParFileUnit,'(I5)') Pol_Active(i)
End Do
Print *, 'Pol_Active=',Pol_Active

! Tb boundary condition above snow surface at vertical and horizontal polarizations (K) (Tsky(Nfreq,2))
Allocate(Tsky(Nfreq,Nangle))
Read(ParFileUnit,*)
Do i=1,Nfreq
  Do j=1,Nangle
    Read(ParFileUnit,'(F10.4)') TskyTemp
    Tsky(i,j)=TskyTemp
   EndDo
End Do
Print *, 'Tsky=', Tsky


! Use prior information, 0 or 1. (UsePrior)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') UsePrior
Print *, 'UsePrior=', UsePrior

! Estimate dZ, -1, 0, 1  (EstimateDz)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') EstimateDz
Print *, 'EstimateDz=', EstimateDz

! Estimate rho, -1, 0, 1  (EstimateRho)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') EstimateRho
Print *, 'EstimateRho=', EstimateRho

! Estimate D, -1, 0, 1  (EstimateD)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') EstimateD
Print *, 'EstimateD=', EstimateD

! Estimate T, -1, 0, 1  (EstimateTsnow)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') EstimateTsnow
Print *, 'EstimateTsnow=', EstimateTsnow

! Estimate Soil Parameter (EstimateSoil)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') EstimateSoil
Print *, 'EstimateSoil=', EstimateSoil

! Estimate P_M, -1,0,1 (EstimateP_M)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') EstimateP_M
Print *, 'EstimateP_M=', EstimateP_M

! Estimate P_Q, -1,0,1 (EstimateP_Q)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') EstimateP_Q
Print *, 'EstimateP_Q=', EstimateP_Q

! Estimate P_SR, -1,0,1 (EstimateP_SR)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') EstimateP_SR
Print *, 'EstimateP_SR=', EstimateP_SR


! Initial Value dZ, -1, 0, 1  (InitialDz)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') InitialDz
Print *, 'InitialDz=', InitialDz

! Initial Value rho, -1, 0, 1  (InitialRho)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') InitialRho
Print *, 'InitialRho=', InitialRho


! Initial Value D, -1, 0, 1  (InitialD)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') InitialD
Print *, 'InitialD=', InitialD

! Initial Value T, -1, 0, 1  (InitialTsnow)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') InitialTsnow
Print *, 'InitialTsnow=', InitialTsnow

! Initial Value S, -1, 0, 1  (InitialSoil)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') InitialSoil
Print *, 'InitialSoil=', InitialSoil


! Constrain rho (ConstrainRho): 0 (No constraint), 1 (surface<bottom), 2 (surface>bottom)Read(1,*)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') ConstrainRho
Print *, 'ConstrainRho=', ConstrainRho

! Constrain rho (ConstrainT): 0 (No constraint), 1 (surface<bottom), 2 (surface>bottom)Read(1,*)
Read(ParFileUnit,*)
Read(ParFileUnit,'(I5)') ConstrainTsnow
Print *, 'ConstrainTsnow=', ConstrainTsnow

! Minimum & maximum limits for layer thickness [m]
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') DzMinLim
Read(ParFileUnit,'(F10.4)') DzMaxLim
Print *, 'DzMinLim=', DzMinLim, '(m)'
Print *, 'DzMaxLim=', DzMaxLim, '(m)'

! Minimum & maximum limits for density [kg/m3]
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') RhoMinLim0
Read(ParFileUnit,'(F10.4)') RhoMaxLim0
Print *, 'RhoMinLim=', RhoMinLim0, '(kg/m3)'
Print *, 'RhoMaxLim=', RhoMaxLim0, '(kg/m3)'


! Minimum & maximum limits for grain diameter [mm]
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') DMinLim
Read(ParFileUnit,'(F10.4)') DMaxLim
Print *, 'DMinLim=', DMinLim, '(mm)'
Print *, 'DMaxLim=', DMaxLim, '(mm)'

! Minimum & maximum limits for snow temperature [degC]
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') TCMinLim
Read(ParFileUnit,'(F10.4)') TCMaxLim

TsnowMaxLim=274.0d0-(TCMinLim+273.15d0) !convert the temperature limit from K to 274-K used in the MCMC
TsnowMinLim=274.0d0-(TCMaxLim+273.15d0) !convert the temperature limit from K to 274-K used in the MCMC
Print *, 'TsnowMinLim=', TsnowMinLim, '(K)'
Print *, 'TsnowMaxLim=', TsnowMaxLim, '(K)'

! Minimum & maximum limits for soil moisture [frac]
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') MvSMinLim
Read(ParFileUnit,'(F10.4)') MvSMaxLim
Print *, 'MvSMinLim=', MvSMinLim, '(frac)'
Print *, 'MvSMaxLim=', MvSMaxLim, '(frac)'

! Minimum & maximum limits for soil rms-height [m]
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') GndSigMinLim
Read(ParFileUnit,'(F10.4)') GndSigMaxLim
Print *, 'GndSigMinLim=', GndSigMinLim, '(m)'
Print *, 'GndSigMaxLim=', GndSigMaxLim, '(m)'

! Minimum & maximum limits for soil temperature [degC]
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') TCMinLim
Read(ParFileUnit,'(F10.4)') TCMaxLim
TsoilMaxLim=279.0d0-(TCMinLim+273.15d0) !convert the temperature limit from K to 279-K used in the MCMC
TsoilMinLim=279.0d0-(TCMaxLim+273.15d0) !convert the temperature limit from K to 279-K used in the MCMC
Print *, 'TsoilMinLim=', TsoilMinLim, '(K)'
Print *, 'TsoilMaxLim=', TsoilMaxLim, '(K)'


! Minimum & Maximum limits for model parameter, P_M
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') P_MMinLim
Read(ParFileUnit,'(F10.4)') P_MMaxLim
Print *, 'P_MMinLim=', P_MMinLim
Print *, 'P_MMaxLim=', P_MMaxLim

! Minimum & Maximum limits for model parameter, P_Q
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') P_QMinLim
Read(ParFileUnit,'(F10.4)') P_QMaxLim
Print *, 'P_QMinLim=', P_QMinLim
Print *, 'P_QMaxLim=', P_QMaxLim


! Minimum & Maximum limits for model parameter, P_SR
Read(ParFileUnit,*)
Read(ParFileUnit,'(F10.4)') P_SRMinLim
Read(ParFileUnit,'(F10.4)') P_SRMaxLim
Print *, 'P_SRMinLim=', P_SRMinLim
Print *, 'P_SRMaxLim=', P_SRMaxLim


Close(ParFileUnit)


End Subroutine ReadRunParams
