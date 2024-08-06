Module CommonVars

Implicit None



! A. Index of input-output files
Integer,Parameter :: ParFileUnit=1,ObsFileUnit=2,TrueFileUnit=3,&
  HyperFileUnit=4,ThetaOutFileUnit=7,AcptFileUnit=8,TbOutFileUnit=9


! B. Basic parameters from RunParams.txt
Integer :: Npits, Niter, Nburn, Nlyr_total, Nfreq, Nangle
Integer :: Np_passive, Np_active, Np_other
Integer :: NsnowVars, NsoilVars, NmodelVars !New variable
Integer :: ScatOpt, ModelOpt
Integer :: UsePrior,EstimateDz,EstimateRho,EstimateTsnow,EstimateD,EstimateSoil
Integer :: EstimateP_M, EstimateP_Q, EstimateP_SR !New variable
Integer :: InitialDz,InitialRho,InitialTsnow,InitialD,InitialSoil
Integer :: ConstrainRho, ConstrainTsnow
Real(8) :: StdTb,StdSigma     !Standard deviation of Observations
!Real(8) :: DzMinLim,DzMaxLim,RhoMinLim,RhoMaxLim,DMinLim,DMaxLim,TsnowMinLim,TsnowMaxLim !Min and Max limit of params
Real(8) :: DzMinLim,DzMaxLim,DMinLim,DMaxLim,TsnowMinLim,TsnowMaxLim !Min and Max limit of params
Real(8) :: MvSMinLim,MvSMaxLim,GndSigMinLim,GndSigMaxLim,TsoilMinLim,TsoilMaxLim
Real(8) :: P_MMinLim,P_MMaxLim,P_QMinLim,P_QMaxLim,P_SRMinLim,P_SRMaxLim
Real(8) :: RhoMinLim0,RhoMaxLim0

Real(8),Dimension(:),Allocatable :: RhoMinLim,RhoMaxLim


Integer,Dimension(:),Allocatable :: Pol_Passive, Pol_Active
Real(8),Dimension(:),Allocatable :: Freq,Angle
Real(8),Dimension(:),Allocatable :: StdOther !Standard deviation of Other Observations. (Np_other).
Real(8),Dimension(:,:),Allocatable :: Tsky   !Downward Tb at snow surface



! C. Parameters from Obs.txt and HyperPar.txt
Real(8),Dimension(:),Allocatable :: DzMu,DzCov,RhoMu,RhoCov,DMu,DCov,TsnowMu,TsnowCov  !Medium and Std of snow parameter priors
Real(8),Dimension(:),Allocatable :: TsoilMu,TsoilCov,MvSMu,MvSCov,GndSigMu,GndSigCov   !Medium and Std of soil parameter priors
Real(8),Dimension(:),Allocatable :: P_MMu,P_MCov,P_QMu,P_QCov,P_SRMu,P_SRCov
Real(8),Dimension(:,:,:),Allocatable :: Obs   !The Observations (Nfreq,Ntheta,Np_passive+Np_active+Np_other).


! D. Intermediate paramters
Integer :: iPit             !The index of current snowpit
Integer :: Ntheta           !Total number of parameters to be iterated (do MCMC running) for each layer plan. (NsnowVars*Nlyr+NsoilVars).
Integer :: NVars            !Total number of variables=NsnowVars+NsoilVars
Integer :: Np               !Np_passive+Np_active+Np_other
Integer :: Nobs             !Non-999 values in Obs (the observations)
Integer :: Nlyr             !Number of layers in snowpack, to be changed for different layer plans
Integer,Dimension(:),Allocatable :: iDz,iRho,iD,iTsnow,iTsoil,iMvS,iGndSig  !Array to save the location of parameters in the Theta Vector
Integer,Dimension(:),Allocatable :: iP_M,iP_Q,iP_SR
Integer,Dimension(:),Allocatable :: Na      !Array to save the acceptance of Each snow&soil parameter (NVars).
Integer,Dimension(:),Allocatable :: ObsI,ObsJ,ObsK   !i,j&k for each non-999 Observations (Nobs)


! E. Intermediate parameters for MCMC computation
Real(8) :: Fu,Fv  !Posterior probability of simulated observations
Real(8) :: PuDz,PuRho,PuD,PuTsnow   !Posterior probability of estimated parameters
Real(8) :: PuTsoil,PuMvS,PuGndSig
Real(8) :: PuP_M, PuP_Q, PuP_SR
Real(8),Dimension(:),Allocatable :: ThetaU, ThetaV    !Estimated variables (Ntheta).V is the updated theta, U is the previous theta.
Real(8),Dimension(:),Allocatable :: JmpDzStd,JmpRhoStd,JmpDStd,JmpTsnowStd  !Jump function of snow (Nlyr).
Real(8),Dimension(:),Allocatable :: JmpTsoilStd,JmpMvSStd,JmpGndSigStd
Real(8),Dimension(:),Allocatable :: JmpP_MStd,JmpP_QStd,JmpP_SRStd
Real(8),Dimension(:),Allocatable :: ObsTemp           !Used in ReadObs.f90 (Np)
Real(8),Dimension(:),Allocatable :: Aux_Soil          !Array to tansfer soil parameters. (7).
Real(8),Dimension(:),Allocatable :: Aux_Model          !Array to tansfer model parameters. (3).
Real(8),Dimension(:),Allocatable :: StdObs            !Std array for observations. (Np).
Real(8),Dimension(:),Allocatable :: Acept             !Save temporary Acept. (NVars).
Real(8),Dimension(:),Allocatable :: AceptGoal         !The acceptance to achieve (NVars).
Real(8),Dimension(:,:,:),Allocatable :: ObsU,ObsV     !Estimated Obs to be compared directly with Obs (Nfreq,Ntheta,Np_passive+Np_active+Np_other)

! Aux_Soil vector speficies 8 auxliliary soil inputs for observation model
! Aux_Soil(1): Soil Temperature [K], 274-Theta(Ntheta-2)
! Aux_Soil(2): Soil Volumbetric water content [frac], Theta(Ntheta-1)
! Aux_Soil(3): Ground roughness [m], Theta(Ntheta)
! Aux_Soil(4): Soil density [kg/m^3], fixed values for MEMLS run but could be revised in SetupMEMLS function
! Aux_Soil(5): Soil sand content [%],fixed values for MEMLS run but could be revised in SetupMEMLS function
! Aux_Soil(6): Soil silt content [%],fixed values for MEMLS run but could be revised in SetupMEMLS function
! Aux_Soil(7): Soil clay content [%],fixed values for MEMLS run but could be revised in SetupMEMLS function
! Aux_Model vector speficies 3 auxliliary observation model parameter
! Aux_Model(1): M
! Aux_Model(2): Q
! Aux_Model(3): SR

! F. Output arrays
Real(8),Dimension(:,:),Allocatable :: ThetaPost !Theta chains. (Ntheta,Niter).
Real(8),Dimension(:,:),Allocatable :: ObsPost    !Simulated observation chains. (Nobs,Niter).



! G. Strings to save filenames
Character*150 fnm(7)




End Module



