classdef MCMCRun4
    
properties


%basic folder info
folder


%MCMC set-up parameters
Npits       %number of snwpits
Niter       %number of iterations
Nburn       %number of burn-in period
NlyrMax     %maximum number of layers
ModelOpt   %observation model
ScatOpt     %scattering coefficient parameterization
UsePrior    %use prior or not
SturmClass  %snow type
Lyrplan     %layerplan option, 1=by algorithm default, 2=choose Nlyr_choose, 3=Two layer at least
Nlyr_choose %constrained number of layers when Lyrplan=2
Run_location %if on osc super computer (=1), or no local desktop(=0)

%basic snowpit information
sp


%sizes
Nf
Nangle  %number of observation angles
Np
Ntheta  %number of variables in theta_post


%% inputs: measurement inputs
active_freq
active_theta
active_pol
active_np

%
passive_freq
passive_theta
passive_pol
passive_np

other_np
other_measurements

%measurements
ObsMCMCIn
Nobs
Obs      %Obs(Run.Nobs,Run.Npits)

%accuracy
stdTb
stdSigma
stdOther
StdObs   %StdObs(Run.Nobs,Run.Npits);
skytb

%combined MCMC freq&theta, including all freq&theta for active and passive
MCMC_freq
MCMC_theta
MCMC_pol  %this is used to save pol of different measurements

%priors, read from hyperpar
DzMu      %DzMu{Nlyr}
DzCov
DMu
DCov
RhoMu
RhoCov
TsnowMu
TsnowCov
TsoilMu
TsoilCov
MvSMu
MvSCov
GndSigMu
GndSigCov
P_MMu
P_MCov
P_QMu
P_QCov
P_SRMu
P_SRCov


%% intermediate(processed) inputs
%converted to Mean and Std
DzMean
DzStd
DMean
DStd
RhoMean
RhoStd
TsnowMean
TsnowStd
TsoilMean
TsoilStd
MvSMean
MvSStd
GndSigMean
GndSigStd
P_MMean
P_MStd
P_QMean
P_QStd
P_SRMean
P_SRStd

%summarized prior
swePr  %priorSWE
sdPr   %priorSd


%% outputs: the general experiment output
theta_post      %theta_post{Nlyr}(Run.Niter, Ntheta(ilyr), Run.Npits)
acceptance
ObsPost         %ObsPost(Run.Niter,Run.Nobs,Run.NlyrMax,Run.Npits)

%chains for all parameters
chains_dz       %chains_dz{Nlyr,Npits}(Niter,Nlyr)
chains_D
chains_rho
chains_Tsnow
chains_Tsoil
chains_mvs
chains_sig
chains_sd
chains_swe
chains_P_M
chains_P_Q
chains_P_SR

%mean of the chains for each layer after Burn-in. Dimensions: {Nlyr,Npits}(Nlyr)
md_dz
md_D
md_rho
md_Tsnow
md_Tsoil
md_mvs
md_sig
md_sd
md_swe
md_P_M
md_P_Q
md_P_SR

%standard derivation of the chains after Burn-in: std{Nlyr,Npits}(Nlyr)
%this standard derivation is not directly taking the std, but the std error
%around the mean above.
std_dz
std_D
std_rho
std_Tsnow
std_Tsoil
std_mvs
std_sig
std_sd
std_swe
std_P_M
std_P_Q
std_P_SR


%% sumarized output
nHat    %number of snow layers
sdHat   %snowdepth
sweHat   %SWE
sdStdHat
sweStdHat  

%
rhoavgHat
TsnowavgHat
DavgHat
rhoavgStdHat %newly-added, bulk density
DavgStdHat   %newly-added, bulk pex
TsnowavgStdHat   %newly-added, bulk T

end %properties



methods    
    %1) constructor
    %2) reads
    %3) analysis
    %4) utilities
    %5) plotting
    %6) writing
    
function Run=Main(Run)
   
   %Set basic running options
   %Run.folder='/Users/office/Downloads/MCMC/MCMCRunData_V3_2/Active_V1';
   %Run.NlyrMax=3;
   %Run.Lyrplan=1;  %%1=use modelSelection to select number of layers; 2=contrained number of layers to Run.Nlyr_choose
   
  
   %% Read inputs
   Run=Run.ReadFilenames;
   Run=Run.ReadRunParams;
   Run=Run.ReadObs;
   Run=Run.ReadHyperPar;
   
   % Read outputs
   Run=Run.ReadAcceptance; 
   Run=Run.ReadTbPost;
   Run=Run.ReadThetaOut;
   
   % Process outputs
   Run=Run.CalcThetaOut;
   Run=Run.ModelSelection;
   
   % Plot outputs
   % input1: 'dz','rho','D','Tsnow','Tsoil','mvs','sig','sd','swe','P_M','P_Q',or,'P_SR'
   % input2: which pit?
   %Run=Run.PlotParamResult('dz',1,symbol);
   %Run=Run.PlotProfileCompare(1,symbol);   
   %Run=Run.PlotObsChains(1,symbol);
end




%% read inputs & outputs
function Run=ReadFilenames(Run)
    
    filename=[Run.folder,'FILENAME.txt'];
    
    if(isempty(Run.sp)==1)
        error('There is no true snowpit information. Check!')
    end
    
end %function ReadFilenames


function Run=ReadRunParams(Run)
    
   filename=[Run.folder,'RunParams.txt'];
   fid=fopen(filename,'r');
   
   fgetl(fid); Run.Npits=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Niter=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Nburn=fscanf(fid,'%f \n',1);   %revised to not read Nburn
   fgetl(fid); Run.NlyrMax=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Nf=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.Nangle=fscanf(fid,'%f \n',1);
   
   fgetl(fid); Run.passive_np=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.active_np=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.other_np=fscanf(fid,'%f \n',1);
   fgetl(fid); fgetl(fid);
   fgetl(fid); fgetl(fid);

   fgetl(fid); Run.ScatOpt=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.ModelOpt=fscanf(fid,'%f \n',1);
   
   fgetl(fid); Run.stdTb=fscanf(fid,'%f \n',1);
   fgetl(fid); Run.stdSigma=fscanf(fid,'%f \n',1);
   fgetl(fid); 
   if(Run.other_np>0)
       Run.StdOther=fscanf(fid,'%f \n',1);
   end
   
   fgetl(fid);
   for iff=1:Run.Nf
       Run.MCMC_freq(iff)=fscanf(fid,'%f \n',1);
   end
   
   fgetl(fid);
   for itt=1:Run.Nangle
       Run.MCMC_theta(iff)=fscanf(fid,'%f \n',1);
   end
   
   fclose(fid);
end %function ReadRunParams



function Run=ReadObs(Run)
    %dimension-1: polarization
    %dimension-2: angle
    %dimension-3: frequency
    %dimension-4: pit

    Run.Np=Run.passive_np + Run.active_np + Run.other_np;
    
    Run.ObsMCMCIn=nan(Run.Np,Run.Nf,Run.Nangle,Run.Npits);
    Run.Obs=nan(Run.Np*Run.Nf*Run.Nangle,Run.Npits);
    Run.StdObs=nan(Run.Np*Run.Nf*Run.Nangle,Run.Npits);
        
    fid=fopen([Run.folder,'tb_obs.txt'],'r');
    
    for ipits=1:Run.Npits
        fgetl(fid);
        icount=0;
        
        for iff=1:Run.Nf
            for itt=1:Run.Nangle
                temp=fscanf(fid,'%f\n',[Run.Np,1])';
                Run.ObsMCMCIn(:,iff,itt,ipits)=temp;
                for ipol=1:length(temp)
                    temp2=temp(ipol);
                    if(isnan(temp2)==0)
                        icount=icount+1;
                        Run.Obs(icount,ipits)=temp2;
                        if(ipol<=Run.passive_np)
                            Run.StdObs(icount,ipits)=Run.stdTb;
                        end
                        
                        if(ipol>Run.passive_np & ipol<=Run.passive_np + Run.active_np)
                            Run.StdObs(icount,ipits)=Run.stdSigma;
                        end
                        
                        if(ipol>Run.passive_np+ Run.active_np)
                            Run.StdObs(icount,ipits)=Run.stdOther;
                        end
                    end
                end
            end
        end
    end
    fclose(fid);
    
    Run.Nobs=icount;
    Run.Obs(icount+1:end,:)=[];   %assuming all snowpits have the same number of measurements
    Run.StdObs(icount+1:end,:)=[];
end %function ReadObs



function Run=ReadHyperPar(Run)
    
    fid=fopen([Run.folder, 'hyperpar.txt'],'r');
    
    for ilyr=1:Run.NlyrMax
        
        %snow parameter
        fgetl(fid);
        
        Run=ReadHyperScan(Run,fid,'Dz',ilyr);
        Run=ReadHyperScan(Run,fid,'D',ilyr);
        Run=ReadHyperScan(Run,fid,'Rho',ilyr);
        Run=ReadHyperScan(Run,fid,'Tsnow',ilyr);
        %this and the following has a single value only.
        Run=ReadHyperScan(Run,fid,'Tsoil',1);
        Run=ReadHyperScan(Run,fid,'MvS',1);
        Run=ReadHyperScan(Run,fid,'GndSig',1);
        Run=ReadHyperScan(Run,fid,'P_M',1);
        Run=ReadHyperScan(Run,fid,'P_Q',1);
        Run=ReadHyperScan(Run,fid,'P_SR',1);
    end
    fclose(fid);  
end


function Run=ReadHyperScan(Run,fid,property,ilyr)
    fgetl(fid);
%     mu=fscanf(fid,'%f \n',ilyr); cov=fscanf(fid,'%f \n',ilyr);
%     [M,V]=lognstat(mu,sqrt(cov));

   %revised, the prior is normal-distrubution instead of lognormal
   %distribution
    M=fscanf(fid,'%f \n',ilyr); V=fscanf(fid,'%f \n',ilyr);
    mu = log((M.^2)./sqrt(V+M.^2));
    cov = log(V./(M.^2)+1); 
    
    eval(['Run.',property,'Mean{ilyr}=M;']) 
    eval(['Run.',property,'Std{ilyr}=sqrt(V);'])
    eval(['Run.',property,'Mu{ilyr}=mu;'])
    eval(['Run.',property,'Cov{ilyr}=cov;'])
end


function Run=ReadAcceptance(Run)
    
   fid=fopen([Run.folder,'acceptance.txt'],'r');

   for ipits=1:Run.Npits
       for ilyr=1:Run.NlyrMax
           fgetl(fid);           
           Run.acceptance.dz(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.rho(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.D(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.Tsnow(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.Tsoil(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.mvs(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.sig(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.P_M(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.P_Q(ipits,ilyr)=fscanf(fid,'%f \n',1);
           Run.acceptance.P_SR(ipits,ilyr)=fscanf(fid,'%f \n',1);
       end
   end
   fclose(fid);
end %function ReadAcceptance



function Run=ReadTbPost(Run)    

    %Nobs, number of real observations; is not following the same
    Run.ObsPost=nan(Run.Niter,Run.Nobs,Run.NlyrMax,Run.Npits);
    
    %read in the Tb output
    fid=fopen([Run.folder,'post_tb.out'],'rb');

    for ipits=1:Run.Npits    
        for ilyr=1:Run.NlyrMax
            for iobs=1:Run.Nobs
                fseek(fid,4,'cof');
                if(Run.Run_location==1)  %need to skip one more row~
                   fseek(fid,4,'cof');
                end
                %A=fread(fid,Run.Niter,'float'); %changed from 4 to 8
                A=fread(fid,Run.Niter,'real*8'); 
                
                Run.ObsPost(1:Run.Niter,iobs,ilyr,ipits)=A;
                
                fseek(fid,4,'cof');
                if(Run.Run_location==1)
                   fseek(fid,4,'cof');
                end
            end
        end      
    end
    fclose(fid);    
end %function ReadTbPost



function Run=ReadThetaOut(Run)
    
    Ntheta=nan(Run.NlyrMax,1);
    for ilyr=1:Run.NlyrMax,
        Ntheta(ilyr)=ilyr*4 + 3 + 3;  %now it is right, to be revised
        Run.theta_post{ilyr}=nan(Run.Niter, Ntheta(ilyr), Run.Npits);
    end

    fid=fopen([Run.folder,'post_theta.out'],'rb');

    for ipits=1:Run.Npits
        for ilyr=1:Run.NlyrMax
            for itheta=1:Ntheta(ilyr)
                
                fseek(fid,4,'cof');
                if(Run.Run_location==1)
                    fseek(fid,4,'cof');
                end 

                %A=fread(fid,Run.Niter,'float');
                A=fread(fid,Run.Niter,'real*8');
                Run.theta_post{ilyr}(1:Run.Niter,itheta,ipits)=A;
                
                fseek(fid,4,'cof');
                if(Run.Run_location==1)
                    fseek(fid,4,'cof');
                end
            end
        end      
    end
    
    fclose(fid);
end%function ReadThetaOut



%% try to deal with the output...
function [iDz,iRho,iD,iT,iTsoil,iMvS,iGndSig,iP_M,iP_Q,iP_SR] = GetPostIndices(Run,n)    
    iDz=1:n;
    iRho=n+1:n*2;
    iD=n*2+1:n*3;
    iT=n*3+1:n*4;
    iTsoil=n*4+1;
    iMvS=n*4+2;
    iGndSig=n*4+3;
    iP_M=n*4+4;
    iP_Q=n*4+5;
    iP_SR=n*4+6;
end %subfunction GetPostIndices


function Run=CalcThetaOut(Run)
    
    for ipits=1:Run.Npits
        for ilyr=1:Run.NlyrMax
            [iDz,iRho,iD,iT,iTsoil,iMvS,iGndSig,iP_M,iP_Q,iP_SR] = Run.GetPostIndices(ilyr);
            
            Run.chains_dz{ilyr,ipits}=Run.theta_post{ilyr}(:,iDz,ipits);
            Run.chains_rho{ilyr,ipits}=Run.theta_post{ilyr}(:,iRho,ipits);
            Run.chains_D{ilyr,ipits}=Run.theta_post{ilyr}(:,iD,ipits);
            Run.chains_Tsnow{ilyr,ipits}=Run.theta_post{ilyr}(:,iT,ipits);     %converted to degC
            Run.chains_Tsoil{ilyr,ipits}=Run.theta_post{ilyr}(:,iTsoil,ipits); %converted to degC
            Run.chains_mvs{ilyr,ipits}=Run.theta_post{ilyr}(:,iMvS,ipits);
            Run.chains_sig{ilyr,ipits}=Run.theta_post{ilyr}(:,iGndSig,ipits);
            Run.chains_P_M{ilyr,ipits}=Run.theta_post{ilyr}(:,iP_M,ipits);
            Run.chains_P_Q{ilyr,ipits}=Run.theta_post{ilyr}(:,iP_Q,ipits);
            Run.chains_P_SR{ilyr,ipits}=Run.theta_post{ilyr}(:,iP_SR,ipits);
            
            if(ilyr>1)
                dims=size(Run.chains_dz{ilyr,ipits});
                if(dims(1)==Run.Niter)
                    ndims=2;
                else if(dims(2)==Run.Niter)
                        ndims=1;
                    else
                        error('Check function CalcThetaOut');
                    end
                end
                Run.chains_sd{ilyr,ipits}=sum(Run.chains_dz{ilyr,ipits},ndims);
                Run.chains_swe{ilyr,ipits}=sum(Run.chains_dz{ilyr,ipits}.*Run.chains_rho{ilyr,ipits},ndims);
            else
                Run.chains_sd{ilyr,ipits}=Run.chains_dz{ilyr,ipits};
                Run.chains_swe{ilyr,ipits}=Run.chains_dz{ilyr,ipits}.*Run.chains_rho{ilyr,ipits};
            end
            
            Run=calc_md(Run,'dz',ilyr,ipits);
            Run=calc_md(Run,'D',ilyr,ipits);
            Run=calc_md(Run,'rho',ilyr,ipits);
            Run=calc_md(Run,'Tsnow',ilyr,ipits);
            Run=calc_md(Run,'Tsoil',ilyr,ipits);
            Run=calc_md(Run,'mvs',ilyr,ipits);
            Run=calc_md(Run,'sig',ilyr,ipits);
            Run=calc_md(Run,'sd',ilyr,ipits);
            Run=calc_md(Run,'swe',ilyr,ipits);
            Run=calc_md(Run,'P_M',ilyr,ipits);
            Run=calc_md(Run,'P_Q',ilyr,ipits);
            Run=calc_md(Run,'P_SR',ilyr,ipits);
        end
    end   
end


function Run=calc_md(Run,property,ilyr,ipits)
    
    %Jinmei, 2017/7/31 revised to mean, check if it is been revised here!
    eval(['data=Run.chains_',property,'{ilyr,ipits};']);
    dims=size(data);
    if(dims(1)==Run.Niter)
        ndims=1;
    else
        ndims=2;
    end
    
    iBurn=Run.Nburn+1:Run.Niter;
    logmean=0; logstd=0;
    if(strcmp(property,'dz')==1 | strcmp(property,'D')==1 | ...
       strcmp(property,'rho')==1 | strcmp(property,'Tsnow')==1)
        ilyrnn=ilyr;
    else
        ilyrnn=1;
    end
    
    for i=1:ilyrnn
        if(ndims==1)
            yy=data(iBurn,i);
        else
            yy=data(i,iBurn);
        end
        

        if(0) %revised for mean, Revised on Apr17, 2017 !!!!!!!!!!!!!!!!
            idx=find(isnan(yy)==1); yy(idx)=[];
            logyy=log(yy);

            sk_yy=skewness(yy); sk_logyy=skewness(logyy);

            if(abs(sk_logyy) < abs(sk_yy) & sk_yy>0)  %sk_yy>0 means, direct mean leans to left
                logmean(i)=exp(mean(logyy));
            else
                logmean(i)=mean(yy);
            end

            dum=(yy-logmean(i)).^2;
            logstd(i)= sqrt(sum(dum)/length(dum-1));
        else
            idx=find(isnan(yy)==1); yy(idx)=[];
            logmean(i)=mean(yy);      %to be revised
%             logmean(i)=median(yy);
%             logmean(i)=mode(yy);
            dum=(yy-logmean(i)).^2;
            logstd(i)= sqrt(sum(dum)/length(dum-1));
        end
    end
    
    eval(['Run.md_',property,'{ilyr,ipits}=logmean;']);
    eval(['Run.std_',property,'{ilyr,ipits}=logstd;']);
end


%% layer result selection
function Run=ModelSelection(Run)   %choose 1,2,3 or N layer plans for each pit?
    
    %
    Run.UsePrior=1;
    
    %Note: this hyper-parameter should be read in from file
    switch Run.SturmClass
        case 'alpine'
            lambda=5;
            disp('Alpine snow: lambda=5');
        otherwise
            lambda=2;      
    end

    iBurn=Run.Nburn+1:Run.Niter;
    J=nan(Run.Npits,Run.NlyrMax);
    
    for ipits=1:Run.Npits  
        Cy=diag(Run.StdObs(:,ipits).^2);
        Y=Run.Obs(:,ipits);
        
        for ilyr=1:Run.NlyrMax
            %Likelihood function. TbPost shape is: [Niter,Nobs,Nlyr,Npits]
            dum=Run.ObsPost(iBurn,:,ilyr,ipits);
            ndims=size(dum);
            if(ndims(1)~=Run.Nobs)
                Yhat=mean(dum)';
            else
                Yhat=mean(dum')';
            end
            fll(ilyr)=sum(log(mvnpdf(Y,Yhat,Cy)));
            
            %Number of layers
            p1ll(ilyr)=log(poisspdf(ilyr,lambda));
            
            %Depth
            dzhat=Run.md_dz{ilyr,ipits}';
            p2ll(ilyr)=log(mvnpdf(log(dzhat),Run.DzMu{ilyr},diag(Run.DzCov{ilyr})));
            
            %Grain size
            Dhat=Run.md_D{ilyr,ipits}';
            p3ll(ilyr)=log(mvnpdf(log(Dhat),Run.DMu{ilyr},diag(Run.DCov{ilyr})));
            
            %Density
            Rhohat=Run.md_rho{ilyr,ipits}';
            p4ll(ilyr)=log(mvnpdf(log(Rhohat),Run.RhoMu{ilyr},diag(Run.RhoCov{ilyr})));   
            
            %Temperature
            That=Run.md_Tsnow{ilyr,ipits}';
            p5ll(ilyr)=log(mvnpdf(log(That),Run.TsnowMu{ilyr},diag(Run.TsnowCov{ilyr}))); 
            
            %does not necessary to add model parameters..
            if Run.UsePrior==1
                J(ipits,ilyr)=fll(ilyr)+p1ll(ilyr)+p2ll(ilyr)+p3ll(ilyr)+p4ll(ilyr)+p5ll(ilyr);
                %disp('from top to bottom: obs, nlyr, thickness, grain size, density, temperature')
                %disp(num2str([fll;p1ll;p2ll;p3ll;p4ll;p5ll]))
            else
                J(ipits,ilyr)=fll(ilyr);
            end
        end 

        switch Run.Lyrplan
            case 1
                [~,Run.nHat(ipits)]=max(J(ipits,:));
            case 2
                Run.nHat(ipits)=Run.Nlyr_choose;
            case 99
                [~,Run.nHat(ipits)]=max(J(ipits,:));
                if(Run.nHat(ipits)==1)
                   Run.nHat(ipits)=2;
                   disp(['change nhat to 2 from 1 for pit:',num2str(ipits)]);
                end
            otherwise
                error('please tell me how to determine number of layers')
        end

        %calculate the retrived SWE and SD
        n=Run.nHat(ipits);
        Run.sdHat(ipits)=Run.md_sd{n,ipits};
        Run.sweHat(ipits)=Run.md_swe{n,ipits};
        Run.sdStdHat(ipits)=Run.std_sd{n,ipits};
        Run.sweStdHat(ipits)=Run.std_swe{n,ipits};

        
        %calculate mean and std of bulk values of D, rho and T
        %using depth-weighted averages
        massm = Run.chains_dz{n,ipits};
        weightm = massm./repmat(sum(massm')',1,n);
        rhom = Run.chains_rho{n,ipits};
        rhoavgm = sum((weightm').*(rhom'));
        Run.rhoavgHat(ipits) = mean(rhoavgm);
        Run.rhoavgStdHat(ipits)= std(rhoavgm);
        
        Dm = Run.chains_D{n,ipits};
        Davgm = sum((weightm').*(Dm'));
        Run.DavgHat(ipits) = mean(Davgm);
        Run.DavgStdHat(ipits) = std(Davgm);
        
        Tm = Run.chains_Tsnow{n,ipits};
        Tavgm = sum((weightm').*(Tm'));
        Run.TsnowavgHat(ipits) = mean(Tavgm);
        Run.TsnowavgStdHat(ipits) = std(Tavgm);
        
        
        %give priors of SWE and SD
        if(n==1)
            Run.sdPr(ipits)=sum(Run.DzMean{n});
            Run.swePr(ipits)=sum(Run.DzMean{n}.*Run.RhoMean{n});
        else
            Run.sdPr(ipits)=(sum(Run.DzMean{n}(2:end))+1) * Run.DzMean{n}(1);
            dz_save=Run.DzMean{n};
            dz_save(2:end)=dz_save(2:end).*dz_save(1);
            Run.swePr(ipits)=sum(dz_save.*Run.RhoMean{n});
        end
        
        %reprocess temperaturs for all parameters
        %chains_Tsnow, chains_Tsoil
        %md_Tsnow, md_Tsoil
        %TsnowavgHat
        for ilyr=1:Run.NlyrMax
            Run.chains_Tsnow{ilyr,ipits}=274-Run.chains_Tsnow{ilyr,ipits}-273.15;
            Run.chains_Tsoil{ilyr,ipits}=279-Run.chains_Tsoil{ilyr,ipits}-273.15;
            Run.md_Tsnow{ilyr,ipits}=274-Run.md_Tsnow{ilyr,ipits}-273.15;
            Run.md_Tsoil{ilyr,ipits}=279-Run.md_Tsoil{ilyr,ipits}-273.15;
        end
        Run.TsnowavgHat=274-Run.TsnowavgHat-273.15;
    end
end         



%% plot the result...
% function PlotParamResult(Run,property,ipits,symbol)
%     
%     MaxLyr=2;
%     X=1:Run.Niter;
%     X0=Run.Nburn;
%     iBurn=Run.Nburn+1:Run.Niter;
%     
%     %plot chains and histgram of chains
%     figure;
%     for ilyr=1:MaxLyr
%         
%         %build legend strings  
%         if(strcmp(property,'dz')==1 | strcmp(property,'D')==1 | ...
%            strcmp(property,'rho')==1 | strcmp(property,'Tsnow')==1)
%                 for i=1:ilyr
%                     lyrstr{i}=['Lyr',num2str(i)];
%                 end
%         else
%             lyrstr={'Lyr1'};
%         end
%         
%         %plot chains
%         subplot(2,MaxLyr,ilyr*2-1)
%         eval(['plot(X,Run.chains_',property,'{ilyr,ipits});'])
%         hold on;
%         a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
%         
%         title([symbol,': Nlyr=' num2str(ilyr)])
%         xlabel('Iteration #');
%         %ylabel(property);
%         switch property
%             case 'dz'
%                 ylabel('dz (m)');
%             case 'D'
%                 ylabel('D (mm)');
%             case 'Tsnow'
%                 ylabel('Tsnow ?degC)');
%             case 'rho'
%                 ylabel('density (kg/m^3)');
%             case 'Tsoil'
%                 ylabel('Tsoil ?degC)');
%             case 'sig'
%                 ylabel('soil roughness (m)');
%             case 'mvs'
%                 ylabel('soil moisture (m^3/m^3)')
%             case 'sd'
%                 ylabel('sd (m)');
%             case 'swe'
%                 ylabel('SWE (mm)');
%             otherwise
%                 ylabel(property);
%         end
%                 
%             
%         legend(lyrstr);
%         
%         
%         %plot histgram
%         subplot(2,MaxLyr,ilyr*2);
%         eval(['hist(Run.chains_',property,'{ilyr,ipits}(iBurn,:),50);'])
%         
%         
%         title([symbol,': Nlyr=' num2str(ilyr)])
%         xlabel(property);
%         ylabel('N');
%         legend(lyrstr);
%     end
% end


function PlotParamResult(Run,property,ipits,symbol)
    
    MaxLyr=2;
    X=1:Run.Niter;
    X0=Run.Nburn;
    iBurn=Run.Nburn+1:Run.Niter;
    
    %plot chains and histgram of chains
    figure;
    for ilyr=1:MaxLyr
        
        %build legend strings  
        if(strcmp(property,'dz')==1 | strcmp(property,'D')==1 | ...
           strcmp(property,'rho')==1 | strcmp(property,'Tsnow')==1)
                for i=1:ilyr
                    lyrstr{i}=['Lyr',num2str(i)];
                    lyrstr2{(i-1)*2+1}=['Lyr',num2str(i)];
                    lyrstr2{i*2}=['Lyr',num2str(i),'-prior'];
                end
        else
            lyrstr={'Lyr1'};
            lyrstr2={'Lyr1','Lyr1-prior'};
        end
        
        %plot chains
        subplot(2,MaxLyr,ilyr*2-1)
        eval(['plot(X,Run.chains_',property,'{ilyr,ipits});'])
        hold on;
        a=axis; plot([X0,X0],[a(3),a(4)],'k--','Linewidth',2)
        
        title([symbol,': Nlyr=' num2str(ilyr)])
        xlabel('Iteration #');
        %ylabel(property);
        switch property
            case 'dz'
                ylabel('dz (m)');
            case 'D'
                ylabel('D (mm)');
            case 'Tsnow'
                ylabel('Tsnow ^oC)');
            case 'rho'
                ylabel('density (kg/m^3)');
            case 'Tsoil'
                ylabel('Tsoil ^oC)');
            case 'sig'
                ylabel('soil roughness (m)');
            case 'mvs'
                ylabel('soil moisture (m^3/m^3)')
            case 'sd'
                ylabel('sd (m)');
            case 'swe'
                ylabel('SWE (mm)');
            otherwise
                ylabel(property);
        end
                
            
        legend(lyrstr);
        
        
        %plot histgram
        subplot(2,MaxLyr,ilyr*2);
        eval(['chains00=Run.chains_',property,'{ilyr,ipits}(iBurn,:);'])
        
        
        %add prior distribution for histogram
        for ilyr2=1:min(ilyr,size(chains00,2))
            %plot histogram
            Is_show_log=0;
            if(Is_show_log==1) %plot log of values
                chains002=log(chains00(:,min(ilyr2,size(chains00,2))));
                [counts0,center0]=hist(chains002,50);

                area00=(counts0(1:end-1)+counts0(2:end))/2.*(center0(2:end)-center0(1:end-1));
                area00_sum=sum(area00);
                if(ily2==1)
                    plot(center0,counts0/area00_sum,'b.-');
                else
                    plot(center0,counts0/area00_sum,'r.-');
                end
                hold on;
                xlabel(['log(',property,')']);
            else %plot values directly
                chains002=chains00(:,min(ilyr2,size(chains00,2)));
                [counts0,center0]=hist(chains002,100);
                area00=(counts0(1:end-1)+counts0(2:end))/2.*(center0(2:end)-center0(1:end-1));
%                 area00_sum=sum(area00);
                area00_sum=length(chains002);
                if(ilyr2==1)
                    plot(center0,counts0/area00_sum,'b.-');
                else
                    plot(center0,counts0/area00_sum,'r.-');
                end
                hold on;
                xlabel(property);
            end


            %legends
            title([symbol,': Nlyr=' num2str(ilyr)])
            ylabel('Probability');
            

            %plot prior;
            Is_show_log=0;
            if(Is_show_log==1)
                %note, this code is not correct for normal distributed values
                chains002_min=min(chains002); %log distribution
                chains002_max=max(chains002);
                chains002_step=(chains002_max-chains002_min)/100;
                chains002_range=(chains002_max-chains002_min);
                YYY=chains002_min:chains002_step:chains002_max;

                switch property
                    case 'dz'
                        Mu0=Run.DzMu{ipits,ilyr}(ilyr2); Cov0=Run.DzCov{ipits,ilyr}(ilyr2);
                        %PrHist=random('Lognormal',Mu0,sqrt(Cov0),2000,1);
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'D'
                        Mu0=Run.DMu{ipits,ilyr}(ilyr2); Cov0=Run.DCov{ipits,ilyr}(ilyr2);
                        %PrHist=random('Lognormal',Mu0,sqrt(Cov0),2000,1);
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'rho'
                        Mu0=Run.RhoMu{ipits,ilyr}(ilyr2); Cov0=Run.RhoCov{ipits,ilyr}(ilyr2);
                        %PrHist=random('Lognormal',Mu0,sqrt(Cov0),2000,1);
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'Tsnow'
                        Mu0=Run.TsnowMu{ipits,ilyr}(ilyr2); Cov0=Run.TsnowCov{ipits,ilyr}(ilyr2);
                        %PrHist=random('Lognormal',Mu0,sqrt(Cov0),2000,1);
                        %PrHist=274-PrHist-273.15;
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'Tsoil'
                        Mu0=Run.TsoilMu{ipits,1}; Cov0=Run.TsoilCov{ipits,1};
                        %PrHist=random('Lognormal',Mu0,sqrt(Cov0),2000,1);
                        %PrHist=279-PrHist-273.15;
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'sig'
                        Mu0=Run.GndSigMu{ipits,1}; Cov0=Run.GndSigCov{ipits,1};
                        %PrHist=random('Normal',Mu0,sqrt(Cov0),2000,1);
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'mvs'
                        Mu0=Run.MvSMu{ipits,1}; Cov0=Run.MvSCov{ipits,1};
                        %PrHist=random('Normal',Mu0,sqrt(Cov0),2000,1);
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'P_M'
                        Mu0=Run.P_MMu{ipits,1}; Cov0=Run.P_MCov{ipits,1};
                        %PrHist=random('Normal',Mu0,sqrt(Cov0),2000,1);
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'P_Q'
                        Mu0=Run.P_QMu{ipits,1}; Cov0=Run.P_QCov{ipits,1};
                        %PrHist=random('Normal',Mu0,sqrt(Cov0),2000,1);
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                    case 'P_SR'
                        Mu0=Run.P_SRMu{ipits,1}; Cov0=Run.P_SRCov{ipits,1};
                        %PrHist=random('Normal',Mu0,sqrt(Cov0),2000,1);
                        PrHist=normpdf(YYY,Mu0,sqrt(Cov0));
                end

                
                area01=(PrHist(1:end-1)+PrHist(2:end))/2.*(YYY(2:end)-YYY(1:end-1));
                area01_sum=sum(area01);
                plot(YYY,PrHist/area01_sum,'ro-');hold on;
            else
                %note, this code is not correct for normal distributed values
%                 chains002_min=min(chains002);
%                 chains002_max=max(chains002);
%                 chains002_step=(chains002_max-chains002_min)/100;
%                 chains002_range=(chains002_max-chains002_min);
%                 YYY=chains002_min:chains002_step:chains002_max;

                switch property
                    case 'dz'
                        if(ilyr2==1)
                            Mu0=Run.DzMean{ipits,ilyr}(ilyr2); Std0=Run.DzStd{ipits,ilyr}(ilyr2);
                        else 
                            Mu0=Run.DzMean{ipits,ilyr}(ilyr2) .* Run.DzMean{ipits,ilyr}(1); 
                            %Std0=Run.DzStd{ipits,ilyr}(ilyr2) .* Run.DzStd{ipits,ilyr}(1); 
                            Std0=Run.DzMean{ipits,ilyr}(ilyr2) .* Run.DzStd{ipits,ilyr}(1) + ...
                                Run.DzStd{ipits,ilyr}(ilyr2) .* Run.DzMean{ipits,ilyr}(1);
                        end
                    case 'D'
                        Mu0=Run.DMean{ipits,ilyr}(ilyr2); Std0=Run.DStd{ipits,ilyr}(ilyr2);
                    case 'rho'
                        Mu0=Run.RhoMean{ipits,ilyr}(ilyr2); Std0=Run.RhoStd{ipits,ilyr}(ilyr2);
                    case 'Tsnow'
                        Mu0=Run.TsnowMean{ipits,ilyr}(ilyr2); Std0=Run.TsnowStd{ipits,ilyr}(ilyr2);
                        Mu0=(274-273.15)-Mu0;
                    case 'Tsoil'
                        Mu0=Run.TsoilMean{ipits,1}; Std0=Run.TsoilStd{ipits,1};
                        Mu0=(279-273.15)-Mu0;
                    case 'sig'
                        Mu0=Run.GndSigMean{ipits,1}; Std0=Run.GndSigStd{ipits,1};
                    case 'mvs'
                        Mu0=Run.MvSMean{ipits,1}; Std0=Run.MvSStd{ipits,1};
                    case 'P_M'
                        Mu0=Run.P_MMean{ipits,1}; Std0=Run.P_MStd{ipits,1};
                    case 'P_Q'
                        Mu0=Run.P_QMean{ipits,1}; Std0=Run.P_QStd{ipits,1};
                    case 'P_SR'
                        Mu0=Run.P_SRMean{ipits,1}; Std0=Run.P_SRStd{ipits,1};
                    case 'swe'
                        mu0000 = Run.DzMean{ipits,ilyr}(1).*Run.RhoMean{ipits,ilyr}(1);
                        cov000 = (Run.DzStd{ipits,ilyr}(1).* Run.RhoStd{ipits,ilyr}(1)).^2 + ...
                            (Run.DzStd{ipits,ilyr}(1).* Run.RhoMean{ipits,ilyr}(1)).^2 + ...
                            (Run.DzMean{ipits,ilyr}(1).* Run.RhoStd{ipits,ilyr}(1)).^2;
                        if(ilyr==1)
                            Mu0= mu0000;
                            Std0=sqrt(cov000);
                        else
                            Mu0 = mu0000.*ilyr;
                            Std0=sqrt(cov000.*ilyr);
                        end
                        
                    case 'sd'
                        mu0000 = Run.DzMean{ipits,ilyr}(1);
                        std0000 = Run.DzStd{ipits,ilyr}(1);
                        if(ilyr==1)
                            Mu0=mu0000;
                            Std0=std0000;
                        else
                            Mu0=mu0000.*ilyr;
                            Std0= sqrt(std0000^2.*ilyr);
                        end
                        
                end
                
%                 PrHist=lognpdf(YYY,Mu0,sqrt(Cov0));
                YYY=Mu0-Std0*4: Std0/20: Mu0+Std0*4;
                PrHist=normpdf(YYY,Mu0,Std0);
                area01=(PrHist(1:end-1)+PrHist(2:end))/2.*(YYY(2:end)-YYY(1:end-1));
%                 area01_sum=sum(area01);
                area01_sum=sum(PrHist);
                plot(YYY,PrHist/area01_sum,'ro-');hold on;
            end
            
            legend(lyrstr2);
        end
    end
end




function Run=PlotChainsLikelihood(Run,property,ilyr,ipits)
    %figure(1), always plot Obs likelihood
    if(strcmp(property,'dz')==1)
        Cy=Run.StdObs(:,ipits);
        Y=Run.Obs(:,ipits);

        mcmc_obsr=nan(Run.Niter,length(Y));

        for i=1:Run.Niter
            Y0=Run.ObsPost(i,:,ilyr,ipits)';

            for j=1:length(Y)            
                mcmc_obsr(i,j)=normpdf(Y0(j),Y(j),Cy(j));
            end
        end

        figure;
        colorstr=jet(length(Y));
        for j=1:length(Y)
            subplot(length(Y),2,(j-1)*2+1);
            plot(Run.ObsPost(:,j,ilyr,ipits)','k.-','color',colorstr(j,:)); hold on;
            plot(Run.Niter,Y(j),'kx','color',colorstr(j,:),'MarkerSize',10)
            
            subplot(length(Y),2,(j-1)*2+2);
            plot([1:Run.Niter]',log(mcmc_obsr(:,j)),'k.-','color',colorstr(j,:)); hold on;
            legendstr{j}=num2str(Y(j));
            title(num2str(Y(j)))
            set(gca,'Ylim',[-5,2]);
        end
        
        figure;
        likeli=nan(Run.Niter,1);
        for i=1:Run.Niter
            temp=mcmc_obsr(i,:);
            temp2=sum(log(temp));
            likeli(i)=temp2;
        end
        plot([1:Run.Niter]',likeli,'b.'); hold on;
        ylabel('sum(log(p(\sigma)))');
        
%         Cy2=diag(Run.StdObs(:,ipits).^2);
%         for i=1:Run.Niter
%             temp=mcmc_obsr(i,:);
%             likeli2(i)=mvnpdf(temp,Y',Cy2);
%         end
%         plot([1:Run.Niter]',likeli,'r.'); hold on;
    end    
    
    
    %figure(2), plot prior probablity
    figure;
    eval(['chains=Run.chains_',property,'{ilyr,ipits};'])
    
    switch property
        case 'dz'
            property='Dz';
        case 'mvs'
            property='MvS';
        case 'rho'
            property='Rho';
        case 'sig'
            property='GndSig'
    end

    eval(['Cy=Run.',property,'Cov{ilyr};'])
    eval(['Y=Run.',property,'Mu{ilyr};'])
    
    mcmc_prob=nan(Run.Niter,ilyr);
    
    for i=1:Run.Niter
        for j=1:ilyr
            Y0=chains(i,j);
            mcmc_prob(i,j)=normpdf(log(Y0),Y(j),Cy(j));
        end
    end

    colorstr=jet(ilyr);
    clear legendstr
    for j=1:ilyr
        plot([1:Run.Niter]',mcmc_prob(:,j),'k.-','color',colorstr(j,:)); hold on;
        legendstr{j}=num2str(j);
    end
    legend(legendstr);

end


function PlotObsChains(Run,ipits,symbol)
    %Obs(Run.Nobs,Run.Npits)
    %ObsPost(Run.Niter,Run.Nobs,Run.NlyrMax,Run.Npits)
    X=1:Run.Niter;
    X0=Run.Nburn;
    
    %build legendstr
    pol_all=[Run.passive_pol,Run.active_pol,'sd'];
    fstr='rgbcmyk';
    degstr={'o','sq','x','+','.','*','d','^','<','>','p','h'};
    polstr={'-','--',':','-.'};
    
    icount=0;    
    for iff=1:Run.Nf
        for itt=1:Run.Nangle
            for ipol=1:Run.Np
                if(Run.ObsMCMCIn(ipol,iff,itt,ipits)<900)
                    icount=icount+1;
                    name=[num2str(Run.MCMC_freq(iff)),'GHz-',num2str(Run.MCMC_theta(itt)),'deg-',pol_all{ipol}];
                                        
                    legendstr{icount}=name;
                    markers{icount}=[fstr(iff),degstr{itt},polstr{ipol}];
                end
            end
        end
    end
    
    
    %plot observations
    figure; hold on;
    for i=1:Run.Nobs
        plot(Run.Niter,Run.Obs(i,ipits),markers{i}(1:2),'MarkerSize',10);
    end
    title(symbol);
    legend(legendstr);
    
    
    %plot chains
    for i=1:Run.Nobs
        plot(X,Run.ObsPost(:,i,Run.nHat(ipits),ipits),[markers{i}(1),markers{i}(3)]);
    end
end



function PlotProfileCompare(Run,ipits,symbol,plotNo)

    %the MCMC result
    n=Run.nHat(ipits);
    dz_post=Run.md_dz{n,ipits}*100;
    rho_post=Run.md_rho{n,ipits};
    D_post=Run.md_D{n,ipits};
    T_post=Run.md_Tsnow{n,ipits};
    
    
    %the measurements
    dz_true=Run.sp.dz'*100;
    rho_true=Run.sp.density';
    T_true=Run.sp.T'-273.15;
    if(Run.ModelOpt==1)
        D_true=Run.sp.pex';
    else
        D_true=Run.sp.dmax';
    end
    
    %plot and compare
    temp=str2num(symbol(4:end));
%     figure(temp);
    figure(plotNo);
    set(gcf,'color','w');
    subplot(1,3,1);
    [h1,h2]=plot_profile(Run,dz_post,rho_post,'ro-');
    [h3,h4]=plot_profile(Run,dz_true,rho_true,'ko-');
    title([symbol,': density']);
    xlabel('density (kg/m^3)'); ylabel('z (cm)');
    legend([h2,h4],{'MCMC','True'});
    
    
    subplot(1,3,2);
    [h1,h2]=plot_profile(Run,dz_post,D_post,'ro-');hold on;
    [h3,h4]=plot_profile(Run,dz_true,D_true,'ko-');
    title([symbol,': micro-structure']);
    if(Run.ModelOpt==1)
        xlabel('p_{ec} (mm)'); 
    else
        xlabel('D_{max} (mm)'); 
    end
    ylabel('z (cm)');
    legend([h2,h4],{'MCMC','True'});
    
    subplot(1,3,3);
    [h1,h2]=plot_profile(Run,dz_post,T_post,'ro-');hold on;
    [h3,h4]=plot_profile(Run,dz_true,T_true,'ko-');
    title([symbol,': temperature']);
    xlabel('T (^oC)'); ylabel('z (cm)');   
    legend([h2,h4],{'MCMC','True'});
end



function [H1,H2]=plot_profile(Run,dz,data,Mark)
    
    %plot MCMC results
    z(1)=dz(1)/2;
    for i=2:length(dz)
        z(i)=z(i-1)+dz(i-1)/2 + dz(i)/2;
    end
    H1=plot(data,z,Mark(1:2),'MarkerSize',8.0); hold on;
    
    %plot limits
    z1=nan(1,length(data));
    zh=nan(1,length(data));
    
    zl(1)=0;
    zh(1)=dz(1);
    for i=2:length(dz)
        zl(i)=zl(i-1)+dz(i-1);
        zh(i)=zh(i-1)+dz(i);
    end
    
    z2=reshape([zl;zh],length(dz)*2,1);
    data2=reshape([data;data],length(dz)*2,1);
    H2=plot(data2,z2,'k-','color',Mark(1),'Linewidth',1.0);
end


end %methods

end %classdef
