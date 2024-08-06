%snowpit structure
%revised to add frost_depth

classdef snowpit

properties
    %basic
    provider
    site
    year
    month
    date
    time
    
    %metereological basic
    Tair   %air temperature in K
    air_pressure  %air pressure in hPa
    rh   %relative humidity in %
    
    %properties
	nlayer
    dz      %in m
    density %in kg/m^3
    T       %in K
    dmax    %in mm
    pex     %in mm
    mv      %m^3/m^3
    
    
    %clustered properties, Y=cell(7)
    Y
    
    %soil properties
    mv_soil          %m^3/m^3
    soilT            %K
    dcL_soil          %soil permittivity at 100MHz
    roughness_fitted %m, fitted using the passive microwave data only~
    mv_soil_fitted
    frost_depth
    
    %simulation parameters
    pex_scaler1      %the factor to rescale Tb, aimed to help closely match Tb at high&low frequencies
    pex_scaler2      %the factor to scale the calibrated pex using Tb, to match radar measurements
    rdv2rv           %the ratio of rdv to rv at different frequencies, used to calibrate q
    
    %simulated Sky temperature using local moisture and air temperature (assuming clear sky condition)
    Tsky
    
    %passive measurements
    %TB is a matrix of (theta,freq)
    theta   %degree
    freq
    tbv     %K
    tbh     %K
    tbv_simu
    tbh_simu
    
    %active measurements
    %sigma is a matrix of (theta_ac,freq_ac)
    theta_ac
    freq_ac
    m
    q
    sigma_vv
    sigma_hh
    sigma_vh
    sigma_hv
    sigma_vv_simu
    sigma_hh_simu
    sigma_vh_simu
    sigma_hv_simu
    
    %statistics
    SD      %m
    SWE     %mm
    avg_density  %kg/m^3
    avg_T        %K
    avg_dmax     %mm
    avg_pex      %mm
    
    %note
    note
end

methods
    function snowpit=summary(snowpit)
        snowpit.SD=sum(snowpit.dz);
        snowpit.SWE=sum(snowpit.dz.*snowpit.density);
        weight=snowpit.dz.*snowpit.density;
        snowpit.avg_density=sum(weight.*snowpit.density)./sum(weight);
        snowpit.avg_T=sum(weight.*snowpit.T)./sum(weight);
        snowpit.avg_dmax=sum(weight.*snowpit.dmax)./sum(weight);
        snowpit.avg_pex=sum(weight.*snowpit.pex)./sum(weight);
        
        Zeros=repmat(0,snowpit.nlayer,1);
        snowpit.Y=cell(7,1);
        
        % inputs to MEMLS, from left to right:
        % layer-number, temp [K], volume fraction of liquid water, density [kg/m3],
        % thickness [cm], Salinity (0 - 0.1) [ppt], expon.corr.length [mm]
        snowpit.Y{1}=[[1:1:snowpit.nlayer]', snowpit.T, Zeros, snowpit.density, snowpit.dz*100, Zeros, snowpit.pex];
        
    end
end

end