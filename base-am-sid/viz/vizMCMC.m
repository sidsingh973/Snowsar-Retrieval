clear all

folder_in='./';

i=1; 
Param=[];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

filename = '/media/ss20/Local/Project2/Ave_SAR/172126_Ave_90m.csv';
Ave = table2cell(readtable(filename));

Ave = cell2mat(Ave);
idx=find(Ave(:,14)==71);
AveG = Ave(idx,:);
AveG(any(isnan(AveG),2),:) = [];
%AveG(1,1)
dt = filename(36:41);
memls = '/media/ss20/Local/Project2/Process_Memls/';
i=1;
for i=1:51
dirn = strcat(memls,dt,'/',strrep(num2str(AveG(i,2)),".","_"),strrep(num2str(AveG(i,1)),".","_"));
disp(i);
test_acceptance = strcat(dirn,'/test/acceptance.txt');
test_tbin = strcat(dirn,'/test/tb_obs.txt');
test_filename = strcat(dirn,'/test/FILENAME.txt');
test_hyperpar = strcat(dirn,'/test/hyperpar.txt');
test_RunParams = strcat(dirn,'/test/RunParams.txt');
test_posttb = strcat(dirn,'/test/post_tb.out');
test_posttheta = strcat(dirn,'/test/post_theta.out');
viz = strcat('/media/ss20/Local/Project2/base-am-main/viz/');

copyfile(test_acceptance,viz);
copyfile(test_filename,viz);
copyfile(test_tbin,viz);
copyfile(test_hyperpar,viz);
copyfile(test_RunParams,viz);
copyfile(test_posttb,viz);
copyfile(test_posttheta,viz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mcmc=MCMCRun4;
% filename=[folder_in,'sd',num2str(i),'/MCMC.mat'];
filename='./MCMC.mat';
load(filename);
%read MCMC results
% mcmc.folder=[folder_in,'sd',num2str(i),'/'];
mcmc.folder='./';
mcmc.Run_location=0;  %skip two 4 bytes, use 1; otherwise, use 0
%     mcmc.Lyrplan=1;
mcmc.Lyrplan=2;mcmc.Nlyr_choose=2;
    
mcmc.SturmClass='taiga';
mcmc=mcmc.Main;
rho =mcmc.md_rho{2};
dz =mcmc.md_dz{2};
Tsnow = mcmc.md_Tsnow{2};
a=[AveG(i,1:5) AveG(i,12:13) rho dz Tsnow]; 


Param = [Param ; a];
%symbol=['pit',num2str(i)]; %set title of the plots
save(strcat(dirn,'/viz/mcmc.mat'),'mcmc');
end

Param = array2table(Param);
Param.Properties.VariableNames(1:13)={'x','y','T2m','AngX','AngKu','VVx','VVku','rho1','rho2','dz1','dz2','Temp1','Temp2'};
writetable(Param, strcat('/media/ss20/Local/Project2/',dt,'_VV_Param_90m.csv'), 'WriteVariableNames', true);

%
%mcmc.PlotProfileCompare(1,symbol,33);
%mcmc.PlotObsChains(1,symbol);
%mcmc.PlotParamResult('dz',1,symbol);
%mcmc.PlotParamResult('rho',1,symbol);
%mcmc.PlotParamResult('Tsoil',1,symbol);
%mcmc.PlotParamResult('sig',1,symbol);
%mcmc.PlotParamResult('P_Q',1,symbol);
