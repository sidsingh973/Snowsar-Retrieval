clc;  clear all;
delete(gcp('nocreate'));
nm = ['20';'40';'60';'80'];
s=1;
parpool(32);

for K=1:1
filename = strcat('/srv/local/data/Project2/SnowSAR_HRRR/181138_SARGRD_30m_',nm(s,:),'.csv');
%Avetab=readtable(filename);
Ave = table2cell(readtable(filename));
Ave = cell2mat(Ave);

dt = filename(39:44);

fl = strcat('/srv/local/data/Project2/GRD/',dt,'/Grd_Param_30m.csv');
soil = table2cell(readtable(fl));
soil = cell2mat(soil);

%fl = strcat('/srv/local/data/Project2/1lyr/',dt,'/Q_',nm(s,:),'/Param_90m.csv');
%MQ =  table2cell(readtable(fl));
%MQ = cell2mat(MQ);

load('/srv/local/data/Project2/SnowSAR_HRRR/181138_9.6 GHz_Aux_GL.mat');

%AveG(1,1)
memls = '/srv/local/data/Project2/Process_Memls/';
fprintf('This message is sent at time %s\n', datestr(now,'HH:MM:SS.FFF'))

start = datestr(now,'HH:MM:SS.FFF');
num_layer=2;
parfor i=1:numel(Ave(:,1))
    try
    dirn = strcat(memls,'/SAR_GRD/',dt,'/2lyr/SWE/30m/','Run1/',strrep(num2str(Ave(i,2)),".","_"),strrep(num2str(Ave(i,1)),".","_"));
    mkdir(dirn);
    copyfile('/srv/local/data/Project2/base-am-main',dirn)
    copyfile('/srv/local/data/Project2/baseam_meta',dirn)

depth = sum(DepthGL(i,:));

swe = SWE(i);
nlayer = numel(nonzeros(DensityGL(i,:)));
rho1 = min(nonzeros(DensityGL(i,:)));
rho2 = max(nonzeros(DensityGL(i,:)));
sg1 = min(nonzeros(DsnowGL(i,:)));
sg2 = max(nonzeros(DsnowGL(i,:)));
    dense = nonzeros(DensityGL(i,:));
    delden=[];
    for j = 1:(numel(dense)-1)
    deldense = (dense(j+1)-dense(j))/dense(j);
    delden = [delden;deldense];
    
    end
[M,nl1] = min(delden);
%nl1 = floor(nlayer/2);
nl2 = nlayer - nl1;

depth1 = sum((DepthGL(i,1:nl1)));
depth2 = sum((DepthGL(i,(nl1+1):nlayer)));
den = mean((DensityGL(i,1:nlayer)));
den1 = sum((DepthGL(i,1:nl1)).*DensityGL(i,1:nl1))/depth1;
den2 = sum((DepthGL(i,(nl1+1):nlayer)).*DensityGL(i,(nl1+1):nlayer))/depth2;
D1 = sum((DepthGL(i,1:nl1)).*DsnowGL(i,1:nl1))/depth1;
D2 = sum((DepthGL(i,(nl1+1):nlayer)).*DsnowGL(i,(nl1+1):nlayer))/depth2;
Temp1 = sum((DepthGL(i,1:nl1)).*TsnowGL(i,1:nl1))/depth1;
Temp2 = sum((DepthGL(i,(nl1+1):nlayer)).*TsnowGL(i,(nl1+1):nlayer))/depth2;

Temp = mean(nonzeros(TsnowGL(i,1:nlayer)));

Ts1 = min(nonzeros(TsnowGL(i,:)));
Ts2 = max(nonzeros(TsnowGL(i,:)));
  
fid = fopen(strcat(dirn,'/src/tb_obs.txt'),'r');
j = 1;
tline = fgetl(fid);
A= [];
A{j} = tline;
while ischar(tline)
    j = j+1;
    tline = fgetl(fid);
    A{j} = tline;
end
fclose(fid);
A{2} = strcat("    ",num2str(Ave(i,7)));
A{3} = strcat("    ",num2str(Ave(i,8)));

fid = fopen(strcat(dirn,'/src/tb_obs.txt'),'w');
for j = 1:numel(A)
    if isstring(A{j+1})
        fprintf(fid,'%s\n', A{j});
    else
        fprintf(fid,'%s\n', A{j});
        break   
    end
end
fclose(fid);


fid = fopen(strcat(dirn,'/src/RunParams.txt'),'r');
j = 1;
tline = fgetl(fid);
A = [];
A{j} = tline;
while ischar(tline)
    j = j+1;
    tline = fgetl(fid);
    A{j} = tline;
end
fclose(fid);

% Change Runparams.txt
A{4} = strcat("    ",num2str(12000));

A{8} = strcat("    ",num2str(num_layer));

A{8} = strcat("    ",num2str(num_layer));
A{36} = strcat("    ",num2str(Ave(i,4)));

rho = swe*1000 / depth;
A{76} = strcat("    ",num2str(0.2*depth));
A{77} = strcat("    ",num2str(0.9*depth));
A{79} = strcat("    ",num2str(0.8*rho1));
A{80} = strcat("    ",num2str(1.2*rho2));

A{82} = strcat("    ",num2str(sg1));
A{83} = strcat("    ",num2str(sg2));

Ts1 = Ts1 - 273.15;
Ts2 = Ts2 - 273.15;

A{85} = strcat("    ",num2str(1.2*Ts1));
A{86} = strcat("    ",num2str(0.8*Ts2));

A{97} = strcat("    ",num2str(0.95*0.11));
A{98} = strcat("    ",num2str(1.05*0.11));

A{100} = strcat("    ",num2str(0.95*0.1));
A{101} = strcat("    ",num2str(1.05*0.1));

A{103} = strcat("    ",num2str(0.95*0.8));
A{104} = strcat("    ",num2str(1.05*0.8));

sthet = soil(i,5) - 0.01*soil(i,5);
srm = soil(i,6) - 0.01*soil(i,6);
A{88} = strcat("    ",num2str(sthet));
A{91} = strcat("    ",num2str(srm));
sthet = soil(i,5) + 0.01*soil(i,5);
srm = soil(i,6) + 0.01*soil(i,6);
A{89} = strcat("    ",num2str(sthet));
A{92} = strcat("    ",num2str(srm));


fid = fopen(strcat(dirn,'/src/RunParams.txt'),'w');
for j = 1:numel(A)
    if ~isnumeric(A{j+1})
        fprintf(fid,'%s\n', A{j});
    else
        fprintf(fid,'%s\n', A{j});
        break
    end
end
fclose(fid);


% Change hyperparams.txt

fid = fopen(strcat(dirn,'/src/hyperpar.txt'),'r');
j = 1;
tline = fgetl(fid);
A = [];
A{j} = tline;
while ischar(tline)
    j = j+1;
    tline = fgetl(fid);
    A{j} = tline;
end

fclose(fid);

A{34} = strcat("  ",num2str(depth1));
A{35} = strcat("  ",num2str(depth2));
A{36} = strcat("  ",num2str(0.1*depth1));
A{37} = strcat("  ",num2str(0.2*depth2));

A{39} = strcat("  ",num2str(D1));
A{40} = strcat("  ",num2str(D2));
A{41} = strcat("  ",num2str(0.2*D1));
A{42} = strcat("  ",num2str(0.2*D2));

A{44} =num2str(den1);
A{45} = num2str(den2);
A{46} =num2str(0.3*den1);
A{47} = num2str(0.3*den2);

Temp1a = 274 - Temp1;
Temp2a = 274 - Temp2;
A{49} = strcat("  ",num2str(Temp1a));
A{50} = strcat("  ",num2str(Temp2a));

A{54} = strcat("  ",num2str(5.0));
A{55} = strcat("  ",num2str(0.5));

A{57} = strcat("  ",num2str(soil(i,5)));
A{58} = strcat("  ",num2str(0.01*soil(i,5)));

A{60} = strcat("  ",num2str(soil(i,6)));
A{61} = strcat("  ",num2str(0.01*soil(i,6)));

A{63} = strcat("  ",num2str(0.11));
A{64} = strcat("  ",num2str(0.001*0.11));

A{66} = strcat("  ",num2str(0.1));
A{67} = strcat("  ",num2str(0.001*0.2));

A{69} = strcat("  ",num2str(0.8));
A{70} = strcat("  ",num2str(0.001*0.8));

fid = fopen(strcat(dirn,'/src/hyperpar.txt'),'w');
for j = 1:numel(A)
    if ~isnumeric(A{j+1})
        fprintf(fid,'%s\n', A{j});
    else
        fprintf(fid,'%s\n', A{j});
        break
    end
end
fclose(fid);
%Here we call base am 
src = strcat(dirn,'/src/');
delete(strcat(dirn,'/src/MetroMEMLS'));
system(strcat("make --directory ",src));
system(strcat(dirn,'/src/MetroMEMLS'));
src_acceptance = strcat(dirn,'/src/acceptance.txt');
src_tbin = strcat(dirn,'/src/tb_obs.txt');
src_filename = strcat(dirn,'/src/FILENAME.txt');
src_hyperpar = strcat(dirn,'/src/hyperpar.txt');
src_RunParams = strcat(dirn,'/src/RunParams.txt');
src_posttb = strcat(dirn,'/src/post_tb.out');
src_posttheta = strcat(dirn,'/src/post_theta.out');
viz = strcat(dirn,'/test/');
copyfile(src_acceptance,viz);
copyfile(src_filename,viz);
copyfile(src_tbin,viz);
copyfile(src_hyperpar,viz);
copyfile(src_RunParams,viz);
copyfile(src_posttb,viz);
copyfile(src_posttheta,viz);
    catch
        disp 'error'
    end
end
s=s+1;
finish = datestr(now,'HH:MM:SS.FFF');
time = (finish-start);
fprintf('This message is sent at time %s\n', datestr(now,'HH:MM:SS.FFF'))
end
