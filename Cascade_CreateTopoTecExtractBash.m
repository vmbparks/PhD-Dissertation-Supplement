%% make server download files

%%% This is constructed to run on Mac/Unix with bash systems, with macports
%%% and sshpass installed.

%% chmod 777 filename

%%
% currently must be run from inside /Volumes/Files/VictoriaFiles/Cascade/
% folder

% Change this per model
timesteps=[1540:20:1560];%11.9e6/5000=2380
ModelFolder='Cascade_N9v55';
LocalSaveFolder='McQ02N9/v55-DefR-v76Orog/topo_tec/';

%LocalSaveFolder='McQ02N9/v54-DefL-v73Orog/topo_tec/';
SaveFolder='/Volumes/Files/VictoriaFiles/Cascade';%uigetdir();

% File type: topo_tec
% which timesteps to download: [1580:20:2000]
% Folder on server (eg username + model_runs + Foldername +
% output/IceCascade/
Serverloc='/esd/esd01/data/vbuford/model_runs/';
CascadeFolder='/output/IceCascade/';
Downloadloc=strcat(Serverloc,ModelFolder,CascadeFolder);
% Folder to save in (on computer) uigetdir
% Bin/bash file
% Folder you will launch from:

% Password file Details (This script does not create this file)
Passfile='Password.txt'; % SAVED IN SaveFolder
Passphrase=strcat({'sshpass -f '},{''''},Passfile,{''''});

%%% SCP phrase
port=6307;
username='vbuford';
IPAdd='134.2.5.40';
scpphrase=strcat({' scp -P '}, num2str(port),{' '}, username, '@',IPAdd,':');

%%% create topo_tec_Extract bash file
% on windows notepad, will need \r\n instead of \n
fidvel=fopen(strcat(SaveFolder,'/topo_tec_extractv55'),'w');
fprintf(fidvel,'#!/bin/bash\n');
fprintf(fidvel,'\n');

for i=1:length(timesteps)
fprintf(fidvel,Passphrase{1});
fprintf(fidvel,scpphrase{1});
fprintf(fidvel,Downloadloc);
if timesteps(i)<10
    fprintf(fidvel,strcat('topo_tec_000',num2str(timesteps(i)),'.dat'));
elseif timesteps(i)<100
    fprintf(fidvel,strcat('topo_tec_00',num2str(timesteps(i)),'.dat'));
elseif timesteps(i)<1000
    fprintf(fidvel,strcat('topo_tec_0',num2str(timesteps(i)),'.dat'));
else
fprintf(fidvel,strcat('topo_tec_',num2str(timesteps(i)),'.dat'));
end
fprintf(fidvel,' ./%s\n', LocalSaveFolder);
end

fclose(fidvel)
disp('Finished')

%% Trying to set it up so it doesn't matter which folder you're in...

% Start in home folder. yahoo!

% Change this per model
timesteps=[1580:20:2380];% 1540:20:2380%11.9e6/5000=2380
ModelFolder='Cascade_N9v54';
LocalSaveFolder='McQ02N9/v54-DefL-v73Orog/topo_tec/';

%LocalSaveFolder='McQ02N9/v54-DefL-v73Orog/topo_tec/';
SaveFolder='/Volumes/Files/VictoriaFiles/Cascade';%uigetdir();

% File type: topo_tec
% which timesteps to download: [1580:20:2000]
% Folder on server (eg username + model_runs + Foldername +
% output/IceCascade/
Serverloc='/esd/esd01/data/vbuford/model_runs/';
CascadeFolder='/output/IceCascade/';
Downloadloc=strcat(Serverloc,ModelFolder,CascadeFolder);
% Folder to save in (on computer) uigetdir
% Bin/bash file
% Folder you will launch from:

% Password file Details (This script does not create this file)
Passfile='/Volumes/Files/VictoriaFiles/Cascade/Password.txt'; % SAVED IN SaveFolder
Passphrase=strcat({'sshpass -f '},{''''},Passfile,{''''});

%%% SCP phrase
port=6307;
username='vbuford';
IPAdd='134.2.5.40';
scpphrase=strcat({' scp -P '}, num2str(port),{' '}, username, '@',IPAdd,':');

%%% create topo_tec_Extract bash file
% on windows notepad, will need \r\n instead of \n
fidvel=fopen(strcat(SaveFolder,'/topo_tec_extract'),'w');
fprintf(fidvel,'#!/bin/bash\n');
fprintf(fidvel,'\n');

for i=1:length(timesteps)
fprintf(fidvel,Passphrase{1});
fprintf(fidvel,scpphrase{1});
fprintf(fidvel,Downloadloc);
if timesteps(i)<10
    fprintf(fidvel,strcat('topo_tec_000',num2str(timesteps(i)),'.dat'));
elseif timesteps(i)<100
    fprintf(fidvel,strcat('topo_tec_00',num2str(timesteps(i)),'.dat'));
elseif timesteps(i)<1000
    fprintf(fidvel,strcat('topo_tec_0',num2str(timesteps(i)),'.dat'));
else
fprintf(fidvel,strcat('topo_tec_',num2str(timesteps(i)),'.dat'));
end
fprintf(fidvel,' %s\n', strcat(SaveFolder,'/',LocalSaveFolder));
end

fclose(fidvel);
disp('Finished')
