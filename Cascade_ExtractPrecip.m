%% Cascade_ExtractPrecip.m
clc; clear;
disp('Locate the folder containing CASCADE output files topo_tec_****.dat.');
OutputFolder = uigetdir('/Volumes/Files/VictoriaFiles/Cascade/PrecipTesting/OrogPrecip/');
% loccs=strfind(string(OutputFolder),'/v');
% OutputDir=dir([OutputFolder(1:loccs)]);
% name = {OutputDir.name}.';
models=[1];
%%
% write for loop to batch plot models version 1:end
% Basically, just add way for MATLAB to switch the folder by finding the
% appropriate folder name...
directoryfolder=dir(OutputFolder);
for j=1:length(models)
    for k=1:length(directoryfolder)
        ifff=strfind(directoryfolder(k).name,strcat('v',num2str(models(j))));
        if ifff==1
            ModelFolder=directoryfolder(k).name;
        end
    end
    CASCADE_Files = dir([directoryfolder(1).folder '/' ModelFolder '/*.dat']);
end
%%
CASCADE_Files = dir([directoryfolder(1).folder '/*.dat']);
%CASCADE_Files = dir([pwd '/*.dat']);
OutputFolder = CASCADE_Files.folder;
%FileNames = {CASCADE_Files.name};
%TotalFiles = length(FileNames);

i=5;
%%% Load CASCADE output files %%%
if (i-1)>=1000
    File = readtable([OutputFolder '/topo_tec_' num2str(i-1) '.dat'], 'HeaderLines', 4);
elseif (i-1)>=100
    File = readtable([OutputFolder '/topo_tec_0' num2str(i-1) '.dat'], 'HeaderLines', 4);
elseif (i-1)>=10
    File = readtable([OutputFolder '/topo_tec_00' num2str(i-1) '.dat'], 'HeaderLines', 4);
else
    File = readtable([OutputFolder '/topo_tec_000' num2str(i-1) '.dat'], 'HeaderLines', 4);
end

OutputVariables={'x_km','y_km','z_km','node','precipitation_mpy',...
    'fluvial_erosion_rate_mpy','diffusion_erosion_rate_mpy',...
    'landslide_erosion_rate_mpy','total_erosion_rate_mpy',...
    'catchment_color','catchment_number',...
    'glacial_erosion_rate_mpy','ice_thickness_m','mass_balance_1py',...
    'total_topography_m','sliding_velocity_mpy','gerode_term_mpy',...
    'rock_contact_km','isostatic_deflection_mpy','slope_mpkm',...
    'totalflexiso_m','constriction','cumulative_erosion_m',...
    'surface_area_km2'}';
File.Properties.VariableNames=OutputVariables';

% Clean up output - Remove NaNs
Nanss = ismissing(File,{NaN}); % creates table size of A, with 1 or 0 for NaN
NanNodes = Nanss(:,4); clear Nanss % pull only the nodes column
File = File(~NanNodes,:); clear NanNodes % remove any row that the node is NaN

% Compute Delaunay triangulation
%Triangulation = delaunay(File.x_km, File.y_km);

xx=min(File.x_km):0.5:max(File.x_km);
yy=min(File.y_km):0.5:max(File.y_km);
[x,y]=meshgrid(xx,yy);%clear xx yy
z=griddata(File.x_km,File.y_km,File.z_km,x,y);
precip_mpy=griddata(File.x_km,File.y_km,File.z_km,File.precipitation_mpy,...
    x,y,z);
%zprec=0:0.5:ZRange(2);
%[zzprec,yyprec]=meshgrid(zprec,yy);
%xprec=max(File.x_km)*ones(size(zzprec));
for k=1:length(y)
    avgprecip(k) = mean(precip_mpy(k,:)) ;
    minprecip(k) = min(precip_mpy(k,:));
    maxprecip(k) = max(precip_mpy(k,:));
end
%avgprecip=smoothdata(avgprecip,'movmean',11); % over 5km, per convo w/ Paul
%avgprecipgrid=ones(size(yyprec)).*avgprecip';

%% Test plot...
set(0,'DefaultAxesFontSize',14)
plot(y(:,1),avgprecip,'b')
hold on
plot(y(:,1),minprecip,'k')
plot(y(:,1),maxprecip,'r')
legend('avg','min','max')
xlabel('along section [km]')
ylabel('precipitation [m]')
hold off
title(ModelFolder)
%%

% some reason still have NaNs created in precip_mpy...
saveas(gcf,strcat(OutputFolder,'/PrecipExtract.png'))
fname=strcat(OutputFolder,'/PrecipExtract.dat');
precipavgminmax=[y(:,1) avgprecip' minprecip' maxprecip'];
% [dist_along_section avg_precip_km min_precip_km max_precip_km]
dlmwrite(fname,precipavgminmax,'delimiter',' ')
disp('Finished')


%%

% tri=delaunay(File.x_km, File.y_km);
% 
% % Elevation Mesh
% trimesh(tri,File.x_km, File.y_km, File.z_km,'FaceColor','flat')
% view([90 0])