%% The following script plots and saves a complete CASCADE model run and optionally creates a video from the topography files
% Victoria Buford
%   edited from Paul Eizenhoefer's version
%% Instructions
% 1. Place script in folder where you want to have the script output (e.g. video and MATLAB figure files) to be saved 
% 2. Run script in MATLAB
% 3. Answer command window questions
% 4. Wait (coffee break?)
% 13, 14, 20:30, 32
%% CASCADE output folder and files
close;clc;clear;
disp('Locate the folder containing CASCADE output files topo_tec_****.dat.');
OutputFolder = uigetdir;
%%

prompt={'X-range [xmin xmax] [km]','Y-range [ymin ymax] [km]',...
    'Z-range [zmin zmax] [km]', 'time per file [ka]',...
    'video (0=no; 1=yes)','Model Name'};

defaultans = {'[0 50]','[0 150]','[0 10]','100','1','CASCADE_Video'};
inp= inputdlg(prompt,'Plot Ranges',1,defaultans);
XRange=str2num(inp{1});
YRange=str2num(inp{2});
ZRange=str2num(inp{3});
TimeStep=str2num(inp{4});
Video=str2num(inp{5});
Model_name=(inp{6});

ZRangePlotColor=[0 7.5];

 PrecMax=2;%max(File.precipitation_mpy);
  

%% Plot Figures and get videoframes 
tic

CASCADE_Files = dir([OutputFolder '/*.dat']);
FileNames = {CASCADE_Files.name};
TotalFiles = length(FileNames);

if Video == 1
    Frames = VideoWriter(strcat(OutputFolder(1:end-8), Model_name), 'MPEG-4');
    Frames.Quality = 100;
    Frames.FrameRate =  1000/TimeStep ;% should equate to 1Ma/s
    open (Frames);
end
f=figure('visible','off'); 
set(f,'Position',[1 1 1300 700]); %1280 704

%XRange=[5 45];YRange=[5 190]; 


for i = 1:TotalFiles

%     for j=1:TotalFiles % Eg.... whoops time out = 5ky
%     i=(j-1)*20+1;
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
    Triangulation = delaunay(File.x_km, File.y_km);
    
     % Plot
%     f=figure;
%     set(f,'Position',[1 1 1300 700])
    % 3D Figure.
    subplot(4,3,[1,5])
    trimesh(Triangulation, File.x_km, File.y_km, File.z_km, 'FaceColor', 'flat')
    daspect([1 1 1])
    xlabel('[km]')
    ylabel('[km]')
    zlabel('[km]')
    ax1=gca;
    c1 = colorbar('Location', 'westoutside');
    xlabel(c1,'Elevation [km]')
     
    %%% Timestep
    % use this one for topo_tec0000 = 0Ma
    % Time = (i-1) * TimeStep / 1e3; % paul's time: model start = 0Ma
    
    %%% Timestep
    % use this one for topo_tec**last = 0Ma
     sz = size(CASCADE_Files);
     TimeOrig = TimeStep * sz(1) - TimeStep;
     if (TimeOrig - (TotalFiles-1) * TimeStep) / 1e3 ==0
         Time = (TimeOrig - (i-1) * TimeStep) / 1e3;
     else
         TimeOrig = TimeStep * sz(1);
         Time = (TimeOrig - (i-1) * TimeStep) / 1e3;
     end
     clear TimeOrig
    %%%%
    
    title(strcat(num2str(Time), ' Ma'), 'FontSize',20);
    demcmap([ZRangePlotColor])
    axis([XRange YRange ZRange])
    light('Position',[0 180 50],'Style','local')
    
    % Precipitation Plot
    subplot(4,3,[3,6])
    trimesh(Triangulation, File.x_km, File.y_km, File.z_km, File.precipitation_mpy,  'FaceColor', 'flat')
    daspect([1 1 1])
    c2 = colorbar('Location', 'westoutside');
    xlabel(c2,'Precipitation [m/yr]')
    xlabel('[km]')
    ylabel('[km]')
    zlabel('[km]')
    ax2=gca;
    pcmap=othercolor('Blues4');
    colormap(ax2,pcmap)
    caxis([0 PrecMax])
    axis([XRange YRange ZRange])
    light('Position',[0 180 50],'Style','local')
    view(0,90)
    
    %%%% X-Sec Plot
    % Precipitation
    subplot(4,3,[7,9])
    xx=min(File.x_km):0.5:max(File.x_km);
    yy=min(File.y_km):0.5:max(File.y_km);
    [x,y]=meshgrid(xx,yy);%clear xx yy
    z=griddata(File.x_km,File.y_km,File.z_km,x,y);
    precip_mpy=griddata(File.x_km,File.y_km,File.z_km,File.precipitation_mpy,...
        x,y,z);
    zprec=0:0.5:ZRange(2);
    [zzprec,yyprec]=meshgrid(zprec,yy);
    caxis([0 PrecMax])
    xprec=max(File.x_km)*ones(size(zzprec));
    for k=1:length(y)
        avgprecip(k) = mean(precip_mpy(k,:)) ;
    end
    avgprecip=smoothdata(avgprecip,'movmean',11); % over 5km, per convo w/ Paul
    avgprecipgrid=ones(size(yyprec)).*avgprecip';
    
    s=surf(xprec,yyprec,zzprec,avgprecipgrid,'EdgeColor','none');
    hold on
    ax3=gca;
    
    %pcmap=othercolor('Blues3');
    colormap(ax3,pcmap)
    caxis([0 PrecMax])
    view([-90 0])
    ax3.Color='none';
    daspect([1 1 1])
    axis([0 50 YRange ZRange])
    ax3.YTick=[];
    ax3.ZTick=[];
    
    % Topography
    subplot(4,3,[10,12])
    trimesh(Triangulation, File.x_km, File.y_km, File.z_km, 'FaceColor', 'flat')
        xlabel('[km]')
    ylabel('[km]')
    zlabel('[km] or [m/y]')
    hold on
    
    plot3(zeros(size(avgprecip)),yy,avgprecip,'b','LineWidth',1)
    daspect([1 1 1])
    ax4=gca;
    demcmap([ZRangePlotColor])
    axis([0 50 YRange ZRange])
    light('Position',[0 180 50],'Style','local')
    view(-90,0)
    
    set(ax4,'Color','None')
    box on
    set(ax4,'Position',[0.13 0.11 0.775 0.1577])
    set(ax3,'Position',[0.13 0.11 0.775 0.1577])
  

    
    % Positions of various axes...
     set(ax1,'Position',[0.1164    0.3141    0.6850    0.5782])
     set(c1,'Position',[0.1297    0.3002    0.0144    0.5941])
     set(ax2,'Position',[0.7781    0.2757    0.1745    0.6157])
     ax1.Color='None';
     ax1.YTick=[0:10:YRange(2)];
     ax2.YTick=[0:10:YRange(2)];
     ax4.YTick=[0:10:YRange(2)];
     

    %savefig([OutputFolder strcat(num2str(Time), ' Ma.fig')])
    %save([OutputFolder strcat(num2str(Time), ' Ma.png')])  
    if Video == 1
        set(gcf, 'Position', [0,0,1300,700])
        Frame = getframe(gcf);
        writeVideo(Frames,Frame);
    end
    
        % for last time step, plot swath topography
    if i==TotalFiles
        % import modern topo (McQ02N)
        filename = '/Volumes/Files/VictoriaFiles/Cascade/McQ02NModernTopo.dat';
        delimiter = ' ';
        formatSpec = '%f%f%f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
        fclose(fileID);
        ModernTopo = table(dataArray{1:end-1}, 'VariableNames', {'x','z_min','z_mean','z_max'});
        clearvars filename delimiter formatSpec fileID dataArray ans;
        plot3(zeros(size(ModernTopo.x)),ModernTopo.x-100,ModernTopo.z_min/1e3,'r','LineWidth',1)
        plot3(zeros(size(ModernTopo.x)),ModernTopo.x-100,ModernTopo.z_mean/1e3,'r','LineWidth',1)
        plot3(zeros(size(ModernTopo.x)),ModernTopo.x-100,ModernTopo.z_max/1e3,'r','LineWidth',1)
        if Video == 1
            set(gcf, 'Position', [0,0,1300,700])
            Frame = getframe(gcf);
            writeVideo(Frames,Frame);
        end
    end

    %close(f);
end
if Video == 1
    close (Frames);
end
close
toc