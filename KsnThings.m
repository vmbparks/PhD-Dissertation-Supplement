%% code to import, process DEM to fill sinks/ create streams etc, then create KSN and topo swaths
% Victoria M Buford Parks
% April 1, 2019

% Requires topo toolbox(, and maybe topographic analysis kit?)
clc;close;clear;

% dem='/Volumes/Files/VictoriaFiles/ArcGIS/DEM/GDEMTan-WGSUTM.tif';
DEM8 = GRIDobj(dem);
DEM8.Z(DEM8.Z>10000)=NaN; 
DEM8.Z(DEM8.Z==0)=NaN;

DEMfill8=fillsinks(DEM8);
%you'll have to get the coordinates to crop with...
%DEMfill8=crop(DEMfill,x5,y5);
DEMfill8=crop(DEMfill);

% %%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%
% If you don't want to crop the DEM yourself, just
% load('DEM8_All.mat')
% %%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%--------------------%%%%%%%%%%%%%%%%%%%%%%

disp([datestr(clock) ' -- DEM 8:FlowObj'])
FDm8=FLOWobj(DEMfill8);%,'dinf');
disp([datestr(clock) ' -- DEM 8:FlowAcc'])
A8=flowacc(FDm8);
disp([datestr(clock) ' -- DEM 8:Streams'])
S8=STREAMobj(FDm8,'minarea',4e3);
disp([datestr(clock) ' -- DEM 8:Hydro Condition'])

% Hydrologically condition
theta_ref=0.5;
segment_length=2000;
zc=mincosthydrocon(S8,DEMfill8,'interp',0.2);
DEMc=GRIDobj(DEMfill8);
DEMc.Z(DEMc.Z==0)=NaN;

DEMc.Z(S8.IXgrid)=zc;
% %%
disp([datestr(clock) ' -- DEM 8:KsnComputation'])
[ksn,ksn_ms]=KSN_Quick(DEMfill8,DEMc,A8,S8,theta_ref,segment_length);

%
disp([datestr(clock) ' -- DEM 8:Export to Shapefile for use in Arc'])
shapewrite(ksn_ms,'ksn_S8-2km.shp')

%%
[xx,yy]=getcoordinates(DEMfill8);
[X,Y]=meshgrid(xx,yy);
disp([datestr(clock) ' -- DEM 8:Gridding Ksn'])
ksn_cell=cell(numel(ksn_ms),1);
for ii=1:numel(ksn_ms)
    ksn_cell{ii}=ones(numel(ksn_ms(ii).X),1)*ksn_ms(ii).ksn;
end
ksn_x=vertcat(ksn_ms.X); ksn_y=vertcat(ksn_ms.Y); ksn_ksn=vertcat(ksn_cell{:});

Fk=scatteredInterpolant(ksn_x,ksn_y,ksn_ksn);
ksn_int=Fk(X,Y);
KSNGrid8=GRIDobj(xx,yy,ksn_int);
%
IDX=GRIDobj(DEM,'logical');
KSNGrid8.Z(IDX.Z)=NaN;

M=IDX;
M.Z(~isnan(DEMfill8.Z))=true;
KSNGrid8=crop(KSNGrid8,M,NaN);


%% % The above Data is stored in DEM8_All.mat
% Rak17/McQ08Swath Then SwathOBJ
width=50;% in km
% 1.0e+06 *
% 
%     0.5666    8.2086
%     0.8113    8.3613
%     0.8113    8.3613
 x=[0.5666;0.8113]*1e6;
 y=[8.2086;8.3613]*1e6;
 disp([datestr(clock) ' -- DEM 8:Swath Topo'])
 swathtopo=SWATHobj(DEMfill8,x,y,'width',width*1e3);
%swathtopo=SWATHobj(DEMfill8,'width',width*1e3);
disp([datestr(clock) ' -- DEM 8:Swath Ksn'])
swathKsn=SWATHobj(KSNGrid8,swathtopo.xy0(1:2,1),swathtopo.xy0(1:2,2),'width',width*1e3);

% %%
figure(1)
imagesc(DEMfill8)
demcmap([0 6500])
colorbar
hold on
plot(swathtopo.X,swathtopo.Y,'r') %shows where swath data is from
plot(swathtopo.xy0(1:2,1),swathtopo.xy0(1:2,2),'k')%Center/sectionline

hold off


% %%
figure(2)
subplot(2,1,1)
plotdz(swathtopo)
axis([0 3e5 0 6500])
title('Swath Topo')

%axis([0 3e5 0 6500])
hold on
%yyaxis right
subplot(2,1,2)
plotdz(swathKsn)
ylabel('Ksn []')
axis([0 3e5 0 1000])
% %%

%% Get Coords of Box
MinX=min(min(swathtopo.X));
ID=find(swathtopo.X==MinX);
MinY=swathtopo.Y(ID);

MinYY=min(min(swathtopo.Y));
ID=find(swathtopo.Y==MinYY);
MinXX=swathtopo.X(ID);

MaxX=max(max(swathtopo.X));
ID=find(swathtopo.X==MaxX);
MaxY=swathtopo.Y(ID);

MaxYY=max(max(swathtopo.Y));
ID=find(swathtopo.Y==MaxYY);
MaxXX=swathtopo.X(ID);
figure(3)
imagesc(DEMfill8)
demcmap([0 6500])
hold on
plot([MinX,MinXX,MaxX,MaxXX,MinX],[MinY,MinYY,MaxY,MaxYY,MinY],'-')%,'MarkerSize',10)
hold off

%%


sz=size(swathtopo.Z);
mean_Sw=zeros(sz(2),1);
median_Sw=zeros(sz(2),1);
mode_Sw=zeros(sz(2),1);
std_Sw=zeros(sz(2),1);
min_Sw=zeros(sz(2),1);
max_Sw=zeros(sz(2),1);

for i=1:sz(2)
    mean_Sw(i)=mean(swathtopo.Z(:,i));
    median_Sw(i)=median(swathtopo.Z(:,i));
    mode_Sw(i)=mode(swathtopo.Z(:,i));
    std_Sw(i)=std(swathtopo.Z(:,i));
    min_Sw(i)=min(swathtopo.Z(:,i));
    max_Sw(i)=max(swathtopo.Z(:,i));
end
% %%
mode_SwAvg=zeros(sz(2),1);
%mean over 3 widths
mode_width=8;
for i=1+mode_width:sz(2)-(1+mode_width)
    %mean_Sw(i)=mean(Perez16.Z(:,i));
    %median_Sw(i)=median(Perez16.Z(:,i));
    datta=swathtopo.Z(:,i-mode_width:i+mode_width);
    mode_SwAvg(i)=mode(datta(:));
    %std_Sw(i)=std(Perez16.Z(:,i));
end
%%
sz=size(swathKsn.Z);
mean_SwKsn=zeros(sz(2),1);
median_SwKsn=zeros(sz(2),1);
mode_SwKsn=zeros(sz(2),1);
std_SwKsn=zeros(sz(2),1);
min_SwKsn=zeros(sz(2),1);
max_SwKsn=zeros(sz(2),1);

for i=1:sz(2)
    mean_SwKsn(i)=mean(swathKsn.Z(:,i));
    median_SwKsn(i)=median(swathKsn.Z(:,i));
    mode_SwKsn(i)=mode(swathKsn.Z(:,i));
    std_SwKsn(i)=std(swathKsn.Z(:,i));
    min_SwKsn(i)=min(swathKsn.Z(:,i));
    max_SwKsn(i)=max(swathKsn.Z(:,i));
end

% %%
mode_SwAvgKsn=zeros(sz(2),1);
%mean over 3 widths
mode_width=8;
for i=1+mode_width:sz(2)-(1+mode_width)
    %mean_Sw(i)=mean(Perez16.Z(:,i));
    %median_Sw(i)=median(Perez16.Z(:,i));
    datta=swathKsn.Z(:,i-mode_width:i+mode_width);
    mode_SwAvgKsn(i)=mode(datta(:));
    %std_Sw(i)=std(Perez16.Z(:,i));
end

%%
f3.Renderer='painters'; % to ensure that the file is editable in AI.
saveas(f3,strcat('/Volumes/Files/AreaMap'),'epsc')

f2.Renderer='painters'; % to ensure that the file is editable in AI.
saveas(f2,strcat('/Volumes/Files/Swaths'),'epsc')



%%
function [ksn,ksn_ms]=KSN_Quick(DEM,DEMc,A,S,theta_ref,segment_length)
	g=gradient(S,DEMc);
	G=GRIDobj(DEM);
	G.Z(S.IXgrid)=g;

	Z_RES=DEMc-DEM;

	ksn=G./(A.*(A.cellsize^2)).^(-theta_ref);
	
	ksn_ms=STREAMobj2mapstruct(S,'seglength',segment_length,'attributes',...
		{'ksn' ksn @mean 'uparea' (A.*(A.cellsize^2)) @mean 'gradient' G @mean 'cut_fill' Z_RES @mean});
end