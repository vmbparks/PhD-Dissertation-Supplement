%% PecubeOutputPlot_batch_DLG v1.6
% Victoria M Buford
% vmb21@pitt.edu
% updated May 11, 2018

% Now with Goodness of Fit!
% AND accurate checks for # of files when not plotting temp &/or Velocity!

% Batch plot all steps fr   om Pecube functionality. 

% Assumptions
    % velocities and structures start at Step0
    % NAMING SCHEME for lines: model_lines_Step***.dat (starting at 0)
                %     vel_new: model_grid_Step***_vel_new.dat (starting at 0)
    % Thermochronometer data must have the following columns:
        % Method, Age, Error, xdist 
        % needs a header, with ^^ in it
        % must be tab sepa44
% Required Information:
    % Folder locations where the following are saved:
        % Structures   (eg Move lines: [ x y z ID ] .dat (space separated)      
        % Ages         (Ages_tec****.dat from Pecube output)    
        % Temperatures (temps_tec***.dat from Pecube output)           
        % Velocities   (**_vel_new**.dat from Pecube output)
        % Save         (folder to save output figures)
    % Velocity Scheme
        % constant, late quick, etc. some txt or str ID for you.
    % Smoothing window: 
        % If you want the age data averaged over a number of points
    % Ao 
        % heat production value (don't actually need this; just labels
        % model more thoroughly)
    % extension  
        % what format you want the figure saved in (eg .png, .epsc)
    % Model       
        % Model Name (will save in plot names)
    % XRange
    	% cross section distance
    % ZRange
    	% Range of depth/elevation to plot
    % Age Range
        % Ages to plot (over which the Tchron Data/ model ran)
    % Which thermochronometers you want to plot
        % MAr, ZFT, ZHe, AFT, and AHe are options
                       
%% 
clear; close all; clc;
%%
set(0,'DefaultUicontrolFontSize',14)
%%

% Define Default Answers
if exist('Answer','var')==1
    LinesDir    = Answer.MoveLinesFolder;
    VelDir      = Answer.VelocityFolder;
    TempDir     = Answer.TempsFolder;
    SaveDir     = Answer.SaveFolder;
    AgesDir     = Answer.AgesFolder;
    Modeldef    = Answer.Model;
    IDdef       = Answer.ID;
    XRdef       = Answer.XRange;
    ZRdef       = Answer.ZRange;
    Agedef      = Answer.AgeRange;
    Figdef      = Answer.extension;
    MArdef      = Answer.tchron{1,2};
    ZFTdef      = Answer.tchron{2,2};
    ZHedef      = Answer.tchron{3,2};
    AFTdef      = Answer.tchron{4,2};
    AHedef      = Answer.tchron{5,2};
    BArdatdef   = Answer.tchron{1,3};
    MArdatdef   = Answer.tchron{2,3};
    ZFTdatdef   = Answer.tchron{3,3};
    ZHedatdef   = Answer.tchron{4,3};
    AFTdatdef   = Answer.tchron{5,3};
    AHedatdef   = Answer.tchron{6,3};
    tmethdef    = Answer.methodt{1};
    winddef     = Answer.methodt{2};
    timedef     = Answer.time{1};
    rangedef    = Answer.time{2};
    leglocdef   = Answer.LegLoc;
    TempContDef = Answer.TempContour;
    TchronDataFileDef = Answer.TchronDataFile;
    VelDef      = Answer.tempvel{1,2};
    TempDef     = Answer.tempvel{2,2};
else 
    LinesDir    = pwd;
    VelDir      = pwd;
    TempDir     = pwd;
    SaveDir     = pwd;
    AgesDir     = pwd;
    Modeldef    = 'Model';
    IDdef       = 'const';
    XRdef       = [0;600];
    ZRdef       = [-50;7.5];
    Agedef      = [0;100];
    Figdef      = '.png';
    BArdef      = true;
    MArdef      = true;
    ZFTdef      = true;
    ZHedef      = true;
    AFTdef      = true;
    AHedef      = true;
    BArdatdef   = true;
    MArdatdef   = true;
    ZFTdatdef   = true;
    ZHedatdef   = true;
    AFTdatdef   = true;
    AHedatdef   = true;
    tmethdef    = 'movmean';
    winddef     = 1;
    timedef     = 'All';
    rangedef    = '18:23';
    leglocdef   = 'right';
    TempContDef = true;
    TchronDataFileDef='/Users/victoriabuford/Box Sync/Research/MoveStuff/SPeru/Tchron/SP1Tchr.txt';%pwd;
    VelDef      = true;
    TempDef     = true;
end
    
Title = 'Pecube Batch Plotting';

%%%% Setting Up the Dialog Window
% Options.WindowStyle = 'modal';
Options.Resize = 'on';
Options.Interpreter = 'tex';
Options.CancelButton = 'on';
Options.ApplyButton = 'off';
Options.ButtonNames = {'Go','Cancel'}; 
Option.Dim = 4; % Horizontal dimension in fields

Prompt = {};
Formats = {};
DefAns = struct([]);

%%% Prompts in the Dialog Window

Prompt(1,:) = {['Select Folders and Parameters.'],[],[]};
Formats(1,1).type = 'text';
Formats(1,1).size = [-1 0];
Formats(1,1).span = [1 2]; % item is 1 field x 4 fields

Prompt(end+1,:) = {'Model Name', 'Model',[]};
Formats(2,1).type = 'edit';
Formats(2,1).format = 'text';
Formats(2,1).size = 200; % automatically assign the height
DefAns(1).Model = Modeldef;

Prompt(end+1,:)={'Model ID','ID',[]};
Formats(2,2).type = 'edit';
Formats(2,2).format = 'text';
Formats(2,2).span=[1 1];
DefAns.ID = IDdef;

Prompt(end+1,:)={'Deformation front on', 'LegLoc',[],};
Formats(2,4).type='edit';
Formats(2,4).format='text';
Formats(2,4).size=[35 25];
DefAns.LegLoc=leglocdef;

Prompt(end+1,:) = {'Move Lines Folder','MoveLinesFolder',[]};
Formats(4,1).type = 'edit';
Formats(4,1).format = 'dir';
Formats(4,1).size = [-1 0];
Formats(4,1).span = [1 3];  % item is 1 field x 3 BAfields
DefAns.MoveLinesFolder = LinesDir;

Prompt(end+1,:)= {'     X Range','XRange',[]};
Formats(4,4).type='edit';
Formats(4,4).format = 'vector';
Formats(4,4).size = [60 50];
DefAns.XRange = XRdef';
%DefAns.XRange = [300 1050];

Prompt(end+1,:) = {'Velocity Folder','VelocityFolder',[]};
Formats(5,1).type = 'edit';
Formats(5,1).format = 'dir';
Formats(5,1).size = [-1 0];
Formats(5,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.VelocityFolder = VelDir;

Prompt(end+1,:)= {'     Z Range','ZRange',[]};
Formats(5,4).type='edit';
Formats(5,4).format = 'vector';
Formats(5,4).size = [60 50];
DefAns.ZRange = ZRdef';
%DefAns.ZRange = [-50 7.5];

Prompt(end+1,:) = {'Ages Folder','AgesFolder',[]};
Formats(6,1).type = 'edit';
Formats(6,1).format = 'dir';
Formats(6,1).size = [-1 0];
Formats(6,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.AgesFolder = AgesDir;

Prompt(end+1,:)= {' Age Range','AgeRange',[]};
Formats(6,4).type='edit';
Formats(6,4).format = 'vector';
Formats(6,4).size = [60 50];
DefAns.AgeRange = Agedef'; %[0 100];

Prompt(end+1,:) = {'Pecube Temperatures Folder','TempsFolder',[]};
Formats(7,1).type = 'edit';
Formats(7,1).format = 'dir';
Formats(7,1).size = [-1 0];
Formats(7,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.TempsFolder = TempDir;

Prompt(end+1,:) = {'Save Folder','SaveFolder',[]};
Formats(8,1).type = 'edit';
Formats(8,1).format = 'dir';
Formats(8,1).size = [-1 0];
Formats(8,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.SaveFolder = SaveDir;

Prompt(end+1,:) = {'Figure Format','extension',[]};
Formats(8,4).format = 'text';
Formats(8,4).size = [45 25];
DefAns.extension = Figdef;

Prompt(end+1,:)={'Tchron Smoothing Method','methodt',[]};
Formats(9,1).type='table';
Formats(9,1).format={{'movmean','movmedian','gaussian','lowess',...
    'rloess','rloess','sgolay'} 'numeric'}; % see matlab help on function smoothdata
Formats(9,1).items={'Method','window size'};
Formats(9,1).size=[200 45];
%Formats(9,1).span=[1,2];
DefAns.methodt={tmethdef winddef}; %tmethdef;

Prompt(end+1,:) = {'Thermochronometers to plot','tchron',[]};
Formats(9,2).type = 'table';
Formats(9,2).format = {'char', 'logical', 'logical'}; 
Formats(9,2).items = {'Row' 'Modeled' 'data T=end'};
Formats(9,2).size = [260 160];
Formats(9,2).span = [3 1];  % item is 3 field x 1 fields
DefAns.tchron = {'BAr' BArdef BArdatdef
                 'MAr' MArdef MArdatdef
                 'ZFT' ZFTdef ZFTdatdef
                 'ZHe' ZHedef ZHedatdef
                 'AFT' AFTdef AFTdatdef
                 'AHe' AHedef AHedatdef};

Prompt(end+1,:)={'Plot?            ','tempvel',[]};
Formats(9,3).type='table'; Formats(9,3).format={'char','logical'};
Formats(9,3).items={'Row' 'Plot'};
Formats(9,3).size=[200 90];
Formats(9,3).span=[1 2];
DefAns.tempvel={'Velocity' VelDef
                'Temperature' TempDef}; 
             
Prompt(end+1,:)={'               Timesteps to plot','time',[]};
Formats(10,1).type='table';
Formats(10,1).format={{'All','Range'},'numeric'};
Formats(10,1).items={'Timesteps','From X:Y'};
Formats(10,1).size=[200 45];
DefAns.time={timedef,rangedef};
% 
Prompt(end+1,:)={'Contour Temperature?' 'TempContour',[]};
Formats(10,4).type='check'; 
%Formats(10,4).format={'char','logical'};
Formats(10,4).size=[200 45];
%Formats(10,4).span=[1 2];
DefAns.TempContour=logical(TempContDef);

MaxTemp=600;

Prompt(end+1,:)={'Tchron data file [txt]','TchronDataFile',[]};
Formats(11,1).type = 'edit';
Formats(11,1).format = 'file';
Formats(11,1).size = [250 20];
% Formats(11,1).span = [1 2];  % item is 1 field x 3 fields
DefAns.TchronDataFile = TchronDataFileDef;

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

%% Answers: 
% Answer.Model              str
% Answer.ID                 str
% Answer.LegLoc             str
% Answer.MoveLinesFolder    str
% Answer.XRange             vector
% Answer.VelocityFolder     str
% Answer.ZRange             vector
% Answer.AgesFolder         str
% Answer.AgeRange           vector
% Answer.TempsFolder        str
% Answer.SaveFolder         str
% Answer.extension          str
% Answer.methodt            1x2 cell
% Answer.tchron             6x2 cell
% Answer.time               1x2 cell
% Answer.tempvel            2x2 cell

%% Interpret Modal Dialog input
Folder_Structures       = Answer.MoveLinesFolder;
Folder_Ages             = Answer.AgesFolder;
Folder_Temperatures     = Answer.TempsFolder;
Folder_Save             = Answer.SaveFolder;
Folder_Velocities       = Answer.VelocityFolder;

window_size = Answer.methodt{2};
method      = Answer.methodt{1};
ID          = Answer.ID;
extension   = Answer.extension;
Model       = Answer.Model;

% Ranges
XMin = Answer.XRange(1);
XMax = Answer.XRange(2);
YMin = Answer.AgeRange(1);
YMax = Answer.AgeRange(2);
ZMin = Answer.ZRange(1);
ZMax = Answer.ZRange(2);

% Thermochronometers
BAr = Answer.tchron{1,2};
MAr = Answer.tchron{2,2};
ZFT = Answer.tchron{3,2};
ZHe = Answer.tchron{4,2};
AFT = Answer.tchron{5,2};
AHe = Answer.tchron{6,2};

% Thermochronometer Data
TchronDatafilename=Answer.TchronDataFile;


% Legend Location
if strcmp(Answer.LegLoc,'right')
    legloc='NorthWest';
elseif strcmp(Answer.LegLoc,'left')
    legloc='NorthEast';
else
    disp 'Legend location must be left or right.'
end

% Plotting temperatures or velocities
if strcmp(Answer.tempvel{1,2},'true')==1 | Answer.tempvel{1,2}==1
    velplot=1;
else
    velplot=0;
end

if strcmp(Answer.tempvel{2,2},'true')==1 | Answer.tempvel{2,2}==1
    tempplot=1;
else
    tempplot=0;
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Import measured tchron data & Canyon info
TchronDataFile='/Users/victoriabuford/Box Sync/Research/MoveStuff/SPeru/Tchron/SP1Tchr-3.txt';
TchronData=readtable(TchronDataFile,'Delimiter','\t');
% filename = TchronDataFile;
% delimiter = '\t';
% startRow = 2;
% endRow = 109;
% formatSpec = '%f%q%q%q%f%f%f%q%f%f%q%f%f%f%q%q%q%[^\n\r]';
% fileID = fopen(filename,'r');
% dataArray = textscan(fileID, formatSpec, endRow-startRow+1, 'Delimiter', delimiter, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines', startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
% fclose(fileID);
% TchronData = table(dataArray{1:end-1}, 'VariableNames', {'MethodNum','Method','Sample_ID','Sample__','S_Lat__WGS','W_Long__WG','Elevation','Fm_age','Age','Error','Study','xdist','indCanyon','indInterfluve','Model_grid','Model_velocity','Model_incisionAge'});
% clearvars filename delimiter startRow endRow formatSpec fileID dataArray ans;
% % %TchronData(109,:)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% indCanyon=1 when samples *are* in the canyon
BArind = strcmp(TchronData.Method,'BAr')-~TchronData.indCanyon;
BArind = BArind>0;
MArind = strcmp(TchronData.Method,'MAr')-~TchronData.indCanyon;
MArind = MArind>0;
ZFTind = strcmp(TchronData.Method,'ZFT');
ZHeind = strcmp(TchronData.Method,'ZHe')-~TchronData.indCanyon;%-~indZr
ZHeind = ZHeind>0;
AFTind = strcmp(TchronData.Method,'AFT')-~TchronData.indCanyon;
AFTind = AFTind>0;
AHeind = strcmp(TchronData.Method,'AHe')-~TchronData.indCanyon;%-~indAp
AHeind = AHeind>0;
% Matrix of thermochronometers sorted by Method
indexCanyon = [BArind MArind ZFTind ZHeind AFTind AHeind];

BArind = strcmp(TchronData.Method,'BAr')-~TchronData.indInterfluve;
BArind = BArind>0;
MArind = strcmp(TchronData.Method,'MAr')-~TchronData.indInterfluve;
MArind = MArind>0;
ZFTind = strcmp(TchronData.Method,'ZFT');
ZHeind = strcmp(TchronData.Method,'ZHe')-~TchronData.indInterfluve;%-~indZr
ZHeind = ZHeind>0;
AFTind = strcmp(TchronData.Method,'AFT')-~TchronData.indInterfluve;
AFTind = AFTind>0;
AHeind = strcmp(TchronData.Method,'AHe')-~TchronData.indInterfluve;%-~indAp
AHeind = AHeind>0;
% Matrix of thermochronometers sorted by Method
indexInterfluve = [BArind MArind ZFTind ZHeind AFTind AHeind];

%%%% New, all samples that are interfluve or canyon
BArind = strcmp(TchronData.Method,'BAr')-indexCanyon(:,1)-indexInterfluve(:,1);
BArind = BArind > 0;
MArind = strcmp(TchronData.Method,'MAr')-indexCanyon(:,2)-indexInterfluve(:,2);
MArind = MArind > 0;
ZFTind = strcmp(TchronData.Method,'ZFT')-indexCanyon(:,3)-indexInterfluve(:,3);
ZFTind = ZFTind > 0;
ZHeind = strcmp(TchronData.Method,'ZHe')-indexCanyon(:,4)-indexInterfluve(:,4);
ZHeind = ZHeind > 0;
AFTind = strcmp(TchronData.Method,'AFT')-indexCanyon(:,5)-indexInterfluve(:,5);
AFTind = AFTind > 0;
AHeind = strcmp(TchronData.Method,'AHe')-indexCanyon(:,6)-indexInterfluve(:,6);
AHeind = AHeind > 0;
% Matrix of thermochronometers sorted by Method
indexMean = [BArind MArind ZFTind ZHeind AFTind AHeind];



%%

mmmm_Id={'mean','canyon','interfluve'};
mmmm_plot={'All','range','range'};
heatpro={'2.0','2.5','3.0'};

tic
% for jkjk=1:length(heatpro) % over heat production
%     IDXages=strfind(Folder_Ages,'ef');
%     Folder_Ages=strcat(Folder_Ages(1:IDXages+8),(heatpro{jkjk}),Folder_Ages(IDXages+12:end));
%     IDXSave=strfind(Folder_Save,'ef');
%     Folder_Save=strcat(Folder_Save(1:IDXSave+8),heatpro{jkjk},Folder_Save(IDXSave+12:end));
%     IDXTemps=strfind(Folder_Temperatures,'ef');
%     Folder_Temperatures=strcat(Folder_Temperatures(1:IDXTemps+8),heatpro{jkjk},Folder_Temperatures(IDXTemps+12:end));
%     
%     IDXId=strfind(ID,'ef');
%     ID=strcat(ID(1:IDXId+8),heatpro{jkjk},ID(IDXId+12:end));


for lllll=[31] %SP4:[10:18] %SP5:[10,12,14,19:29] % for velocities
   % Folder_Temperatures=strcat(Folder_Temperatures(1:42),num2str(lllll),Folder_Temperatures(45:end));
    %Folder_Velocities=strcat(Folder_Velocities(1:end-10),num2str(lllll),Folder_Velocities(end-7:end));
    Folder_Ages=strcat(Folder_Ages(1:42),num2str(lllll),Folder_Ages(45:end));
    Folder_Save=strcat(Folder_Save(1:42),num2str(lllll),Folder_Save(45:end));
    ID=strcat(ID(1:2),num2str(lllll),ID(5:end));
    
    Answer.time{1}='All';
    
    for mmm=1:3 % e.g. canyon, interfluve, etc
        if mmm==1
            Folder_Ages=strrep(Folder_Ages,mmmm_Id{mmm+1},mmmm_Id{mmm});
            Folder_Ages=strrep(Folder_Ages,mmmm_Id{mmm+2},mmmm_Id{mmm});
            Folder_Save=strrep(Folder_Save,mmmm_Id{mmm+1},mmmm_Id{mmm});
            Folder_Save=strrep(Folder_Save,mmmm_Id{mmm+2},mmmm_Id{mmm});
            ID = strrep(ID,mmmm_Id{mmm+1},mmmm_Id{mmm})
            ID = strrep(ID,mmmm_Id{mmm+2},mmmm_Id{mmm})
            Answer.time{1}='All';
        elseif mmm>1
            Folder_Ages=strrep(Folder_Ages,mmmm_Id{mmm-1},mmmm_Id{mmm})
            Folder_Save=strrep(Folder_Save,mmmm_Id{mmm-1},mmmm_Id{mmm});
            ID = strrep(ID,mmmm_Id{mmm-1},mmmm_Id{mmm})
            Answer.time{1}='Range';
        end
       
        % %% Determine file names and number of files / plots
        % Pecube ages
        Files_Ages = dir([Folder_Ages '/*.dat']);
        FileNames_Ages = {Files_Ages.name};
        TotalFiles_Ages = length(FileNames_Ages);
%        TotalFiles=TotalFiles_Ages;

        TotalFiles = TotalFiles_Ages;
        i=TotalFiles;
        disp('Processing')
        % tic
        for i = TotalFiles % 1:TotalFiles
            
            % Import and Clean Up Age data
            if i>=1000
                File_Age = importdata([Folder_Ages '/Ages_tec' num2str(i) '.dat'], ' ', 4);
            elseif i>=100
                File_Age = importdata([Folder_Ages '/Ages_tec0' num2str(i) '.dat'], ' ', 4);
            elseif i>=10
                File_Age = importdata([Folder_Ages '/Ages_tec00' num2str(i) '.dat'], ' ', 4);
            else
                File_Age = importdata([Folder_Ages '/Ages_tec000' num2str(i) '.dat'], ' ', 4);
            end
            File_Age.data((File_Age.data(:,4) > XMax),:) = [];
            File_Age.data((File_Age.data(:,4) < XMin),:) = [];
            File_Age.data((File_Age.data(:,5) ~= 2),:) = [];
            

            
            A=size(File_Age.data);
            for j=1:A(2)
                File_Age.data(isnan(File_Age.data(:,j)),:) = [];
            end
            
            % Sort Ages for line plot
            File_Age.data = sortrows(File_Age.data,4); % x.
            
            % Compute moving mean
            Age_Avg=[zeros(size(File_Age.data))];
            % id ; x ; y ; real x ; real y ; real z ; AHe ; AFT ; ZHe ; ZFT ; MAr;
            Age_Avg(:,4) = File_Age.data(:,4);
            %window_size=7;
            for k=7:A(2)
                Age_Avg(:,k)=smoothdata(File_Age.data(:,k),method,window_size);
            end
            %Age_Avg=File_Age.data;
            

            legdata={'BAr data','MAr data','ZFT data','ZHe data','AFT data','AHe data'}; %legend

            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% GOF %% Calc for each canyon/interfluve
        if strcmp (mmmm_Id{mmm},'canyon')==1
            index=indexCanyon;
        elseif strcmp (mmmm_Id{mmm},'interfluve')==1
            index=indexInterfluve;
        else
            index=indexMean;
        end
        [d,e]=size(index);
        if i==TotalFiles
            %%% Goodness Of Fit
            GOF_version=1.1;
            samplelocll=zeros(d,e);
            ModelVal=zeros(d,e);
            Matches=zeros(d,e);
            
            ModelError=1; % in [Ma]
            for mm=2 % to calc mult xsec km errors
                XsecError=mm; % in [km]
                for g=1:e % AFT then ZHe
                    if g==6 %AHe
                        ii=7;
                    elseif g==5 %AFT
                        ii=8;
                    elseif  g==4 %ZHe
                        ii=9; % which index value tchronometer I want
                    elseif g==3 % ZFT
                        ii=10;
                    elseif g==2 %MAr
                        ii=11;
                    elseif g==1 %BAr
                        ii=12;
                    end
                    
                    % for each sample in table
                    for jj=1:d
                        % if that sample exists in this thermochronometer
                        if index(jj,g)==1
                            Matches(jj,g)=3;
                            %%%%%% Variables
                            %%%% Measured Tchron Value
                            %%% TchronData.Age(jj)
                            
                            %%%% Measured Tchron Error
                            %%% TchronData.Error(jj)
                            
                            %%%% Sample Location (Adjusted to Pecube Grid)
                            %%% TchronData.xdist(jj)
                            sampleloc(jj,g)=round(TchronData.xdist(jj)*2)/2;
                            
                            % Loop over X-error in cross-section
                            % Without looping over 1km wide window, N10const had 28.57%
                            % fit. With 1km error, 42.86% fit.
                            ll=[sampleloc(jj,g)-(XsecError): 0.5: sampleloc(jj,g)+(XsecError)];
                            
                            for lll= 1:length(ll)
                                samplelocll(jj,g)=ll(lll);
                                
                                % Find Sample Location in MODELLED
                                IndexSampleLocMod=find(Age_Avg(:,4)==samplelocll(jj,g));
                                % Get Modelled Value at that location
                                ModelVal(jj,g)=Age_Avg(IndexSampleLocMod,ii);
                                
                                % Compare Modelled to Measured
                                %%%% To accurately compare, we have 4 cases to test:
                                % Case 1. Model + Model Error within measured Age
                                % Case 2. Model - Model Error within measured Age
                                % Case 3. Measured + Measured error within modelled Age
                                % Case 4. Measured - Measured error within modelled Age
                                
                                % Case 1. Model Upper Error within measured Age
                                if (ModelVal(jj,g) + ModelError)< (TchronData.Age(jj) + TchronData.Error(jj)) ...
                                        && ...
                                        (ModelVal(jj,g) + ModelError)> (TchronData.Age(jj) - TchronData.Error(jj))
                                    Matches(jj,g)=1; % it does match
                                    
                                    % Case 2. Model Lower Error within measured Age
                                elseif (ModelVal(jj,g) - ModelError)< (TchronData.Age(jj) + TchronData.Error(jj)) ...
                                        && ...
                                        (ModelVal(jj,g) - ModelError)> (TchronData.Age(jj) - TchronData.Error(jj))
                                    Matches(jj,g)=1; % it does match
                                    
                                    % Case 3. Measured Upper error within modelled Age
                                elseif (TchronData.Age(jj) + TchronData.Error(jj)) < (ModelVal(jj,g) + ModelError) ...
                                        && ...
                                        (TchronData.Age(jj) + TchronData.Error(jj)) >(ModelVal(jj,g) - ModelError)
                                    Matches(jj,g)=1; % it does match
                                    
                                    % Case 4. Measured Lower error within modelled Age
                                elseif (TchronData.Age(jj) - TchronData.Error(jj)) < (ModelVal(jj,g) + ModelError) ...
                                        && ...
                                        (TchronData.Age(jj) - TchronData.Error(jj)) >(ModelVal(jj,g) - ModelError)
                                    Matches(jj,g)=1; % it does match
                                    % No else statement needed bc line above defines
                                    % Matches(jj,g)=3 %Not Matching and this loop only
                                    % overwrites that if the sample is within the error
                                end
                            end
                        end
                    end
                end
                
                % Percent Fit
                YesMatches=[];
                NoMatches=[];
                YesMatches=find(Matches==1);
                NoMatches=find(Matches==3);
                PercFit=length(YesMatches)/(length(NoMatches)+length(YesMatches))
                
                YesMatchesTchron=[];
                nsamples=[];
                for kk=1:e
                    YesM=length((find(Matches(:,kk)==1)));
                    NoM=length((find(Matches(:,kk)==3)));
                    nsamples(kk)=YesM+NoM;
                    YesMatchesTchron(kk)=YesM/(NoM+YesM);
                end
                
                %%% Print goodness of fit to bulk file
                %     % Model_Name	Model_ID	Model_Error	Sample_Loc_Error	%Fit
                %     GOFLoc='/Volumes/Files/VictoriaFiles/Pecube/GoodnessOfFit.txt';
                %     fidGOF=fopen(GOFLoc,'at');
                %     fprintf(fidGOF,'\n%10s\t%18s\t%0.1f\t%16.1f\t%7.4f',Model,ID,ModelError,XsecError,PercFit);
                %     fclose(fidGOF);
                
                
                %%%% Print goodness of fit + matches + sample info to individual file
                GOFLoc_Ind=strcat(Folder_Save(1:end-9),'/GoodnessOfFit_',Model,'_',ID,'_',num2str(mm),'.txt');
                fidGOF_ind=fopen(GOFLoc_Ind,'w');
                fprintf(fidGOF_ind,'Goodness of Fit Computation; Version %1.1f\n',GOF_version);
                fprintf(fidGOF_ind,'Model_Name\tModel_ID\t\t\tModel_Error\tSample_Loc_Error\tPercFit\n');
                fprintf(fidGOF_ind,'%s\t\t%s\t\t%0.1f\t\t%f\t%f\n',Model,ID,ModelError,XsecError,PercFit);
                fprintf(fidGOF_ind,'\n');
                fprintf(fidGOF_ind,'Summary \tBAr\tMAr\tZFT\tZHe\tAFT\tAHe\n');
                fprintf(fidGOF_ind,'PercFit \t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n',YesMatchesTchron');
                %fprintf(fidGOF_ind,'\n');
                fprintf(fidGOF_ind,'No. samples \t%0.0f\t%0.0f\t%0.0f\t%0.0f\t%0.0f\t%0.0f\n',nsamples');
                fprintf(fidGOF_ind,'\nMatches [n=%.0f; e.g. No=3; Yes=1]\n',(length(NoMatches)+length(YesMatches)));
                fprintf(fidGOF_ind,'\n');
                fprintf(fidGOF_ind,'Sample \t%s \t%s \t%s \t%s \t%s \t%s \n',legdata{:});
                for tt=1:length(Matches)
                    fprintf(fidGOF_ind,'%10s \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \t%.0f \n',TchronData.Sample_ID{tt}, Matches(tt,:));
                    %fprintf(fidGOF_ind,'%s \t\t%.0f \t\t%.0f \t\t%.0f \t\t%.0f \t\t%.0f \n',num2str(tt), Matches(tt,:));
                end
                fclose(fidGOF_ind);
                

                MatchSum=sum(Matches');
                      
                GOFMaster='/Volumes/Files/VictoriaFiles/Pecube/SP5/GOFMaster2.txt';
                fidGOF_Master=fopen(GOFMaster,'a');
%                 fidGOF_Master=fopen(GOFMaster,'w');
%                 
%                 fprintf(fidGOF_Master,'Goodness of Fit Computation; Version %1.1f\n',GOF_version);
%                 
%                 fmt=['%s ' repmat(' %s',1,numel(TchronData.Sample_ID)) '\n'];
%                 fprintf(fidGOF_Master,fmt,'Model',TchronData.Sample_ID{:})
%                 
                
                fmt=['%s =' repmat(' %1.0f',1,numel(MatchSum)) '\n'];
                fprintf(fidGOF_Master,fmt,strcat(Model,'_',ID),MatchSum)
                fclose(fidGOF_Master)

                %fprintf(fidGOF_Master,'\t%.0f')
          
%                 
%                 for tt=1:length(MatchSum)
%                     fprintf(fidGOF_Master,'%10s \t%.0f \n',TchronData.Sample_ID{tt}, MatchSum(tt));
%                     %fprintf(fidGOF_ind,'%s \t\t%.0f \t\t%.0f \t\t%.0f \t\t%.0f \t\t%.0f \n',num2str(tt), Matches(tt,:));
%                 end
%                 
            end
            %%%%% End GOF
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        
    end
% end
end
% disp('Finished')
t=toc;
disp('Finished; time elapsed: ')
disp(datestr(datenum(0,0,0,0,0,t),'HH:MM:SS'))
%%
load handel; sound(y,2.75*Fs) % full list of sounds help audiovideo

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Goodness Of Fit %%%%% & Final Plot %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if i==TotalFiles
   % XMax=1035;XMin=625; MAr=false;AHe=false;ZFT=false;legloc='NorthEast';
   %%% Plot ages %%%
    f=figure('visible','off');  % to make loop run faster: comment out for testing.
    subplot(2,1,1); hold on
    leg={};    
% MAr
    if MAr==1
    plot(Age_Avg(:,4), Age_Avg(:,11),...
        'o-',...
        'MarkerEdgeColor',[0.6350 0.0780 0.1840],...
        'MarkerFaceColor', [0.6350 0.0780 0.1840],...
        'MarkerSize',4);
    leg(end+1)={['MAr']};
    end
% ZFT
    if ZFT==1
    plot(Age_Avg(:,4), Age_Avg(:,10),...
        'd-',...
        'color',[0.4940, 0.1840, 0.5560],...
        'MarkerFaceColor', [0.4940, 0.1840, 0.5560],...
        'MarkerSize',4);
    leg(end+1)={['ZFT']};
    end
% ZHe
    if ZHe==1
    plot(Age_Avg(:,4), Age_Avg(:,9),...
        'v-',...        
        'color',[0.8500, 0.3250, 0.0980],...
        'MarkerFaceColor', [0.8500, 0.3250, 0.0980],...
        'MarkerSize',4);
    leg(end+1)={['ZHe']};
    end
% AFT
    if AFT==1
    plot(Age_Avg(:,4), Age_Avg(:,8),...
        's-',...
        'color',[0, 0.4470, 0.7410],...
        'MarkerFaceColor', [0, 0.4470, 0.7410],...
        'MarkerSize',4);
    leg(end+1)={['AFT']};
    end
% AHe
    if AHe==1
    plot(Age_Avg(:,4), Age_Avg(:,7),...
        'p-',...
        'color',[0.4660, 0.6740, 0.1880],...
        'MarkerFaceColor', [0.4660, 0.6740, 0.1880],...
        'MarkerSize',4);
    leg(end+1)={['AHe']};
    end
% Formatting
    axis([XMin, XMax, YMin, YMax]);
    set(gca,'XMinorTick','on','YMinorTick','on');
    grid on;
    xlabel('Distance [km]');
    ylabel('Age [Ma]');
    title(strcat(Model, {' '}, ID,{' '},'Timestep',{' '}, num2str(i-1)), 'Interpreter', 'none');
    box on
% Plot Tchron Data
        if max(cell2mat(Answer.tchron(:,3)))==1
            % Plot the thermochronometer data
            for g=1:e % for each thermochronometer
                if Answer.tchron{g,3} % if user wants to plot this Tchron
                    if max(index(:,g))==1 % if the tchron data exists
                        errorbar(TchronData.xdist(index(:,g)),...
                            TchronData.Age((index(:,g))),...
                            TchronData.Error((index(:,g))),...
                            colordata{g},'MarkerFaceColor',markcol(g,:),...
                            'MarkerEdgeColor','k')
                        leg(end+1)={[legdata{g}]};
                    end
                end
            end
        end 
    % Legend         
    if      length(leg)==1
                legend(leg{1},'Location',legloc)
    elseif  length(leg)==2
                legend(leg{1},leg{2},'Location',legloc)
    elseif  length(leg)==3
                legend(leg{1},leg{2},leg{3},'Location',legloc)
    elseif  length(leg)==4
                legend(leg{1},leg{2},leg{3},leg{4},'Location',legloc)
    elseif  length(leg)==5
                legend(leg{1},leg{2},leg{3},leg{4},leg{5},'Location',legloc)
    elseif  length(leg)==6
                legend(leg{1},leg{2},leg{3},leg{4},leg{5},...
                    leg{6},'Location',legloc)
    elseif  length(leg)==7
                legend(leg{1},leg{2},leg{3},leg{4},leg{5},...
                    leg{6},leg{7},'Location',legloc)    
    elseif  length(leg)==8
                legend(leg{1},leg{2},leg{3},leg{4},leg{5},...
                    leg{6},leg{7},leg{8},'Location',legloc)          
    elseif  length(leg)==9
                legend(leg{1},leg{2},leg{3},leg{4},leg{5},...
                    leg{6},leg{7},leg{8},leg{9},'Location',legloc)
    elseif  length(leg)==10
                legend(leg{1},leg{2},leg{3},leg{4},leg{5},...
                    leg{6},leg{7},leg{8},leg{9},leg{10},'Location',legloc)
    end
    hold off
    subplot(2,1,2); 
    hold on
    %%% Plot temperature field %%%
    if tempplot==1
        [x, z] = meshgrid(XMin:XMax, ZMin:ZMax);
        Temperature = griddata(File_Temperatures.data(:,1),...
            File_Temperatures.data(:,3), File_Temperatures.data(:,5), x, z);
        surf(x, z, Temperature, 'EdgeColor', 'none');
        colormap jet
        c = colorbar('Location', 'southoutside');
        caxis([0 800])
        xlabel(c, 'Temperature [�C]','FontSize',14);
        velcol={'white'};
    else
        velcol={[0.25 0.5 0.85]};
    end
    
    %%% Plot structures %%%
    n = 1;
    while n < length(File_Structures.data(:,4))
        LineID = File_Structures.data(n,4);
        LineID_Indices = find(File_Structures.data(:,4) == LineID);
        plot3(File_Structures.data(LineID_Indices, 1), File_Structures.data(LineID_Indices, 3),...
            (zeros(length(File_Structures.data(LineID_Indices, 3)), 1) + 1e3),...
            'color', 'black',...
            'LineWidth', 1);
        hold on
        n = LineID_Indices(length(LineID_Indices)) + 1;
    end

    %%% Plot velocity field
    if velplot==1
        QuiverPlotDensity = 5;
        quiver3(File_Velocities.data(1:QuiverPlotDensity:end,2), ...
            File_Velocities.data(1:QuiverPlotDensity:end,4),...
            File_Velocities.data(1:QuiverPlotDensity:end,4)+1000,...
            File_Velocities.data(1:QuiverPlotDensity:end,5),...
            File_Velocities.data(1:QuiverPlotDensity:end,7), ...
            File_Velocities.data(1:QuiverPlotDensity:end,6),...
            0.5,...
            'LineWidth', 0.5,...
            'Color',velcol{1},...
            'MaxHeadSize', 0.010,...
            'AutoScaleFactor', 0.5);
    end
    if tempplot==0 && velplot==1
        set(gca,'Color',[0.90 0.90 0.90])
    end
    
    % Plot modifications
    axis([XMin, XMax, ZMin, ZMax]);
    set(gcf,'Position',[1 5 1280 700]) % fullscreen for Macbook pro 13" retina
    daspect([1 1 1])  % set subplot aspect ratio: no vertical exaggeration
    set(gca,'XMinorTick','on','YMinorTick','on');
    grid on;
    xlabel('Distance [km]');
    ylabel('Elevation [km]');
    box on
    hold off

    FileName = char(strcat(Model,'_',ID,...
        '_Timestep_', num2str(i-1),'WData',extension));
    saveas(gcf,fullfile(Folder_Save, FileName)) %saves figure
    close % closes figure
    
    
    %%% Goodness Of Fit
    GOF_version=1.0;
    samplelocll=zeros(d,g);
    ModelVal=zeros(d,g);
    Matches=zeros(d,g);
    
    ModelError=1; % in [Ma]
    for mm=1:2
        XsecError=mm; % in [km]
    for g=3:4 % AFT then ZHe
        if g==5 %AHe
            ii=7;
        elseif  g==4 %AFT
            ii=8; % which index value tchronometer I want
        elseif g==3 % ZHe
            ii=9;
        elseif g==2 %ZFT
            ii=10
        elseif g==1 %MAr
            ii=11
        end
        
        % for each sample in table
        for jj=1:d
            % if that sample exists in this thermochronometer
            if index(jj,g)==1
                Matches(jj,g)=3;
                %%%%%% Variables
                %%%% Measured Tchron Value
                %%% TchronData.Age(jj)
                
                %%%% Measured Tchron Error
                %%% TchronData.Error(jj)
                
                %%%% Sample Location (Adjusted to Pecube Grid)
                %%% TchronData.xdist(jj)
                sampleloc(jj,g)=round(TchronData.xdist(jj)*2)/2;
                
                % Loop over X-error in cross-section
                % Without looping over 1km wide window, N10const had 28.57%
                % fit. With 1km error, 42.86% fit.
                ll=[sampleloc(jj,g)-(XsecError): 0.5: sampleloc(jj,g)+(XsecError)];
                
                for lll= 1:length(ll)
                    samplelocll(jj,g)=ll(lll);
                    
                    % Find Sample Location in MODELLED
                    IndexSampleLocMod=find(Age_Avg(:,4)==samplelocll(jj,g));
                    % Get Modelled Value at that location
                    ModelVal(jj,g)=Age_Avg(IndexSampleLocMod,ii);
                    
                    % Compare Modelled to Measured
                    %%%% To accurately compare, we have 4 cases to test:
                    % Case 1. Model + Model Error within measured Age
                    % Case 2. Model - Model Error within measured Age
                    % Case 3. Measured + Measured error within modelled Age
                    % Case 4. Measured - Measured error within modelled Age
                    
                    % Case 1. Model Upper Error within measured Age
                    if (ModelVal(jj,g) + ModelError)< (TchronData.Age(jj) + TchronData.Error(jj)) ...
                            && ...
                            (ModelVal(jj,g) + ModelError)> (TchronData.Age(jj) - TchronData.Error(jj))
                        Matches(jj,g)=1; % it does match
                        
                        % Case 2. Model Lower Error within measured Age
                    elseif (ModelVal(jj,g) - ModelError)< (TchronData.Age(jj) + TchronData.Error(jj)) ...
                            && ...
                            (ModelVal(jj,g) - ModelError)> (TchronData.Age(jj) - TchronData.Error(jj))
                        Matches(jj,g)=1; % it does match
                        
                        % Case 3. Measured Upper error within modelled Age
                    elseif (TchronData.Age(jj) + TchronData.Error(jj)) < (ModelVal(jj,g) + ModelError) ...
                            && ...
                            (TchronData.Age(jj) + TchronData.Error(jj)) >(ModelVal(jj,g) - ModelError)
                        Matches(jj,g)=1; % it does match
                        
                        % Case 4. Measured Lower error within modelled Age
                    elseif (TchronData.Age(jj) - TchronData.Error(jj)) < (ModelVal(jj,g) + ModelError) ...
                            && ...
                            (TchronData.Age(jj) - TchronData.Error(jj)) >(ModelVal(jj,g) - ModelError)
                        Matches(jj,g)=1; % it does match
                        % No else statement needed bc line above defines
                        % Matches(jj,g)=3 %Not Matching and this loop only
                        % overwrites that if the sample is within the error
                    end
                end
                %
            end
        end
    end
    
    % Percent Fit
    YesMatches=find(Matches==1);
    NoMatches=find(Matches==3);
    PercFit=length(YesMatches)/(length(NoMatches)+length(YesMatches))
    
    YesMatchesTchron=[];
    nsamples=[];
    for kk=1:e
        YesM=length((find(Matches(:,kk)==1)));
        NoM=length((find(Matches(:,kk)==3)));
        nsamples(kk)=YesM+NoM;
        YesMatchesTchron(kk)=YesM/(NoM+YesM);
    end
    
    %%% Print goodness of fit to bulk file
%     % Model_Name	Model_ID	Model_Error	Sample_Loc_Error	%Fit
%     GOFLoc='/Volumes/Files/VictoriaFiles/Pecube/GoodnessOfFit.txt';
%     fidGOF=fopen(GOFLoc,'at');
%     fprintf(fidGOF,'\n%10s\t%18s\t%0.1f\t%16.1f\t%7.4f',Model,ID,ModelError,XsecError,PercFit);
%     fclose(fidGOF);
    
   %%%% Print goodness of fit + matches + sample info to individual file
    GOFLoc_Ind=strcat(Folder_Save,'/GoodnessOfFit_',Model,'_',ID,'_',num2str(mm),'.txt');
    fidGOF_ind=fopen(GOFLoc_Ind,'w');
    fprintf(fidGOF_ind,'Goodness of Fit Computation; Version %1.1f\n',GOF_version);
    fprintf(fidGOF_ind,'Model_Name\tModel_ID\tModel_Error\tSample_Loc_Error\tPercFit\n');
    fprintf(fidGOF_ind,'%s\t\t%s\t\t%0.1f\t\t%f\t\t%f\n',Model,ID,ModelError,XsecError,PercFit);
    fprintf(fidGOF_ind,'\n')
    fprintf(fidGOF_ind,'Summary \tMAr\tZFT\tZHe\tAFT\tZHe\n');
    fprintf(fidGOF_ind,'PercFit \t%0.4f\t%0.4f\t%0.4f\t%0.4f\t%0.4f\n',YesMatchesTchron');
    fprintf(fidGOF_ind,'No. samples\t%0.0f\t%0.0f\t%0.0f\t%0.0f\t%0.0f\n',nsamples');
    fprintf(fidGOF_ind,'\nMatches [n=%.0f; e.g. No=3; Yes=1]\n',(length(NoMatches)+length(YesMatches)));
    fprintf(fidGOF_ind,'Sample \t%s \t%s \t%s \t%s \t%s \n',legdata{:});
    for tt=1:length(Matches)
        fprintf(fidGOF_ind,'%s \t\t%.0f \t\t%.0f \t\t%.0f \t\t%.0f \t\t%.0f \n',TchronData.sample_num{tt}, Matches(tt,:));
        %fprintf(fidGOF_ind,'%s \t\t%.0f \t\t%.0f \t\t%.0f \t\t%.0f \t\t%.0f \n',num2str(tt), Matches(tt,:));
    end
    fclose(fidGOF_ind);
    
    end
    
end
% XMax=1035;XMin=625; MAr=false;AHe=false;ZFT=false;legloc='NorthEast';a=39;b=39;
%% Uncomment this if you want it to alert you after plotting finishes
% only when running full script
 load handel; sound(y,1.75*Fs) % full list of sounds help audiovideo
%% Comments on alterations
% changes:

% v1.5
% added option to not plot temperatures or velocities

% v1.4
% fixed errors in thermochronometer data plotting
% fixed problems with vertical exaggeration

% v1.3
% added ability to plot thermochronometer data
% made cross sections have no vertical exaggeration
% 261: removed Ages from y other than xsec.

% v1.2
% Change vel scheme & Ao to text ID
% changed so that default is automatically the last run.
% allow user to select which type of moving average : use smoothdata
% add which time steps to plot

%%%% Changes still to do

%%%General Code
% redo how we import files: would prefer to select files and just plot them
% rather than search for name
    % # of files sees hidden/deleted temp files...
    

% resizing the subplot (2,1,2) window (line 372: set(gcf, 'Position'...)
%      resizes all of the fonts except for the colorbar (axis and label)

% Make it so that daspect doesn't change width, only height.
% add K Ar/Ar and B Ar/Ar order==(ZHe, K Ar/Ar, B Ar/Ar, ZFT)

% line 279: remove zeros also in the ages?
% add some way of determining goodness of fit between data & model?
    %  interpolate model to data locations (interp2)
    %  perhaps in the interpolation, average nearest ~3 model points

%%%GUI
% Make table with tchron smoothing options column with method wider.
% add help popups?
% tick boxes: what if I only want to plot ages and lines?
% ability to toggle off temperature or velocities. 
    % change quiver color if temperatures are ticked off.

% edit folder by hand without popups (use popup button)
% Xrange, ZRange, sliders... ?/]
% at group meeting: mention the box folder.

% New > App > Guide > Blank GUI...