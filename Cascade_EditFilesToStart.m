%% Cascade_EditFilesToStart
% written by Victoria Buford; vmb21@pitt.edu

% Edit MOVE grids, topos, and lines in order to input into CASCADE
% CASCADE analyzes from x=0 to x=ModelSize [kms]
% Orographic Precipitation (as of version 1.8) only works correctly with
% the deformation front at x=ModelSize (eg, on the right,)
% August 15, 2018

% Note: Trim Grid is currently off/not available, until CASCADE won't crash
% if UIDs "appear" in subsequent grids.

%clear; close all; clc;
%%
set(0,'DefaultUicontrolFontSize',14)
%%
% Locate the folder w/ grid %
% Locate the folder w/ topo %
% Locate the folder w/ movelines %
% Locate the folder to save new lines %
% How big do you want the model space
% Where is the deformation front in your model?
% Is the deformation front on the Right or Left?
% Do you want to delete beyond the Cascade model space?
clear;
disp('run')
% Define Default Answers

if exist('Answer','var')==1
    LinesDir    = Answer.MoveLinesFolder;
    GridDir     = Answer.GridFolder;
    TopoDir     = Answer.TopoFolder;
    SaveDir     = Answer.SaveFolder;
    ModelSize   = Answer.ModelSize;
    leglocdef   = Answer.LegLoc;
    leglocintdef= Answer.LegLocInt;
    TrimDef     = Answer.Trim;
else
    LinesDir    = pwd;
    GridDir     = pwd;
    TopoDir     = pwd;
    SaveDir     = pwd;
    ModelSize   = 200;
    leglocdef   = 'right';
    leglocintdef= 500;
    TrimDef     = 'no';
end 
    
    
Title = 'Cascade Shift and Rotate';

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
%Formats(1,1).size = [-1 0];
%Formats(1,1).span = [1 2]; % item is 1 field x 4 fields

Prompt(end+1,:) = {'Model Size','ModelSize',[]};
Formats(1,2).type='edit';
Formats(1,2).format='integer';
Formats(1,2).size=[35 25];
DefAns(1).ModelSize=ModelSize;

Prompt(end+1,:) = {'Trim Grid to Model','Trim',[]};
Formats(1,3).type='edit';
Formats(1,3).format='text';
Formats(1,3).size=[35 25];
DefAns.Trim=TrimDef;

Prompt(end+1,:)={'Deformation front on', 'LegLoc',[],};
Formats(2,1).type='edit';
Formats(2,1).format='text';
Formats(2,1).size=[35 25];
DefAns.LegLoc=leglocdef;

Prompt(end+1,:)={'Location', 'LegLocInt',[],};
Formats(2,3).type='edit';
Formats(2,3).format='integer';
Formats(2,3).size=[50 25];
DefAns.LegLocInt=leglocintdef;

Prompt(end+1,:) = {'Grid Folder','GridFolder',[]};
Formats(3,1).type = 'edit';
Formats(3,1).format = 'dir';
Formats(3,1).size = [-1 0];
Formats(3,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.GridFolder = GridDir;

Prompt(end+1,:) = {'Topo Folder','TopoFolder',[]};
Formats(4,1).type = 'edit';
Formats(4,1).format = 'dir';
Formats(4,1).size = [-1 0];
Formats(4,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.TopoFolder = TopoDir;

Prompt(end+1,:) = {'Move Lines Folder','MoveLinesFolder',[]};
Formats(5,1).type = 'edit';
Formats(5,1).format = 'dir';
Formats(5,1).size = [-1 0];
Formats(5,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.MoveLinesFolder = LinesDir;

Prompt(end+1,:) = {'Save Folder','SaveFolder',[]};
Formats(6,1).type = 'edit';
Formats(6,1).format = 'dir';
Formats(6,1).size = [-1 0];
Formats(6,1).span = [1 3];  % item is 1 field x 3 fields
DefAns.SaveFolder = SaveDir;

[Answer,Cancelled] = inputsdlg(Prompt,Title,Formats,DefAns,Options);

%% Interpret Model Dialog
% Answer.ModelSize          integer
% Answer.Trim               str
% Answer.LegLoc             str     'right' or 'left'
% Answer.LegLocInt          integer
% Answer.MoveLinesFolder    str     dir
% Answer.GridFolder         str     dir
% Answer.TopoFolder         str     dir
% Answer.SaveFolder         str     dir

%CASCADE_Files = dir([OutputFolder '/*.dat']);

buffer=5;
% Determine rotation & Model Shift
if strcmp(Answer.LegLoc,'right')
    rot     = 1; % Don't rotate
    amount  = Answer.LegLocInt + buffer - Answer.ModelSize;
else 
    rot     = -1; % Rotate it
    amount  = -Answer.LegLocInt + buffer - Answer.ModelSize;
end

% Define Directories
Line_Files  = dir([Answer.MoveLinesFolder '/*.dat']);
Topo_Files  = dir([Answer.TopoFolder '/*.dat']);
Grid_Files  = dir([Answer.GridFolder '/*.dat']);
Save_Folder = Answer.SaveFolder;

% Trim Model?
if strcmp(Answer.Trim,'yes')
    Trim=1;
else
    Trim=0;
end

% Loop Sizes
a=1; 
b=length(Grid_Files);
bb=length(Topo_Files);
bbb=length(Line_Files);


%% Edit Grid Files
tic
for i=a:b
    filename = strcat(Grid_Files(i).folder,'/',Grid_Files(i).name);
    savename = strcat(Save_Folder,'/',Grid_Files(i).name);
    delimiter = ' ';
    formatSpec = '%f%f%f%f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
    fclose(fileID);
    
    % Create output variable
    grid = [dataArray{1:end-1}];
    grid2 = grid;
    
    %Clear temporary variables
    clearvars  formatSpec fileID dataArray;

    % Edit data
    if rot == 1 % don't rotate; shift
        grid2(:,1)=(grid2(:,1)-(amount))*rot;
    elseif rot == -1 % rotate, then shift
        grid2(:,1)=(grid2(:,1)*rot)-(amount);
    end
    
%     % Trim Grid
%     idx=find(grid2(:,1)<0 || grid2(:,1)>Answer.ModelSize);
%     idxx=find(grid2(:,1)>Answer.ModelSize);
%     idx=[idx;idxx];
%     grid2(idx,:)=[]; % Remove rows outside of the model space
    
    % plot to verify correct offset
    figure
    subplot(2,1,1)
    plot(grid(:,1),grid(:,3),'o');
    axis([-400 1200 -100 100])
    subplot(2,1,2)
    plot(grid2(:,1),grid2(:,3),'o');
    axis([-400 1200 -100 100])
    title(Grid_Files(i).name)
    hold off
    %pause
     %export editted grid
    dlmwrite(savename,grid2,'delimiter',delimiter);
    
end
toc
%% Edit topo files
clear grid grid2
close all
%%%%%%% To edit/ add a  bit of topography..
% plot(topo.x,topo.z) % original
% zplu=linspace(0,1,length(topo.x));
% hold on
% topot=topo;
% topot.z=topot.z+zplu'; % edited...
% plot(topot.x,topot.z)
% axis([0 200 0 6])
% hold off
% %
% savename='/Users/victoriabuford/Box Sync/Research/MoveStuff/Adam/ES/Cascade/NBFTp5_31_t-Cascade.dat';
% 
% dlmwrite(savename,table2array(topot),'delimiter',' ');
%%%%%%% original edit topo files continues here
tic
for i=a:bb
    filename = strcat(Topo_Files(i).folder,'/',Topo_Files(i).name);
    savename = strcat(Save_Folder,'/',Topo_Files(i).name);    
    delimiter = ' ';
    formatSpec = '%8f%8f%7f%f%[^\n\r]';
    % Open the text file.
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
    
    % Close the text file.
    fclose(fileID);
    % Create output variable
    topo = [dataArray{1:end-1}];
    topo2 = topo;
    % Clear temporary variables
    clearvars  formatSpec fileID dataArray ans;
    
    % Edit data
    if rot == 1 % don't rotate; shift
        topo2(:,1)=(topo2(:,1)-(amount))*rot;
    elseif rot == -1 % rotate, then shift
        topo2(:,1)=(topo2(:,1)*rot)-(amount);
    end
    
    % export editted topo lines
    dlmwrite(savename,topo2,'delimiter',delimiter);

    %plot to verify
    figure
    subplot(2,1,1)
    plot(topo(:,1),topo(:,3));
    axis([-600 1200 -50 50])
    subplot(2,1,2)
    plot(topo2(:,1),topo2(:,3));
    axis([-600 1200 -50 50])
    title(Topo_Files(i).name)
    hold off
end
toc
%% Edit MOVE lines.
clear topo topo2
close all
%%%%% DON'T RUN THIS YET!! Doesn't rotate and whatnot correctly.
tic
for i=a:bbb
    filename = strcat(Line_Files(i).folder,'/',Line_Files(i).name);
    savename = strcat(Save_Folder,'/',Line_Files(i).name);    
    % Import move lines
    File_Structures = importdata([filename], ' ', 1);
    
    % Create output variable
    File_Structures2=File_Structures;
    
    % edit Movelines
    File_Structures2.data(:,1)=(File_Structures2.data(:,1)-amount)*rot;
    
    % Edit data
    if rot == 1 % don't rotate; shift
        File_Structures2.data(:,1)=(File_Structures2.data(:,1)-(amount))*rot;
    elseif rot == -1 % rotate, then shift
        File_Structures2.data(:,1)=(File_Structures2.data(:,1)*rot)-(amount);
    end
    
    % export editted move lines
    dlmwrite(savename,File_Structures2.data,'delimiter',delimiter)
    
    %%% Plot structures to verify%%%
    figure
    subplot(2,1,1)
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
        title(strcat('Original: Step ',num2str(i)))
        axis([-300 1025 -60 20])
        daspect([1 1 1])
        view(2)
    end
    
    subplot(2,1,2)
    n = 1;
    while n < length(File_Structures2.data(:,4))
        LineID = File_Structures2.data(n,4);
        LineID_Indices = find(File_Structures2.data(:,4) == LineID);
        plot3(File_Structures2.data(LineID_Indices, 1), File_Structures2.data(LineID_Indices, 3),...
            (zeros(length(File_Structures2.data(LineID_Indices, 3)), 1) + 1e3),...
            'color', 'black',...
            'LineWidth', 1);
        hold on
        n = LineID_Indices(length(LineID_Indices)) + 1;
        title(strcat('Editted: Step ',num2str(i)))
        axis([-300 1025 -60 20])
        view(2)
        daspect([1 1 1])
    end
    
    %    subplot(2,1,2)
    %     n = 1;
    %     while n < length(Movelines2(:,4))
    %         LineID = Movelines2(n,4);
    %         LineID_Indices = find(Movelines2(:,4) == LineID);
    %         plot3(Movelines2(LineID_Indices, 1), Movelines2(LineID_Indices, 3),...
    %             (zeros(length(Movelines2(LineID_Indices, 3)), 1) + 1e3),...
    %             'color', 'black',...
    %             'LineWidth', 1);
    %         hold on
    %         n = LineID_Indices(length(LineID_Indices)) + 1;
    %         title(strcat('editted : Step ',num2str(i)))
    %     end
    
    
end
toc

