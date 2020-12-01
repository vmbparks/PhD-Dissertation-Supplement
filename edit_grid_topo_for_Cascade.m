%% Shift Grid and Topo lines for CASCADE models.
clear;clc;close;

%% Setup
amount=(1022+5); %(amount of shift; currently moves from right to left)
amount=amount-200;

rot=1; % to swap east and west = -1; to just shift =1
gridloc='/Volumes/Files/VictoriaFiles/Cascade/McQ02N9/gridtopolines_v40-on/grid/McQ02N9_grid_Step';
topoloc='/Volumes/Files/VictoriaFiles/Cascade/McQ02N9/gridtopolines_v40-on/topo/McQ02N9_topo_Step';
moveloc='/Volumes/Files/VictoriaFiles/Cascade/McQ02N9/gridtopolines_v40-on/lines/McQ02N9_lines_Step';
a=26; % which step to start on
b=38; % which step to end on

%% Edit Grid Files
tic
for i=a:b
    filename = strcat(gridloc,num2str(i),'.dat');
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
    grid2(:,1)=(grid2(:,1)-(amount))*rot;
    %export editted grid
    dlmwrite(filename,grid2,'delimiter',delimiter);
    
    % plot to verify correct offset
    figure
    subplot(2,1,1)
    plot(grid(:,1),grid(:,3),'o');
    axis([-400 1200 -100 100])
    subplot(2,1,2)
    plot(grid2(:,1),grid2(:,3),'o');
    axis([-400 1200 -100 100])
    title(num2str(i))
    hold off
end
toc

%% Edit topo files
tic
for i=a:b
    filename = strcat(topoloc,num2str(i),'.dat');
    delimiter = ' ';
    formatSpec = '%8f%8f%7f%f%[^\n\r]';
    % Open the text file.
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN,  'ReturnOnError', false);
    
    % Close the text file.
    fclose(fileID);
    % Create output variable
    McQ02N2Cascadetopo = [dataArray{1:end-1}];
    McQ02N2Cascadetopo2 = McQ02N2Cascadetopo;
    % Clear temporary variables
    clearvars  formatSpec fileID dataArray ans;
    
    % edit topo
    McQ02N2Cascadetopo2(:,1)=(McQ02N2Cascadetopo2(:,1)-(amount))*rot;
    % export editted topo lines
    dlmwrite(filename,McQ02N2Cascadetopo2,'delimiter',delimiter);

    %plot to verify
    figure
    subplot(2,1,1)
    plot(McQ02N2Cascadetopo(:,1),McQ02N2Cascadetopo(:,3));
    axis([-600 1200 -50 50])
    subplot(2,1,2)
    plot(McQ02N2Cascadetopo2(:,1),McQ02N2Cascadetopo2(:,3));
    axis([-600 1200 -50 50])
    hold off
end
toc
%% Edit MOVE lines.
tic
for i=a:b
    filename =  strcat(moveloc,num2str(i),'.dat');
    % Import move lines
    File_Structures = importdata([filename], ' ', 1);
    
    % Create output variable
    File_Structures2=File_Structures;
    
    % edit Movelines
    File_Structures2.data(:,1)=(File_Structures2.data(:,1)-amount)*rot;
    
    % export editted move lines
    dlmwrite(filename,File_Structures2.data,'delimiter',delimiter)
    
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
%% Impose initial gradient for topo