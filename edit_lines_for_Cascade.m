%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Volumes/Files/McQ02N2/Cascade/MoveLines/McQ02N2_lines_Step26.dat
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2018/03/05 10:53:12

%% Initialize variables.
filename = '/Volumes/Files/McQ02N2/Cascade/MoveLines/McQ02N2_lines_Step26.dat';
delimiter = ' ';
startRow = 2;

%Format for each line of text:
%   column1: double (%f)
%	column2: double (%f)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%f%f%f%f%[^\n\r]';

%Open the text file.
fileID = fopen(filename,'r');

% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true, 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
%%
%Close the text file.
fclose(fileID);
%
%Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%Create output variable
McQ02N2lines = [dataArray{1:end-1}];
%Clear temporary variables
clearvars startRow formatSpec fileID dataArray ans;

% Create output variable
McQ02N2lines = [dataArray{1:end-1}];
%%
McQ02N2lines2 = McQ02N2lines;

 %Clear temporary variables
%clearvars  formatSpec fileID ;

% ^^ Import Data
% Edit data

McQ02N2lines2(:,1)=McQ02N2lines2(:,1)-(1025-167);
%% test

subplot(2,1,1)
n = 1;
    while n < length(McQ02N2lines(:,4))
        LineID = McQ02N2lines(n,4);
        LineID_Indices = find(McQ02N2lines(:,4) == LineID);
        plot3(McQ02N2lines(LineID_Indices, 1), McQ02N2lines(LineID_Indices, 3),...
            (zeros(length(McQ02N2lines(LineID_Indices, 3)), 1) + 1e3),...
            'color', 'black',...
            'LineWidth', 1);
        hold on
        n = LineID_Indices(length(LineID_Indices)) + 1;
    end
    hold off
%plot(McQ02N2lines(:,1),McQ02N2lines(:,3),'o');
axis([-400 1200 -100 100])
%%

subplot(2,1,2)
n = 1;
    while n < length(McQ02N2lines2(:,4))
        LineID = McQ02N2lines2(n,4);
        LineID_Indices = find(McQ02N2lines2(:,4) == LineID);
        plot3(McQ02N2lines2(LineID_Indices, 1), McQ02N2lines2(LineID_Indices, 3),...
            (zeros(length(McQ02N2lines2(LineID_Indices, 3)), 1) + 1e3),...
            'color', 'black',...
            'LineWidth', 1);
        hold on
        n = LineID_Indices(length(LineID_Indices)) + 1;
    end
    hold off
%plot(McQ02N2lines2(:,1),McQ02N2lines2(:,3),'o');
axis([-400 1200 -100 100])
daspect([1 1 1])
hold off

%% Clean up lines
l=find(McQ02N2lines2(:,1)>134.3&McQ02N2lines2(:,1)<134.5)
m=find(McQ02N2lines2(l,3)>4.5)
lineee=McQ02N2lines2(m,4)
zzzzz=find(McQ02N2lines2(:,4)==lineee);
McQ02N2lines2(zzzzz,:)=[];


%% Export and rewrite data
dlmwrite(filename,McQ02N2lines2,'delimiter',delimiter);
