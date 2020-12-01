%% Instructions and notes
% Script that extracts user-defined topography swaths and normalised river
% steepness values.
% Simply follow the prompts in the command window.
% Note: Imported DEMs need to be in WGS 1984 (for now...)
%
% Modified by Victoria M. Buford Parks (Dec 2020) from 
% Paul R. Eizenhöfer, Eberhard-Karls-Universität Tübingen
% paul-reinhold.eizenhoefer@uni-tuebingen.de
% January 2020

% October 2020: 
%   - added local relief output
%   - added precipitation output (requires TRMM precipitation file)
%   - cleaned up plots

%% Initiate script
tic
clear

% Analyses options
UserInput = 'DEM Type? (GeoTiff = 0, CASCADE = 1) ';
DEMInput = input(UserInput);

if DEMInput == 1
    disp('Locate CASCADE output file topo_tec_****.dat.');
    [File_CASCADE, Path_CASCADE] = uigetfile('topo_tec_****.dat', 'Select CASCADE output file.');
    
    CASCADE_File = readtable([Path_CASCADE File_CASCADE], 'HeaderLines', 4);
    
    OutputVariables = {'x_km','y_km','z_km','node','precipitation_mpy',...
        'fluvial_erosion_rate_mpy','diffusion_erosion_rate_mpy',...
        'landslide_erosion_rate_mpy','total_erosion_rate_mpy',...
        'catchment_color','catchment_number',...
        'glacial_erosion_rate_mpy','ice_thickness_m','mass_balance_1py',...
        'total_topography_m','sliding_velocity_mpy','gerode_term_mpy',...
        'rock_contact_km','isostatic_deflection_mpy','slope_mpkm',...
        'totalflexiso_m','constriction','cumulative_erosion_m',...
        'surface_area_km2'}';
    CASCADE_File.Properties.VariableNames = OutputVariables';
    
    % Clean up output - Remove NaNs
    Nanss = ismissing(CASCADE_File,{NaN}); % creates table size of A, with 1 or 0 for NaN
    NanNodes = Nanss(:,4); clear Nanss % pull only the nodes column
    CASCADE_File = CASCADE_File(~NanNodes,:); clear NanNodes % remove any row that the node is NaN
    
    % Interpolation to regular grid
    
    UserInput = ('Choose DEM resolution in [m]: ');
    DEMResolution = input(UserInput);
    
    UserInput = ('Choose DEM width in [m] (should correspond to CASCADE model width): ');
    DEMWidth = input(UserInput);
    
    UserInput = ('Choose DEM length in [m] (should correspond to CASCADE model length): ');
    DEMLength = input(UserInput);
    
    [x_DEM,y_DEM] = meshgrid(0:DEMResolution:DEMWidth, 0:DEMResolution:DEMLength);
    z_DEM = griddata(CASCADE_File.x_km * 1e3, CASCADE_File.y_km * 1e3, CASCADE_File.z_km * 1e3, x_DEM, y_DEM, 'natural');
    Erosion_DEM = ...
        griddata(CASCADE_File.x_km * 1e3, CASCADE_File.y_km * 1e3, CASCADE_File.total_erosion_rate_mpy, x_DEM, y_DEM, 'natural');
    Precipitation_DEM = ...
        griddata(CASCADE_File.x_km * 1e3, CASCADE_File.y_km * 1e3, CASCADE_File.precipitation_mpy, x_DEM, y_DEM, 'natural');
    
    DEM = GRIDobj(x_DEM, y_DEM, z_DEM);
    DEM_Erosion = GRIDobj(x_DEM, y_DEM, Erosion_DEM);
    DEM_Precipitation = GRIDobj(x_DEM, y_DEM, Precipitation_DEM);
    
    % Define Swath Geometry
    UserInput = 'Enter [longitude latitude] of profile start in [m]: ';
    ProfileStartInput = input(UserInput);
    
    UserInput = 'Enter [longitude latitude] of profile end in [m]: ';
    ProfileEndInput = input(UserInput);
    
else
    disp('Select DEM for geomorphic analysis...')
    [File_DEM, Path_DEM] = uigetfile('****.tif', 'Locate DEM');
    
    UserInput = 'DEM resolution in [m]? ';
    DEM_Resolution = input(UserInput);
    
    disp('DEM is loading...')
    
    DEM = GRIDobj([Path_DEM File_DEM]);
    DEM = reproject2utm(DEM, DEM_Resolution);
    
    disp('Done.')
    
    disp('Select TRMM file for precipitation swath...')
    [File_Precipitation, Path_Precipitation] = uigetfile('****.tif', 'Locate TRMM file');
    
    UserInput = 'Precipitation data resolution in [m]? ';
    Precipitation_Resolution = input(UserInput);
    
    UserInput = 'Enter UTM Zone [integer]: ';
    UTM_Zone = input(UserInput);

    UserInput = 'Northern or southern hemisphere (N = 0; S = 1)? ';
    HemisphereCheck = input(UserInput);
    
    disp('Precipitation data are loading...')
    
    PrecipitationData = GRIDobj([Path_Precipitation File_Precipitation]);
    
    disp('Done.')
    
    % Define Swath Geometry
    UserInput = 'Enter [longitude latitude] of profile start in [decimal degrees]: ';
    ProfileStartInput = input(UserInput);
    
    UserInput = 'Enter [longitude latitude] of profile end in [decimal degrees]: ';
    ProfileEndInput = input(UserInput);
end

UserInput = 'Swath profile resolution? (in [m]) ';
ProfileResolution = input(UserInput);

UserInput = 'Width of analysis swath? (in [m]) ';
ProfileWidth = input(UserInput);

%% TOPOGRAPHY SWATH

disp('Extracting topography and precipitation swaths...')

if DEMInput == 1
    ProfileStart = ProfileStartInput;
    ProfileEnd = ProfileEndInput;
else
    [ProfileStart(1), ProfileStart(2)] = ll2utm(ProfileStartInput(2), ProfileStartInput(1));
    [ProfileEnd(1), ProfileEnd(2)] = ll2utm(ProfileEndInput(2), ProfileEndInput(1));
end

ProfileLength = sqrt((ProfileStart(1) - ProfileEnd(1))^2 + (ProfileStart(2) - ProfileEnd(2))^2);
ProfileAngle = asin(abs(ProfileStart(2) - ProfileEnd(2)) / ProfileLength); % Note: angle in radians

% Extract Swath Topography

% Define temporary profile end points
% Swath to the left
for i = 1:round(ProfileWidth / ProfileResolution / 2)
    ProfileStart(i+1,1) = ProfileStart(1,1) - i * cos(pi/2 - ProfileAngle) * ProfileResolution;
    ProfileStart(i+1,2) = ProfileStart(1,2) + i * sin(pi/2 - ProfileAngle) * ProfileResolution;
    ProfileEnd(i+1,1) = ProfileEnd(1,1) - i * cos(pi/2 - ProfileAngle) * ProfileResolution;
    ProfileEnd(i+1,2) = ProfileEnd(1,2) + i * sin(pi/2 - ProfileAngle) * ProfileResolution;
end

% Swath to the right
for j = 1:round(ProfileWidth / ProfileResolution / 2)
    ProfileStart(i+j+1,1) = ProfileStart(1,1) + j * sin(ProfileAngle) * ProfileResolution;
    ProfileStart(i+j+1,2) = ProfileStart(1,2) - j * cos(ProfileAngle) * ProfileResolution;
    ProfileEnd(i+j+1,1) = ProfileEnd(1,1) + j * sin(ProfileAngle) * ProfileResolution;
    ProfileEnd(i+j+1,2) = ProfileEnd(1,2) - j * cos(ProfileAngle) * ProfileResolution;
end

% Extract profiles
% Extract topographic profiles
if DEMInput == 0
    if (UTM_Zone > 30 && HemisphereCheck == 0)
        PrecipitationData = crop(PrecipitationData, [0 0 180.02 180.02], [0 36.02 0 36.02]);
    elseif (UTM_Zone > 30 && HemisphereCheck == 1)
        PrecipitationData = crop(PrecipitationData, [0 0 180.02 180.02], [0 -36.02 0 -36.02]);
    elseif (UTM_Zone < 30 && HemisphereCheck == 0)
        PrecipitationData = crop(PrecipitationData, [0 0 -180.02 -180.02], [0 36.02 0 36.02]);
    elseif (UTM_Zone < 30 && HemisphereCheck == 1)
        PrecipitationData = crop(PrecipitationData, [0 0 -180.02 -180.02], [0 -36.02 0 -36.02]);
    end
    PrecipitationData = reproject2utm(PrecipitationData, Precipitation_Resolution);
end

for k = 1:(j+i+1)
    [temp_d,temp_z,~,~] = ...
        demprofile(DEM,round(ProfileLength/ProfileResolution),...
        [ProfileStart(k,1) ProfileEnd(k,1)], [ProfileStart(k,2) ProfileEnd(k,2)]);
    if DEMInput == 1
        [~,temp_Erosion,~,~] = ...
            demprofile(DEM_Erosion, round(ProfileLength/ProfileResolution),...
            [ProfileStart(k,1) ProfileEnd(k,1)], [ProfileStart(k,2) ProfileEnd(k,2)]);
        [~,temp_Precipitation,~,~] = ...
            demprofile(DEM_Precipitation, round(ProfileLength/ProfileResolution),...
            [ProfileStart(k,1) ProfileEnd(k,1)], [ProfileStart(k,2) ProfileEnd(k,2)]);
        ProfileErosion(k,:) = temp_Erosion;
        ProfilePrecipitation(k,:) = temp_Precipitation;
    else
        [~,temp_Precipitation,~,~] = ...
            demprofile(PrecipitationData, round(ProfileLength/ProfileResolution),...
            [ProfileStart(k,1) ProfileEnd(k,1)], [ProfileStart(k,2) ProfileEnd(k,2)]);
        ProfilePrecipitation(k,:) = temp_Precipitation;
    end
    ProfileDistance = temp_d;
    ProfileElevation(k,:) = temp_z;
   
end
% Calculate minimum, maximum and mean topographies and precipitation
for k = 1:length(ProfileDistance)
    MaximumTopography(k) = max(ProfileElevation(:,k));
    MeanTopography(k) = nanmean(ProfileElevation(:,k));
    MinimumTopography(k) = min(ProfileElevation(:,k));
    
    if DEMInput == 1
        MeanErosion(k) = nanmean(ProfileErosion(:,k));
        MeanPrecipitation(k) = nanmean(ProfilePrecipitation(:,k));
    else
        MaximumPrecipitation(k) = max(ProfilePrecipitation(:,k)/1e3);
        MeanPrecipitation(k) = nanmean(ProfilePrecipitation(:,k)/1e3);
        MinimumPrecipitation(k) = min(ProfilePrecipitation(:,k)/1e3);
    end
end

Topography = [ProfileDistance MeanTopography' MaximumTopography' MinimumTopography'];

if DEMInput == 1
    Erosion = [ProfileDistance MeanErosion'];
    Precipitation = [ProfileDistance MeanPrecipitation'];
else
    Precipitation = [ProfileDistance MeanPrecipitation' MaximumPrecipitation' MinimumPrecipitation'];
end

disp('Done.')

%% KSN SWATH
disp('Performing local gradient and Ksn analysis...')

% DEM flow routing and sink filling
FD = FLOWobj(DEM,'preprocess','carve');
DEM = imposemin(FD,DEM);
DEM = fillsinks(DEM);
FDm = FLOWobj(DEM,'dinf');
% FDm = Flowobj(DEM, 'multi');
A = flowacc(FDm);

% Calculate stream network
S = STREAMobj(FDm, 'unit', 'mapunits', 'minarea', 1e7); % Define minimum drainage area here

% Calculate Ksn values and drainage area for all streams
G = gradient8(DEM);
DRAINAGE_AREA = A.*(A.cellsize^2);
KSN = G./DRAINAGE_AREA.^-0.45; % Define stream concavity index here

% Extract streams
[x,y,ksn,localgradient,drainage_area] = STREAMobj2XY(S,KSN,G,DRAINAGE_AREA);
Streams = [x, y, ksn, localgradient, drainage_area];

% Filter data (NaN's and ksn = 0 will be removed)
Streams(any(isnan(Streams),2),:) = [];
Streams(Streams(:,3) == 0,:) = [];

% Rotation
RotationAngle = -rad2deg(ProfileAngle);
Rotation = rotz(RotationAngle); % insert rotation angle here (counterclockwise)
RotationCentre = [ProfileStart(1,1); ProfileStart(1,2); 0]; % insert rotation centre here, ideally beginning of section line

% All streams
Streams_rotated = (Rotation * ([Streams(:,1)'; Streams(:,2)'; Streams(:,3)'] - RotationCentre)) + RotationCentre;
Streams_rotated = [Streams_rotated; Streams(:,4:5)'];
disp('Done.')

% Extract ksn values
% Area of interest
x_range = [RotationCentre(1); RotationCentre(1) + ProfileLength]; % define length of section starting at the rotation centre
y_range = [RotationCentre(2) - ProfileWidth;  RotationCentre(2) + ProfileWidth]; % define area across section

% All streams
Streams_filtered = Streams_rotated(:,...
    Streams_rotated(1,:) > x_range(1) & Streams_rotated(1,:) < x_range(2) &...
    Streams_rotated(2,:) > y_range(1) & Streams_rotated(2,:) < y_range(2));
KsnValues = bin(Streams_filtered(1,:), Streams_filtered(3,:), round(ProfileLength/ProfileResolution));
LocalGradient = bin(Streams_filtered(1,:), Streams_filtered(4,:), round(ProfileLength/ProfileResolution));% Define number of bins
% Streams_binned_SmoothMedian = smoothdata(Streams_binned(:,2), 'movmedian'); % Define width of moving average window

%% Plots
% Topography and Ksn
fig = figure;
set(fig,'defaultAxesColorOrder',[[0 0 0]; [0 0 0]]);

if DEMInput == 1
    subplot(3,4,1)
else
    subplot(3,3,1)
end
imageschs(DEM,[],'colormap',[1 1 1],'colorbar',false,'ticklabel','nice');

hold on

plot([ProfileStart(1,1); ProfileEnd(1,1)], [ProfileStart(1,2); ProfileEnd(1,2)], 'r:', 'LineWidth', 2)
plot([ProfileStart(i+1,1); ProfileEnd(i+1,1)], [ProfileStart(i+1,2); ProfileEnd(i+1,2)], 'r', 'LineWidth', 2)
plot([ProfileStart(i+j+1,1); ProfileEnd(i+j+1,1)], [ProfileStart(i+j+1,2); ProfileEnd(i+j+1,2)], 'r', 'LineWidth', 2)
plot(ProfileStart(:,1), ProfileStart(:,2), 'r', 'LineWidth', 2)
plot(ProfileEnd(:,1), ProfileEnd(:,2), 'r', 'LineWidth', 2)

xlabel('Longitude [m]')
ylabel('Latitude [m]')
title({'River Steepness';...
    ['[' num2str(round(ProfileStartInput(1), 2)) ' ' num2str(round(ProfileStartInput(2), 2))...
    '] to [' num2str(round(ProfileEndInput(1), 2)) ' ' num2str(round(ProfileEndInput(2), 2)) ']']});

plotc(S,KSN); 
h1 = colorbar; hcb_title1 = get(h1, 'Title'); set(hcb_title1, 'String', 'K_{sn}');
colormap(h1, jet); 
caxis([0 1000])

% Topography and Local Relief
if DEMInput == 1
    subplot(3,4,2)
else
    subplot(3,3,2)
end
imageschs(DEM,G,'colormap','parula','ticklabel','nice','colorbarylabel','Slope [ ]','caxis',[0 1])

hold on

plot([ProfileStart(1,1); ProfileEnd(1,1)], [ProfileStart(1,2); ProfileEnd(1,2)], 'r:', 'LineWidth', 2)
plot([ProfileStart(i+1,1); ProfileEnd(i+1,1)], [ProfileStart(i+1,2); ProfileEnd(i+1,2)], 'r', 'LineWidth', 2)
plot([ProfileStart(i+j+1,1); ProfileEnd(i+j+1,1)], [ProfileStart(i+j+1,2); ProfileEnd(i+j+1,2)], 'r', 'LineWidth', 2)
plot(ProfileStart(:,1), ProfileStart(:,2), 'r', 'LineWidth', 2)
plot(ProfileEnd(:,1), ProfileEnd(:,2), 'r', 'LineWidth', 2)

xlabel('Longitude [m]')
ylabel('Latitude [m]')
title({'Local Relief';...
    ['[' num2str(round(ProfileStartInput(1), 2)) ' ' num2str(round(ProfileStartInput(2), 2))...
    '] to [' num2str(round(ProfileEndInput(1), 2)) ' ' num2str(round(ProfileEndInput(2), 2)) ']']});


h2 = colorbar; hcb_title2 = get(h2, 'Title'); set(hcb_title2, 'String', 'Slope [ ]');
colormap(parula); 

% Precipitation
if DEMInput == 1
    subplot(3,4,3)
    imageschs(DEM,DEM_Precipitation,'colormap','cool','ticklabel','nice','colorbarylabel','Precipitation [m/yr]')
    
    hold on
    
    plot([ProfileStart(1,1); ProfileEnd(1,1)], [ProfileStart(1,2); ProfileEnd(1,2)], 'r:', 'LineWidth', 2)
    plot([ProfileStart(i+1,1); ProfileEnd(i+1,1)], [ProfileStart(i+1,2); ProfileEnd(i+1,2)], 'r', 'LineWidth', 2)
    plot([ProfileStart(i+j+1,1); ProfileEnd(i+j+1,1)], [ProfileStart(i+j+1,2); ProfileEnd(i+j+1,2)], 'r', 'LineWidth', 2)
    plot(ProfileStart(:,1), ProfileStart(:,2), 'r', 'LineWidth', 2)
    plot(ProfileEnd(:,1), ProfileEnd(:,2), 'r', 'LineWidth', 2)
    
    xlabel('Long. [m]')
    ylabel('Lat. [m]')
    title({'Precipitation';...
        ['[' num2str(round(ProfileStartInput(1), 2)) ' ' num2str(round(ProfileStartInput(2), 2))...
        '] to [' num2str(round(ProfileEndInput(1), 2)) ' ' num2str(round(ProfileEndInput(2), 2)) ']']});
else
    subplot(3,3,3)
    
    imagesc(PrecipitationData)
    
    hold on
    
    plot([ProfileStart(1,1); ProfileEnd(1,1)], [ProfileStart(1,2); ProfileEnd(1,2)], 'r:', 'LineWidth', 2)
    plot([ProfileStart(i+1,1); ProfileEnd(i+1,1)], [ProfileStart(i+1,2); ProfileEnd(i+1,2)], 'r', 'LineWidth', 2)
    plot([ProfileStart(i+j+1,1); ProfileEnd(i+j+1,1)], [ProfileStart(i+j+1,2); ProfileEnd(i+j+1,2)], 'r', 'LineWidth', 2)
    plot(ProfileStart(:,1), ProfileStart(:,2), 'r', 'LineWidth', 2)
    plot(ProfileEnd(:,1), ProfileEnd(:,2), 'r', 'LineWidth', 2)
    
    xlabel('Longitude [m]')
    ylabel('Latitude [m]')
    title({'TRMM Precipitation';...
        ['[' num2str(round(ProfileStartInput(1), 2)) ' ' num2str(round(ProfileStartInput(2), 2))...
        '] to [' num2str(round(ProfileEndInput(1), 2)) ' ' num2str(round(ProfileEndInput(2), 2)) ']']});
    xlim(xlim(h1))
    ylim(ylim(h1))
end

caxis([0 5])
h3 = colorbar; hcb_title3 = get(h3, 'Title'); set(hcb_title3, 'String', 'Precipitation [m/yr]');
colormap(h3, cool); 



% Erosion (LEMs only, unless erosion maps exist at some point...)
if DEMInput == 1
    subplot(3,4,4)
    imageschs(DEM,DEM_Erosion,'ticklabel','nice','colorbarylabel','Erosion [m/yr]')
    
    hold on
    
    plot([ProfileStart(1,1); ProfileEnd(1,1)], [ProfileStart(1,2); ProfileEnd(1,2)], 'r:', 'LineWidth', 2)
    plot([ProfileStart(i+1,1); ProfileEnd(i+1,1)], [ProfileStart(i+1,2); ProfileEnd(i+1,2)], 'r', 'LineWidth', 2)
    plot([ProfileStart(i+j+1,1); ProfileEnd(i+j+1,1)], [ProfileStart(i+j+1,2); ProfileEnd(i+j+1,2)], 'r', 'LineWidth', 2)
    plot(ProfileStart(:,1), ProfileStart(:,2), 'r', 'LineWidth', 2)
    plot(ProfileEnd(:,1), ProfileEnd(:,2), 'r', 'LineWidth', 2)
    
    xlabel('Long. [m]')
    ylabel('Lat. [m]')
    title({'Erosion';...
        ['[' num2str(round(ProfileStartInput(1), 2)) ' ' num2str(round(ProfileStartInput(2), 2))...
        '] to [' num2str(round(ProfileEndInput(1), 2)) ' ' num2str(round(ProfileEndInput(2), 2)) ']']});
    h4 = colorbar; hcb_title4 = get(h4, 'Title'); set(hcb_title4, 'String', 'Erosion [m/yr]');
    %colormap(h4, turbo);
    caxis([-5e6 0])
end

% Ksn and Topography Swaths
if DEMInput == 1
    subplot(3,4,[5,8])
else
    subplot(3,3,[4,6])
end
x_temp = [Topography(:,1)', fliplr(Topography(:,1)')];
TopographySwath = [Topography(:,3)', fliplr(Topography(:,4)')];
fill(x_temp, TopographySwath, [0.86, 0.86, 0.86]);
hold on;
plot(Topography(:,1)', Topography(:,3)','r');
plot(Topography(:,1)', Topography(:,2)', 'k');
plot(Topography(:,1)', Topography(:,4)','b');
xlabel('Distance along Swath [m]')
ylabel('Elevation [m]')
title({'Topography and River Steepness';...
    ['[' num2str(round(ProfileStartInput(1), 2)) ' ' num2str(round(ProfileStartInput(2), 2))...
    '] to [' num2str(round(ProfileEndInput(1), 2)) ' ' num2str(round(ProfileEndInput(2), 2)) ']']});
box on
grid on

yyaxis right
plot(KsnValues(:,1) - RotationCentre(1), KsnValues(:,2), 'm-');

ylabel('Median K_{sn}')
legend([num2str(ProfileWidth) ' m Swath'], 'Max. Topography', 'Mean Topography', 'Min. Topography',...
    'Median K_{sn}')
ax=gca;
ax.XDir='reverse'

% Precipitation and Local Relief
if DEMInput == 1
    subplot(3,4,[9,12])
    
    plot(Precipitation(:,1)', Precipitation(:,2)','b');
    hold on
    plot(Erosion(:,1)', -Erosion(:,2)', 'r', 'LineWidth', 2)
    ylabel('Precipitation/Erosion [m/yr]')
    box on
    grid on
    
    yyaxis right
    
    plot(LocalGradient(:,1) - RotationCentre(1), LocalGradient(:,2), 'g-');
    xlabel('Distance along Swath [m]')
    ylabel('Slope [ ]')
    title({'Local Gradient and Precipitation';...
        ['[' num2str(round(ProfileStartInput(1), 2)) ' ' num2str(round(ProfileStartInput(2), 2))...
        '] to [' num2str(round(ProfileEndInput(1), 2)) ' ' num2str(round(ProfileEndInput(2), 2)) ']']});
    legend('Precipitation', 'Erosion', 'Local Relief')
    ax=gca;
    ax.XDir='reverse';
else
    subplot(3,3,[7,9])
    
    x_temp = [Precipitation(:,1)', fliplr(Precipitation(:,1)')];
    PrecipitationSwath = [Precipitation(:,3)', fliplr(Precipitation(:,4)')];
    fill(x_temp, PrecipitationSwath, [0.86, 0.86, 0.86]);
    hold on;
    plot(Precipitation(:,1)', Precipitation(:,3)','r');
    plot(Precipitation(:,1)', Precipitation(:,2)', 'k');
    plot(Precipitation(:,1)', Precipitation(:,4)','b');
    ylabel('Precipitation [m/yr]')
    box on
    grid on
    
    yyaxis right
    plot(LocalGradient(:,1) - RotationCentre(1), LocalGradient(:,2), 'g-');
    xlabel('Distance along Swath [m]')
    ylabel('Slope [ ]')
    title({'Local Gradient and Precipitation';...
        ['[' num2str(round(ProfileStartInput(1), 2)) ' ' num2str(round(ProfileStartInput(2), 2))...
        '] to [' num2str(round(ProfileEndInput(1), 2)) ' ' num2str(round(ProfileEndInput(2), 2)) ']']});
    legend([num2str(ProfileWidth) ' m Swath'], 'Max. Precipitation', 'Mean Precipitation', 'Min. Precipitation',...
        'Local Relief')
end

%% Clean up
% CASCADE file input
clear File_CASCADE File_DEM OutputVariables Path_CASCADE Path_DEM UserInput x_DEM y_DEM z_DEM
clear Erosion_DEM Precipitation_DEM localgradient

% Topography section
clear temp_d temp_x temp_y temp_z UserInout ProfileLatitude ProfileLongitude ProfileAngle Path_DEM i j k KsnOutput
clear HemisphereCheck File_DEM ProfileDistance ProfileElevation x_temp ProfileStart ProfileEnd TopographyOutput
clear MaximumTopography MinimumTopography MeanTopography TopographySwath UserInput
clear ProfileErosion ProfilePrecipitation temp_Erosion temp_Precipitation MeanErosion MeanPrecipitation

% Ksn section
clear Rotation A drainage_area FD h h1 h2 h3 h4 hcb_title1 hcb_title2 hcb_title3 hcb_title4
clear ksn KsnSwath RotationAngle Streams
clear Streams_rotated Streams_filtered x y x_range y_range fig

toc