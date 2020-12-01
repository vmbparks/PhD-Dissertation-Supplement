%% Creating Pecube.in, velocity.txt, and run_pecube.sbatch files
% Victoria Buford; May 10 2018
% vmb21@pitt.edu

% Have to change Pecubein.txt to Pecube.in outside of MATLAB
% - In Windows, Use Notepad ++ (NOT NOTEPAD); open Pecubein.txt:
% save as, Pecube.in
% - on Mac, right click on Pecubein.txt, click ?Get Info?. In the Name
% & Extension: blank, change the name to Pecube.in; when a dialog
% pops up, click ?Use .in?
%clc;clear;close;
%% Change per Move Model
% In user input now...
% ageinitiation=100; % should not be == ages(1)
% agemodelstart=50;

Xmin=0; % bottom left of start model
Xmax= 600;
%n=21; % number of deformation steps
% Xmin=350; % bottom left of start model
% Xmax= 1130;
% n=38; % number of deformation steps
email='vmb21@pitt.edu';

% don't normally change, but we can:
modelsizem=[500 1000]; % pecube pt grid spacing x and y;

% Which thermochronometers do you want....
% (a) AHe age (Default kinetics)
% (b) AHe age (20 um grain size)
% (c) AHe age (40 um grain size)
% (d) AHe age (70 um grain size)
% (e) AHe age (Low radiation damage)
% (f) AHe age (Moderate radiation damage)
% (g) AHe age (high radiation damage)
% (h) ZHe age
% (i) AFT age
% (j) ZFT age, D0 = 0.001_8, energy = 208.32768_8, grain_size = 3.158_8, geometry_factor = 55.0_8
% (k) Muscovite Ar/Ar age, D0 = 4.0e-8_8,  energy = 180.0_8, grain_size = 750.0_8, geometry_factor = 8.65_8
% (l) Ar in K-feldspar, D0 = 5.6, Ea = 120.0
% (m) Ar in Biotite, D0 = 160.0, Ea = 211.0
% (n) Ar in Muscovite, D0 = 13.2, Ea = 183.0
% (o) Ar in Hornblende, D0 = 14.0, Ea = 176.0
% (p) Apatite U-Th / Pb, D0 = 2.0e-8_8, energy = 230.0_8, grain_size = 50.0_8, geometry_factor = 55.0_8
% (q) Biotite, D0 = 2.0e-13_8, energy = 105.0_8, grain_size = 500.0_8, geometry_factor = 8.65_8
% (r) Ruite U-Pb, D0 = 1.6e-10_8, energy = 243.0_8, grain_size = 250.0_8, geometry_factor = 55.0_8
% (s) Titanite U-Pb, D0 = 1.1e-4_8, energy = 331.0_8, grain_size = 500.0_8, geometry_factor = 55.0_8
% (t) Zircon U-Pb, D0 = 7.8e-3_8, energy = 544.0_8, grain_size = 50.0_8, geometry_factor = 55.0_8
% (u) Titanite U-Th / He, D0 = 5.9e-3, energy = 188, grain_size = 250
% (v) ZHe, low damage D0z = 4.6
% (w) ZHe, med damage D0z = 0.3

AHea=1; %AHe
AHeb=0;
AHec=0;
AHed=0;
AHee=0;
AHef=0;
AHeg=0;
Zheh=1; % ZHe
AFTi=1; % AFT
ZFTj=1; % ZFT
MArArk=1; % MAr
ArKsparl=0;
ArBiom=1; %BAr
ArMuscn=0;
ArHorno=0;
ApUThp=0;
Bioq=0;
RuiUPbr=0;
TitUPbs=0;
ZirUPbt=0;
TitUThu=0;
ZHelowv=0;
ZHemedw=0;
ZHehighx=0;

ages=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change per Pecube Model
modelname=input('ModelName (for Velocity_file_name): ','s');%'McQ02N9_VelVar5';
abbrevmodelname=input('Abbrev. Model name (8 char display in queue): ','s');%'PN9VV5a3ef12'; % only 8 characters display on server list
NodeAssign=input('Node to assign: 133,134,135,136,\n esd_node1, esd_node2, arduino, \n argand, or 0 for no assignment: ','s');%'0'; % Node to assign: '38', '39','40' or  '0' for no assignment
%%
ef=input('e-folding depth [km]: ');% 12;
ao=input('thermal heat production[muW/m3] (ao): ');% 3.0;
atmlapse=input('atmospheric lapse rate [C/km](4.16): ');% 5.3;
% Navarro-Serrano_etal_2020- S. Peru - 4.16
%%
%%%% Velocities!
VV6=[42.18	40.83	39.45	38.27	37.03	35.90	34.43	33.07 ...
    32.11	31.14	30.18	28.80	27.21	25.83	24.61	23.13	21.75 ...
    20.53	19.41	18.13	16.60	15.09	14.06	12.86	12.01	11.39 ...
    10.76	10.15	9.49	8.74	7.94	6.81	5.75	4.33	3.11 ...
    1.83	0.83	0.00];
ageinitiation=input('Age of thermal initiation: ');
% agemodelstart=input('Age of Step 0: ');
% x=input('Enter ages in brackets [] for Steps 1:end: ');
% sz=size(x);
% if sz(1)>sz(2)
%     ages=[ageinitiation;agemodelstart;x];
% else
%     ages=[ageinitiation;agemodelstart;x'];
% end

x=input('Enter ages in brackets [] for Steps 0:end: ');
sz=size(x);
if sz(1)>sz(2)
    ages=[ageinitiation;x];
else
    ages=[ageinitiation;x'];
end

%%
% get grid names
disp('Select Grid Files')
[gridfilename,gridfolder]=uigetfile('*.dat','Multiselect','on');
%%
disp('Select Topography Files')
[topofilename,topofolder]=uigetfile('*.dat','Multiselect','on');
%%
disp('Select Save Folder');
SaveFolder=uigetdir;
SaveFolder=strcat(SaveFolder,'/');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% modelname=input('ModelName (for Velocity_file_name): ','s');%'McQ02N9_VelVar5';
% abbrevmodelname=input('Abbrev. Model name (8 char display in queue): ','s');
% VV6=[42.18	40.83	39.45	38.27	37.03	35.90	34.43	33.07 ...
%     32.11	31.14	30.18	28.80	27.21	25.83	24.61	23.13	21.75 ...
%     20.53	19.41	18.13	16.60	15.09	14.06	12.86	12.01	11.39 ...
%     10.76	10.15	9.49	8.74	7.94	6.81	5.75	4.33	3.11 ...
%     1.83	0.83	0.00];
% ageinitiation=input('Age of thermal initiation: ');
% agemodelstart=input('Age of Step 0: ');
% x=input('Enter ages in brackets [] for Steps 1:end: ');
% sz=size(x);
% if sz(1)>sz(2)
%     ages=[ageinitiation;agemodelstart;x];
% else
%     ages=[ageinitiation;agemodelstart;x'];
% end
% disp('Save Folder')
% SaveFolder=uigetdir;
% SaveFolder=strcat(SaveFolder,'/');
%% Model Setup
velname=strcat(modelname,'.txt');
modelsizepts=[(Xmax-Xmin)*(1000/modelsizem(1)) 5];
%%%%%%%%%% Thermal Parameters %%%%%%%%%%%%%%%%%%%%
% % (a) Model thickness (km)
% % (b) Number of z-node planes/layers in the z direction (integer)
% % NOTE: If this value is zero, Pecube will automatically define the z-node plane
% %       distribution such that the elements have a 1:1 (x/y to z) aspect ratio
% %       down to 5 km below the surface, 3:1 down to 15 km below the surface and
% %       ~9:1 down to the model base.
% % (c) Thermal conductivity (W/m K)
% % (d) Specific heat capacity (J/kg K)
% %   *NOTE: diffusivity is now caluclated in Pecube, rather than defined here*
% % (e) Crustal density (kg/m^3)
% % (f) Mantle density (kg/m^3)
% % ** SECOND LINE **
% % (g) Temperature at the base of the model (degrees C)
% % (h) Temperature at z=0 (degrees C)
% %   If lapse=0 this will be the surface temperature everywhere
% % (i) Atmospheric lapse rate (degrees C/km)
% %   NOTE: Positve lapse rate => decreasing T with elevation
% %         Negative lapse rate => increasing T with elevation
% % (j) Crustal volumetric heat production (uW/m^3)
% % (k) e-folding depth of crustal heat production (km)
% %   NOTE: Crustal heat production is constant at the given value for all nodes
% %     above sea level and decreases exponentially below msl. Also, if efold=0,
% %     then crustal heat production will be constant everywhere
% % (l) Mantle volumetric heat production (uW/m^3)
% %   NOTE: mantle HP not yet implemented - does nothing
% %         Also, mantle heat production is assumed to be constant
% % (m) Shear heating
% %  Set brittle shear heating constant below
% %  1 = on
% %  0 = off
% % (n) Shear heating constant (unitless)
% %   Scales shear heating within the brittle realm.
% %   Implemented in same form as used by F. Herman (02/08)
% %   1 = Full (unscaled) brittle shear heating
% %   0 = No brittle shear heating
therml1=[110 220 2.5 800 2500 3300];
therml2=[1300 23.0 atmlapse ao ef 0.01 0 0];

%%%%%%%%%% Thermochronometers to plot %%%%%%%%%%%%%%%%%%%%

% % % (a) AHe age (Default kinetics)
% AHea=1;
% % % (b) AHe age (20 um grain size)
% AHeb=0;
% % % (c) AHe age (40 um grain size)
% AHec=0;
% % % (d) AHe age (70 um grain size)
% AHed=0;
% % % (e) AHe age (Low radiation damage)
% AHee=0;
% % % (f) AHe age (Moderate radiation damage)
% AHef=0;
% % % (g) AHe age (high radiation damage)
% AHeg=0;
% % % (h) ZHe age
% Zheh=1;
% % % (i) AFT age
% AFTi=1;
% % % (j) ZFT age, D0 = 0.001_8, energy = 208.32768_8, grain_size = 3.158_8, geometry_factor = 55.0_8
% ZFTj=1;
% % % (k) Muscovite Ar/Ar age, D0 = 4.0e-8_8,  energy = 180.0_8, grain_size = 750.0_8, geometry_factor = 8.65_8
% MArArk=1;
% % % (l) Ar in K-feldspar, D0 = 5.6, Ea = 120.0
% ArKsparl=0;
% % % (m) Ar in Biotite, D0 = 160.0, Ea = 211.0
% ArBiom=0;
% % % (n) Ar in Muscovite, D0 = 13.2, Ea = 183.0
% ArMuscn=0;
% % % (o) Ar in Hornblende, D0 = 14.0, Ea = 176.0
% ArHorno=0;
% % % (p) Apatite U-Th / Pb, D0 = 2.0e-8_8, energy = 230.0_8, grain_size = 50.0_8, geometry_factor = 55.0_8
% ApUThp=0;
% % % (q) Biotite, D0 = 2.0e-13_8, energy = 105.0_8, grain_size = 500.0_8, geometry_factor = 8.65_8
% Bioq=0;
% % % (r) Ruite U-Pb, D0 = 1.6e-10_8, energy = 243.0_8, grain_size = 250.0_8, geometry_factor = 55.0_8
% RuiUPbr=0;
% % % (s) Titanite U-Pb, D0 = 1.1e-4_8, energy = 331.0_8, grain_size = 500.0_8, geometry_factor = 55.0_8
% TitUPbs=0;
% % % (t) Zircon U-Pb, D0 = 7.8e-3_8, energy = 544.0_8, grain_size = 50.0_8, geometry_factor = 55.0_8
% ZirUPbt=0;
% % % (u) Titanite U-Th / He, D0 = 5.9e-3, energy = 188, grain_size = 250
% TitUThu=0;
% % % (v) ZHe, low damage D0z = 4.6
% ZHelowv=0;
% % % (w) ZHe, med damage D0z = 0.3
% ZHemedw=0;
% % % (x) ZHe, high damage D0z = 4.6E+05
% ZHehighx=0;

tchrons=[AHea AHeb AHec AHed AHee AHef AHeg Zheh AFTi ZFTj MArArk ...
    ArKsparl ArBiom ArMuscn ArHorno ApUThp Bioq RuiUPbr TitUPbs ZirUPbt ...
    TitUThu ZHelowv ZHemedw ZHehighx];

%%%%%%%%%% check size %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear TotalFiles
if length(topofilename)==length(gridfilename) && length(topofilename)==length(ages)-1
    TotalFiles=length(topofilename);
    n=length(topofilename)-1;
    fprintf('Good to go\n')
else
    errordlg('Topo, Ages, and Grid lengths do not match!')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% Create Pecube files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% create velocity.txt file
% on windows notepad, will need \r\n instead of \n
fidvel=fopen(strcat(SaveFolder,velname),'w');
nrows=TotalFiles;
fprintf(fidvel,'10\n');
fprintf(fidvel,'5\n');
fprintf(fidvel,'10\n');
fprintf(fidvel,'1.0\n');
for row=2:nrows+1
    fprintf(fidvel,'%2.2f %s\n', ages(row),gridfilename{row-1});
end
fclose(fidvel);

%%%%%%%%%%% Create txts to copy into Pecube.in %%%%%%%%%%
%% Pecube.in complete.
% - In Windows, Use Notepad ++ (NOT NOTEPAD); open Pecubein.txt:
% save as, Pecube.in
% - on Mac, right click on Pecubein.txt, click ?Get Info?. In the Name
% & Extension: blank, change the name to Pecube.in; when a dialog
% pops up, click ?Use .in?

% Pecube params we change
% input 1, 2,4,6,7?,9, 10, 12a, 14(thermal), 17 (thermochronometers),
% 22 (name of velocity file.

% Need to change name to pecube.in
%topofilelist=strcat(modelname,'_Pecubein.txt');
fidcomp=fopen(strcat(SaveFolder,'Pecubein.txt'),'w');
nrows=TotalFiles;

%%%% pecube.in writeup
fprintf(fidcomp,'$ *********************************************************************************\n');
fprintf(fidcomp,'$ *** Pecube-D\n');
fprintf(fidcomp,'$ *****\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Input file for running Pecube-D.  University of Tuebingen, Germany (18 May, 2015).\n');
fprintf(fidcomp,'$ Report bugs noticed with this version to Todd Ehlers (todd.ehlers@uni-tuebingen.de).\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ This version of Pecube is based on the distrubtion by Jean Braun.  It has been\n');
fprintf(fidcomp,'$ Modified substatially to account for\n');
fprintf(fidcomp,'$ a.  Calculation of predicted ages for different thermo- geochronometer systems.\n');
fprintf(fidcomp,'$ b.  Many different options for user defined velocity input fields (e.g. McQuarrie and\n');
fprintf(fidcomp,'$     Ehlers, 2015).\n');
fprintf(fidcomp,'$ c.  Different output format options.\n');
fprintf(fidcomp,'$ d.  calculation of detrital cooling ages for user defined sample points on the topog.\n');
fprintf(fidcomp,'$      (e.g. Whipp et al., 2009)\n');
fprintf(fidcomp,'$ e.  Iterative Inversion of cooling ages for topographic change scenarios (e.g. Olen et al. 2012)\n');
fprintf(fidcomp,'$ f.  Monte Carlo Inversion of cooling ages to identify the denudation histories that\n');
fprintf(fidcomp,'$     that can produce observed ages.  (e.g. Thiede & Ehlers, 2013)\n');
fprintf(fidcomp,'$ g.  Coupling with CASCADE or IceCascade (e.g Yanites and Ehlers, 2013 EPSL)\n');
fprintf(fidcomp,'$ and numerous other significant changes to age prediction, heat production, shear\n');
fprintf(fidcomp,'$ heating, kinematics, and thermal history output\n');
fprintf(fidcomp,'$ are highlighted in the readme folders in  /docs folder. Different options have also been added\n');
fprintf(fidcomp,'$ for simulation other kinematic fields including an ellipsoidal exhumation field, as well as\n');
fprintf(fidcomp,'$ a coupling to 2D Move restoration files. Note - use of the 2D Move output for velocity fields\n');
fprintf(fidcomp,'$ also requires the program velocity.py.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Significant program changes have been made to this version with thanks to Willi Kappler,\n');
fprintf(fidcomp,'$ David Whipp, and Chris Spath.  If you use this program for publications the references that\n');
fprintf(fidcomp,'$ describe the methods used are:.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Braun , J., 2003. Pecube: A new finite element code to solve the 3D heat\n');
fprintf(fidcomp,'$  transport equation including the effects of a time-varying, finite\n');
fprintf(fidcomp,'$  amplitude surface topography.  Computers and Geosciences, v.29, pp.787-794.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Braun , J., 2002. Quantifying the effect of recent relief changes on age-elevation\n');
fprintf(fidcomp,'$  relationships.  Earth and Planetary Science Letters, v.200, pp.331-343.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Reference to use concerning program changes made by T. Ehlers group:\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Olen, S., Ehlers, T.A., Densmore, M.S., 2012, Limits to reconstructing paleotopography\n');
fprintf(fidcomp,'$ from ther- mochronometer data, J. Geophysical Res ??? Earth Surface, v. 117,\n');
fprintf(fidcomp,'$ doi:10.1029/2011/ JF001985\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Whipp, D.M. Jr., Ehlers, T.A., Braun, J., Spath, C.D., 2009, Effects of exhumation\n');
fprintf(fidcomp,'$ kinematics and topo- graphic evolution on detrital thermochronometer data,\n');
fprintf(fidcomp,'$ J. Geophysical Res. ??? Earth Surface, V. 114, F04021, doi:10.1029/2008JF001195.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Thiede, R.C., Ehlers, T.A., 2013, Large spatial and temporal variations in Himalayan\n');
fprintf(fidcomp,'$ denudation, Earth and Planetary Science Letters, 371-372, pp. 278-293.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ McQuarrie, N., and Ehlers, T.A., 2015 Influence of thrust belt geometry\n');
fprintf(fidcomp,'$ and shortening rate on thermochronometer cooling ages: Insights from the\n');
fprintf(fidcomp,'$ Bhutan Himalaya, Tectonics. 34, doi:10.1002/2014TC003783.\n');
fprintf(fidcomp,'$ Related programs to this distribution:\n');
fprintf(fidcomp,'$       * Bivar - this program takes topography output from cascade and Ice and formats them for input\n');
fprintf(fidcomp,'$       into Pecube-D\n');
fprintf(fidcomp,'$   * Cascade and Ice (Univ. Tuebingen modified versions, see Yanites and Ehlers, 2012 In review).\n');
fprintf(fidcomp,'$       These programs can be run prior to\n');
fprintf(fidcomp,'$       to Pecube-D to provide input topographies.\n');
fprintf(fidcomp,'$       * velocity.py - used to take 2D Move output and create transient velocity fields\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ NOTE - not all of the above files may be in this distribution.  We''re in the process of\n');
fprintf(fidcomp,'$ assembling a more complete distribution.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ You can add as many comment lines as you wish as long as they start with a\n');
fprintf(fidcomp,'$ dollar sign\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ *** HOW TO RUN / EXECUTE THIS PROGRAM ***\n');
fprintf(fidcomp,'$ This version of Pecube is compiled with MPI so that it can run on multicore machines.  This is not\n');
fprintf(fidcomp,'$ normally needed if you are doing a single simulation.  However, if you are using the\n');
fprintf(fidcomp,'$ Monte Carlo or Genetic search algorithm then it can run multiple jobs at once.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ The implications of this are that:\n');
fprintf(fidcomp,'$ 1. You have to have a version of MPI installed on your machine.  We are using OpenMPI, but others will likely\n');
fprintf(fidcomp,'$ work.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ 2. Compile the program using scons.  On our system this is done with scons --use-mpi, and then scons -c to clean\n');
fprintf(fidcomp,'$ out the .o files.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ 3. To run the program you need to have the pecube.in file in the same directory as the executable, and then\n');
fprintf(fidcomp,'$ start the job with mpi run.  For example, on our system we do:\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$  mpirun -n pecube\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Where -n is the number of cores you want to use.  If you are only running 1 job, then there is likely no speed\n');
fprintf(fidcomp,'$ difference if N=1 or N=8.  So, for 1 job, you write\n');
fprintf(fidcomp,'$ mpirun -n 1 pecube.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ *********************************************************************************\n');

fprintf(fidcomp,'\n');
fprintf(fidcomp,'\n');
fprintf(fidcomp,'$ Set mode of pecube operation\n');
fprintf(fidcomp,'$ Valid options:\n');
fprintf(fidcomp,'$ normal_mode: use this option if you just want to run pecube once\n');
fprintf(fidcomp,'$ error_iteration: use this option if you want to run a sequence of pecube simulation\n');
fprintf(fidcomp,'$       one after the other, in order to optimize the erosion rates iteravely\n');
fprintf(fidcomp,'$ monte_carlo: use this option if you want to run thousands of pecube simulation\n');
fprintf(fidcomp,'$       (some of them in parallel) in roder to optimize erosion rates using\n');
fprintf(fidcomp,'$       monte carlo randomisazion\n');
fprintf(fidcomp,'pecube_run_mode: normal_mode\n');
fprintf(fidcomp,'\n');
fprintf(fidcomp,'$ If run mode is error iteration:\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ Maximum number of error iteration\n');
fprintf(fidcomp,'error_iter_limit: 15\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ Error tolerance to exit iteration before limit is reached\n');
fprintf(fidcomp,'$ error_iter_tolerance: 1.0\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ Flag whether observables should be created from pecube or not\n');
fprintf(fidcomp,'$ This is currently not working\n');
fprintf(fidcomp,'$ error_iter_create_observables: off\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ If run mode is monte carlo:\n');
fprintf(fidcomp,'$\n');

fprintf(fidcomp,'$ The maximum value for the randomized erosion rates\n');
fprintf(fidcomp,'mc_max_erosion_rate: 4.0\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ The total number of simulations\n');
fprintf(fidcomp,'$mc_num_of_simulations: 10000\n');
fprintf(fidcomp,'mc_num_of_simulations: 20000\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ The tolerance for chi squared\n');
fprintf(fidcomp,'mc_tolerance_chi_squared: 4.0\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ Flag whether to check the minimum threshold and correct it\n');
fprintf(fidcomp,'mc_check_min_threshold: on\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ The actual minimum threshold value\n');
fprintf(fidcomp,'mc_min_threshold_factor: 0.20\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ This factor sets the scale factor for the jitter that is added to\n');
fprintf(fidcomp,'$ each erosion rate values in order to randomize them\n');
fprintf(fidcomp,'$ Lower values means less jitter and higher chance to get stuck in\n');
fprintf(fidcomp,'$ a local minimum\n');
fprintf(fidcomp,'$ Higher values means more jitter and higher chance to miss an optimal value\n');
fprintf(fidcomp,'$ but also avoids to get stuck in a local minimum\n');
fprintf(fidcomp,'$ Useful ranges: 0.01 to 0.1\n');
fprintf(fidcomp,'$ If this value is set to >= 1.0, then evolutionary / genetic algorithm is disabled\n');
fprintf(fidcomp,'$mc_random_jitter_factor: 0.01\n');
fprintf(fidcomp,'mc_random_jitter_factor: 0.05\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ The input file containing the ages and errors\n');
fprintf(fidcomp,'$ This must be a text file with semicolon separated columns\n');
fprintf(fidcomp,'mc_csv_input_file: mhsxls.csv\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ IMPORTANT: This marks the end of the pecube mode configuration\n');
fprintf(fidcomp,'pecube_end_run_mode\n');

fprintf(fidcomp,'\n');
fprintf(fidcomp,'$ *********************************************************************************\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 1) Name of the run (also the name of the folder in which the solution is stored)\n');
fprintf(fidcomp,'$ NOTE: You might need to create this folder manually before running Pecube.\n');
fprintf(fidcomp,'$ AKA: #1:output folder \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'output/Pecube-D\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 2) Number of topography files to be loaded (should be number of time steps+1)\n');
fprintf(fidcomp,'$ 0 = No topography file will be loaded\n');
fprintf(fidcomp,'$ 1 = The same topo file will be used for all time steps. (note relief can still\n');
fprintf(fidcomp,'$     change as specified below). Note the number of time steps, or steps in the\n');
fprintf(fidcomp,'$     tectonomorphic scenario is defined later.\n');
fprintf(fidcomp,'$ >1 = A new topo file will be loaded for each timestep\n');
fprintf(fidcomp,'$ If fewer topo files are listed then the number of time steps, the model will\n');
fprintf(fidcomp,'$   use the last topo file for the subsequent/remaining time steps.\n');
fprintf(fidcomp,'$ When multiple topography files are loaded Pecube will exponentially morph\n');
fprintf(fidcomp,'$   between the two topographies over the given time step.  tau, specified below\n');
fprintf(fidcomp,'$   determines the exponential rate of change in the topography.\n');
fprintf(fidcomp,'$ AKA: #2: # of topo files \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'%d\n',n+2);
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 3) Flag for topography input\n');
fprintf(fidcomp,'$ 1 = User will list all topography file names below, one file on each new line\n');
fprintf(fidcomp,'$ 2 = User specifies file prefix (see Input 4) and Pecube will load all\n');
fprintf(fidcomp,'$     files with that prefix plus a 4 digit number after it.\n');
fprintf(fidcomp,'$ 3 = All listed topography files are exported from 2d move\n');
fprintf(fidcomp,'$ AKA: #3: Topo files from move? \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'3\n');
fprintf(fidcomp,'\n');


fprintf(fidcomp,'$ (Input 4) Name of the topo file used\n');
fprintf(fidcomp,'$ "Nil" = Topography is assumed to be flat for that time step\n');
fprintf(fidcomp,'$ Otherwise the file should contain nx by ny elevation points (see below)\n');
fprintf(fidcomp,'$    defining the topography in meters\n');
fprintf(fidcomp,'$ Note that the evolution of this topography (in amplitude and elevation offset)\n');
fprintf(fidcomp,'$    can change at each time step, as specified below in Input 12.\n');
fprintf(fidcomp,'$ If multiple topography files are being loaded with a user defined filename\n');
fprintf(fidcomp,'$    prefix filename (e.g., option 2 above) then format is:\n');
fprintf(fidcomp,'$    prefix = "topo_input" (or another user-defined name)\n');
fprintf(fidcomp,'$    which will load files topo_input0000.dat, topo_input0001.dat, etc.\n');
fprintf(fidcomp,'$ Note: If detrital age calculation for a Cascade mesh is specified (value of 1\n');
fprintf(fidcomp,'$       for Input 18) then the topo files must be named topo_pecube_0000.dat,\n');
fprintf(fidcomp,'$       topo_pecube_0001.dat, etc or the prefix topo_pecube_ if automating it\n');
fprintf(fidcomp,'$ AKA: #4: Topo files \n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'%s\n', topofilename{1});
for row=1:nrows;
    fprintf(fidcomp,'%s\n', topofilename{row});
end
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 5) Coordinate system flag for Pecube input\n');
fprintf(fidcomp,'$ 1 = Degrees\n');
fprintf(fidcomp,'$ 2 = UTM (meters)\n');
fprintf(fidcomp,'$ AKA: #5: Coord syst flag \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'2 \n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 6) Number of points (nx, ny) in the longitude and latitude directions\n');
fprintf(fidcomp,'$   of the topography file being loaded.\n');
fprintf(fidcomp,'$ Note: The shell script make_topo.sh will output this information to the screen\n');
fprintf(fidcomp,'$   if you are using this to create your topo files from ArcGIS grids\n');
fprintf(fidcomp,'$ AKA: #6: model size [points] \n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'%d %d \n',modelsizepts);
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 7) Spacing of longitude and latitude points (in degrees or meters) in\n');
fprintf(fidcomp,'$   the topography input file.\n');
fprintf(fidcomp,'$ Note: The shell script make_topo.sh will also output this information if you\n');
fprintf(fidcomp,'$   use the script to export an ArcGIS DEM grid to Pecube format. Units of the\n');
fprintf(fidcomp,'$   values below should agree with what is specified in Input 5.\n');
fprintf(fidcomp,'$ AKA: #7: model size [m]\n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'%d %d \n',modelsizem);
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 8) Skipping factor (nskip) for points in the topo input file\n');
fprintf(fidcomp,'$ 1 = All points of the topography are used\n');
fprintf(fidcomp,'$ 2 = Every second point is used, etc.\n');
fprintf(fidcomp,'$ Note: nx, ny AND nskip define the resolution of the finite element mesh in the\n');
fprintf(fidcomp,'$   horizontal directions\n');
fprintf(fidcomp,'$ AKA: #8: skipping factor \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'1\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 9) Geographic location for the origin (bottom left corner) of the\n');
fprintf(fidcomp,'$   Pecube grid.\n');
fprintf(fidcomp,'$ Specify the longitude and latitude (in degrees or meters) of the bottom left\n');
fprintf(fidcomp,'$   corner of the topography file. Units must match above units.\n');
fprintf(fidcomp,'$ NOTE: a) You can set this value to be 0,0 for synthetic topography, or it\n');
fprintf(fidcomp,'$   can be 85670 (utm x), 983443 (utm y) or 109.756 (degrees long), 42.235\n');
fprintf(fidcomp,'$   (degrees lat) if you want Pecube to georeference the grid to your gegraphic\n');
fprintf(fidcomp,'$   area of study.\n');
fprintf(fidcomp,'$ NOTE: b) If you are using a DEM to generate the topography you want to specify\n');
fprintf(fidcomp,'$   an offset below that is 1/2 the topo file spacing specified in Input 7.\n');
fprintf(fidcomp,'$   (e.g., (DEM reolution / 2))\n');
fprintf(fidcomp,'$ AKA: #9: model start loc [m] \n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'%d %d \n', [Xmin*1000, 0]);
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 10) Number of time steps in the tectonomorphic scenario for your\n');
fprintf(fidcomp,'$   simulation\n');
fprintf(fidcomp,'$ An integer number (>= 1) is required. The value should be 1 less than the\n');
fprintf(fidcomp,'$   number of time step inputs defined in Input 12 below.\n');
fprintf(fidcomp,'$ Examples: a value of 1 will require two input lines for Input 12 below (a line\n');
fprintf(fidcomp,'$   for the starting time condition and one for the final time step condition).\n');
fprintf(fidcomp,'$   A value of 2 below will require 3 lines in Input 12 below. In this case, the\n');
fprintf(fidcomp,'$   first line would be the starting time condition, the second line would be\n');
fprintf(fidcomp,'$   the condition at some intermediate time, and final (third) line would be the\n');
fprintf(fidcomp,'$   final model condition.\n');
fprintf(fidcomp,'$ AKA: #10: # of tect steps (#def+1) \n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'%d \n', n+1);
fprintf(fidcomp,'\n');
% 
fprintf(fidcomp,'$ (Input 11) Erosional time scale (tau, in My) for topographic change\n');
fprintf(fidcomp,'$   This input allows the user to have non-linear morphing of topography with\n');
fprintf(fidcomp,'$   time.  A large value (e.g., 1000) will generate essentially linear changes\n');
fprintf(fidcomp,'$   between the input topography files.  Effectively, this is the e-folding time\n');
fprintf(fidcomp,'$   for the topographic evolution, a.k.a. the exponential decay rate of\n');
fprintf(fidcomp,'$   topography.\n');
fprintf(fidcomp,'$ AKA: #11: Erosional Time Scale \n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'1000\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 12a) Definition of the tectonomorphic time steps\n');
fprintf(fidcomp,'$ NOTE: The number of lines should be 1 greater than the value specified in\n');
fprintf(fidcomp,'$   Input 10.\n');
fprintf(fidcomp,'$ Each line formatted as follows:\n');
fprintf(fidcomp,'$ (a) Time (in My in the past)\n');
fprintf(fidcomp,'$   NOTES: (i) The first time step (first line) calculates a steady state\n');
fprintf(fidcomp,'$     thermal solution with the prescribed parameters\n');
fprintf(fidcomp,'$   (ii) Any transient features will occur between the previous listed time step\n');
fprintf(fidcomp,'$     line and the current time step line. For example, for a model with 3 time\n');
fprintf(fidcomp,'$     steps, a 50% decrease in topographic relief and change in the velocity\n');
fprintf(fidcomp,'$     field desired in the final time step would be listed on the last two\n');
fprintf(fidcomp,'$     lines, where the desired final relief and velocity field over that time\n');
fprintf(fidcomp,'$     are listed on the final time step line\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (b) Amplification factor for relief change\n');
fprintf(fidcomp,'$   1 = static topography\n');
fprintf(fidcomp,'$   2 = 200% increase in relief over this time step\n');
fprintf(fidcomp,'$   0.5 = 50% decrease in relief over this time step\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (c) Vertical offset factor (in km) for static topography elevation shifts\n');
fprintf(fidcomp,'$   during simulation.\n');
fprintf(fidcomp,'$   0 = No shift in surface elevations\n');
fprintf(fidcomp,'$   2 = Increase in all surface elevations by 2 km over this time step\n');
fprintf(fidcomp,'$   Why would you use this?  Well, if relief is 2 km, with a mean of 1 km, and\n');
fprintf(fidcomp,'$     relief is decreasing by 50% then if you specify a value of 0.5 (km) here\n');
fprintf(fidcomp,'$     it would shift your mean elevation such that it would remain at 1 km.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (d) Flag for output of time-temperature histories\n');
fprintf(fidcomp,'$   Enabling this will output temperature, time, x, y and z positions for all\n');
fprintf(fidcomp,'$     surface points at each step where listed.\n');
fprintf(fidcomp,'$   0 = No output of thermal history at the time step\n');
fprintf(fidcomp,'$   1 = Output of thermal history at the time step\n');
fprintf(fidcomp,'$   Note: Because the first time step is a steady state calculation, there is no\n');
fprintf(fidcomp,'$     thermal history available for the first time step.\n');
fprintf(fidcomp,'$   If the entire thermal history is wanted for surface points at t = 0Ma, then\n');
fprintf(fidcomp,'$     the user should set this flag for thermal output at the last time step\n');
fprintf(fidcomp,'$     specified below\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (e) Kinematic field flag (details of kinematic field specified in subsequent\n');
fprintf(fidcomp,'$   inputs)\n');
fprintf(fidcomp,'$   1 = vertical movement (erosion only)\n');
fprintf(fidcomp,'$   2 = uniform diagonal movement\n');
fprintf(fidcomp,'$   3 = listric fault\n');
fprintf(fidcomp,'$   4 = New Nepal thrust belt model for rotated model (Whipp testing - 10/07)\n');
fprintf(fidcomp,'$   5 = Parabolic uplift field (S. Olen testing - 06/10)\n');
fprintf(fidcomp,'$   6 = Ellipsoid uplift field (M. Schmiddunser - 07/11)\n');
fprintf(fidcomp,'$       An inner and outer ellipse and the uplift rates for the three corresponding\n');
fprintf(fidcomp,'$       areas are specified. The uplift rate between the inner and outer ellipse will\n');
fprintf(fidcomp,'$       be interpolated between the inner and the outer uplift rate\n');
fprintf(fidcomp,'$   7 = As 6, but uses a different function for uplift rate calculation that produces\n');
fprintf(fidcomp,'$       a smoothed uplift profile\n');
fprintf(fidcomp,'$   8 = velocity file names\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (f) Details of kinematics (Peklet)\n');
fprintf(fidcomp,'$   If e=1, value here is the erosion rate (mm/yr)\n');
fprintf(fidcomp,'$   If e=2, value here is the magnitude of the velocity vector at which material\n');
fprintf(fidcomp,'$     is moving laterally (mm/yr) in an Eulerian reference frame\n');
fprintf(fidcomp,'$   If e=3, value here is the maximum slip velocity on fault (mm/yr)\n');
fprintf(fidcomp,'$   If e=4, enter 1, values for velocities are computed within code and scaled\n');
fprintf(fidcomp,'$     by 1 here\n');
fprintf(fidcomp,'$   If e=5, enter 1, values for velocities are based on maximum and minimum\n');
fprintf(fidcomp,'$           velocities defined at (k) and (l)\n');
fprintf(fidcomp,'$   If e=6 or e=7, enter value for uplift rate inside the inner ellipse (mm/yr)\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ ADDITIONAL OPTIONAL PARAMETERS (depending on kinematic field used)\n');
fprintf(fidcomp,'$ (g)\n');
fprintf(fidcomp,'$   If e=1, enter 0, no additional input required\n');
fprintf(fidcomp,'$   If e=2, enter fault dip angle theta (degrees). This is the angle from\n');
fprintf(fidcomp,'$     horizontal (positive down) defining the dip of the velocity vector.\n');
fprintf(fidcomp,'$   If e=3, enter the longitude or utm x position of one endpoint of the listric\n');
fprintf(fidcomp,'$     fault trace.\n');
fprintf(fidcomp,'$   Note: If you think of traveling along a line that starts at the first point\n');
fprintf(fidcomp,'$     and ends at the second, the fault would dip off to the left of that line\n');
fprintf(fidcomp,'$   Note:  This value must be entered in km (S. Olen, 06/21/2010)\n');
fprintf(fidcomp,'$   If e=4, enter the horizontal convergence rate (mm/yr) across the Main\n');
fprintf(fidcomp,'$     Frontal Thrust. Note: Fault geometries are hard coded in Pecube\n');
fprintf(fidcomp,'$   If e=5, enter the x-value for the lower right point of the parabola axis (km)\n');
fprintf(fidcomp,'$   If e=6 or e=7, enter the x-value for the first focus (in km)\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (h)\n');
fprintf(fidcomp,'$   If e=1, enter 0, no additional input required\n');
fprintf(fidcomp,'$   If e=2, enter angle phi (degrees), the azimuth of the velocity vector in the\n');
fprintf(fidcomp,'$     x-y plane.\n');
fprintf(fidcomp,'$   If e=3, enter the latitude or utm y position of the first endpoint of the\n');
fprintf(fidcomp,'$     listric fault trace in item (g) above.\n');
fprintf(fidcomp,'$   Note:  This value must be entered in km (S. Olen, 06/21/2010)\n');
fprintf(fidcomp,'$   If e=4, enter the horizontal convergence rate (mm/yr) across the Main\n');
fprintf(fidcomp,'$     Boundary Thrust\n');
fprintf(fidcomp,'$   If e=5, enter the y-value for the lower right point of the parabola axis (km)\n');
fprintf(fidcomp,'$   If e=6 or e=7, enter the y-value for the first focus (in km)\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (i)\n');
fprintf(fidcomp,'$   If e=1 or e=2, enter 0, no additional input required\n');
fprintf(fidcomp,'$   If e=3, enter the longitude or utm x of the second end point of the listric\n');
fprintf(fidcomp,'$     fault.\n');
fprintf(fidcomp,'$   Note:  This value must be entered in km (S. Olen, 06/21/2010)\n');
fprintf(fidcomp,'$   If e=4, enter the horizontal convergence rate (mm/yr) across the Main\n');
fprintf(fidcomp,'$     Central Thrust\n');
fprintf(fidcomp,'$   If e=5, enter the x-value for the upper left point of the parabola axis (km)\n');
fprintf(fidcomp,'$   If e=6 or e=7, enter the x-value for the second focus (in km)\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (j)\n');
fprintf(fidcomp,'$   If e=1 or 2, enter 0, no additional input required\n');
fprintf(fidcomp,'$   If e=3, enter the latitutde or utm y of the second endpoint of the listric\n');
fprintf(fidcomp,'$     fault.\n');
fprintf(fidcomp,'$   Note:  This value must be entered in km (S. Olen, 06/21/2010)\n');
fprintf(fidcomp,'$   If e=4, enter the horizontal extension rate (mm/yr) across the South Tibetan\n');
fprintf(fidcomp,'$     Detachment\n');
fprintf(fidcomp,'$   If e=5, enter the y-value for the upper left point of the parabola axis (km)\n');
fprintf(fidcomp,'$   If e=6 or e=7, enter the y-value for the second focus (in km)\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (k)\n');
fprintf(fidcomp,'$   If e=1 or 2, enter 0, no additional input required\n');
fprintf(fidcomp,'$   If e=3, enter the soling depth (km) of the fault. Note: Fault has an\n');
fprintf(fidcomp,'$     exponential shape.\n');
fprintf(fidcomp,'$   If e=4, enter 0 or 1 for whether or not you want underplating in the Sub-\n');
fprintf(fidcomp,'$     Himalaya during this time step (0=no; 1=yes)\n');
fprintf(fidcomp,'$   If e=5, enter the minimum velocity (mmyr-1)\n');
fprintf(fidcomp,'$   If e=6 or e=7, enter semi-major axis of inner elipse (in km)\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (l)\n');
fprintf(fidcomp,'$   If e=1 or 2, enter 0, no additional input required\n');
fprintf(fidcomp,'$   If e=3, enter surface dip angle of the fault in degrees\n');
fprintf(fidcomp,'$   If e=4, enter 0 or 1 for whether or not you want underplating in the Lesser\n');
fprintf(fidcomp,'$     Himalaya during this time step (0=no; 1=yes)\n');
fprintf(fidcomp,'$   If e=5, enter the maximum velocity (mmyr-1)\n');
fprintf(fidcomp,'$   If e=6 or e=7, enter semi-major axis of outer elipse (in km)\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ (m)\n');
fprintf(fidcomp,'$   If e=1,2,4,5, enter 0., no additonal input required\n');
fprintf(fidcomp,'$   if e=3, enter additional uplift rate (mm/yr, z velocity)\n');
fprintf(fidcomp,'$   If e=6 or e=7, enter uplift rate (mm/yr) outside of outer ellipse.\n');
fprintf(fidcomp,'$       The uplift rate between inner and outer ellipse will be interpolated\n');
fprintf(fidcomp,'$       between the two values$ (Input 12 b) If e=8 in Input 12 above, this defines the min and max\n');
fprintf(fidcomp,'$ AKA: #12a: Tect Time Steps \n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'$(a)   (b) (c) (d) (e) (f) (g) (h) (i) (j) (k) (l) (m)\n');
params=[1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.];
%agelist=strcat(modelname,'_age.txt');
%nrows=TotalFiles;
formatspec='%2.2f \t %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f  %1.0f  %1.0f  %1.0f\n';
for row=1:nrows+1;
    fprintf(fidcomp,formatspec, ages(row),params);
end

fprintf(fidcomp,'\n');
fprintf(fidcomp,'$ (Input 12 b) If e=8 in Input 12 above, this defines the min and max\n');
fprintf(fidcomp,'$ AKA: #12b: allowed velocities \n'); % we don't change
fprintf(fidcomp,'$ allowed velocity values in [mm/years]\n');
fprintf(fidcomp,'\n');
fprintf(fidcomp,'$ vx_min vx_max\n');
fprintf(fidcomp,'-100.0 100.0\n');
fprintf(fidcomp,'\n');
fprintf(fidcomp,'$ vy_min vy_max\n');
fprintf(fidcomp,'-100.0 100.0\n');
fprintf(fidcomp,'\n');
fprintf(fidcomp,'$ vz_min vz_max\n');
fprintf(fidcomp,'-100.0 100.0\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 13) Isostasy (NOTE: All values listed on one line)\n');
fprintf(fidcomp,'$  Remove input (b) and (c), which are included in input 14... necessary\n');
fprintf(fidcomp,'$  for isostacy to turn on (S. Olen, 1/20/2010)\n');
fprintf(fidcomp,'$ (a) Flag for isostasy\n');
fprintf(fidcomp,'$   1 = isostasy on\n');
fprintf(fidcomp,'$   0 = isostasy off\n');
fprintf(fidcomp,'$ (b) Crustal density (kg/m^3)  ***Remove***\n');
fprintf(fidcomp,'$ (c) Mantle density (kg/m^3)  ***Remove***\n');
fprintf(fidcomp,'$ (d) Young modulus (Pa)\n');
fprintf(fidcomp,'$ (e) Poisson''s ratio\n');
fprintf(fidcomp,'$ (f) Elastic plate thickness (*km*)\n');
fprintf(fidcomp,'$ (g) Size of the FFT grid for elastic rebound calculations (typically 1024 1024\n');
fprintf(fidcomp,'$   but must be a power of 2)\n');
fprintf(fidcomp,'$ AKA: #13: Isostasy \n') ;% we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'0 1.d11 0.25 15. 1024 1024\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 14) Thermal model input parameters\n');
fprintf(fidcomp,'$ NOTE: Values on two lines\n');
fprintf(fidcomp,'$ Note: Pecube currently assumes homogeneous medium.\n');
fprintf(fidcomp,'$ ** FIRST LINE **\n');
fprintf(fidcomp,'$ (a) Model thickness (km)\n');
fprintf(fidcomp,'$ (b) Number of z-node planes/layers in the z direction (integer)\n');
fprintf(fidcomp,'$ NOTE: If this value is zero, Pecube will automatically define the z-node plane\n');
fprintf(fidcomp,'$       distribution such that the elements have a 1:1 (x/y to z) aspect ratio\n');
fprintf(fidcomp,'$       down to 5 km below the surface, 3:1 down to 15 km below the surface and\n');
fprintf(fidcomp,'$       ~9:1 down to the model base.\n');
fprintf(fidcomp,'$ (c) Thermal conductivity (W/m K)\n');
fprintf(fidcomp,'$ (d) Specific heat capacity (J/kg K)\n');
fprintf(fidcomp,'$   *NOTE: diffusivity is now caluclated in Pecube, rather than defined here*\n');
fprintf(fidcomp,'$ (e) Crustal density (kg/m^3)\n');
fprintf(fidcomp,'$ (f) Mantle density (kg/m^3)\n');
fprintf(fidcomp,'$ ** SECOND LINE **\n');
fprintf(fidcomp,'$ (g) Temperature at the base of the model (degrees C)\n');
fprintf(fidcomp,'$ (h) Temperature at z=0 (degrees C)\n');
fprintf(fidcomp,'$   If lapse=0 this will be the surface temperature everywhere\n');
fprintf(fidcomp,'$ (i) Atmospheric lapse rate (degrees C/km)\n');
fprintf(fidcomp,'$   NOTE: Positve lapse rate => decreasing T with elevation\n');
fprintf(fidcomp,'$         Negative lapse rate => increasing T with elevation\n');
fprintf(fidcomp,'$ (j) Crustal volumetric heat production (uW/m^3)\n');
fprintf(fidcomp,'$ (k) e-folding depth of crustal heat production (km)\n');
fprintf(fidcomp,'$   NOTE: Crustal heat production is constant at the given value for all nodes\n');
fprintf(fidcomp,'$     above sea level and decreases exponentially below msl. Also, if efold=0,\n');
fprintf(fidcomp,'$     then crustal heat production will be constant everywhere\n');
fprintf(fidcomp,'$ (l) Mantle volumetric heat production (uW/m^3)\n');
fprintf(fidcomp,'$   NOTE: mantle HP not yet implemented - does nothing\n');
fprintf(fidcomp,'$         Also, mantle heat production is assumed to be constant\n');
fprintf(fidcomp,'$ (m) Shear heating\n');
fprintf(fidcomp,'$  Set brittle shear heating constant below\n');
fprintf(fidcomp,'$  1 = on\n');
fprintf(fidcomp,'$  0 = off\n');
fprintf(fidcomp,'$ (n) Shear heating constant (unitless)\n');
fprintf(fidcomp,'$   Scales shear heating within the brittle realm.\n');
fprintf(fidcomp,'$   Implemented in same form as used by F. Herman (02/08)\n');
fprintf(fidcomp,'$   1 = Full (unscaled) brittle shear heating\n');
fprintf(fidcomp,'$   0 = No brittle shear heating\n');
fprintf(fidcomp,'$ AKA: #14: Thermal Params \n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'%d %d %2.1f %d %d %d \n',therml1);
fprintf(fidcomp,'%d %2.2f %2.2f %1.2f %2.1f %.3f %d %d \n',therml2);
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 15) Thermal model input parameters for Nepal model geometry\n');
fprintf(fidcomp,'$   **NOTE** This is not used unless the geometry flag above is set to 4.\n');
fprintf(fidcomp,'$   On each line, there are five values. The values are for the Indian Shield,\n');
fprintf(fidcomp,'$   Sub-Himalaya, Lesser Himalaya, Greater Himalaya and Tethyan Himalaya\n');
fprintf(fidcomp,'$   Line 1: Volumetric heat production (uW/m^3)\n');
fprintf(fidcomp,'$   Line 2: Thermal conductivity (W/m K)\n');
fprintf(fidcomp,'$   Line 3: Rock density (kg/m^3)\n');
fprintf(fidcomp,'$   Line 4: Specific heat capacity (J/kg K)\n');
fprintf(fidcomp,'$ AKA: #15:Nepal Thermal \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'0.8 0.8 0.8 1.9 0.8 \n');
fprintf(fidcomp,'2.75 2.75 2.75 2.75 2.75 \n');
fprintf(fidcomp,'2700. 2700. 2700. 2700. 2700. \n');
fprintf(fidcomp,'1000. 1000. 1000. 1000. 1000. \n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 16) Option to read in thermochron data and compare to predicted ages\n');
fprintf(fidcomp,'$ First specify the number of data files for comparison\n');
fprintf(fidcomp,'$ 0 = No data file(s) will be read in\n');
fprintf(fidcomp,'$ For each file name that is specified, the file format should be as follows:\n');
fprintf(fidcomp,'$ First line in file = number of samples (and lines) in rest of file.\n');
fprintf(fidcomp,'$ Each line after that is for an individual sample and should contain (space\n');
fprintf(fidcomp,'$   separated)\n');
fprintf(fidcomp,'$ (a) Sample longitude or utm x (value in km)\n');
fprintf(fidcomp,'$ (b) Sample latitude or utm y (value in km)\n');
fprintf(fidcomp,'$ (c) Sample elevation\n');
fprintf(fidcomp,'$ (d) Flag for type of AHe age to predict:\n');
fprintf(fidcomp,'$   1=Default diffusion kinetics; 2-4=Use grain size of 20, 40 or 70 um, resp.\n');
fprintf(fidcomp,'$   5-7=Use low, moderate or high eU (radiation damage) values resp.\n');
fprintf(fidcomp,'$   (Schuster et al., 2006)\n');
fprintf(fidcomp,'$   NOTE: These values can be modifed in the Mad_He.f90 subroutine\n');
fprintf(fidcomp,'$   Comments in that subroutine further explain the differences above\n');
fprintf(fidcomp,'$ (e) AHe age (Ma), negative age if non-existant\n');
fprintf(fidcomp,'$ (f) AHe age error, 1s.d. (Ma), use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (g) AFT age (Ma), negative age if non-existant\n');
fprintf(fidcomp,'$ (h) AFT age error, 1s.d. (Ma), use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (i) ZHe age (Ma), negative age if non-existant\n');
fprintf(fidcomp,'$ (j) ZHe age error, 1s.d. (Ma), use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (k) ZFT age (Ma), negative age if non-existant\n');
fprintf(fidcomp,'$ (l) ZFT age error, 1s.d. (Ma), use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (m) MAr age (Ma), negative age if non-existant\n');
fprintf(fidcomp,'$ (n) MAr age error, 1s.d. (Ma), use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (o) Ar in K-feldspar, negative age if non-existant\n');
fprintf(fidcomp,'$ (p) Ar in K-feldspar, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (q) Ar in Biotite, negative age if non-existant\n');
fprintf(fidcomp,'$ (r) Ar in Biotite, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (s) Ar in Muscovite, negative age if non-existant\n');
fprintf(fidcomp,'$ (t) Ar in Muscovite, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (u) Ar in Hornblende, negative age if non-existant\n');
fprintf(fidcomp,'$ (v) Ar in Hornblende, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (w) Apatite U-Th / Pb, negative age if non-existant\n');
fprintf(fidcomp,'$ (x) Apatite U-Th / Pb, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (y) Biotite, negative age if non-existant\n');
fprintf(fidcomp,'$ (z) Biotite, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (a1) Ruite U-Pb, negative age if non-existant\n');
fprintf(fidcomp,'$ (b1) Ruite U-Pb, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (c1) Titanite U-Pb, negative age if non-existant\n');
fprintf(fidcomp,'$ (d1) Titanite U-Pb, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (e1) Zircon U-Pb, negative age if non-existant\n');
fprintf(fidcomp,'$ (f1) Zircon U-Pb, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (g1) Titanite U-Th, negative age if non-existant\n');
fprintf(fidcomp,'$ (h1) Titanite U-Th, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (i1) ZHe, low damage, negative age if non-existant\n');
fprintf(fidcomp,'$ (j1) ZHe, low damage, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (k1) ZHe, med damage, negative age if non-existant\n');
fprintf(fidcomp,'$ (l1) ZHe, med damage, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (m1) ZHe, high damage, negative age if non-existant\n');
fprintf(fidcomp,'$ (n1) ZHe, high damage, use 0 if previous value is negative\n');
fprintf(fidcomp,'$ (o1) Sample ID\n');
fprintf(fidcomp,'$ AKA: #16:Read Tchron? \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'0 \n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ cstmt_ahe_errinc.dat\n');
fprintf(fidcomp,'$ cstmt_aft.dat\n');

fprintf(fidcomp,'$ NOTE: Input 17 and 18 are used for detrital age calculation.  Instructions\n');
fprintf(fidcomp,'$ on use are in 00README_detrital_ages file if user needs more information.\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 17) Flags for which ages to output\n');
fprintf(fidcomp,'$ 0 = Does not calculate or output predicted ages for this system\n');
fprintf(fidcomp,'$ 1 = Calculates and outputs specified system''s ages\n');
fprintf(fidcomp,'$ NOTE: See Mad_He.f90 subroutine to modify the predicted AHe ages below\n');
fprintf(fidcomp,'$ (a) AHe age (Default kinetics)\n');
fprintf(fidcomp,'$ (b) AHe age (20 um grain size)\n');
fprintf(fidcomp,'$ (c) AHe age (40 um grain size)\n');
fprintf(fidcomp,'$ (d) AHe age (70 um grain size)\n');
fprintf(fidcomp,'$ (e) AHe age (Low radiation damage)\n');
fprintf(fidcomp,'$ (f) AHe age (Moderate radiation damage)\n');
fprintf(fidcomp,'$ (g) AHe age (high radiation damage)\n');
fprintf(fidcomp,'$ (h) ZHe age\n');
fprintf(fidcomp,'$ (i) AFT age\n');
fprintf(fidcomp,'$ (j) ZFT age, D0 = 0.001_8, energy = 208.32768_8, grain_size = 3.158_8, geometry_factor = 55.0_8\n');
fprintf(fidcomp,'$ (k) Muscovite Ar/Ar age, D0 = 4.0e-8_8,  energy = 180.0_8, grain_size = 750.0_8, geometry_factor = 8.65_8\n');
fprintf(fidcomp,'$ (l) Ar in K-feldspar, D0 = 5.6, Ea = 120.0\n');
fprintf(fidcomp,'$ (m) Ar in Biotite, D0 = 160.0, Ea = 211.0\n');
fprintf(fidcomp,'$ (n) Ar in Muscovite, D0 = 13.2, Ea = 183.0\n');
fprintf(fidcomp,'$ (o) Ar in Hornblende, D0 = 14.0, Ea = 176.0\n');
fprintf(fidcomp,'$ (p) Apatite U-Th / Pb, D0 = 2.0e-8_8, energy = 230.0_8, grain_size = 50.0_8, geometry_factor = 55.0_8\n');
fprintf(fidcomp,'$ (q) Biotite, D0 = 2.0e-13_8, energy = 105.0_8, grain_size = 500.0_8, geometry_factor = 8.65_8\n');
fprintf(fidcomp,'$ (r) Ruite U-Pb, D0 = 1.6e-10_8, energy = 243.0_8, grain_size = 250.0_8, geometry_factor = 55.0_8\n');
fprintf(fidcomp,'$ (s) Titanite U-Pb, D0 = 1.1e-4_8, energy = 331.0_8, grain_size = 500.0_8, geometry_factor = 55.0_8\n');
fprintf(fidcomp,'$ (t) Zircon U-Pb, D0 = 7.8e-3_8, energy = 544.0_8, grain_size = 50.0_8, geometry_factor = 55.0_8\n');
fprintf(fidcomp,'$ (u) Titanite U-Th / He, D0 = 5.9e-3, energy = 188, grain_size = 250\n');
fprintf(fidcomp,'$ (v) ZHe, low damage D0z = 4.6\n');
fprintf(fidcomp,'$ (w) ZHe, med damage D0z = 0.3\n');
fprintf(fidcomp,'$ (x) ZHe, high damage D0z = 4.6E+05\n');
fprintf(fidcomp,'$ AKA: #17: which thermochronometers \n'); % we DO change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'$(a)(b)(c)(d)(e)(f)(g)(h)(i)(j)(k)(l)(m)(n)(o)(p)(q)(r)(s)(t)(u)(v)(w)(x)\n');
formatspec2='%1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f %1.0f \n';
fprintf(fidcomp,formatspec2,tchrons);
fprintf(fidcomp,'\n');
 
fprintf(fidcomp,'$ (Input 18) Flag to calculate detrital age distributions for catchments\n');
fprintf(fidcomp,'$ First line:  0 = no detrital calculation; 1 = detrital calcuation\n');
fprintf(fidcomp,'$ Note:  If a series of CASCADE topographies were loaded in then\n');
fprintf(fidcomp,'$ set Input 18 of the Pecube.in file to be 1, Pecube will output the detrital\n');
fprintf(fidcomp,'$ ages of every cascade catchment at every timestep.  These files will be created\n');
fprintf(fidcomp,'$ in the ''catchments'' folder within your output run directory.  It will create the\n');
fprintf(fidcomp,'$ ''catchments'' folder if it does not exist there.  The files will be named as\n');
fprintf(fidcomp,'$ ''Timestep_0001_Catchment_0001.dat'' and so on for all catchments and timesteps.\n');
fprintf(fidcomp,'$ All the program needs to run properly is to have the tecplot formatted cascade\n');
fprintf(fidcomp,'$ output files for every timestep (eg. ''topo_tec_0001.dat'') in\n');
fprintf(fidcomp,'$ the ''output/Cascade'' directory\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Specifying user defined basins in a file causes Pecube to open the file and read it.\n');
fprintf(fidcomp,'$ The file should contain lines with the following syntax:\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ <xpos> <ypos> <age type> <Nil or filename> <yes/Yes or no/No>\n');
fprintf(fidcomp,'$ Ex: 64.8495 128.5548 1 Nil No\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Where xpos is the x value of the basin outlet and ypos is the y value.  The age type is a number 1-11 with the\n');
fprintf(fidcomp,'$ following coding:\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ 1 = Apatite Helium Age - Farley, 2000\n');
fprintf(fidcomp,'$ 2 = Apatite Helium Age - Small grain size\n');
fprintf(fidcomp,'$ 3 = Apatite Helium Age - Medium grain size\n');
fprintf(fidcomp,'$ 4 = Apatite Helium Age - Large grain size\n');
fprintf(fidcomp,'$ 5 = Apatite Helium Age - Low radiation damage\n');
fprintf(fidcomp,'$ 6 = Apatite Helium Age - Medium radiation damage\n');
fprintf(fidcomp,'$ 7 = Apatite Helium Age - High radiation damage\n');
fprintf(fidcomp,'$ 8 = Apatite Fission Track Age\n');
fprintf(fidcomp,'$ 9 = Zircon Helium Age\n');
fprintf(fidcomp,'$ 10 = Zircon Fission Track Age\n');
fprintf(fidcomp,'$ 11 = Muscovite Age\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ The next entry for each basin line can be ''Nil'' or a filename.  If ''Nil'' is specified, then Pecube uses\n');
fprintf(fidcomp,'$ the TAPES-Grid method of finding upstream points of the basin outlet.  Then, it will write out the x, y,\n');
fprintf(fidcomp,'$ and z positions along with the age data for the upstream points into the main run output directory with\n');
fprintf(fidcomp,'$ the naming convention of ''Timestep_0001_Basin_X_64.8495_Y_128.5548.dat''.  Also, the PDF for each basin is\n');
fprintf(fidcomp,'$ written to the folder ''pdf_data'' within the run output directory with the naming convention of\n');
fprintf(fidcomp,'$ ''Timestep_0001_Basin_X_64.8495_Y_128.5548_Agetype_01_pdf.dat''. If a filename is specified as the entry\n');
fprintf(fidcomp,'$ on the line, then Pecube will open this file and read in each line of data with the format:\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ <age> <error>\n');
fprintf(fidcomp,'$ Ex: 35.756469 1.07269407\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Where the age and error are absolute values. The PDF for each of these basins is created in the ''pdf_data''\n');
fprintf(fidcomp,'$ folder in the run output directory with the naming convention of ''Basin_X_64.8495_Y_128.5548_Agetype_01_pdf.dat''.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ The final entry on each line is a ''yes'' or ''no'' on whether the user wants to run the Monte Carlo test for that specified\n');
fprintf(fidcomp,'$ basin.  For the Monte Carlo test to run properly there must be two basins with the same outlet point (x and y positions) and\n');
fprintf(fidcomp,'$ EACH must have a ''yes'' to run the monte carlo routine. Also, the age types MUST be the same or the Monte Carlo test is not run.\n');
fprintf(fidcomp,'$ An example of correct syntax for the Monte Carlo routine to run properly is as follows:\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ 64.8495 128.5548 7 Nil Yes\n');
fprintf(fidcomp,'$ 64.8495 128.5548 7 datafile.txt Yes\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Please note that if any or all of these criteria for running the Monte Carlo test are not met, then the program simple does not\n');
fprintf(fidcomp,'$ run the comparison (skips it) and does not output anything for it.\n');
fprintf(fidcomp,'$\n');
fprintf(fidcomp,'$ Note: The subdirectories in the run output directory where Pecube writes most of these files will be\n');
fprintf(fidcomp,'$       automatically created by Pecube if they do not exist already.\n');
fprintf(fidcomp,'\n');
fprintf(fidcomp,'$pdf_tester_for_data.txt\n');
fprintf(fidcomp,'$ AKA: #18: Detrital Age Flag \n'); % we don't change
fprintf(fidcomp,'0\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 19) Minimum number of nodes for a catchment to be output\n');
fprintf(fidcomp,'$ This is a threshold value of nodes that a catchment needs to have in order for an\n');
fprintf(fidcomp,'$   output file to be written for that catchment at that timestep.\n');
fprintf(fidcomp,'$ AKA: #19: # of nodes \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'100\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 20) Name of directory where Cascade tecplot formatted output files are located\n');
fprintf(fidcomp,'$ to be read in by Pecube for PDF calculation\n');
fprintf(fidcomp,'$ Note: This only matters if a ''1'' is selected above (Input 18) for use of Cascade catchments\n');
fprintf(fidcomp,'$ The program will disregard whatever is here if ''0'' or a filename is specified\n');
fprintf(fidcomp,'$ AKA: #20: Cascade ouput folder \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'output/Cascade \n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 21) Name of the temperature file if needed. Otherwise Nil\n');
fprintf(fidcomp,'$ AKA: #21: temperature file \n'); % we don't change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'Nil\n');
fprintf(fidcomp,'\n');

fprintf(fidcomp,'$ (Input 22) If e=8 in Input 12 above, filename for velocity file. Otherwise Nil\n');
fprintf(fidcomp,'$ AKA: #22: Velocity.txt file \n');% we change
fprintf(fidcomp,'\n');
fprintf(fidcomp,'%s\n',velname);
fclose(fidcomp);
%% SBatch file

fidsbatch=fopen(strcat(SaveFolder,'run_pecube.txt'),'w');

fprintf(fidsbatch,'#!/bin/bash -l \n');
fprintf(fidsbatch,'## Run script for pecube monte carlo on esd slurm \n');
fprintf(fidsbatch,' \n');
fprintf(fidsbatch,' \n');
fprintf(fidsbatch,'## General configuration options \n');
fprintf(fidsbatch,'#SBATCH -J %s \n',abbrevmodelname);
fprintf(fidsbatch,'#SBATCH -e Pecube_Error%%j \n');
fprintf(fidsbatch,'#SBATCH -o Pecube_Screen%%j \n');
fprintf(fidsbatch,'#SBATCH --mail-user=%s \n',email);
fprintf(fidsbatch,'#SBATCH --mail-type=ALL \n');
fprintf(fidsbatch,' \n');
fprintf(fidsbatch,' \n');
fprintf(fidsbatch,'## Machine and CPU configuration \n');
fprintf(fidsbatch,'## Number of tasks per job: \n');
fprintf(fidsbatch,'#SBATCH -n 1 \n');
fprintf(fidsbatch,'## Number of nodes: \n');
fprintf(fidsbatch,'#SBATCH -N 1 \n');
% if strcmp(NodeAssign,'0');
%         fprintf(fidsbatch,' \n');
% else
%   fprintf(fidsbatch,'#SBATCH -w u-005-s0%s \n',NodeAssign);
%   fprintf(fidsbatch,' \n');
% end
if strcmp(NodeAssign,'133')||strcmp(NodeAssign,'134') || strcmp(NodeAssign,'135')|| strcmp(NodeAssign,'136')
    fprintf(fidsbatch,'#SBATCH -w u-017-s%s \n',NodeAssign);
    %fprintf(fidsbatch,'#SBATCH -w u-005-s0%s \n',NodeAssign);
    fprintf(fidsbatch,' \n');
elseif strcmp(NodeAssign,'esd_node1')||strcmp(NodeAssign,'esd_node2')||strcmp(NodeAssign,'arduino')
    fprintf(fidsbatch,'#SBATCH -w %s \n',NodeAssign);
    fprintf(fidsbatch,' \n');
else
     fprintf(fidsbatch,' \n');
end

fprintf(fidsbatch,'module load pecube \n');
fprintf(fidsbatch,' \n');

fprintf(fidsbatch,'pecube \n');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% If you just want to copy and paste into pecube.in.. %%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% list of topo
% topofilelist=strcat(modelname,'_topo.txt');
% fidtopo=fopen(strcat(velfolder,topofilelist),'w');
% nrows=TotalFiles;
% 
% fprintf(fidtopo,'% Input 4: Topo files \n'); % we DO change
% fprintf(fidtopo,'%s\n', topofilename{1});
% for row=1:nrows;
%     fprintf(fidtopo,'%s\n', topofilename{row});
% end
% fprintf(fidtopo,'\n');
% fclose(fidtopo)
% 
% %% list of tect ages
% fidtect = fopen(strcat(velfolder,modelname,'_tectlist.txt')
% fprintf(fidtect,'% Input 12a: Tect Time Steps \n'); % we DO change
% params=[1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.];
% agelist=strcat(modelname,'_age.txt');
% nrows=TotalFiles;
% formatspec='%2.2f \t %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f\n';
% fprintf(fidtect,formatspec, ageinitiation,params);
% for row=1:nrows;
%     fprintf(fidtect,formatspec, ages(row),params);
% end
% fprintf(fidtect,'\n');
% fclose(fidtect)

% %% list of ages and params
% params=[1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.];
% agelist=strcat(modelname,'_age.txt');
% fidages=fopen(strcat(velfolder,agelist),'w');
% 
% nrows=TotalFiles;
% formatspec='%2.2f \t %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f   %1.0f\n';
% 
% fprintf(fidages,formatspec, ageinitiation,params);
% 
% for row=1:nrows
%     fprintf(fidages,formatspec, ages(row),params);
% end
% fclose(fidages);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%ORIGINAL PECUBE CODE/WRITEUP%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % $ *********************************************************************************
% % $ *** Pecube-D
% % $ *****
% % $
% % $ Input file for running Pecube-D.  University of Tuebingen, Germany (18 May, 2015).
% % $ Report bugs noticed with this version to Todd Ehlers (todd.ehlers@uni-tuebingen.de).
% % $
% % $ This version of Pecube is based on the distrubtion by Jean Braun.  It has been
% % $ Modified substatially to account for
% % $ a.  Calculation of predicted ages for different thermo- geochronometer systems.
% % $ b.  Many different options for user defined velocity input fields (e.g. McQuarrie and
% % $     Ehlers, 2015).
% % $ c.  Different output format options.
% % $ d.  calculation of detrital cooling ages for user defined sample points on the topog.
% % $      (e.g. Whipp et al., 2009)
% % $ e.  Iterative Inversion of cooling ages for topographic change scenarios (e.g. Olen et al. 2012)
% % $ f.  Monte Carlo Inversion of cooling ages to identify the denudation histories that
% % $     that can produce observed ages.  (e.g. Thiede & Ehlers, 2013)
% % $ g.  Coupling with CASCADE or IceCascade (e.g Yanites and Ehlers, 2013 EPSL)
% % $ and numerous other significant changes to age prediction, heat production, shear
% % $ heating, kinematics, and thermal history output
% % $ are highlighted in the readme folders in  /docs folder. Different options have also been added
% % $ for simulation other kinematic fields including an ellipsoidal exhumation field, as well as
% % $ a coupling to 2D Move restoration files. Note - use of the 2D Move output for velocity fields
% % $ also requires the program velocity.py.
% % $
% % $ Significant program changes have been made to this version with thanks to Willi Kappler,
% % $ David Whipp, and Chris Spath.  If you use this program for publications the references that
% % $ describe the methods used are:.
% % $
% % $ Braun , J., 2003. Pecube: A new finite element code to solve the 3D heat
% % $  transport equation including the effects of a time-varying, finite
% % $  amplitude surface topography.  Computers and Geosciences, v.29, pp.787-794.
% % $
% % $ Braun , J., 2002. Quantifying the effect of recent relief changes on age-elevation
% % $  relationships.  Earth and Planetary Science Letters, v.200, pp.331-343.
% % $
% % $ Reference to use concerning program changes made by T. Ehlers group:
% % $
% % $ Olen, S., Ehlers, T.A., Densmore, M.S., 2012, Limits to reconstructing paleotopography
% % $ from ther- mochronometer data, J. Geophysical Res ??? Earth Surface, v. 117,
% % $ doi:10.1029/2011/ JF001985
% % $
% % $ Whipp, D.M. Jr., Ehlers, T.A., Braun, J., Spath, C.D., 2009, Effects of exhumation
% % $ kinematics and topo- graphic evolution on detrital thermochronometer data,
% % $ J. Geophysical Res. ??? Earth Surface, V. 114, F04021, doi:10.1029/2008JF001195.
% % $
% % $ Thiede, R.C., Ehlers, T.A., 2013, Large spatial and temporal variations in Himalayan
% % $ denudation, Earth and Planetary Science Letters, 371-372, pp. 278-293.
% % $
% % $ McQuarrie, N., and Ehlers, T.A., 2015 Influence of thrust belt geometry
% % $ and shortening rate on thermochronometer cooling ages: Insights from the
% % $ Bhutan Himalaya, Tectonics. 34, doi:10.1002/2014TC003783.
% % $ Related programs to this distribution:
% % $       * Bivar - this program takes topography output from cascade and Ice and formats them for input
% % $       into Pecube-D
% % $   * Cascade and Ice (Univ. Tuebingen modified versions, see Yanites and Ehlers, 2012 In review).
% % $       These programs can be run prior to
% % $       to Pecube-D to provide input topographies.
% % $       * velocity.py - used to take 2D Move output and create transient velocity fields
% % $
% % $ NOTE - not all of the above files may be in this distribution.  We're in the process of
% % $ assembling a more complete distribution.
% % $
% % $ You can add as many comment lines as you wish as long as they start with a
% % $ dollar sign
% % $
% % $ *** HOW TO RUN / EXECUTE THIS PROGRAM ***
% % $ This version of Pecube is compiled with MPI so that it can run on multicore machines.  This is not
% % $ normally needed if you are doing a single simulation.  However, if you are using the
% % $ Monte Carlo or Genetic search algorithm then it can run multiple jobs at once.
% % $
% % $ The implications of this are that:
% % $ 1. You have to have a version of MPI installed on your machine.  We are using OpenMPI, but others will likely
% % $ work.
% % $
% % $ 2. Compile the program using scons.  On our system this is done with scons --use-mpi, and then scons -c to clean
% % $ out the .o files.
% % $
% % $ 3. To run the program you need to have the pecube.in file in the same directory as the executable, and then
% % $ start the job with mpi run.  For example, on our system we do:
% % $
% % $  mpirun -n pecube
% % $
% % $ Where -n is the number of cores you want to use.  If you are only running 1 job, then there is likely no speed
% % $ difference if N=1 or N=8.  So, for 1 job, you write
% % $ mpirun -n 1 pecube.
% % $
% % $ *********************************************************************************
% % 
% % 
% % $ Set mode of pecube operation
% % $ Valid options:
% % $ normal_mode: use this option if you just want to run pecube once
% % $ error_iteration: use this option if you want to run a sequence of pecube simulation
% % $       one after the other, in order to optimize the erosion rates iteravely
% % $ monte_carlo: use this option if you want to run thousands of pecube simulation
% % $       (some of them in parallel) in roder to optimize erosion rates using
% % $       monte carlo randomisazion
% % pecube_run_mode: normal_mode
% % 
% % $ If run mode is error iteration:
% % $
% % $ Maximum number of error iteration
% % error_iter_limit: 15
% % 
% % $ Error tolerance to exit iteration before limit is reached
% % error_iter_tolerance: 1.0
% % 
% % $ Flag wether observables should be created from pecube or not
% % $ This is currently not working
% % error_iter_create_observables: off
% % 
% % $ If run mode is monte carlo:
% % $
% % $ The maximum value for the randomized erosion rates
% % mc_max_erosion_rate: 4.0
% % 
% % $ The total number of simulations
% % $mc_num_of_simulations: 10000
% % mc_num_of_simulations: 20000
% % 
% % $ The tolerance for chi squared
% % mc_tolerance_chi_squared: 4.0
% % 
% % $ Flag wether to check the minimum threshold and correct it
% % mc_check_min_threshold: on
% % 
% % $ The actual minimum threshold value
% % mc_min_threshold_factor: 0.20
% % 
% % $ This factor sets the scale factor for the jitter that is added to
% % $ each erosion rate values in order to randomize them
% % $ Lower values means less jitter and higher chance to get stuck in
% % $ a local minimum
% % $ Higher values means more jitter and higher chance to miss an optimal value
% % $ but also avoids to get stuck in a local minimum
% % $ Useful ranges: 0.01 to 0.1
% % $ If this value is set to >= 1.0, then evolutionary / genetic algorithm is disabled
% % $mc_random_jitter_factor: 0.01
% % mc_random_jitter_factor: 0.05
% % 
% % $ The input file containing the ages and errors
% % $ This must be a text file with semicolon separated columns
% % mc_csv_input_file: mhsxls.csv
% % 
% % $ IMPORTANT: This marks the end of the pecube mode configuration
% % pecube_end_run_mode
% % 
% % $ (Input 1) Name of the run (also the name of the folder in which the solution is stored)
% % $ NOTE: You might need to create this folder manually before running Pecube.
% % 
% % output/Pecube-D
% % 
% % $ (Input 2) Number of topography files to be loaded (should be number of time steps+1)
% % $ 0 = No topography file will be loaded
% % $ 1 = The same topo file will be used for all time steps. (note relief can still
% % $     change as specified below). Note the number of time steps, or steps in the
% % $     tectonomorphic scenario is defined later.
% % $ >1 = A new topo file will be loaded for each timestep
% % $ If fewer topo files are listed then the number of time steps, the model will
% % $   use the last topo file for the subsequent/remaining time steps.
% % $ When multiple topography files are loaded Pecube will exponentially morph
% % $   between the two topographies over the given time step.  tau, specified below
% % $   determines the exponential rate of change in the topography.
% % 
% % 30
% % 
% % $ (Input 3) Flag for topography input
% % $ 1 = User will list all topography file names below, one file on each new line
% % $ 2 = User specifies file prefix (see Input 4) and Pecube will load all
% % $     files with that prefix plus a 4 digit number after it.
% % $ 3 = All listed topography files are exported from 2d move
% % $
% % 
% % 3
% % 
% % $ (Input 4) Name of the topo file used
% % $ "Nil" = Topography is assumed to be flat for that time step
% % $ Otherwise the file should contain nx by ny elevation points (see below)
% % $    defining the topography in meters
% % $ Note that the evolution of this topography (in amplitude and elevation offset)
% % $    can change at each time step, as specified below in Input 12.
% % $ If multiple topography files are being loaded with a user defined filename
% % $    prefix filename (e.g., option 2 above) then format is:
% % $    prefix = "topo_input" (or another user-defined name)
% % $    which will load files topo_input0000.dat, topo_input0001.dat, etc.
% % $ Note: If detrital age calculation for a Cascade mesh is specified (value of 1
% % $       for Input 18) then the topo files must be named topo_pecube_0000.dat,
% % $       topo_pecube_0001.dat, etc or the prefix topo_pecube_ if automating it
% % 
% % McQ02N2_topo_Step0.dat
% % McQ02N2_topo_Step0.dat
% % McQ02N2_topo_Step1.dat
% % McQ02N2_topo_Step2.dat
% % McQ02N2_topo_Step3.dat
% % McQ02N2_topo_Step4.dat
% % McQ02N2_topo_Step5.dat
% % McQ02N2_topo_Step6.dat
% % McQ02N2_topo_Step7.dat
% % McQ02N2_topo_Step8.dat
% % McQ02N2_topo_Step9.dat
% % McQ02N2_topo_Step10.dat
% % McQ02N2_topo_Step11.dat
% % McQ02N2_topo_Step12.dat
% % McQ02N2_topo_Step13.dat
% % McQ02N2_topo_Step14.dat
% % McQ02N2_topo_Step15.dat
% % McQ02N2_topo_Step16.dat
% % McQ02N2_topo_Step17.dat
% % McQ02N2_topo_Step18.dat
% % McQ02N2_topo_Step19.dat
% % McQ02N2_topo_Step20.dat
% % McQ02N2_topo_Step21.dat
% % McQ02N2_topo_Step22.dat
% % McQ02N3_topo_Step23.dat
% % McQ02N4_topo_Step24.dat
% % McQ02N4_topo_Step25.dat
% % McQ02N4_topo_Step26.dat
% % McQ02N4_topo_Step27.dat
% % McQ02N4_topo_Step28.dat
% % 
% % $ (Input 5) Coordinate system flag for Pecube input
% % $ 1 = Degrees
% % $ 2 = UTM (meters)
% % 
% % 2
% % 
% % $ (Input 6) Number of points (nx, ny) in the longitude and latitude directions
% % $   of the topography file being loaded.
% % $ Note: The shell script make_topo.sh will output this information to the screen
% % $   if you are using this to create your topo files from ArcGIS grids
% % 
% % 1560 5
% % 
% % $ (Input 7) Spacing of longitude and latitude points (in degrees or meters) in
% % $   the topography input file.
% % $ Note: The shell script make_topo.sh will also output this information if you
% % $   use the script to export an ArcGIS DEM grid to Pecube format. Units of the
% % $   values below should agree with what is specified in Input 5.
% % 
% % 500 1000
% % 
% % $ (Input 8) Skipping factor (nskip) for points in the topo input file
% % $ 1 = All points of the topography are used
% % $ 2 = Every second point is used, etc.
% % $ Note: nx, ny AND nskip define the resolution of the finite element mesh in the
% % $   horizontal directions
% % 
% % 1
% % 
% % $ (Input 9) Geographic location for the origin (bottom left corner) of the
% % $   Pecube grid.
% % $ Specify the longitude and latitude (in degrees or meters) of the bottom left
% % $   corner of the topography file. Units must match above units.
% % $ NOTE: a) You can set this value to be 0,0 for synthetic topography, or it
% % $   can be 85670 (utm x), 983443 (utm y) or 109.756 (degrees long), 42.235
% % $   (degrees lat) if you want Pecube to georeference the grid to your gegraphic
% % $   area of study.
% % $ NOTE: b) If you are using a DEM to generate the topography you want to specify
% % $   an offset below that is 1/2 the topo file spacing specified in Input 7.
% % $   (e.g., (DEM reolution / 2))
% % 
% % 315000 0.0
% % 
% % $ (Input 10) Number of time steps in the tectonomorphic scenario for your
% % $   simulation
% % $ An integer number (>= 1) is required. The value should be 1 less than the
% % $   number of time step inputs defined in Input 12 below.
% % $ Examples: a value of 1 will require two input lines for Input 12 below (a line
% % $   for the starting time condition and one for the final time step condition).
% % $   A value of 2 below will require 3 lines in Input 12 below. In this case, the
% % $   first line would be the starting time condition, the second line would be
% % $   the condition at some intermediate time, and final (third) line would be the
% % $   final model condition.
% % 
% % 29
% % 
% % $ (Input 11) Erosional time scale (tau, in My) for topographic change
% % $   This input allows the user to have non-linear morphing of topography with
% % $   time.  A large value (e.g., 1000) will generate essentially linear changes
% % $   between the input topography files.  Effectively, this is the e-folding time
% % $   for the topographic evolution, a.k.a. the exponential decay rate of
% % $   topography.
% % 
% % 1000
% % 
% % $ (Input 12a) Definition of the tectonomorphic time steps
% % $ NOTE: The number of lines should be 1 greater than the value specified in
% % $   Input 10.
% % $ Each line formatted as follows:
% % $ (a) Time (in My in the past)
% % $   NOTES: (i) The first time step (first line) calculates a steady state
% % $     thermal solution with the prescribed parameters
% % $   (ii) Any transient features will occur between the previous listed time step
% % $     line and the current time step line. For example, for a model with 3 time
% % $     steps, a 50% decrease in topographic relief and change in the velocity
% % $     field desired in the final time step would be listed on the last two
% % $     lines, where the desired final relief and velocity field over that time
% % $     are listed on the final time step line
% % $
% % $ (b) Amplification factor for relief change
% % $   1 = static topography
% % $   2 = 200% increase in relief over this time step
% % $   0.5 = 50% decrease in relief over this time step
% % $
% % $ (c) Vertical offset factor (in km) for static topography elevation shifts
% % $   during simulation.
% % $   0 = No shift in surface elevations
% % $   2 = Increase in all surface elevations by 2 km over this time step
% % $   Why would you use this?  Well, if relief is 2 km, with a mean of 1 km, and
% % $     relief is decreasing by 50% then if you specify a value of 0.5 (km) here
% % $     it would shift your mean elevation such that it would remain at 1 km.
% % $
% % $ (d) Flag for output of time-temperature histories
% % $   Enabling this will output temperature, time, x, y and z positions for all
% % $     surface points at each step where listed.
% % $   0 = No output of thermal history at the time step
% % $   1 = Output of thermal history at the time step
% % $   Note: Because the first time step is a steady state calculation, there is no
% % $     thermal history available for the first time step.
% % $   If the entire thermal history is wanted for surface points at t = 0Ma, then
% % $     the user should set this flag for thermal output at the last time step
% % $     specified below
% % $
% % $ (e) Kinematic field flag (details of kinematic field specified in subsequent
% % $   inputs)
% % $   1 = vertical movement (erosion only)
% % $   2 = uniform diagonal movement
% % $   3 = listric fault
% % $   4 = New Nepal thrust belt model for rotated model (Whipp testing - 10/07)
% % $   5 = Parabolic uplift field (S. Olen testing - 06/10)
% % $   6 = Ellipsoid uplift field (M. Schmiddunser - 07/11)
% % $       An inner and outer ellipse and the uplift rates for the three corresponding
% % $       areas are specified. The uplift rate between the inner and outer ellipse will
% % $       be interpolated between the inner and the outer uplift rate
% % $   7 = As 6, but uses a different function for uplift rate calculation that produces
% % $       a smoothed uplift profile
% % $   8 = velocity file names
% % $
% % $ (f) Details of kinematics (Peklet)
% % $   If e=1, value here is the erosion rate (mm/yr)
% % $   If e=2, value here is the magnitude of the velocity vector at which material
% % $     is moving laterally (mm/yr) in an Eulerian reference frame
% % $   If e=3, value here is the maximum slip velocity on fault (mm/yr)
% % $   If e=4, enter 1, values for velocities are computed within code and scaled
% % $     by 1 here
% % $   If e=5, enter 1, values for velocities are based on maximum and minimum
% % $           velocities defined at (k) and (l)
% % $   If e=6 or e=7, enter value for uplift rate inside the inner ellipse (mm/yr)
% % $
% % $ ADDITIONAL OPTIONAL PARAMETERS (depending on kinematic field used)
% % $ (g)
% % $   If e=1, enter 0, no additional input required
% % $   If e=2, enter fault dip angle theta (degrees). This is the angle from
% % $     horizontal (positive down) defining the dip of the velocity vector.
% % $   If e=3, enter the longitude or utm x position of one endpoint of the listric
% % $     fault trace.
% % $   Note: If you think of traveling along a line that starts at the first point
% % $     and ends at the second, the fault would dip off to the left of that line
% % $   Note:  This value must be entered in km (S. Olen, 06/21/2010)
% % $   If e=4, enter the horizontal convergence rate (mm/yr) across the Main
% % $     Frontal Thrust. Note: Fault geometries are hard coded in Pecube
% % $   If e=5, enter the x-value for the lower right point of the parabola axis (km)
% % $   If e=6 or e=7, enter the x-value for the first focus (in km)
% % $
% % $ (h)
% % $   If e=1, enter 0, no additional input required
% % $   If e=2, enter angle phi (degrees), the azimuth of the velocity vector in the
% % $     x-y plane.
% % $   If e=3, enter the latitude or utm y position of the first endpoint of the
% % $     listric fault trace in item (g) above.
% % $   Note:  This value must be entered in km (S. Olen, 06/21/2010)
% % $   If e=4, enter the horizontal convergence rate (mm/yr) across the Main
% % $     Boundary Thrust
% % $   If e=5, enter the y-value for the lower right point of the parabola axis (km)
% % $   If e=6 or e=7, enter the y-value for the first focus (in km)
% % $
% % $ (i)
% % $   If e=1 or e=2, enter 0, no additional input required
% % $   If e=3, enter the longitude or utm x of the second end point of the listric
% % $     fault.
% % $   Note:  This value must be entered in km (S. Olen, 06/21/2010)
% % $   If e=4, enter the horizontal convergence rate (mm/yr) across the Main
% % $     Central Thrust
% % $   If e=5, enter the x-value for the upper left point of the parabola axis (km)
% % $   If e=6 or e=7, enter the x-value for the second focus (in km)
% % $
% % $ (j)
% % $   If e=1 or 2, enter 0, no additional input required
% % $   If e=3, enter the latitutde or utm y of the second endpoint of the listric
% % $     fault.
% % $   Note:  This value must be entered in km (S. Olen, 06/21/2010)
% % $   If e=4, enter the horizontal extension rate (mm/yr) across the South Tibetan
% % $     Detachment
% % $   If e=5, enter the y-value for the upper left point of the parabola axis (km)
% % $   If e=6 or e=7, enter the y-value for the second focus (in km)
% % $
% % $ (k)
% % $   If e=1 or 2, enter 0, no additional input required
% % $   If e=3, enter the soling depth (km) of the fault. Note: Fault has an
% % $     exponential shape.
% % $   If e=4, enter 0 or 1 for whether or not you want underplating in the Sub-
% % $     Himalaya during this time step (0=no; 1=yes)
% % $   If e=5, enter the minimum velocity (mmyr-1)
% % $   If e=6 or e=7, enter semi-major axis of inner elipse (in km)
% % $
% % $ (l)
% % $   If e=1 or 2, enter 0, no additional input required
% % $   If e=3, enter surface dip angle of the fault in degrees
% % $   If e=4, enter 0 or 1 for whether or not you want underplating in the Lesser
% % $     Himalaya during this time step (0=no; 1=yes)
% % $   If e=5, enter the maximum velocity (mmyr-1)
% % $   If e=6 or e=7, enter semi-major axis of outer elipse (in km)
% % $
% % $ (m)
% % $   If e=1,2,4,5, enter 0., no additonal input required
% % $   if e=3, enter additional uplift rate (mm/yr, z velocity)
% % $   If e=6 or e=7, enter uplift rate (mm/yr) outside of outer ellipse.
% % $       The uplift rate between inner and outer ellipse will be interpolated
% % $       between the two values$ (Input 12 b) If e=8 in Input 12 above, this defines the min and max
% % 
% % 
% % $(a)   (b) (c) (d) (e) (f) (g) (h) (i) (j) (k) (l) (m)
% % 100.   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 50.0   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 49.87  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 47.82  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 45.55  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 43.53  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 41.61  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 39.85  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 37.76  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 35.84  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 33.53  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 31.85  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 29.95  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 27.97  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 25.95  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 24.03  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 22.78  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 20.98  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 19.19  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 16.69  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 14.19  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 11.69  1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 8.94   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 6.86   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 4.38   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 3.63   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 3.00   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 2.34   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 1.27   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 0.00   1.  0.  1.  8  1.0 0.  0.  0.  0.  0.  0.  0.  0.  0.  0.
% % 
% % 
% % 
% % $ (Input 12 b) If e=8 in Input 12 above, this defines the min and max
% % 
% % $ allowed velocity values in [mm/years]
% % $ vx_min vx_max
% % -100.0 100.0
% % 
% % $ vy_min vy_max
% % -100.0 100.0
% % 
% % $ vz_min vz_max
% % -100.0 100.0
% % 
% % $ (Input 13) Isostasy (NOTE: All values listed on one line)
% % $  Remove input (b) and (c), which are included in input 14... necessary
% % $  for isostacy to turn on (S. Olen, 1/20/2010)
% % $ (a) Flag for isostasy
% % $   1 = isostasy on
% % $   0 = isostasy off
% % $ (b) Crustal density (kg/m^3)  ***Remove***
% % $ (c) Mantle density (kg/m^3)  ***Remove***
% % $ (d) Young modulus (Pa)
% % $ (e) Poisson's ratio
% % $ (f) Elastic plate thickness (*km*)
% % $ (g) Size of the FFT grid for elastic rebound calculations (typically 1024 1024
% % $   but must be a power of 2)
% % 
% % 0 1.d11 0.25 15. 1024 1024
% % 
% % $ (Input 14) Thermal model input parameters
% % $ NOTE: Values on two lines
% % $ Note: Pecube currently assumes homogeneous medium.
% % $ ** FIRST LINE **
% % $ (a) Model thickness (km)
% % $ (b) Number of z-node planes/layers in the z direction (integer)
% % $ NOTE: If this value is zero, Pecube will automatically define the z-node plane
% % $       distribution such that the elements have a 1:1 (x/y to z) aspect ratio
% % $       down to 5 km below the surface, 3:1 down to 15 km below the surface and
% % $       ~9:1 down to the model base.
% % $ (c) Thermal conductivity (W/m K)
% % $ (d) Specific heat capacity (J/kg K)
% % $   *NOTE: diffusivity is now caluclated in Pecube, rather than defined here*
% % $ (e) Crustal density (kg/m^3)
% % $ (f) Mantle density (kg/m^3)
% % $ ** SECOND LINE **
% % $ (g) Temperature at the base of the model (degrees C)
% % $ (h) Temperature at z=0 (degrees C)
% % $   If lapse=0 this will be the surface temperature everywhere
% % $ (i) Atmospheric lapse rate (degrees C/km)
% % $   NOTE: Positve lapse rate => decreasing T with elevation
% % $         Negative lapse rate => increasing T with elevation
% % $ (j) Crustal volumetric heat production (uW/m^3)
% % $ (k) e-folding depth of crustal heat production (km)
% % $   NOTE: Crustal heat production is constant at the given value for all nodes
% % $     above sea level and decreases exponentially below msl. Also, if efold=0,
% % $     then crustal heat production will be constant everywhere
% % $ (l) Mantle volumetric heat production (uW/m^3)
% % $   NOTE: mantle HP not yet implemented - does nothing
% % $         Also, mantle heat production is assumed to be constant
% % $ (m) Shear heating
% % $  Set brittle shear heating constant below
% % $  1 = on
% % $  0 = off
% % $ (n) Shear heating constant (unitless)
% % $   Scales shear heating within the brittle realm.
% % $   Implemented in same form as used by F. Herman (02/08)
% % $   1 = Full (unscaled) brittle shear heating
% % $   0 = No brittle shear heating
% % 
% % 110. 220 2.5 800. 2500. 3300.
% % 1300. 23.0 5.3 3.2 12.0 0.01 0 0
% % 
% % $ (Input 15) Thermal model input parameters for Nepal model geometry
% % $   **NOTE** This is not used unless the geometry flag above is set to 4.
% % $   On each line, there are five values. The values are for the Indian Shield,
% % $   Sub-Himalaya, Lesser Himalaya, Greater Himalaya and Tethyan Himalaya
% % $   Line 1: Volumetric heat production (uW/m^3)
% % $   Line 2: Thermal conductivity (W/m K)
% % $   Line 3: Rock density (kg/m^3)
% % $   Line 4: Specific heat capacity (J/kg K)
% % 
% % 0.8 0.8 0.8 1.9 0.8
% % 2.75 2.75 2.75 2.75 2.75
% % 2700. 2700. 2700. 2700. 2700.
% % 1000. 1000. 1000. 1000. 1000.
% % 
% % $ (Input 16) Option to read in thermochron data and compare to predicted ages
% % $ First specify the number of data files for comparison
% % $ 0 = No data file(s) will be read in
% % $ For each file name that is specified, the file format should be as follows:
% % $ First line in file = number of samples (and lines) in rest of file.
% % $ Each line after that is for an individual sample and should contain (space
% % $   separated)
% % $ (a) Sample longitude or utm x (value in km)
% % $ (b) Sample latitude or utm y (value in km)
% % $ (c) Sample elevation
% % $ (d) Flag for type of AHe age to predict:
% % $   1=Default diffusion kinetics; 2-4=Use grain size of 20, 40 or 70 um, resp.
% % $   5-7=Use low, moderate or high eU (radiation damage) values resp.
% % $   (Schuster et al., 2006)
% % $   NOTE: These values can be modifed in the Mad_He.f90 subroutine
% % $   Comments in that subroutine further explain the differences above
% % $ (e) AHe age (Ma), negative age if non-existant
% % $ (f) AHe age error, 1s.d. (Ma), use 0 if previous value is negative
% % $ (g) AFT age (Ma), negative age if non-existant
% % $ (h) AFT age error, 1s.d. (Ma), use 0 if previous value is negative
% % $ (i) ZHe age (Ma), negative age if non-existant
% % $ (j) ZHe age error, 1s.d. (Ma), use 0 if previous value is negative
% % $ (k) ZFT age (Ma), negative age if non-existant
% % $ (l) ZFT age error, 1s.d. (Ma), use 0 if previous value is negative
% % $ (m) MAr age (Ma), negative age if non-existant
% % $ (n) MAr age error, 1s.d. (Ma), use 0 if previous value is negative
% % $ (o) Ar in K-feldspar, negative age if non-existant
% % $ (p) Ar in K-feldspar, use 0 if previous value is negative
% % $ (q) Ar in Biotite, negative age if non-existant
% % $ (r) Ar in Biotite, use 0 if previous value is negative
% % $ (s) Ar in Muscovite, negative age if non-existant
% % $ (t) Ar in Muscovite, use 0 if previous value is negative
% % $ (u) Ar in Hornblende, negative age if non-existant
% % $ (v) Ar in Hornblende, use 0 if previous value is negative
% % $ (w) Apatite U-Th / Pb, negative age if non-existant
% % $ (x) Apatite U-Th / Pb, use 0 if previous value is negative
% % $ (y) Biotite, negative age if non-existant
% % $ (z) Biotite, use 0 if previous value is negative
% % $ (a1) Ruite U-Pb, negative age if non-existant
% % $ (b1) Ruite U-Pb, use 0 if previous value is negative
% % $ (c1) Titanite U-Pb, negative age if non-existant
% % $ (d1) Titanite U-Pb, use 0 if previous value is negative
% % $ (e1) Zircon U-Pb, negative age if non-existant
% % $ (f1) Zircon U-Pb, use 0 if previous value is negative
% % $ (g1) Titanite U-Th, negative age if non-existant
% % $ (h1) Titanite U-Th, use 0 if previous value is negative
% % $ (i1) ZHe, low damage, negative age if non-existant
% % $ (j1) ZHe, low damage, use 0 if previous value is negative
% % $ (k1) ZHe, med damage, negative age if non-existant
% % $ (l1) ZHe, med damage, use 0 if previous value is negative
% % $ (m1) ZHe, high damage, negative age if non-existant
% % $ (n1) ZHe, high damage, use 0 if previous value is negative
% % $ (o1) Sample ID
% % 0
% % 
% % $ cstmt_ahe_errinc.dat
% % $ cstmt_aft.dat
% % 
% % 
% % $ NOTE: Input 17 and 18 are used for detrital age calculation.  Instructions
% % $ on use are in 00README_detrital_ages file if user needs more information.
% % 
% % $ (Input 17) Flags for which ages to output
% % $ 0 = Does not calculate or output predicted ages for this system
% % $ 1 = Calculates and outputs specified system's ages
% % $ NOTE: See Mad_He.f90 subroutine to modify the predicted AHe ages below
% % $ (a) AHe age (Default kinetics)
% % $ (b) AHe age (20 um grain size)
% % $ (c) AHe age (40 um grain size)
% % $ (d) AHe age (70 um grain size)
% % $ (e) AHe age (Low radiation damage)
% % $ (f) AHe age (Moderate radiation damage)
% % $ (g) AHe age (high radiation damage)
% % $ (h) ZHe age
% % $ (i) AFT age
% % $ (j) ZFT age, D0 = 0.001_8, energy = 208.32768_8, grain_size = 3.158_8, geometry_factor = 55.0_8
% % $ (k) Muscovite Ar/Ar age, D0 = 4.0e-8_8,  energy = 180.0_8, grain_size = 750.0_8, geometry_factor = 8.65_8
% % $ (l) Ar in K-feldspar, D0 = 5.6, Ea = 120.0
% % $ (m) Ar in Biotite, D0 = 160.0, Ea = 211.0
% % $ (n) Ar in Muscovite, D0 = 13.2, Ea = 183.0
% % $ (o) Ar in Hornblende, D0 = 14.0, Ea = 176.0
% % $ (p) Apatite U-Th / Pb, D0 = 2.0e-8_8, energy = 230.0_8, grain_size = 50.0_8, geometry_factor = 55.0_8
% % $ (q) Biotite, D0 = 2.0e-13_8, energy = 105.0_8, grain_size = 500.0_8, geometry_factor = 8.65_8
% % $ (r) Ruite U-Pb, D0 = 1.6e-10_8, energy = 243.0_8, grain_size = 250.0_8, geometry_factor = 55.0_8
% % $ (s) Titanite U-Pb, D0 = 1.1e-4_8, energy = 331.0_8, grain_size = 500.0_8, geometry_factor = 55.0_8
% % $ (t) Zircon U-Pb, D0 = 7.8e-3_8, energy = 544.0_8, grain_size = 50.0_8, geometry_factor = 55.0_8
% % $ (u) Titanite U-Th / He, D0 = 5.9e-3, energy = 188, grain_size = 250
% % $ (v) ZHe, low damage D0z = 4.6
% % $ (w) ZHe, med damage D0z = 0.3
% % $ (x) ZHe, high damage D0z = 4.6E+05
% % 
% % 
% % $(a)(b)(c)(d)(e)(f)(g)(h)(i)(j)(k)(l)(m)(n)(o)(p)(q)(r)(s)(t)(u)(v)(w)(x)
% % 1  0  0  0  0  0  0  1  1  1  1  0  0  0  0  0  0  0  0  0  0  0  0  0
% % 
% % $ (Input 18) Flag to calculate detrital age distributions for catchments
% % $ First line:  0 = no detrital calculation; 1 = detrital calcuation
% % $ Note:  If a series of CASCADE topographies were loaded in then
% % $ set Input 18 of the Pecube.in file to be 1, Pecube will output the detrital
% % $ ages of every cascade catchment at every timestep.  These files will be created
% % $ in the 'catchments' folder within your output run directory.  It will create the
% % $ 'catchments' folder if it does not exist there.  The files will be named as
% % $ 'Timestep_0001_Catchment_0001.dat' and so on for all catchments and timesteps.
% % $ All the program needs to run properly is to have the tecplot formatted cascade
% % $ output files for every timestep (eg. 'topo_tec_0001.dat') in
% % $ the 'output/Cascade' directory
% % $
% % $ Specifying user defined basins in a file causes Pecube to open the file and read it.
% % $ The file should contain lines with the following syntax:
% % $
% % $ <xpos> <ypos> <age type> <Nil or filename> <yes/Yes or no/No>
% % $ Ex: 64.8495 128.5548 1 Nil No
% % $
% % $ Where xpos is the x value of the basin outlet and ypos is the y value.  The age type is a number 1-11 with the
% % $ following coding:
% % $
% % $ 1 = Apatite Helium Age - Farley, 2000
% % $ 2 = Apatite Helium Age - Small grain size
% % $ 3 = Apatite Helium Age - Medium grain size
% % $ 4 = Apatite Helium Age - Large grain size
% % $ 5 = Apatite Helium Age - Low radiation damage
% % $ 6 = Apatite Helium Age - Medium radiation damage
% % $ 7 = Apatite Helium Age - High radiation damage
% % $ 8 = Apatite Fission Track Age
% % $ 9 = Zircon Helium Age
% % $ 10 = Zircon Fission Track Age
% % $ 11 = Muscovite Age
% % $
% % $ The next entry for each basin line can be 'Nil' or a filename.  If 'Nil' is specified, then Pecube uses
% % $ the TAPES-Grid method of finding upstream points of the basin outlet.  Then, it will write out the x, y,
% % $ and z positions along with the age data for the upstream points into the main run output directory with
% % $ the naming convention of 'Timestep_0001_Basin_X_64.8495_Y_128.5548.dat'.  Also, the PDF for each basin is
% % $ written to the folder 'pdf_data' within the run output directory with the naming convention of
% % $ 'Timestep_0001_Basin_X_64.8495_Y_128.5548_Agetype_01_pdf.dat'. If a filename is specified as the entry
% % $ on the line, then Pecube will open this file and read in each line of data with the format:
% % $
% % $ <age> <error>
% % $ Ex: 35.756469 1.07269407
% % $
% % $ Where the age and error are absolute values. The PDF for each of these basins is created in the 'pdf_data'
% % $ folder in the run output directory with the naming convention of 'Basin_X_64.8495_Y_128.5548_Agetype_01_pdf.dat'.
% % $
% % $ The final entry on each line is a 'yes' or 'no' on whether the user wants to run the Monte Carlo test for that specified
% % $ basin.  For the Monte Carlo test to run properly there must be two basins with the same outlet point (x and y positions) and
% % $ EACH must have a 'yes' to run the monte carlo routine. Also, the age types MUST be the same or the Monte Carlo test is not run.
% % $ An example of correct syntax for the Monte Carlo routine to run properly is as follows:
% % $
% % $ 64.8495 128.5548 7 Nil Yes
% % $ 64.8495 128.5548 7 datafile.txt Yes
% % $
% % $ Please note that if any or all of these criteria for running the Monte Carlo test are not met, then the program simple does not
% % $ run the comparison (skips it) and does not output anything for it.
% % $
% % $ Note: The subdirectories in the run output directory where Pecube writes most of these files will be
% % $       automatically created by Pecube if they do not exist already.
% % 
% % $pdf_tester_for_data.txt
% % 0
% % 
% % $ (Input 19) Minimum number of nodes for a catchment to be output
% % $ This is a threshold value of nodes that a catchment needs to have in order for an
% % $   output file to be written for that catchment at that timestep.
% % 
% % 100
% % 
% % $ (Input 20) Name of directory where Cascade tecplot formatted output files are located
% % $ to be read in by Pecube for PDF calculation
% % $ Note: This only matters if a '1' is selected above (Input 18) for use of Cascade catchments
% % $ The program will disregard whatever is here if '0' or a filename is specified
% % 
% % output/Cascade
% % 
% % $ (Input 21) Name of the temperature file if needed. Otherwise Nil
% % 
% % Nil
% % 
% % $ (Input 22) If e=8 in Input 12 above, filename for velocity file. Otherwise Nil
% % 
% % McQ02N4_VV5.txt


