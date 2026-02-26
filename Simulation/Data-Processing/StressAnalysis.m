%% Stress analysis from DEM simulations (ESyS-Particle)
%
% This script post-processes stress tensor data obtained from DEM
% simulations performed with ESyS-Particle. It extracts the full stress
% tensor at cell centers, computes pressure, and analyzes the evolution
% of internal stresses during cyclic compression–decompression tests.
%
% The script:
%  - Reads cell-wise stress tensors and boundary forces.
%  - Computes pressure and stress components over time.
%  - Selects representative "sensor" cells at different heights.
%  - Filters and smooths stress signals.
%  - Compares internal stresses with external confining force.
%  - Computes cycle-averaged stress responses.
%
% Output figures are used to analyze hysteresis and stress transmission
% in granular assemblies under quasi-static loading.
%
% Author: Noelia Olivera Rodríguez
% Requirements: MATLAB R2021a or newer


clear all
close all
sym='os^dvos^dv';

%% Initialization of variables
% Temporal storage for stress tensor components, pressure, and time steps.
steps = [];
points = [];
pressuresTB = [];
sigmaxxTB = [];
sigmayyTB = [];
sigmazzTB = [];
sigmaxyTB = [];

% Simulation and physical parameters
BoxSide = 0.4;           % Box side length [m]
BoxArea = BoxSide^2;    % Cross-sectional area [m^2]
dt = 3e-7;              % DEM timestep [s]
density = 3e3;          % Particle density [kg/m^3]
rmin = 0.003;            % Minimum particle radius [m]
rmax = 0.007;            % Maximum particle radius [m]

%% Data extraction
% Directory containing stress tensor files and boundary force data.
directory   = 'C:\Users\Noeli\Documents\Quasi-static\Simulation\Data\Stress'; % path to the Stress folder

% Reaction force at the bottom plate (divided by 2 due to ESyS output convention)
FloorData = readmatrix(fullfile(directory, "floorForce.dat"));

% Cell center positions
FloorPosition = readmatrix(fullfile(directory, "floorPosition.dat"));

% List of stress tensor files
files = {dir(directory).name};
% Loop over stress tensor files


for i = 3:length(files)-2
    data = readmatrix(fullfile(directory, files{i}));

    % Extract timestep index from filename
    steps = [steps, str2double(extract(convertCharsToStrings(files{i}), digitsPattern))];

    % Stress tensor components (divided by 2 due to ESyS convention)
    sigma_xx = data(:,4);
    sigma_xy = data(:,5);
    sigma_xz = data(:,6);
    sigma_yy = data(:,7);
    sigma_yz = data(:,8);
    sigma_zz = data(:,9);

    % Local pressure definition
    pressure = -(sigma_xx + sigma_yy + sigma_zz) / 3;

    % Time series storage (sign convention adjusted)
    sigmaxxTB = [sigmaxxTB, -sigma_xx];
    sigmayyTB = [sigmayyTB, -sigma_yy];
    sigmazzTB = [sigmazzTB, -sigma_zz];
    sigmaxyTB = [sigmaxyTB, -sigma_xy];
    pressuresTB = [pressuresTB, pressure];
end


%% Cell center extraction

x = data(:,1);
y = data(:,2);
z = data(:,3);

points = [x, y, z];

%% Time ordering
% Time vector
times = steps .* dt;

% Sort data chronologically (files are alphabetically ordered)
[~, order] = sort(steps);
sigmaxxTB = sigmaxxTB(:,order);
sigmayyTB = sigmayyTB(:,order);
sigmazzTB = sigmazzTB(:,order);
sigmaxyTB = sigmaxyTB(:,order);
pressuresTB = pressuresTB(:,order);
times = sort(times);

% Confining pressure computed from the normal force at the bottom wall
ConfiningPressure = FloorData(1:length(steps),2) ./ BoxSide^2 /1000;


%% Reference plots of sensor locations

% Y-coordinates defining the top, middle and bottom planes
TOP    = 0.3639;
MIDDLE = 0.2536;
BOTTOM = 0.033083;

% Particles located at the bottom plane (Y = BOTTOM)
places = find(abs(points(:,1) - 0.3667) < 1e-2);
placesxmin = places;
newpoints  = points(places,:);
%%
% 3D scatter plot of all sensor positions
figure
scatter3(newpoints(:,1), newpoints(:,2), newpoints(:,3),10,places, 'filled');
xlabel('x', 'FontSize', 11, 'Interpreter', 'latex');
ylabel('y', 'FontSize', 11, 'Interpreter', 'latex');
zlabel('z', 'FontSize', 11, 'Interpreter', 'latex');
title('All points', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
box(gca,'on');
hold(gca,'off');

%% Comparison of Vertical Force from Internal Stress and Base
places = find(abs(points(:,2)-  0.0331) < 1e-2);
placesymin = places;

N = length(places);                
A = BoxSide^2 / N; 
forcemean = sum(sigmayyTB(places,:) * A);

figure
plot(times, forcemean, 'LineWidth',2)
hold on
plot(times, FloorData(1:length(times),2), 'LineWidth',2)
legend('Net Force in Section','Wall Force','Interpreter', 'latex')
title('Consistency check: Section vs Wall Force', 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex', 'FontSize', 16);
xlabel('Time (s)', 'FontSize', 16, 'Interpreter', 'latex')
ylabel('Fuerza (N)', 'FontSize', 16, 'Interpreter', 'latex')
%% Places Finder;
results = find( abs(points(:,1) - 0.1889) < 1e-2  & abs(points(:,2) - 0.1875) < 2.5e-1 &  abs(points(:,3) - 0.1889) < 1e-2);
%% Pressure vs force for selected sensors

% Stress signals from six manually selected sensors are extracted,
% converted to kPa, and plotted against the confining force.
% Loading and unloading branches are separated using posmax.

%% Selection of virtual stress sensors

% Six cells are manually selected to represent bottom and side regions
% of the sample and act as virtual stress sensors.

S1Bottom = 3268;
S2Bottom = 1972;
S3Bottom = 2620;
S1Side   =  2735;
S2Side   = 2771;
S3Side   =  2807;

% Sensor positions
S1BottomPos = points(S1Bottom,:);
S2BottomPos = points(S2Bottom,:);
S3BottomPos = points(S3Bottom,:);
S1SidePos   = points(S1Side,:);
S2SidePos   = points(S2Side,:);
S3SidePos   = points(S3Side,:);

% Stress signals extracted at sensor locations and converted to kPa
stress1 = sigmayyTB(S1Bottom,:) / 1000;
stress2 = sigmayyTB(S2Bottom,:) / 1000;
stress3 = sigmayyTB(S3Bottom,:) / 1000;
stress4 = sigmaxxTB(S1Side,:)   / 1000;
stress5 = sigmaxxTB(S2Side,:)   / 1000;
stress6 = sigmaxxTB(S3Side,:)   / 1000;

% Time subdivision corresponding to loading–unloading cycles
timediv = 0:0.3:1.8;

% Indices of time subdivision points in the time vector
timedivindex = round(linspace(1, length(times), 9));

% Matrix containing stress signals from all sensors
stress_matrix = [stress1; stress2; stress3; stress4; stress5; stress6];

%% Signal cleaning and smoothing

% Local outlier removal based on deviation from a moving average,
% followed by moving-average smoothing of the stress signals.
remove_outliers_local = @(x,w,thresh) ...
    (abs(x - movmean(x,w)) < thresh * movstd(x,w));

win    = 3;   % smoothing window size (number of points)
thresh = 5;   % threshold in units of standard deviation

for i = 1:6
    mask = remove_outliers_local(stress_matrix(i,:), win, thresh);
    stress_matrix(i, ~mask) = NaN;
    stress_matrix(i,:) = movmean(stress_matrix(i,:), win, 'omitnan');
end

% Store initially smoothed signals
stress_matrix_ini = stress_matrix;

% Additional smoothing
stress_matrix = zeros(size(stress_matrix_ini));
for i = 1:6
    stress_matrix(i, :) = smoothdata(stress_matrix_ini(i,:), 'movmean', 5);
end

% Reassign to individual variables
stress1 = stress_matrix(1,:);
stress2 = stress_matrix(2,:);
stress3 = stress_matrix(3,:);
stress4 = stress_matrix(4,:);
stress5 = stress_matrix(5,:);
stress6 = stress_matrix(6,:);

%% Pressure vs force: individual sensors and cycle-averaged response

colc = ["#7EA67E", "#807EA6", "#F79797", "#A84275", "#A84242", "#A87542", ...
        "#7EA67E", "#807EA6", "#F79797", "#A84275", "#A84242", "#A87542"];
cold = ["#5F8A5F", "#5E5C8A", "#D16B6B", "#7E2F5A", "#7E2F2F", "#7E5A2F",...
        "#5F8A5F", "#5E5C8A", "#D16B6B", "#7E2F5A", "#7E2F2F", "#7E5A2F"];
figure

for i = 1:size(stress_matrix,1)
    % Full compression–decompression cycles (light colors),
    % split into loading and unloading branches using timedivindex.
    st = stress_matrix(i,:);
    %Loading (4 cycles)
    plot(ConfiningPressure(timedivindex(1):timedivindex(2)), ...
    st(timedivindex(1):timedivindex(2)), ...
    'Marker','.', 'LineStyle','none', 'Color',cold(i));
    hold on
    plot(ConfiningPressure(timedivindex(3):timedivindex(4)), ...
    st(timedivindex(3):timedivindex(4)), ...
    'Marker','.', 'LineStyle','none', 'Color',cold(i));
    
    plot(ConfiningPressure(timedivindex(5):timedivindex(6)), ...
    st(timedivindex(5):timedivindex(6)), ...
    'Marker','.', 'LineStyle','none', 'Color',cold(i));
    
    plot(ConfiningPressure(timedivindex(7):timedivindex(8)), ...
    st(timedivindex(7):timedivindex(8)), ...
    'Marker','.', 'LineStyle','none', 'Color',cold(i));
    
    %Unloading (4 cycles)
    plot(ConfiningPressure(timedivindex(2):timedivindex(3)), ...
    st(timedivindex(2):timedivindex(3)), ...
    'Marker','.', 'LineStyle','none', 'Color',cold(i));
    
    plot(ConfiningPressure(timedivindex(4):timedivindex(5)), ...
    st(timedivindex(4):timedivindex(5)), ...
    'Marker','.', 'LineStyle','none', 'Color',cold(i));
    
    plot(ConfiningPressure(timedivindex(6):timedivindex(7)), ...
    st(timedivindex(6):timedivindex(7)), ...
    'Marker','.', 'LineStyle','none', 'Color',cold(i));
    
    plot(ConfiningPressure(timedivindex(8):timedivindex(9)), ...
    st(timedivindex(8):timedivindex(9)), ...
    'Marker','.', 'LineStyle','none', 'Color',cold(i));
end

% Initialization of cycle-averaged matrices
ini = zeros(6,120);
fin = zeros(6,119);

%% Cycle averaging

hold on
% Stress signals are averaged over equivalent loading and unloading
% stages across multiple compression–decompression cycles.
for i = 1:6

    tomean = [stress_matrix(i, timedivindex(1):timedivindex(2));
              stress_matrix(i, timedivindex(3):timedivindex(4));
              stress_matrix(i, timedivindex(5):timedivindex(6));
              stress_matrix(i, timedivindex(7):timedivindex(8))];

    ini(i,:) = mean(tomean, 'omitnan');

    plot(ConfiningPressure(timedivindex(1):timedivindex(2)), ...
         ini(i,:), 'LineWidth',5, 'LineStyle','-', 'Color',cold(i));

    len = min([timedivindex(3)-timedivindex(2), ...
               timedivindex(5)-timedivindex(4), ...
               timedivindex(7)-timedivindex(6)]);

    tomean = [stress_matrix(i, timedivindex(2):(timedivindex(2)+len-1));
              stress_matrix(i, timedivindex(4):(timedivindex(4)+len-1));
              stress_matrix(i, timedivindex(6):(timedivindex(6)+len-1));
              stress_matrix(i, timedivindex(8):(timedivindex(8)+len-1))];

    fin(i,:) = mean(tomean, 'omitnan');

    plot(ConfiningPressure(timedivindex(2):(timedivindex(2)+len-1)), ...
         fin(i,:), 'LineWidth',5, 'LineStyle','--', 'Color',cold(i));
end

ylabel('Internal Stress (kPa)','Interpreter','latex');
xlabel('External Force (kN)','Interpreter','latex');
legend('y=0.0993','','y=0.1213','','y=0.1434','','y=0.1654','','y=0.1875','','y=0.2095')
set(gca,'TickLabelInterpreter','latex','FontSize',14);

%% Stress–pressure relation averaged over sensors

c_comp = [0.08,0.51,0.51];   % green
c_decomp = [0.61,0.61,0.61];    % grey

c_comp_p =  [0.54 0.75 0.75];   %green
c_decomp_p = [0.8 0.8 0.8];    % grey

figure


for i = 1:size(stress_matrix,1)
    st = stress_matrix(i,:);
    %Loading (4 cycles)
    plot(ConfiningPressure(timedivindex(1):timedivindex(2)), ...
    st(timedivindex(1):timedivindex(2)), ...
    'Marker','.', 'LineStyle','none', 'Color',c_comp_p);
    hold on
    plot(ConfiningPressure(timedivindex(3):timedivindex(4)), ...
    st(timedivindex(3):timedivindex(4)), ...
    'Marker','.', 'LineStyle','none', 'Color',c_comp_p);
    
    plot(ConfiningPressure(timedivindex(5):timedivindex(6)), ...
    st(timedivindex(5):timedivindex(6)), ...
    'Marker','.', 'LineStyle','none', 'Color',c_comp_p);
    
    plot(ConfiningPressure(timedivindex(7):timedivindex(8)), ...
    st(timedivindex(7):timedivindex(8)), ...
    'Marker','.', 'LineStyle','none', 'Color',c_comp_p);
    
    %Unloading (4 cycles)
    plot(ConfiningPressure(timedivindex(2):timedivindex(3)), ...
    st(timedivindex(2):timedivindex(3)), ...
    'Marker','.', 'LineStyle','none', 'Color',c_decomp_p );
    
    plot(ConfiningPressure(timedivindex(4):timedivindex(5)), ...
    st(timedivindex(4):timedivindex(5)), ...
    'Marker','.', 'LineStyle','none', 'Color',c_decomp_p );
    
    plot(ConfiningPressure(timedivindex(6):timedivindex(7)), ...
    st(timedivindex(6):timedivindex(7)), ...
    'Marker','.', 'LineStyle','none', 'Color',c_decomp_p );
    
    plot(ConfiningPressure(timedivindex(8):timedivindex(9)), ...
    st(timedivindex(8):timedivindex(9)), ...
    'Marker','.', 'LineStyle','none', 'Color',c_decomp_p );
end

plot(ConfiningPressure(timedivindex(1):timedivindex(2)), ...
     mean(ini(1:3,:), 'omitnan'), 'LineWidth',5, 'LineStyle','-','Color',c_comp); hold on
plot(ConfiningPressure(timedivindex(1):timedivindex(2)), ...
     mean(ini(4:6,:), 'omitnan'), 'LineWidth',5, 'LineStyle','--','Color',c_comp);

plot(ConfiningPressure(timedivindex(2):(timedivindex(2)+len-1)), ...
     mean(fin(1:3,:), 'omitnan'), 'LineWidth',5, 'LineStyle','-','Color',c_decomp);
plot(ConfiningPressure(timedivindex(2):(timedivindex(2)+len-1)), ...
     mean(fin(4:6,:), 'omitnan'), 'LineWidth',5, 'LineStyle','--','Color',c_decomp);

ylabel('Internal Stress (kPa)','Interpreter','latex');
xlabel('Confining Pressure (kPa)','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',14);
box(gca,'on');
hold(gca,'off');
annotation('arrow', [0.63 0.7367], [0.27 0.29], 'Color', c_comp, 'LineWidth', 3,'LineStyle','--');
annotation('arrow', [0.6400 0.7527], [0.44 0.51], 'Color', c_comp, 'LineWidth', 3);
annotation('arrow', [0.89 0.7867], [0.71,0.6525], 'Color', c_decomp, 'LineWidth', 3);
annotation('arrow', [0.8707 0.7690], [ 0.4745 0.4342], 'Color', c_decomp, 'LineWidth', 3,'LineStyle','--');

text(0.85,0.15,'$\sigma_{yy}$','Units','normalized','Interpreter','latex','FontSize',16)
text(0.85,0.80,'$\sigma_{xx}$','Units','normalized','Interpreter','latex','FontSize',16)

%% Consistency check: internal stress vs bottom wall force
ForceY = FloorData(:,2);

figure
plot(times, ForceY(1:length(times)), 'LineWidth',5);
xlabel('Time (s)','Interpreter','latex');
ylabel('Internal Stress (kPa)','Interpreter','latex');
set(gca,'TickLabelInterpreter','latex','FontSize',14);
box(gca,'on');
hold(gca,'off');
