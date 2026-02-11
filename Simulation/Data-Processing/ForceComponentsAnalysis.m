%%  Analysis of force components vs confining pressure
%  Compression / decompression cycles
%
%  This script:
%  - Loads averaged and maximum force components by region
%  - Computes confining pressure from floor force
%  - Identifies compression and decompression branches
%  - Computes cycle-averaged force–pressure curves
%  - Plots mean and maximum force components for top, middle
%    and bottom regions

% Author: Noelia Olivera Rodríguez
% Requirements: MATLAB R2021a or newer

clear; clc; close all;

%% 1. Data loading
% Table containing force statistics per timestep and region
scriptDir = fileparts(mfilename('fullpath'));
directory   = fullfile(scriptDir, '..', 'Data\Forces\xyz_Fuerzas_componentes_entrecentroyesquina_corrected.csv');
region = 'Medium'; %Center, Corner, Medium
T = readtable(directory);

% Extract relevant columns
timestep   = T.timestep;
zona       = string(T.zona);

Fy_mean    = T.Fy_mean;
Fx_mean    = T.Fx_mean;
Fz_mean    = T.Fz_mean;

Fx_max     = T.Fx_max;
Fy_max     = T.Fy_max;
Fz_max     = T.Fz_max;

% If paths or variables change, modify only this section


%% Floor force and confining pressure

% Load floor position and force (second column only)
FloorPos   = readmatrix('C:\Users\Noeli\Camila\3CiclosBajo\floorPosition.dat');
FloorForce = readmatrix('C:\Users\Noeli\Camila\3CiclosBajo\floorForce.dat');

FloorPos   = FloorPos(:,2);
FloorForce = FloorForce(:,2);

% Confining pressure [kPa]
% Force converted to kN and divided by floor area (0.16 m^2)
pressure = (FloorForce/1000) / 0.16;


%% Compression and descompression identification
% Position increment
dPos = diff(FloorPos);

is_comp   = false(size(FloorPos));
is_decomp = false(size(FloorPos));

% Compression: increasing floor position
% Decompression: decreasing floor position
is_comp(2:end)   = dPos > 0;
is_decomp(2:end) = dPos < 0;


%% Cycle averaging function
function [x_common, y_c, y_d] = average_cycle(x, y, comp, decomp)
% Averages force–pressure cycles for compression and decompression
% separately and interpolates them onto a common pressure axis.

x = x(:);
y = y(:);

% --- Compression branch ---
xc = x(comp);
yc = y(comp);

% --- Decompression branch ---
xd = x(decomp);
yd = y(decomp);

% Safety check
if isempty(xc) || isempty(xd)
    x_common = [];
    y_c = [];
    y_d = [];
    return
end

% Sort by pressure
[xc, ic] = sort(xc); yc = yc(ic);
[xd, id] = sort(xd); yd = yd(id);

% Average repeated pressure values
[xc_u, ~, icu] = unique(xc);
yc_m = accumarray(icu, yc, [], @mean);

[xd_u, ~, idu] = unique(xd);
yd_m = accumarray(idu, yd, [], @mean);

% Common pressure interval
xmin = max(min(xc_u), min(xd_u));
xmax = min(max(xc_u), max(xd_u));

if xmin >= xmax
    x_common = [];
    y_c = [];
    y_d = [];
    return
end

% Interpolation grid
x_common = linspace(xmin, xmax, 200);

y_c = interp1(xc_u, yc_m, x_common);
y_d = interp1(xd_u, yd_m, x_common);
end


%% Plotting parameters
zonas = ["top","middle","bottom"];

% Colors
c_comp = [0.08,0.51,0.51];   % green
c_decomp = [0.61,0.61,0.61];    % grey

c_comp_p =  [0.54 0.75 0.75];   %green
c_decomp_p = [0.8 0.8 0.8];    % grey


%% Mean Fy

figure
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for i = 1:3
    nexttile; hold on

    mask = zona == zonas(i);

    Fy = Fy_mean(mask);
    P  = pressure(1:length(Fy));

    comp   = is_comp(1:length(Fy));
    decomp = is_decomp(1:length(Fy));

    [p, Fy_c, Fy_d] = average_cycle(P, Fy, comp, decomp);

    Fy_c = smoothdata(Fy_c, 'sgolay', 21);
    Fy_d = smoothdata(Fy_d, 'sgolay', 21);

    plot(P(comp),   Fy(comp),   '.', 'MarkerSize',6,'Color',c_comp_p)
    plot(P(decomp), Fy(decomp), '.', 'MarkerSize',6,'Color',c_decomp_p)
    plot(p, Fy_c, '-',  'Color', c_comp,   'LineWidth',4)
    plot(p, Fy_d, '--', 'Color', c_decomp, 'LineWidth',4)

    title(upper(zonas(i)), 'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',20)
    box(gca,'on');
    hold(gca,'off');
end
set(gcf,'Position',[100 100 820 633])
xlabel(t,'Confining pressure (kPa)', 'Interpreter','latex','FontSize',20)
ylabel(t,'$F_x^{mean}$', 'Interpreter','latex','FontSize',20)
exportgraphics(gcf,'..\Figures\'+ string(region) +'\meanfx.png','Resolution',300)

%% Mean Fx
figure
t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for i = 1:3
    nexttile; hold on

    mask = zona == zonas(i);
    Fpar = Fx_mean(mask);

    P  = pressure(1:length(Fpar));
    comp   = is_comp(1:length(Fpar));
    decomp = is_decomp(1:length(Fpar));

    [p, Fc, Fd] = average_cycle(P, Fpar, comp, decomp);

    Fc = smoothdata(Fc, 'sgolay', 21);
    Fd = smoothdata(Fd, 'sgolay', 21);

    plot(P(comp),Fpar(comp),'.','MarkerSize',6,'Color',c_comp_p)
    plot(P(decomp),Fpar(decomp),'.','MarkerSize',6,'Color',c_decomp_p)
    plot(p,Fc,'-','Color',c_comp,'LineWidth',4)
    plot(p,Fd,'--','Color',c_decomp,'LineWidth',4)

    title(upper(zonas(i)),'Interpreter','latex')
    set(gca,'TickLabelInterpreter','latex','FontSize',20)
    box(gca,'on');
    hold(gca,'off');
end
    set(gcf,'Position',[100 100 820 633])
    xlabel(t,'Confining pressure (kPa)','Interpreter','latex','FontSize',20)
    ylabel(t,'$F_y^{mean}$','Interpreter','latex','FontSize',20)
    exportgraphics(gcf,'..\Figures\'+ string(region) +'\meanfy.png','Resolution',300)
%% Maximum force components

ForceMax = {Fx_max, Fy_max, Fz_max};
labels   = {'$F_y^{max}$','$F_x^{max}$'};

for k = 1:2
    figure
    t = tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    for i = 1:3
        nexttile; hold on

        mask = zona == zonas(i);
        Fpar = ForceMax{k}(mask);

        P  = pressure(1:length(Fpar));
        comp   = is_comp(1:length(Fpar));
        decomp = is_decomp(1:length(Fpar));

        [p, Fc, Fd] = average_cycle(P, Fpar, comp, decomp);

        Fc = smoothdata(Fc, 'sgolay', 21);
        Fd = smoothdata(Fd, 'sgolay', 21);

        plot(P(comp),Fpar(comp),'.','MarkerSize',6,'Color',c_comp_p)
        plot(P(decomp),Fpar(decomp),'.','MarkerSize',6,'Color',c_decomp_p)
        plot(p,Fc,'-','Color',c_comp,'LineWidth',4)
        plot(p,Fd,'--','Color',c_decomp,'LineWidth',4)

        title(upper(zonas(i)),'Interpreter','latex')
        set(gca,'TickLabelInterpreter','latex','FontSize',20)
        box(gca,'on');
        hold(gca,'off');
    end
    set(gcf,'Position',[100 100 820 633])
    xlabel(t,'Confining Pressure (kPa)','Interpreter','latex','FontSize',20)
    ylabel(t,labels{k},'Interpreter','latex','FontSize',20)
    exportgraphics(gcf,'..\Figures\'+ string(region) +'\'+string(labels{k}) + ".png",'Resolution',300)
end
