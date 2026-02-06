clear; clc; close all;


%% =========================================================
%  CARGAR DATOS DE FUERZA EN PAREDES
%% =========================================================

directory = 'C:\Users\Noeli\Camila\3CiclosBajo';

FloorPosition = readmatrix(fullfile(directory,"floorPosition.dat"));
FloorPosition = FloorPosition(:,2);

FloorForce = readmatrix(fullfile(directory,"floorForce.dat"));
FloorForce = FloorForce(:,2);

roofForce = readmatrix(fullfile(directory,"roofForce.dat"));
roofForce = roofForce(:,2);

xForce = readmatrix(fullfile(directory,"x+WallForce.dat"));
xForce = xForce(:,1);

xForcemin = readmatrix(fullfile(directory,"x-WallForce.dat"));
xForcemin = xForcemin(:,1);

zForce = readmatrix(fullfile(directory,"z+WallForce.dat"));
zForce = zForce(:,3);

zForcemin = readmatrix(fullfile(directory,"z-WallForce.dat"));
zForcemin = zForcemin(:,3);

times = linspace(0,0.5,length(FloorForce));

%% =========================================================
%  STRAIN
%% =========================================================

L0 = 0.4;
strain = (L0 - FloorPosition) / L0;
strain_pct = 100 - 100*strain;

%% =========================================================
%  COMPRESIÓN / DESCOMPRESIÓN
%% =========================================================

dPos = diff(FloorPosition);
idx_comp   = find(dPos > 0) + 1;
idx_decomp = find(dPos < 0) + 1;

%% =========================================================
%  FUNCIÓN DE PROMEDIO POR CICLO
%% =========================================================

function [x_common, mean_comp, mean_decomp] = ...
         average_cycle(x, signal, idx_comp, idx_decomp, Npts)

if nargin < 5
    Npts = 300;
end

% --- COMPRESIÓN ---
xC = x(idx_comp);
fC = signal(idx_comp);
[xC, ord] = sort(xC);
fC = fC(ord);
[xC_u,~,ic] = unique(xC);
fC_m = accumarray(ic,fC,[],@mean);

% --- DESCOMPRESIÓN ---
xD = x(idx_decomp);
fD = signal(idx_decomp);
[xD, ord] = sort(xD);
fD = fD(ord);
[xD_u,~,id] = unique(xD);
fD_m = accumarray(id,fD,[],@mean);

% --- EJE COMÚN ---
x_common = linspace( ...
    max(min(xC_u),min(xD_u)), ...
    min(max(xC_u),max(xD_u)), ...
    Npts);

% --- INTERPOLACIÓN ---
mean_comp   = interp1(xC_u,fC_m,x_common,'linear');
mean_decomp = interp1(xD_u,fD_m,x_common,'linear');
end

%% =========================================================
%  STRESS
%% =========================================================

stress_zz = (abs(zForce)/1000)/0.16;
[strain_m, zz_comp, zz_decomp] = ...
    average_cycle(strain_pct, stress_zz, idx_comp, idx_decomp);

stress_yy = (abs(FloorForce)/1000)/0.16;
[~, yy_comp, yy_decomp] = ...
    average_cycle(strain_pct, stress_yy, idx_comp, idx_decomp);

%% =========================================================
%  ESTILO GRÁFICO (FIGURA DE REFERENCIA)
%% =========================================================

set(groot,'defaultAxesFontSize',16)
set(groot,'defaultTextFontSize',16)
set(groot,'defaultAxesLineWidth',1.2)
set(groot,'defaultAxesTickDir','out')
set(groot,'defaultAxesBox','on')

% colores EXACTOS del ejemplo
c_comp       = [0.08, 0.51, 0.51];     % verde oscuro
c_decomp     = [0.61, 0.61, 0.61];     % gris oscuro

c_comp_p     = [0.54, 0.75, 0.75];     % verde claro (puntos)
c_decomp_p   = [0.80, 0.80, 0.80];     % gris claro  (puntos)

%% =========================================================
%  STRESS vs STRAIN (ESTILO FINAL)
%% =========================================================

figure; hold on

% --- puntos (trayectorias individuales) ---
%plot(stress_zz(idx_comp), strain_pct(idx_comp),'Marker','.', 'LineStyle','none', 'Color',c_comp_p)

%plot(stress_zz(idx_decomp), strain_pct(idx_decomp),  'Marker','.', 'LineStyle','none', 'Color',c_decomp_p)

plot(stress_yy(idx_comp),strain_pct(idx_comp), ...
     'Marker','.', 'LineStyle','none', 'Color',c_comp_p)

plot( stress_yy(idx_decomp),strain_pct(idx_decomp), ...
     'Marker','.', 'LineStyle','none', 'Color',c_decomp_p)

% --- promedios ---
%plot( zz_comp, strain_m,  '--',  'Color',c_comp,   'LineWidth',4)
plot( yy_comp,strain_m,   '-', 'Color',c_comp,   'LineWidth',4)

%plot( zz_decomp,strain_m, '--',  'Color',c_decomp, 'LineWidth',4)
plot( yy_decomp,strain_m, '-', 'Color',c_decomp, 'LineWidth',4)

annotation('arrow', [0.5 0.6], [0.5 0.6], 'Color', c_comp, 'LineWidth', 3);
annotation('arrow', [0.6 0.5],  [0.88 0.8], 'Color', c_decomp, 'LineWidth', 3);

ylabel('$\mathrm{Relative\ Strain}\ (\%)$','Interpreter','latex')
xlabel('$\mathrm{Confining\ Pressure}\ (\mathrm{kPa})$','Interpreter','latex')
box(gca,'on');
hold(gca,'off');
set(gca,'FontSize',14,'TickLabelInterpreter','latex')
grid off
