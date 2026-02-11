clear all
close all

%% estimate the dissipated energy in every odeometric cycle for all the present datasets
% path changes may change lines 27, 63, 125, 183

addpath 'C:\Users\Camila Sedofeito\Desktop\Dropbox\PhD\Quasi-staticData'

cd = {'c31-Jan-2020_Mov_1\',...
    'c03-Feb-2020_Mov_1\',...
    'c12-Jul-2020_Mov_1 (1)\',...
    'c12-Aug-2020_Mov_2 (1)\',...
    'c10-02_cycles\',...
    'c01-03_Data\',...
    };

cmap = linspace(0.1, 0.8, size(cd,2));   % 0.7 = gris claro, 0 = negro
colors = repmat(cmap(:), 1, 3);  
%%
months = containers.Map({'Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', ...
    'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'}, 1:12);

fechas = NaT(1, length(cd));
folderNameCleans = strings(1, length(cd));

for i = 1:length(cd)
    [~, folderName] = fileparts(cd{i}(1:end-1));
    
    tokens = regexp(folderName, 'c?(\d+)-([A-Za-z]+)-(\d{4})', 'tokens');
    if ~isempty(tokens)
        parts = tokens{1};
        day = str2double(parts{1});
        if isKey(months, parts{2})
            month = months(parts{2});
        else
            warning('Mes no reconocido: %s', parts{2});
            month = 1; % fallback a enero
        end
        year = str2double(parts{3});
        fechas(i) = datetime(year, month, day);
    else
        year = 2021;
        fechas(i) = datetime(year, month, day);
        
    end
    
    
end

% Calcular diferencias de días entre fechas
diferencia = seconds(fechas - fechas(1))/60/60/24;
%%
count_i_anterior = 0;
lens=0;


Cycle_global = [];
Energy_global = [];
cycle_offset = 0;


for i = 1:6
 readyFile = [cd{i}, 'ReadyQuasiStaticData_', cd{i}(2:end-1), '.mat'];
 
 load(readyFile);
 
 smoothedPressure = movmean(ExternalPressure, 200);
 
 stops = round(stops(:,1));
 stops = stops(:);
 lens=lens+length(stops)/2

 nCycles = floor(length(stops)/2);
E = zeros(nCycles,1);

for c = 1:nCycles

    if c == 1
        % primer ciclo: desde el inicio
        i_comp_start = 1;
    else
        i_comp_start = stops(2*c-2);
    end

    i_comp_end  = stops(2*c-1);
    i_rel_start = stops(2*c-1);
    i_rel_end   = stops(2*c);

    % --- Compress
    eps_c = Strain(i_comp_start:i_comp_end);
    sig_c = ExternalPressure(i_comp_start:i_comp_end);

    % --- Relax
    eps_r = Strain(i_rel_start:i_rel_end);
    sig_r = ExternalPressure(i_rel_start:i_rel_end);

    [eps_c, ic] = sort(eps_c); sig_c = sig_c(ic);
    [eps_r, ir] = sort(eps_r); sig_r = sig_r(ir);

    Wc = trapz(eps_c, sig_c);
    Wr = trapz(eps_r, sig_r);

    E(c) = Wc - Wr;
end


nCycles = numel(E);

Cycle_global  = [Cycle_global;  (cycle_offset+1:cycle_offset+nCycles)'];
Energy_global = [Energy_global; E(:)];

cycle_offset = cycle_offset + nCycles;

end


colors = repmat(cmap(:), 1, 3); 
cmap = linspace(0.1, 0.8, size(cd,2));   % 0.7 = gris claro, 0 = negro
figure(2)
clf
hold on

lens = 0; 
for i = 1:6
    readyFile = [cd{i}, 'ReadyQuasiStaticData_', cd{i}(2:end-1), '.mat'];
    load(readyFile);
  
    stops = round(stops(:,1));
    stops = stops(:);

    num_ciclos_actual = round(length(stops)/2-1);
    
    current_indices = lens +1: (lens + num_ciclos_actual);
    
    % Plot
    plot(Cycle_global(current_indices), ...
         Energy_global(current_indices), ...
         'Color', colors(i,:), ...
         'LineStyle', 'none', ...
         'Marker', 's', ...
         'MarkerFaceColor', 'none', ... 
         'LineWidth', 1.5);
    
    hold on;
    
    lens = lens+round(length(stops)/2-1);
end
xlabel('Cycle number','Interpreter','latex')
ylabel('Dissipated energy  (kPa)','Interpreter','latex')
box on
ylim([0 200])
a=gca;
a.FontSize=16;
a.TickLabelInterpreter = 'latex';
box on

ax1 = gca;
pos = ax1.Position;

ax2 = axes('Position', pos, ...
    'XAxisLocation', 'top', ...
    'YAxisLocation', 'right', ...
    'Color', 'none', ...
    'XColor', 'k', ...
    'YColor', 'none', ...
    'FontSize', 16);
set(ax2, 'XTick', [1 9/2-1 18/2-1 29/2-1 42/2-1 57/2-1]);
set(ax2,'XTickLabel', {'0','3','163','194','376','395'});
xlabel(ax2, 'Day','Interpreter','latex')  % Cambiá esto si querés otra cosa arriba
a=gca;
a.FontSize=16;
a.TickLabelInterpreter = 'latex';
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';
  
    
figure(2)
ax2.XLim = ax1.XLim;

%%
i=3;
readyFile = [cd{i}, 'ReadyQuasiStaticData_', cd{i}(2:end-1), '.mat'];
load(readyFile);
figure(100),clf,h1=plot(ExternalPressure(stops(4)+250:15:stops(6)+100),Strain(stops(4)+250:15:stops(6)+100),'k--','LineWidth',1.5),
hold on
xlabel('External Pressure (kPa)','Interpreter','latex')
ylabel('Strain ($\%$)','Interpreter','latex')
xlabel('Confining pressure $p$ (kPa)','Interpreter','latex')
ylabel('Strain ($\%$)','Interpreter','latex')
a=gca;
a.FontSize=16;
a.TickLabelInterpreter = 'latex';
box on
dataX = ExternalPressure(stops(4)+250:15:stops(6)+100);
dataY = Strain(stops(4)+250:15:stops(6)+100);

hold on
hFill = fill(dataX, dataY, [0 0 0], 'EdgeColor', 'none', 'FaceAlpha', 0.1);

h1 = plot(dataX, dataY, 'k-', 'LineWidth', 2);
a=gca;
a.FontSize=16;
