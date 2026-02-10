clear all
close all

%% cummulative strain through all the comprression paths for each dataset
% changing path may change lines 26, 57, 168

addpath 'C:\Users\Came\Dropbox\PhD\Quasi-staticData'

cd = {'c31-Jan-2020_Mov_1\',...
    'c03-Feb-2020_Mov_1\',...
    'c12-Jul-2020_Mov_1 (1)\',...
    'c12-Aug-2020_Mov_2 (1)\',... 
    'c10-02_cycles\',...
    'c01-03_Data\',...
};
%%
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

diferencia = seconds(fechas - fechas(1))/60/60/24;

%% plot results
figure(10)
clf
count_i_anterior = 0;
lens=0;
for i = 1:6
    readyFile = [cd{i}, 'ReadyQuasiStaticData_', cd{i}(2:end-1), '.mat'];
    
    load(readyFile);

    smoothedPressure = movmean(ExternalPressure, 200);
    p100=find(abs(diff(smoothedPressure<100)));
    p70=find(abs(diff(smoothedPressure<70)));
    
    
    
        figure(i*10000)
        plot(ExternalPressure,Strain)
        hold on
        
        figure
        plot(ExternalPressure./max(ExternalPressure))
        hold on
        plot(Strain./max(Strain),'r')

    
    stops = round(stops(:,1));
    stops = stops(:);
    lens=lens+length(stops)
    

    cumm_strain = diff([Strain(stops(1:2:end))']);
    cumm_strain3 = diff([Strain((p100(2:2:end)))']);
    cumm_strain4 = diff([Strain((p70(2:2:end)))']);
    
    datos_it = [cumm_strain(:), cumm_strain3(:), cumm_strain4(:)];
    
    promcumm = mean(datos_it, 2);
    stdcumm  = std(datos_it, 0, 2); % '0' normalización estándar (N-1)
    
    if i==6
        figure(1)
        hold on
        
        x_values = linspace(1, numel(cumm_strain), numel(cumm_strain)) + count_i_anterior;
        
       
        errorbar(x_values, promcumm, stdcumm, ...
            'Color', colors(i,:), ...
            'LineStyle', 'none', ...
            'Marker', 's', ...
            'MarkerFaceColor', colors(i,:), ... 
            'LineWidth', 1.5, ...
            'CapSize', 10); 
        
        count_i_anterior = count_i_anterior + numel(cumm_strain);
        
    else
       
        figure(1)
        hold on
        x_values = linspace(1, numel(cumm_strain), numel(cumm_strain)) + count_i_anterior;
        
        errorbar(x_values, promcumm, stdcumm, ...
            'Color', colors(i,:), ...
            'LineStyle', 'none', ...
            'Marker', 's', ...
            'MarkerFaceColor', colors(i,:), ...
            'LineWidth', 1.5, ...
            'CapSize', 10); 
        count_i_anterior = count_i_anterior + numel(cumm_strain);
        
    end
    
    
    figure(i*100)
    plot( cumm_strain ,'ko-')
    hold on
    plot( cumm_strain3 ,'bo-')
    
    
end

figure(10)
figure(1)
xlabel('Cycle number','Interpreter','latex')
ylabel('Cummulative strain, $\Delta \epsilon$, ($\%$)','Interpreter','latex')
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
count_i_anterior = count_i_anterior + numel(cumm_strain);

figure(10)
ax2.XLim = ax1.XLim;
%% inset
i=3;
readyFile = [cd{i}, 'ReadyQuasiStaticData_', cd{i}(2:end-1), '.mat'];
load(readyFile);
figure(100),clf,h1=plot(ExternalPressure(stops(4)+250:15:stops(6)+100),Strain(stops(4)+250:15:stops(6)+100),'k-'),
hold on
h2=plot(ExternalPressure(stops(6)+100:15:stops(8)),Strain(stops(6)+100:15:stops(8)),'k--');
xlabel('External Pressure (kPa)','Interpreter','latex')
ylabel('Strain ($\%$)','Interpreter','latex')
xlabel('Confining pressure $p$ (kPa)','Interpreter','latex')
ylabel('Strain ($\%$)','Interpreter','latex')
a=gca;
a.FontSize=14;
a.TickLabelInterpreter = 'latex';
box on

xline(196.5, '-', 'LineWidth', 2, 'Color', 'k');
xline(70, '-', 'LineWidth', 2, 'Color', 'k');
xline(100, '-', 'LineWidth', 2, 'Color', 'k');

annotation('textbox',...
    [0.669565217391306 0.630289532293986 0.0660869565217381 0.0601336302895324],...
    'String','$\Delta\epsilon$',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox', [0.2 0.2 0.1 0.1], ...
    'String', char(9651), ...   % △ white up triangle
    'FontSize', 15, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'LineStyle', 'none');

annotation('textbox', [0.2 0.3 0.1 0.1], ...
    'String', char(9633), ...   % □ hollow square
    'FontSize', 15, ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center');

annotation('textbox', [0.2 0.4 0.1 0.1], ...
    'String', char(9675), ...   % ○ hollow circle
    'FontSize', 15, ...
    'LineStyle', 'none', ...
    'HorizontalAlignment', 'center');


