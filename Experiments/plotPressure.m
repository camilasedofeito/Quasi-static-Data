date={'31-Jan-2020_Mov_1',...
    '03-Feb-2020_Mov_1',...    
    '19-Feb-2020_Mov_1',...
    '12-Jul-2020_Mov_1 (1)',...
    '12-Aug-2020_Mov_1_1 (1)',...
    '10-02_cycles',...
    '01-03_Data'};

DataName='ReadyQuasiStaticData';

for i=1:length(date)-2
    load([DataName,'_',date{i},'.mat']);
    
    % External pressure 
    load(['dt_',date{i}(1:11),'.mat'])
    time=0:dt:(length(ExternalPressure)-1).*dt;
    figure(i),plot(time,ExternalPressure)
    xlabel('Tiempo (s)')
    ylabel('Stress (kPa)')
end


for i=length(date)-1:length(date)
    load([DataName,'_',date{i},'.mat']);
    
    % External pressure 
    load(['dt_',date{i}(1:5),'.mat'])
    time=0:dt:(length(ExternalPressure)-1).*dt;
    figure(i),plot(time,ExternalPressure)
    xlabel('Tiempo (s)')
    ylabel('Stress (kPa)')
end
