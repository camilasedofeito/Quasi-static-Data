clear all
close all

%% load data and plot only once to get indices
figname='Glass beads';
load([figname 'Glob'])

ContPresure=nanmean(DataTotGlob(:,4:6),2)';
Npick=size(Indices,1)-1

% Npick=12
% figure
% plot(InputPresureGlob)
% X=ginput(Npick);
% Indices=round(X);
Indup=[];
for kk=3:2:Npick
    Indup=[Indup Indices(kk,1):Indices(kk+1,1)];
end
Inddwn=[];
for kk=2:2:Npick
    Inddwn=[Inddwn Indices(kk,1):Indices(kk+1,1)];
end
%% split data up
r=0.037/2;
PresureBotup=nanmean(DataTotGlob(Indup,4:6)/1000,2)';
Presuresidup=nanmean(DataTotGlob(Indup,1:3)/1000,2)';
PresureBotsigup=std(DataTotGlob(Indup,4:6)/1000,[],2)';
Presuresidsigup=std(DataTotGlob(Indup,1:3)/1000,[],2)';
InputPresup=InputPresureGlob(Indup).*(pi*r.^2)/(0.5^2).*1e5/1000;

Indup=[];
for kk=7:2:Npick
    Indup=[Indup Indices(kk,1):Indices(kk+1,1)];
end
%% velocity, pressures
InputPres=InputPresureGlob(Indup(1):end).*(pi*r.^2)/(0.5^2).*1e5/1000;
VeldltGlobup=VeldltGlob(Indup(1):end);

distk=3.5;
Pasok=1;
Pmax=max(InputPresureGlob).*(pi*r.^2)/(0.5^2).*1e5/1000-5;
Pmaxglob=170;
Pintmax=max(max(DataTotGlob/1000));
veckk=5:Pasok:Pmax;
for kk=1:length(veckk)
    InputPresResup(kk)=nanmean(InputPresup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
    PresureBotResup(kk)=nanmean(PresureBotup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));    
    PresuresidResup(kk)=nanmean(Presuresidup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
    InputPresRes(kk)=nanmean(InputPres(InputPres>=(veckk(kk)-distk) & InputPres<=(veckk(kk)+distk)));
    VeldltGlobResup(kk)=nanmean(VeldltGlobup(InputPres>=(veckk(kk)-distk) & InputPres<=(veckk(kk)+distk)));
    VeldltGlobResstdup(kk)=std(VeldltGlobup(InputPres>=(veckk(kk)-distk) & InputPres<=(veckk(kk)+distk)));
    PresureBotsigResup(kk)=mean(PresureBotsigup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
    PresuresidsigResup(kk)=mean(Presuresidsigup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
end

%%

r=0.037/2;
PresureBotup=mean(DataTotGlob(Indup,4:6)/1000,2)';
PresureBotsigup=std(DataTotGlob(Indup,4:6)/1000,[],2)';
Presuresidup=mean(DataTotGlob(Indup,1:3)/1000,2)';
Presuresidsigup=std(DataTotGlob(Indup,1:3)/1000,[],2)';
InputPresup=InputPresureGlob(Indup).*(pi*r.^2)/(0.5^2).*1e5/1000;
VeldltGlobup=VeldltGlob(Indup);
PresureBotdwn=mean(DataTotGlob(Inddwn,4:6)/1000,2)';
PresureBotsigdwn=std(DataTotGlob(Inddwn,4:6)/1000,[],2)';
Presuresiddwn=mean(DataTotGlob(Inddwn,1:3)/1000,2)';
Presuresidsigdwn=std(DataTotGlob(Inddwn,1:3)/1000,[],2)';
InputPresdwn=InputPresureGlob(Inddwn).*(pi*r.^2)/(0.5^2).*1e5/1000;
VeldltGlobdwn=VeldltGlob(Inddwn);


InputPres=InputPresureGlob(Indup(1):end).*(pi*r.^2)/(0.5^2).*1e5/1000;

distk=3.5;
Pasok=1;
Pmax=max(InputPresureGlob).*(pi*r.^2)/(0.5^2).*1e5/1000-5;
Pintmax=max(max(DataTotGlob/1000));
veckk=1:Pasok:Pmax;
for kk=1:length(veckk)
    InputPresResup(kk)=mean(InputPresup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
    PresureBotResup(kk)=mean(PresureBotup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));    
    PresureBotsigResup(kk)=mean(PresureBotsigup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
    PresuresidResup(kk)=mean(Presuresidup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
    PresuresidsigResup(kk)=mean(Presuresidsigup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
    VeldltGlobResup(kk)=mean(VeldltGlobup(InputPresup>=(veckk(kk)-distk) & InputPresup<=(veckk(kk)+distk)));
    
    InputPresResdwn(kk)=mean(InputPresdwn(InputPresdwn>=(veckk(kk)-distk) & InputPresdwn<=(veckk(kk)+distk)));
    PresureBotResdwn(kk)=mean(PresureBotdwn(InputPresdwn>=(veckk(kk)-distk) & InputPresdwn<=(veckk(kk)+distk)));    
    PresureBotsigResdwn(kk)=mean(PresureBotsigdwn(InputPresdwn>=(veckk(kk)-distk) & InputPresdwn<=(veckk(kk)+distk)));
    PresuresidResdwn(kk)=mean(Presuresiddwn(InputPresdwn>=(veckk(kk)-distk) & InputPresdwn<=(veckk(kk)+distk)));
    PresuresidsigResdwn(kk)=mean(Presuresidsigdwn(InputPresdwn>=(veckk(kk)-distk) & InputPresdwn<=(veckk(kk)+distk)));
    VeldltGlobResdwn(kk)=mean(VeldltGlobdwn(InputPresdwn>=(veckk(kk)-distk) & InputPresdwn<=(veckk(kk)+distk)));
end
%% data prep


kk=1
Pmin=15
Veckk=Indices(kk,1):Indices(kk+1,1);
indfix(kk)=find(InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000>Pmin,1);
Fixstrain=(DisplacementGlob(Veckk(indfix(kk)))-DisplacementGlob(Veckk(1)))*1E-3/0.5*100;
Pressure4StrainTotal=[];
StrainTotal=[];
VltTotal=[];
for kk=1:2:Npick
    if kk==13
        Veckk=Indices(kk,1):Indices(kk+1,1);
        indfix(kk)=length(Veckk);
        diffstrain=(DisplacementGlob(Veckk(indfix(kk))))*1E-3/0.5*100-Fixstrain;
    else
        Veckk=Indices(kk,1):Indices(kk+1,1);
        indfix(kk)=find(InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000>Pmin,1);
        diffstrain=(DisplacementGlob(Veckk(indfix(kk))))*1E-3/0.5*100-Fixstrain;
        %figure(Nfig)
        %hold on
        %       plot(DisplacementGlob(Veckk)*1E-3/0.5*100-diffstrain,InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000,'k')
    end
    if kk>5
        Veckk=(Indices(kk,1)-2):Indices(kk+1,1);
    else
        Veckk=Indices(kk,1):Indices(kk+1,1);
    end  
    %plot(InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000,DisplacementGlob(Veckk)*1E-3/0.5*100-diffstrain,'.','Color',[0.6,0.6,0.6])
    Pressure4Strain=InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000;
    Strain4Av=DisplacementGlob(Veckk)*1E-3/0.5*100-diffstrain;
    Vltloc=VeldltGlob(Veckk);
    veckk=1:Pasok:Pmax;
    for jj=1:length(veckk)
        InputPresReStrain(jj,(kk+1)/2)=nanmean(Pressure4Strain(Pressure4Strain>=(veckk(jj)-distk) & Pressure4Strain<=(veckk(jj)+distk)));
        StrainRes(jj,(kk+1)/2)=nanmean(Strain4Av(Pressure4Strain>=(veckk(jj)-distk) & Pressure4Strain<=(veckk(jj)+distk)));
        Vltpercycle(jj,(kk+1)/2)=nanmean(Vltloc(Pressure4Strain>=(veckk(jj)-distk) & Pressure4Strain<=(veckk(jj)+distk)));
    end
    Pressure4StrainTotal=[Pressure4StrainTotal Pressure4Strain];
    StrainTotal=[StrainTotal Strain4Av];
    VltTotal=[VltTotal Vltloc];
end
NP50=find(InputPresResup>50,1)
[coefssid,errssid] = errormc(InputPresResup(NP50:end),PresuresidResup(NP50:end),1);
[coefsbot,errsbot] = errormc(InputPresResup(NP50:end),PresureBotResup(NP50:end),1);
NP50s=find(InputPresReStrain>50,1)
[coefstrain,errsstrain] = errormc(nanmean(InputPresReStrain(NP50s:end-2,:),2),nanmean(StrainRes(NP50s:end-2,:),2),1);

kk=1
Pmin=15
Veckk=Indices(kk,1):Indices(kk+1,1);
Fixstrainup=(DisplacementGlob(Veckk(indfix(kk)))-DisplacementGlob(Veckk(1)))*1E-3/0.5*100;
Fixstraindwn=polyval(coefstrain,(Pmax-12));
Pressure4StrainTotalAll=[];
StrainTotalAll=[];
VltTotalAll=[];
for kk=1:1:Npick
    if kk==13
        Veckk=Indices(kk,1):Indices(kk+1,1);
        indfix(kk)=length(Veckk);
        diffstrain=(DisplacementGlob(Veckk(indfix(kk))))*1E-3/0.5*100-Fixstrain;
    else
        Veckk=Indices(kk,1):Indices(kk+1,1); 
        if round((kk+1)/2)==(kk+1)/2
            indfix(kk)=find(InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000>Pmin,1);       
            diffstraini=(DisplacementGlob(Veckk(indfix(kk))))*1E-3/0.5*100-Fixstrainup
            indfix(kk)=find(InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000>(Pmax-12),1);       
            diffstrainp=(DisplacementGlob(Veckk(indfix(kk))))*1E-3/0.5*100-Fixstraindwn
        end
        %figure(Nfig)
        %hold on
        %       plot(DisplacementGlob(Veckk)*1E-3/0.5*100-diffstrain,InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000,'k')
    end
    if kk>5
        Veckk=(Indices(kk,1)-2):Indices(kk+1,1);
    else
        Veckk=Indices(kk,1):Indices(kk+1,1);
    end  
    %plot(InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000,DisplacementGlob(Veckk)*1E-3/0.5*100-diffstrain,'.','Color',[0.6,0.6,0.6])
    Pressure4Strain=InputPresureGlob(Veckk).*(pi*r.^2)/(0.5^2).*1e5/1000;
    if round((kk+1)/2)==(kk+1)/2
        Strain4Av=DisplacementGlob(Veckk)*1E-3/0.5*100-diffstraini;
    else
        Strain4Av=DisplacementGlob(Veckk)*1E-3/0.5*100-diffstrainp;
    end
    Vltloc=VeldltGlob(Veckk);
    veckk=1:Pasok:Pmax;
    for jj=1:length(veckk)
        InputPresReStrainAll(jj,kk)=nanmean(Pressure4Strain(Pressure4Strain>=(veckk(jj)-distk) & Pressure4Strain<=(veckk(jj)+distk)));
        StrainResAll(jj,kk)=nanmean(Strain4Av(Pressure4Strain>=(veckk(jj)-distk) & Pressure4Strain<=(veckk(jj)+distk)));
        VltpercycleAll(jj,kk)=nanmean(Vltloc(Pressure4Strain>=(veckk(jj)-distk) & Pressure4Strain<=(veckk(jj)+distk)));
    end
    Pressure4StrainTotalAll=[Pressure4StrainTotalAll Pressure4Strain];
    StrainTotalAll=[StrainTotalAll Strain4Av];
    VltTotalAll=[VltTotalAll Vltloc];
%     Nfig=42;
%     figure(Nfig)
%     close
%     figure(Nfig)
%     W=4;
%     h=plot(Pressure4StrainTotalAll,StrainTotalAll,'.','Color',[0.6,0.6,0.6]);
%     pause
end
%% plot todo
strain_raw = DisplacementGlob * 1e-3 / 0.5 * 100;

pressure_raw_cycle = InputPresureGlob* (pi*r.^2)/(0.5^2).*1e5/1000;
strain_raw_cycle   = strain_raw;
figure
plot(pressure_raw_cycle,strain_raw_cycle,'o')

ii=1;
strain_ini=strain_raw_cycle(Indices(ii,1));
figure

h1=plot(pressure_raw_cycle(Indices(ii,1):end),strain_raw_cycle(Indices(ii,1):end)-strain_ini,'s','Color','#b8b8b8','MarkerFaceColor','#b8b8b8');
hold on
slopes = zeros(length(Indices),1);
for ii=1:2:length(Indices)-1
    if ii==length(Indices)-1
        h2=plot(pressure_raw_cycle(Indices(ii,1):Indices(ii+1,1)),strain_raw_cycle(Indices(ii):Indices(ii+1))-strain_ini,'s','Color','#54a1a1','MarkerFaceColor','#54a1a1');
        i1 = Indices(ii,1);
        i2 = Indices(ii+1,1);
        y= pressure_raw_cycle(i1:i2);
        x  = strain_raw_cycle(i1:i2)./100;
        
        mask = y > 50;
        x = x(mask);
        y = y(mask);
        
        p = polyfit(x,y,  1);
        slopes(ii) = p(1);       
        
    else
        plot(pressure_raw_cycle(Indices(ii,1):Indices(ii+1,1)),strain_raw_cycle(Indices(ii):Indices(ii+1))-strain_ini,'s','Color','#54a1a1','MarkerFaceColor','#54a1a1')
        i1 = Indices(ii,1);
        i2 = Indices(ii+1,1);
        y = pressure_raw_cycle(i1:i2);
        x= strain_raw_cycle(i1:i2)./100;
        
        mask = y > 50;
        x = x(mask);
        y = y(mask); 
         
        p = polyfit(x,y,  1);
        slopes(ii) = p(1); 
    end
end
ylim([0 3.5])
%legend([h2 h1], {'Compression','Relaxation'})
xlabel('Confining pressure $p$ (kPa)','Interpreter','latex')
ylabel('Strain ($\%$)','Interpreter','latex')
a=gca;
a.FontSize=16;
a.TickLabelInterpreter = 'latex';
xlim([0 165])
xticks([0 50 100 150])

pbaspect([1.8 3 1]);

ax_in = axes('Position',[0.5 0.18 0.30 0.25]);  % [x y width height] in figure-normalized units
box on; hold on;
plot(1:length(slopes(slopes>0)), slopes(slopes>0)./1e3,'s','Color','#54a1a1','MarkerFaceColor','#54a1a1');
xlabel('Cycle','Interpreter','latex','FontSize',15)
ylabel('M (MPa)','Interpreter','latex','FontSize',15)
set(ax_in,'FontSize',10,'TickLabelInterpreter','latex')
xticks([1 2 3 4 5])
ylim([5 20])
xlim([1 5])

a=gca;
a.FontSize=14;
a.TickLabelInterpreter = 'latex';
annotation('textbox',...
    [0.782609195402299 0.50958904109589 0.0464137931034483 0.0356164383561643],...
    'String','$\epsilon_1^r$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');

annotation('textbox',...
[0.778298850574713 0.604109589041096 0.0464137931034483 0.0356164383561643],...
    'String','$\epsilon_2^r$',...
    'Interpreter','latex',...
    'FontSize',15,...
    'FitBoxToText','off',...
    'EdgeColor','none');


%saveas(gcf, [figname,'_strainstress.fig'],'fig');
%% plot sacando continua
Nfig=51

figure(Nfig)
hold off
W=4
h=plot(Pressure4StrainTotalAll,StrainTotalAll,'s','MarkerSize',1,'Color',[0.80,0.80,0.80],'MarkerFaceColor',[0.80,0.80,0.80])
hold on
h=plot(Pressure4StrainTotal,StrainTotal,'s','MarkerSize',1,'Color',[0.54,0.75,0.75],'MarkerFaceColor',[0.54,0.75,0.75])
set(h,'Linewidth',W);
h=plot(nanmean(InputPresReStrain,2),nanmean(StrainRes,2),'-','Color','#54a1a1','MarkerFaceColor','#54a1a1')
set(h,'Linewidth',W);
h=plot(nanmean(InputPresReStrain,2),nanmean(StrainResAll(:,2:2:end),2),'-','Color','#b8b8b8','MarkerFaceColor','#b8b8b8')
set(h,'Linewidth',W);
xlim([0 Pmaxglob])
ylim([0 1.5])
ylabel('Relative Strain (\%)','Interpreter','LaTex')
xlabel('Confining pressure $p$ (kPa)','Interpreter','LaTex')
a=gca;
a.FontSize=16;
a.TickLabelInterpreter = 'latex';
%text(20, 1.4,figname,'fontsize',16,'Interpreter','LaTex')

%% plot up
C1 = [0.3294 0.6314 0.6314];   % start (#54a1a1)
C2 = [0.09,0.33,0.33];         % end
N = size(StrainResAll,2);
colors = zeros(N,3);
for i = 1:N
    t = (i-1)/(N-1);          % goes from 0 to 1
    colors(i,:) = (1-t)*C1 + t*C2;
end


figure
hold on
baseColor = [0.3294 0.6314 0.6314];
for ii=1%:2:size(StrainResAll,2)
    shadeFactor = ii/(5-1); 
    col = baseColor * shadeFactor;
	plot(nanmean(InputPresReStrainAll(:,ii),2),nanmean(StrainResAll(:,ii),2),'s','MarkerSize',5,'Color',colors(ii,:),'MarkerFaceColor',colors(ii,:))
end
colormap(colors);
cb = colorbar;
cb.Ticks = linspace(0,1,N);
cb.TickLabels = 1:N;
ylabel('Relative Strain (\%)','Interpreter','LaTex')
xlabel('Confining pressure $p$ (kPa)','Interpreter','LaTex')
a=gca;
a.FontSize=16;
a.TickLabelInterpreter = 'latex';
ylim([0 2])
box on

%% internal stresses


Pmax=170
figure(11)
hold off
plot(InputPresureGlob(Inddwn).*(pi*r.^2)/(0.5^2).*1e5/1000,DataTotGlob(Inddwn,:)/1000,'.','Color',[0.80,0.80,0.80])
hold on
plot(InputPresureGlob(Indup).*(pi*r.^2)/(0.5^2).*1e5/1000,DataTotGlob(Indup,:)/1000,'.','Color',[0.54,0.75,0.75])
h=plot(InputPresResdwn,PresureBotResdwn,'-','Color','#b8b8b8')
W=4
set(h,'Linewidth',W);
h=plot(InputPresResdwn,PresuresidResdwn,'--','Color','#b8b8b8')
set(h,'Linewidth',W);
h=plot(InputPresResup,PresureBotResup,'-','Color','#54a1a1')
set(h,'Linewidth',W);
h=plot(InputPresResup,PresuresidResup,'--','Color','#54a1a1')
set(h,'Linewidth',W);
xlabel('Confining pressure $p$(kPa)','Interpreter','LaTex')
ylabel('Internal Stress (kPa)','Interpreter','LaTex')
xlim([0 Pmax+15])
ylim([0 Pintmax])
text(Pmax+1,200,'\sigma_{xx}','Fontsize',16)
text(Pmax+1,65,'\sigma_{yy}','Fontsize',16)
annotation(gcf,'arrow',[0.715467806841046 0.814059356136821],...
    [0.508971988795518 0.59467787114846],'Color','#54a1a1','LineWidth',4,'Headwidth',15);
annotation(gcf,'arrow',[0.753244466800804 0.63556338028169],...
    [0.851848739495797 0.724705882352941],'Color','#b8b8b8','LineWidth',4,'Headwidth',15);
annotation(gcf,'arrow',[0.849999999999999 0.698214285714285],...
    [0.411904761904764 0.400000000000003],'Color','#b8b8b8','LineWidth',4,'LineStyle','--','Headwidth',15);
annotation(gcf,'arrow',[0.719039235412473 0.844642857142857],...
    [0.257039215686276 0.285714285714286],'Color','#54a1a1','LineWidth',4,...
    'LineStyle','--','Headwidth',15);
a=gca;
a.FontSize=16;
a.TickLabelInterpreter = 'latex';
