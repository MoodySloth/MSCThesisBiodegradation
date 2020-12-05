% ReproducibilityPISOX
clear all
close all
clc

%% Load and clean data
[Acc, Vessel] = xlsread('Soilblanks_Reprod.xlsx');
[PISOX] = xlsread('PisoxAcc.xlsx');
[Tpis] = xlsread('Tpis.xlsx');

Tpis = Tpis(:,:)/17.5;

Exp2 = Acc(:,[4:6])/17;
Exp3 = Acc(:,[7:9])/15;
Thes = Acc(:,[10:13])/17.5;
PISOX = PISOX(:,:)/17.5;
Vessel = Vessel(1,:);

Exp2Corr = filloutliers(Exp2,'linear','movmedian',50);
Exp3Corr = filloutliers(Exp3,'linear','movmedian',50);

% Remove NAN
Exp2Corr = Exp2Corr(1:609,:);
ThesCorr = Thes(1:1072,:);

% Fill gaps Exp3
nanx = isnan(Exp3Corr);
t = 1:numel(Exp3Corr);
Exp3Corr(nanx) = interp1(t(~nanx),Exp3Corr(~nanx),t(nanx));

figure(1)
[p1B] = plot(Exp2Corr,'b','LineWidth',3)
hold on
[p2B] = plot(Exp3Corr,'r','LineWidth',3)
hold on
[p3B] = plot(ThesCorr,'g','LineWidth',3)
% hold on
% plot(PISOX,'m')

legend([p1B(1) p2B(1) p3B(1)],{'Experiment 2','Experiment 3','Thesis'})
grid minor
ax = gca;
set(gca,'FontSize',30)
title('Soil Blank Vessels')
yticks(0:0.1:1.3);
xlabel('Time (days)')
ylabel('Accumulated CO2 (mg/g soil)')
names = {'0','5','10','15','20','25','30','35','40','45','50','55','60'};
set(gca,'xtick',[0:120:1320],'xticklabel',names);

%% Mean + std
[Exp2Mean,Exp2Std]=accmeanstd(Exp2Corr);
[Exp3Mean,Exp3Std]=accmeanstd(Exp3Corr);
[ThesMean,ThesStd]=accmeanstd(ThesCorr);

figure(2)
[p1]=stdplot(1:609,Exp2Mean,Exp2Std,'b','b','b',0.3)

[p2]=stdplot(1:1289,Exp3Mean,Exp3Std,'r','r','r',0.15)

[p3]=stdplot(1:1072,ThesMean,ThesStd,'g','g','g',0.15)

legend([p1 p2 p3],{'Experiment 2','Experiment 3','Thesis'})
grid minor
ax = gca;
set(gca,'FontSize',30)
title('Mean and Standard Deviation of Soil Blanks') 
yticks(0:0.1:1.3);
ylim([0 1.3]);
xlim([0 1320]);
xlabel('Time (days)')
ylabel('Accumulated CO2 (mg/g soil)')
names = {'0','5','10','15','20','25','30','35','40','45','50','55','60'};
set(gca,'xtick',[0:120:1320],'xticklabel',names);

% max values
MaxExp2 = max(Exp2Std)
MaxExp3 = max(Exp3Std)
MaxThes = max(ThesStd)
%% PISOX reproducibility mean + std
figure(3)
% PISOX with blank from exp 3
Exp3Pis = Exp3Mean(1:1072,:); % adjust length blank exp 3 to PISOX
[Rpis3,Bpis3,Mpis3,Spis3] = repbiodegsm(Exp3Pis,PISOX,Tpis);
Xp = 1:1072;
[p2] = stdplot(Xp,Mpis3,Spis3,'r','r','r',0.3);

% PISOX with blank from Thesis
[RpisT,BpisT,MpisT,SpisT] = repbiodegsm(ThesMean,PISOX,Tpis);
[p3] = stdplot(Xp,MpisT,SpisT,'g','g','g',0.3);

% PISOX with blank from exp 2
PISOX_Exp2 = PISOX(1:609,:); % adjust length PISOX to blank exp 2
[Rpis2,Bpis2,Mpis2,Spis2] = repbiodegsm(Exp2Mean,PISOX_Exp2,Tpis);
Xp = 1:609;
[p1] = stdplot(Xp,Mpis2,Spis2,'b','b','b',0.3);

legend([p1 p2 p3],{'PISOX Exp2','PISOX Exp3','PISOX Thesis'}) % change legend and p's accordingly, and scale
grid minor
ax = gca;
set(gca,'FontSize',30)
title('Reproducibility') 
yticks(0:10:100);
ylim([-2 100]);
%% Statistics
Nans = NaN(463,1);
Mpis2Long = [Mpis2;Nans];
MeanPISOXCom = [Mpis2Long Mpis3 MpisT]';
MeanPISOXAll = nanmean(MeanPISOXCom)';
PISOXstd = nanstd(MeanPISOXCom)';

figure(4)
p4 = stdplot(Xp,MeanPISOXAll,PISOXstd,'g','g','g',0.3);
legend([p4],{'Combined PISOX'},'FontSize',20)
grid minor
ax = gca;
set(gca,'FontSize',30)
title('Biodegradation mean with std of 3 different blanks') 
yticks(-10:10:100);
ylim([-10 100]);

MaxStd = max(PISOXstd)

%%
Vmean = nanvar(MeanPISOXCom);
figure(5)
plot(Vmean,'LineWidth',3)
MaxVar = max(Vmean)
grid minor
ax = gca;
xlabel('Time (days)')
ylabel('Variance in Biodegradation (%)')
set(gca,'FontSize',30)
title('Variance of the Mean Between Experiments') 
xlim([0 1080]);
names = {'0','5','10','15','20','25','30','35','40','45'};
set(gca,'xtick',[0:120:1080],'xticklabel',names);
