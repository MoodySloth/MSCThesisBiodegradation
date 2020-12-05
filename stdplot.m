function [p] = stdplot(Xp,Mean,Std,MColor,FColor,EColor,Transparency)
% Function to plot mean and std
% p = plot output for legend
% Xp = matrix of the number of horizontal points
% Mean = matrix of the mean values (has to be same length/size as Xp)
% Std = matrix of the std values (has to be same length/size as Xp)
% MColor = color of the mean line (can be RGB value)
% FColor = color of the filled area (can be RGB value)
% EColor = color of the edge of filled area (can be RGB value)
Upper = Mean + Std;  % vertical
Lower = Mean - Std;

U = Upper'; % horizontal
L = Lower';

jbfill(Xp,U,L,FColor,EColor,0,Transparency);
hold on
p = plot(Mean,'Color',MColor,'LineWidth',3);
hold on
xlabel('Time (days)')
ylabel('Biodegradation %')
yticks(0:10:100);
ylim([0 100]);
xlim([0 1080]);
names = {'0','5','10','15','20','25','30','35','40','45','50'};
set(gca,'xtick',[0:120:1080],'xticklabel',names);
grid on
end