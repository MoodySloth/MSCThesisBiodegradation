%% %%% Seawater - Conductivity %%% %%
%% Clean start
clear all
close all
clc
%% Load data %%
% Load seawater Conductivity
% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 19);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:S1169";

% Specify column names and types
opts.VariableNames = ["Time", "V7", "V36", "V50", "V61", "V65", "V22", "V42", "V59", "V54", "V62", "V63", "V4", "V8", "V35", "V49", "V64", "V91", "V92"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Time", "InputFormat", "");

% Import the data
Cond = readtable("/Users/Debora/Documents/My_Documents/Uva_Master/Master_Thesis/My_Data/Seawater_Cond.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts
%% Load C Ratio 
opts = spreadsheetImportOptions("NumVariables", 10);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:J2";

% Specify column names and types
opts.VariableNames = ["PISOX", "Hexanediol", "Isosorbide", "OxalicAcid", "NAOX", "NaAc", "Cellulose", "FDCA", "MonomerMix", "PET"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
Cratio = readtable("/Users/Debora/Documents/My_Documents/Uva_Master/Master_Thesis/My_Data/C_ratio.xlsx", opts, "UseExcel", false);

% To array
CratioA = table2array(Cratio); 

% Clear temporary variables
clear opts
%% Load Mass Seawater
opts = spreadsheetImportOptions("NumVariables", 14);

% Specify sheet and range
opts.Sheet = "Seawater mass";
opts.DataRange = "A2:N3";

% Specify column names and types
opts.VariableNames = ["Cellulose20u_seawater_2", "Cellulose20u_seawater_1", "Cellulose20u_seawater_3", "SodiumAc_seawater_3", "SodiumAc_seawater_2", "SodiumAc_seawater_1", "PISOX25HD_seawater_1", "PISOX25HD_seawater_2", "PISOX25HD_seawater_3", "PISOX25HD_seawater_4", "PISOX25HD_seawater_5", "PISOX1", "PISOX2", "PISOX3"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
SWMass = readtable("/Users/Debora/Documents/My_Documents/Uva_Master/Master_Thesis/My_Data/Mass_materials.xlsx", opts, "UseExcel", false);

% To array and mg
SWMassA = table2array(SWMass);
SWMassA = [SWMassA(1,:) ; SWMassA(2,:) * 1000];

% Clear temporary variables
clear opts
%% Clean data
Vessel = fieldnames(Cond);
Vessel = Vessel(2:19,:)';  % selects vessel numbers

CondT = Cond(:,2:19);    % selects data
CondA = table2array(CondT); % transform into array

% Remove data from 67h to 101h
CondA(67:101,1:16) = NaN;
% Replace values with interpolation
figure(1)
plot(CondA)

nanx = isnan(CondA);
t = 1:numel(CondA);
CondA(nanx) = interp1(t(~nanx),CondA(~nanx),t(nanx));

figure(2)
plot(CondA)

% Fix Sodium Acetate (added 4 days later)
A = NaN(96,3);
SdAc = [CondA(97:1168,9:11);A];
CondAll = [CondA(:,1:8) SdAc(:,:) CondA(:,12:18)];
%% Remove outliers
CondAllO = filloutliers(CondA,'linear','movmedian',50); 

diff = CondAll - CondAllO;
figure(3)  % a few vessels before and after outlier removal
plot(1:1168,CondAll(:,[1:5]),'b',1:1168,CondAllO(:,[1:5]),'g')
legend('Blue: original', 'Green: outliers removed')
%% Captured CO2
% Set correct matrix size
rows = 1168;
col = 18;
%% Captured CO2 - Create Correct Ct0 matrix
% Find max value for Ct0 in first 24 hours
s = zeros(rows,col);
s(1,:) = 1;  % first row filled with ones
CT0 = zeros(rows,col);

for i = 1:rows;
    for j = 1:col;
        if s(i,j) == 1;
          CT0(i,j) = max(CondAllO(i:i+24,j)); % matrix in which the first row is filled with starting ct0 values
        end
    end
end

% Find max value for Ct0 for 24 hours after new KOH added
for i = 2:rows;
    for j = 1:col;
        if CondAllO(i,j) - CondAllO(i-1,j) >= 4.5;  % difference in cond
            CT0(i,j) = max(CondAllO(i:i+24,j)); % matrix with the max value in the location of diff
        end
    end
end

% Fill zeros with actual Ct0 values
for i = 2:rows;
    for j = 1:col;
        if CT0(i,j) == 0;
            CT0(i,j) = CT0(i-1,j); % if element is zero, replace it for the element before it
        end
    end
end

% Fix starting peaks in first 24 hours 
RangeFval = zeros(rows,col);
RangeFval(1:24,:) = CondAllO(1:24,:); % selects only the first 24h
MaxFloc = zeros(rows,col);         
for i = 1:rows;
    for j = 1:col;
        if RangeFval(i,j) == CT0(i,j);  % if in the first 24h, a value is equal to the max,
            MaxFloc(i,j) = 1;           % then that position will be filled with a 1 in MaxFloc
        end
    end
end

BoundariesF = MaxFloc + s;  % s has the first row filled with ones
                            % BoundariesF ends in the location of the max

for i = 1:1144;         % 1168 minus 24 to not exceed array bounds
    for j = 1:col;
        if any(BoundariesF(i+1:i+24,j) == 1); % if any value in a 24h range is 1,
            BoundariesF(i,j) = 1;             % then that position will be filled with a 1
        end
    end
end

% Replace CT0 in interval from new koh to max with original conductivity for
% the first 24h of incubation
for i = 1:rows;
    for j = 1:col;
        if BoundariesF(i,j) == 1;
            CT0(i,j) = CondAllO(i,j); % assign ct0 values as themselves, instead of 
        end                           % using the max value, to avoid little peaks in Acc CO2
    end
end

% Fix peak values after new KOH
% Find location new KOH
Loc = zeros(rows,col);

for i = 2:rows;
    for j = 1:col;
        if CondAllO(i,j) - CondAllO(i-1,j) >= 4.5;
    Loc(i,j) = 1; % if there is a diff of 5mS, fills that location with 1 in the Loc matrix
        end
    end
end

% Create Range of 24 hours after new KOH
Range = zeros(rows,col);

for i = 1:rows;
    for j = 1:col;
        if Loc(i,j) == 1;
         Range(i:i+24,j) = 1;
        end
    end
end

% Show values in 24 hour range after new KOH
Rangeval = zeros(rows,col);

for i = 1:rows;
    for j = 1:col;
        if Range(i,j) == 1;
            Rangeval(i,j) = CondAllO(i,j);
        end
    end
end

% Replace 0 with NaN
for i = 1:rows;
    for j = 1:col;
        if Rangeval(i,j) == 0;
            Rangeval(i,j) = NaN;
        end
    end
end

% Location of max value
Maxloc = zeros(rows,col);
for i = 1:rows;
    for j = 1:col;
        if Rangeval(i,j) == CT0(i,j);
            Maxloc(i,j) = 1;
        end
    end
end

% Location of boundaries: diff of 5mS until max value 
Boundaries = Maxloc + Loc;

% Fill with NaNs 
for i = 1:rows;
    for j = 1:col;
        if Range(i,j) == 0;        % outside the range of 24h after new koh
            Boundaries(i,j) = NaN; % fill those positions with nan
        end
    end
end

NaNLoc = isnan(Boundaries);

for i = 1:1144;         % rows minus 24 to not exceed array bounds
    for j = 1:col;
        if any(Boundaries(i+1:i+24,j) == 1) & NaNLoc(i,j) == 0; % if any value in a 24h range is 1 and not a nan,
            Boundaries(i,j) = 1;   % then that position will be filled with a 1
        end
    end
end
        
% Replace CT0 in range from new koh to max with original conductivity
for i = 1:rows;
    for j = 1:col;
        if Boundaries(i,j) == 1;       
            CT0(i,j) = CondAllO(i,j); % assign ct0 values as themselves, instead of
        end                            % using the max value, to avoid little peaks in Acc CO2
    end
end

% Calculate Captured CO2
A = 219;
CO2All = A*((CT0-CondAllO)./CT0);

figure(4) % a few vessels
plot(CO2All(:,[1 2 3 8 10 15]))
legend(Vessel)
%% Accumulated CO2 - Find range to sum and find add value
% Find value to add and store in correct location
AddValLoc = zeros(rows,col); % location of the value to be added
for i = 2:rows;
    for j = 1:col;
        if Loc(i,j) == 1;    % in the locations equal to 1
            AddValLoc(i,j) = CO2All(i-1,j); % gets the CO2 value from one position 
                                            % before the change of KOH (diff of 5mS)
                                            % and store it where Loc equals 1                                  
        end
    end
end

% Replace zeros after the CO2 value in AddValLoc with AddVal
AddVal = AddValLoc;
for i = 2:rows;
    for j = 1:col;
        if AddVal(i-1,j) > 0 ;  % if previous value is higher than zero
            AddVal(i,j) = AddVal(i-1,j);  % fills that position with previous value
        end
    end
end
%% Accumulatted CO2 - add value
AccCO2All = CO2All + AddVal;

figure(3) % all
plot(AccCO2All)
legend(Vessel)

Acc = AccCO2All;
%% Blanks
Blanks = Acc(:,1:5);

% Validity of blanks
% Applying "Blanks must be within (ACTUALLY UNTIL) 30% of the mean at the plateau level of the test
% material, Pisox"

figure(5)
plot(Acc(:,[12:16]))  % pisox
MP = mean(Acc([600:1168],[12:16])); % choosing 600h as start of plateau phase

% the limit is 30% of MP
Lim = 0.3*MP;
figure(6)
plot(Acc(:,[1:5])) 
legend(Vessel)  % vessels 36,50,65 , except at 668h

% Choosing 4.8832 (from Lim) as the limit, the vessels 36,50,65 are valid by looking at the graph (except at 668h)
% but vessels 50 and 65 are negative from 0h to 65h (cond increases), so
% vessel 36 will be chosen as the blank for that period and for the rest the
% 3 of them (as the mean)
V36 = Acc([1:101],2); % only vessel 36 for beginning
ValidVs = Acc([102:1168],[2 3 5]); % the 3 vessels: 36,50,65
MValidVs = mean(ValidVs,2);   % mean of vessels 36,50,65 for the remaining time

Blank = [V36; MValidVs];

figure(7)
plot(Acc(:,1:5),'LineWidth',3) 
legend(Vessel)
title('Seawater Blanks')
xlabel('Time (hours)')
ylabel('Accumulated CO2 (mg)') 
ax = gca;
set(gca,'FontSize',30)
grid minor
%% Calculate theoretical carbon and percentage of biodegradation 
%Th = (44/12)*m*w  m is mass of material in mg, w is carbon content as mass fraction, 
% 44 is the molecular mass of carbon dioxide, and 12 is the atomic mass of carbon.

%Bio = (Resp/Th)*100 is for each test flask
%% Theoretical CO2
Tcel = (44/12)*SWMassA(2,1:2)*CratioA(:,8);
Tsod = (44/12)*SWMassA(2,4:6)*CratioA(:,7);
Tpis = (44/12)*SWMassA(2,7:11)*CratioA(:,1);
%% Cellulose
Xp = 1:1168; % Standard length of horizontal axis/time

% Extract Cellulose
Cel = Acc(:,6:7); % without col 8 (v59) bc it's too different than rest

% Respiration & Biodegradation
[Rcel,Bcel,Mcel,Scel] = repbiodeg(Blank,Cel,Tcel);           % Without filter
[Rcel,BcelSm,McelSm,ScelSm] = repbiodegsm(Blank,Cel,Tcel);   % With filter

figure(8)   % Figure 8 from here out with no Smoothing
[p1] = stdplot(Xp,Mcel,Scel,'r','r','r',0.3);

figure(9)   % Figure 9 from here out with Smoothing
[p1sm] = stdplot(Xp,McelSm,ScelSm,'r','r','r',0.3);
%% Sodium acetate
Xp = 1:1072; % Adjust horizontal points for shorter time
% Extract Sodium Acetate
Sod = Acc(1:1072,9:11); % shorter interval
BlankSh = Blank(1:1072); % Adjust Blank to correct time/length

% Respiration & Biodegradation
[Rsod,Bsod,Msod,Ssod] = repbiodeg(BlankSh,Sod,Tsod);           % Without filter
[Rsod,BsodSm,MsodSm,SsodSm] = repbiodegsm(BlankSh,Sod,Tsod);   % With filter

figure(8)
[p2] = stdplot(Xp,Msod,Ssod,'g','g','g',0.3);

figure(9)
[p2sm] = stdplot(Xp,MsodSm,SsodSm,'g','g','g',0.3);
%% Pisox
Xp = 1:1168; % Standard length of horizontal axis/time
% Extract Pisox
Pis = Acc(:,12:16); 
TpisEx = Tpis(:,1:5); % Adjust T for exclusion

% Respiration & Biodegradation
[Rpis,Bpis,Mpis,Spis] = repbiodeg(Blank,Pis,TpisEx);           % Without filter
[Rpis,BpisSm,MpisSm,SpisSm] = repbiodegsm(Blank,Pis,TpisEx);   % With filter

figure(8)
[p3] = stdplot(Xp,Mpis,Spis,'m','m','m',0.3);
yticks(0:1:8);
ylim([0 8]);
legend([p1 p2 p3],{'Cellulose','Sodium Acetate','Pisox'},'FontSize',20)
xlim([0 1200]);
names = {'0','5','10','15','20','25','30','35','40','45','50'};
set(gca,'xtick',[0:120:1200],'xticklabel',names);
grid minor

figure(9)
[p3sm] = stdplot(Xp,MpisSm,SpisSm,'m','m','m',0.3);
yticks(0:5:60);
ylim([0 60]);
legend([p1sm p2sm p3sm],{'Cellulose','Sodium Acetate','PISOX'},'FontSize',20)
xlim([0 1200]);
set(gca,'xtick',[0:120:1200],'xticklabel',names);
grid minor
ax = gca;
set(gca,'FontSize',30)
title('Biodegradation in Seawater') 
%%
% % PISOX in MiliQ water
% Water = Acc(:,17:18); 
% figure (7)
% plot(Water)

