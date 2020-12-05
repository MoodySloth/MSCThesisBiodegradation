%% %%% Interface - Conductivity %%% %%
%% Clean start
clear all
close all
clc
%% Load data %%
% Load seawater Conductivity
% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 22);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:V1169";

% Specify column names and types
opts.VariableNames = ["Time", "V1", "V10", "V41", "V56", "V58", "V3", "V17", "V40", "V52", "V55", "V57", "V9", "V46", "V48", "V53", "V60", "V2", "V18", "V34", "V47", "V51"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Time", "InputFormat", "");

% Import the data
Cond = readtable("/Users/Debora/Documents/My_Documents/Uva_Master/Master_Thesis/My_Data/Interface_Cond.xlsx", opts, "UseExcel", false);

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
%% Load Mass Interface
% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 16);

% Specify sheet and range
opts.Sheet = "Interface mass";
opts.DataRange = "A2:P3";

% Specify column names and types
opts.VariableNames = ["Cellulose20u_sedimentseawater_3", "Cellulose20u_sedimentseawater_1", "Cellulose20u_sedimentseawater_2", "SodiumAc_sedimentseawater_2", "SodiumAc_sedimentseawater_1", "SodiumAc_sedimentseawater_3", "PISOX25HD_sedimentseawater_1", "PISOX25HD_sedimentseawater_3", "PISOX25HD_sedimentseawater_2", "PISOX25HD_sedimentseawater_4", "PISOX25HD_sedimentseawater_5", "PISOX25HDlowconcentration_sedimentseawater_5", "PISOX25HDlowconcentration_sedimentseawater_4", "PISOX25HDlowconcentration_sedimentseawater_3", "PISOX25HDlowconcentration_sedimentseawater_2", "PISOX25HDlowconcentration_sedimentseawater_1"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
IMass = readtable("/Users/Debora/Documents/My_Documents/Uva_Master/Master_Thesis/My_Data/Mass_materials.xlsx", opts, "UseExcel", false);

% To array and mg
IMassA = table2array(IMass);
IMassA = [IMassA(1,:) ; IMassA(2,:) * 1000];

% Clear temporary variables
clear opts
%% Clean data 
Vessel = fieldnames(Cond);
Vessel = Vessel(2:22,:)';  % selects vessel numbers

CondT = Cond(:,2:22);    % selects data
CondA = table2array(CondT); % transform into array
plot(CondA)
ylim([0 200])
title('Conductivity (mS)')

% Remove data from 67h to 101h
CondA(67:101,:) = NaN;
% Replace values with interpolation
figure(1)
plot(CondA)

nanx = isnan(CondA);
t = 1:numel(CondA);
CondA(nanx) = interp1(t(~nanx),CondA(~nanx),t(nanx));

figure(2)
plot(CondA)

% Fix Sodium Acetate and PISOX Low (added 4 days later)
A = NaN(96,3);
B = NaN(96,5);
SdAc = [CondA(97:1168,9:11);A];
PLow = [CondA(97:1168,17:21);B];
CondAll = [CondA(:,1:8) SdAc(:,:) CondA(:,12:16) PLow(:,:)];
%% Remove outliers
CondAllO = filloutliers(CondA,'linear','movmedian',50); 

figure(3)
diff = CondA - CondAllO;
figure(3)  % a few vessels before and after outlier removal
plot(1:1168,CondA(:,[1:5]),'b',1:1168,CondAllO(:,[1:5]),'g')
legend('Blue: original', 'Green: outliers removed')
%% Captured CO2
% Set correct matrix size
rows = 1168;
col = 21;
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
        if CondAllO(i,j) - CondAllO(i-1,j) >= 10;  % difference in cond
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
        if CondAllO(i,j) - CondAllO(i-1,j) >= 10;
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

% Location of boundaries: diff of 10mS until max value 
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

figure(2) % a few vessels
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
figure(4)
plot(Acc(:,1:5)) 
legend(Vessel)
plot(Acc(:,1:5),'LineWidth',3) 
legend(Vessel)
title('Interface Blanks')
xlabel('Time (hours)')
ylabel('Accumulated CO2 (mg)') 
ax = gca;
set(gca,'FontSize',30)
grid minor

Mbla = mean(Blanks,2);   
Sbla = std(Blanks,0,2);

Blank = Mbla; % mean of the blanks
%% Calculate theoretical carbon and percentage of biodegradation 
%Th = (44/12)*m*w  m is mass of material in mg, w is carbon content as mass fraction, 
% 44 is the molecular mass of carbon dioxide, and 12 is the atomic mass of carbon.

%Bio = (Resp/Th)*100 is for each test flask
%% Theoretical CO2
Tcel = (44/12)*IMassA(2,1:3)*CratioA(:,7);
Tsod = (44/12)*IMassA(2,4:6)*CratioA(:,6);
Tpis = (44/12)*IMassA(2,7:11)*CratioA(:,1);
Tpisl = (44/12)*IMassA(2,12:16)*CratioA(:,1);
%% Cellulose
Xp = 1:1168; % Standard length of horizontal axis/time

% Extract Cellulose
Cel = Acc(:,6:8);

% Respiration & Biodegradation
[Rcel,Bcel,Mcel,Scel] = repbiodeg(Blank,Cel,Tcel);           % Without filter
[Rcel,BcelSm,McelSm,ScelSm] = repbiodegsm(Blank,Cel,Tcel);   % With filter

figure(5)   % Figure 5 from here out with no Smoothing
[p1] = stdplot(Xp,Mcel,Scel,'r','r','r',0.3);

figure(6)   % Figure 6 from here out with Smoothing
[p1sm] = stdplot(Xp,McelSm,ScelSm,'r','r','r',0.3);
%% Sodium acetate
Xp = 1:1072; % Adjust horizontal points for shorter time
% Extract Sodium Acetate
Sod = Acc(1:1072,9:11);
BlankSh = Blank(1:1072); % Adjust Blank to correct time/length

% Respiration & Biodegradation
[Rsod,Bsod,Msod,Ssod] = repbiodeg(BlankSh,Sod,Tsod);           % Without filter
[Rsod,BsodSm,MsodSm,SsodSm] = repbiodegsm(BlankSh,Sod,Tsod);   % With filter

figure(5)
[p2] = stdplot(Xp,Msod,Ssod,'g','g','g',0.3);

figure(6)
[p2sm] = stdplot(Xp,MsodSm,SsodSm,'g','g','g',0.3);
%% PISOX
Xp = 1:1168; % Standard length of horizontal axis/time
% Extract Pisox
Pis = Acc(:,12:16);

% Respiration & Biodegradation
[Rpis,Bpis,Mpis,Spis] = repbiodeg(Blank,Pis,Tpis);           % Without filter
[Rpis,BpisSm,MpisSm,SpisSm] = repbiodegsm(Blank,Pis,Tpis);   % With filter

figure(5)
[p3] = stdplot(Xp,Mpis,Spis,'m','m','m',0.3);

figure(6)
[p3sm] = stdplot(Xp,MpisSm,SpisSm,'m','m','m',0.3); 
%% PISOX Low
Xp = 1:1072; % Adjust horizontal points for shorter time
% Extract Pisox Low
Pisl = Acc(1:1072,17:21);
BlankSh = Blank(1:1072); % Adjust Blank to correct time/length

% Respiration & Biodegradation
[Rpisl,Bpisl,Mpisl,Spisl] = repbiodeg(BlankSh,Pisl,Tpisl);           % Without filter
[Rpisl,BpislSm,MpislSm,SpislSm] = repbiodegsm(BlankSh,Pisl,Tpisl);   % With filter

figure(5)
[p4] = stdplot(Xp,Mpisl,Spisl,'c','c','c',0.3);
yticks(0:5:70);
ylim([0 70]);
legend([p1 p2 p3 p4],{'Cellulose','Sodium Acetate','Pisox','Pisox low'},'FontSize',20)
xlim([0 1200]);
names = {'0','5','10','15','20','25','30','35','40','45','50'};
set(gca,'xtick',[0:120:1200],'xticklabel',names);
grid minor

figure(6)
[p4sm] = stdplot(Xp,MpislSm,SpislSm,'c','c','c',0.3);
yticks(0:5:60);
ylim([0 60]);
legend([p1sm p2sm p3sm p4sm],{'Cellulose','Sodium Acetate','PISOX','PISOX Low'},'FontSize',20)
xlim([0 1200]);
set(gca,'xtick',[0:120:1200],'xticklabel',names);
grid minor
ax = gca;
set(gca,'FontSize',30)
title('Biodegradation at the seawater-sediment interface') 

