%% %%% Soil - Conductivity %%% %%
%% Clean start
clear all
close all
clc
%% Load data %%
% Load soil Conductivity
% Setup the Import Options and import the data
opts = spreadsheetImportOptions("NumVariables", 44);

% Specify sheet and range
opts.Sheet = "Sheet1";
opts.DataRange = "A2:AR1073";

% Specify column names and types
opts.VariableNames = ["Time", "V5", "V11", "V27", "V43", "V77", "V86", "V87", "V88", "V13", "V20", "V33", "V28", "V69", "V79", "V6", "V15", "V31", "V25", "V30", "V71", "V24", "V67", "V68", "V89", "V90", "V16", "V21", "V78", "V23", "V45", "V66", "V14", "V38", "V72", "V26", "V39", "V44", "V80", "V81", "V82", "V83", "V84", "V85"];
opts.VariableTypes = ["datetime", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify variable properties
opts = setvaropts(opts, "Time", "InputFormat", "");

% Import the data
Cond = readtable("/Users/Debora/Documents/My_Documents/Uva_Master/Master_Thesis/My_Data/Soil_Cond.xlsx", opts, "UseExcel", false);

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
%% Load Mass Soil
opts = spreadsheetImportOptions("NumVariables", 35);

% Specify sheet and range
opts.Sheet = "Soil mass";
opts.DataRange = "A2:AI4";

% Specify column names and types
opts.VariableNames = ["Cellulose20u1", "Cellulose20u2", "Cellulose20u3", "Hexanediollowconcentration_1", "Hexanediollowconcentration_2", "Hexanediollowconcentration_3", "IsosorbidethenISO_1", "IsosorbidethenISO_2", "IsosorbidethenISO_3", "IsosorbidethenPISOX25HD_1", "IsosorbidethenPISOX25HD_2", "IsosorbidethenPISOX25HD_3", "Oxalicacidlowconcentration_1", "Oxalicacidlowconcentration_2", "Oxalicacidlowconcentration_3", "PETfibers1", "PETfibers2", "PETfibers3", "PISOX25HDmonomermixture_1", "PISOX25HDmonomermixture_2", "PISOX25HDmonomermixture_3", "PSIXO25HDthenISO_1", "PSIXO25HDthenISO_2", "PSIXO25HDthenISO_3", "PSIXO25HDthenPISOX25HD_1", "PSIXO25HDthenPISOX25HD_2", "PSIXO25HDthenPISOX25HD_3", "FDCA_1", "FDCA_2", "FDCA_3", "NaOxalate_1", "NaOxalate_2", "NaOxalate_3", "Oxalicacidlowconcentration_4", "Oxalicacidlowconcentration_5"];
opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Import the data
Mass = readtable("/Users/Debora/Documents/My_Documents/Uva_Master/Master_Thesis/My_Data/Mass_materials.xlsx", opts, "UseExcel", false);

% To array and mg
MassA = table2array(Mass);
MassA(2,:) = MassA(2,:) * 1000;
MassA(3,:) = MassA(3,:) * 1000;

% Clear temporary variables
clear opts
%% Clean data, separate 80-90 and 79
Vessel = fieldnames(Cond);
Vessel = Vessel(2:44,:)';

CondT = Cond(:,2:44);
CondA = table2array(CondT);
plot(CondA(:,12:14)) % v79 is a bad data vessel, removed later on
plot(CondA(:,6:7)) % v86 and v87 are a straight line with a square peak

Cond8090 = CondA(:,[8 24 25 38:42]);     % Extract 80-90 because they start later and remove 85,86,87 (bad vessels)
Vessel8090 = Vessel(:,[8 24 25 38:42]); % Correct vessel order for 80-90

CondN = CondA(:,[1:5 9:23 26:37]);     % Extract normal vessels 
VesselN = Vessel(:,[1:5 9:23 26:37]); % Correct vessel order for normal vessels

CondB = CondA(:,[6 7 43]) ;        % Extract bad vessels (85,86,87)
VesselB = Vessel(:,[6 7 43]) ;     % Correct vessel order for bad vessels
%% Fix time 80-90 and create complete matrix for all correct vessels
A = NaN(238,8);
Cond8090C = [Cond8090(239:1072,:);A]; % fill the matrix with nan values

B = NaN(215,1);
Cond79 = [CondA(216:1072,14);B];

CondAll = [CondN(:,[1:5]) Cond8090C(:,1) CondN(:,[6:10]) Cond79(:,1) CondN(:,[12:20]) Cond8090C(:,[2 3])...
    CondN(:,[21:32]) Cond8090C(:,[4:8])];   % combine CondN and Cond8090C correct
VesselAll = Vessel(:,[1:5 8:42]);
%% Remove outliers
CondAllO = filloutliers(CondAll,'linear','movmedian',50); 
diff = CondAll - CondAllO;
figure(1)  % a few vessels before and after outlier removal
plot(1:1072,CondAll(:,[1:5]),'b',1:1072,CondAllO(:,[1:5]),'g')
legend('Blue: original', 'Green: outliers removed')
%% Captured CO2
% Set correct matrix size
rows = 1072;
col = 40;
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
        if CondAllO(i,j) - CondAllO(i-1,j) >= 5;  % difference of 5mS in cond
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

for i = 1:1048;         % 1072 minus 24 to not exceed array bounds
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
        if CondAllO(i,j) - CondAllO(i-1,j) >= 5;
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

for i = 1:1048;         % 1072 minus 24 to not exceed array bounds
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
plot(CO2All(:,[1 2 3 39 40]))
legend(VesselAll)
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

% V83 (col39) V84 (col40) have two new koh and need to be corrected
% Put back original 2nd KOH value
for i = 1:rows;
    for j = 1:col;
        if AddVal(i,j) < AddValLoc(i,j);
            AddVal(i,j) = AddValLoc(i,j) + AddVal(i,j);
        end
    end
end

% Replace values after 2nd KOH value with correct value
for i = 2:rows;
    for j = 1:col;
        if AddVal(i-1,j) > AddVal(i,j);
            AddVal(i,j) = AddVal(i-1,j);
        end
    end
end
%% Accumulated CO2 - add value
Acc = CO2All + AddVal;

figure(3) % all
plot(Acc)
legend(VesselAll)
%% Validity of Blanks
figure(6)
x = plot(Acc(:,[1:6]))  
legend(VesselAll)
title('Soil Blanks')
xlabel('Time (h)')
ylabel('Accum. CO2 (mg)')
xlim([0 1100])
ylim([0 26])
set(x,'LineWidth',3)
grid on
grid minor
ax = gca;
ax.FontSize = 30;

Blank = mean(Acc(:,[1:4]),2);  % only blanks 5, 11, 27, 43 selected (77 had mold, 88 is shorter)
figure(7)
plot(Blank)
%close all % change this line to see all graphs

% Accumulated CO2 mean & std of Blanks
Xp = 1:1072;
Bla = Acc(:,1:4);
[ABlaM,ABlaS] = accmeanstd(Bla);
figure(8) % Figure 8 from here out is Accumulated CO2
[pBA] = stdplot(Xp,ABlaM,ABlaS,[0.8500, 0.3250, 0.0980],[0.8500, 0.3250, 0.0980],[0.8500, 0.3250, 0.0980],0.15); 
%% Calculate theoretical carbon and percentage of biodegradation
%Th = (44/12)*m*w  m is mass of material in mg, w is carbon content as mass fraction, 
% 44 is the molecular mass of carbon dioxide, and 12 is the atomic mass of carbon.

%Bio = (Resp/Th)*100 is for each test flask
%% Theoretical CO2
Tcel = (44/12)*MassA(2,1:3)*CratioA(:,7); 
Thex = (44/12)*MassA(2,4:6)*CratioA(:,2);
Toxa = (44/12)*MassA(2,[13:15 34 35])*CratioA(:,4);   
Tpet = (44/12)*MassA(2,16:18)*CratioA(:,10);        
Tmix = (44/12)*MassA(2,19:21)*CratioA(:,9);             
Tpis = (44/12)*MassA(2,22:27)*CratioA(:,1);  
Tfdca = (44/12)*MassA(2,28:30)*CratioA(:,8);
Tnaox = (44/12)*MassA(2,31:32)*CratioA(:,5);

dlmwrite('Tpis.txt',Tpis,'delimiter','\t','newline','pc')

% ThCO2 for double groups
Tiso = (44/12)*MassA(2,7:9)*CratioA(:,3);
Tiso2 = (Tiso) + (44/12)*MassA(3,7:9)*CratioA(:,3);
Tpis2 = (Tiso) + (44/12)*MassA(3,10:12)*CratioA(:,1);
%% Cellulose
Xp = 1:1072; % Standard length of horizontal axis/time

% Extract Cellulose
Cel = Acc(:,7:9); 

% Accumulated CO2 mean & std
% [ACelM,ACelS] = accmeanstd(Cel);
% figure(8) % Figure 8 from here out is Accumulated CO2
% [p1A] = stdplot(Xp,ACelM,ACelS,'r','r','r',0.15);

% Respiration & Biodegradation
[Rcel,Bcel,Mcel,Scel] = repbiodeg(Blank,Cel,Tcel);           % Without filter
[Rcel,BcelSm,McelSm,ScelSm] = repbiodegsm(Blank,Cel,Tcel);   % With filter

figure(9)   % Figure 9 from here out with no Smoothing
[p1] = stdplot(Xp,Mcel,Scel,'r','r','r',0.3);

figure(10)   % Figure 10 from here out with Smoothing (Cell, Hex, Oxa, Naox, Pis)
[p1sm10] = stdplot(Xp,McelSm,ScelSm,'r','r','r',0.3);

figure(11)   % Figure 10 from here out with Smoothing (Cell, Pet, Mix, Pis)
[p1sm11] = stdplot(Xp,McelSm,ScelSm,'r','r','r',0.3);
%% Hexanediol Low
% Extract Hexanediol Low
Hex = Acc(:,10:11); % Exclude vessel 79 (column 12) due to shorter time
ThexEx = Thex(:,1:2); % Adjust T for exclusion

% Accumulated CO2 mean & std
% [AHexM,AHexS] = accmeanstd(Hex);
% figure(8) % Figure 8 from here out is Accumulated CO2
% [p2A] = stdplot(Xp,AHexM,AHexS,'g','g','g',0.15);

% Respiration & Biodegradation
[Rhex,Bhex,Mhex,Shex] = repbiodeg(Blank,Hex,ThexEx);           % Without filter
[Rhex,BhexSm,MhexSm,ShexSm] = repbiodegsm(Blank,Hex,ThexEx);   % With filter

figure(9)
[p2] = stdplot(Xp,Mhex,Shex,'g','g','g',0.3);

figure(10)
[p2sm] = stdplot(Xp,MhexSm,ShexSm,'g','g','g',0.3);
%% ISO then ISO - Iso double
% Extract ISO then ISO
Isod = Acc(:,13:15);

% Accumulated CO2 mean & std
[AIsodM,AIsodS] = accmeanstd(Isod);
figure(8) % Figure 8 from here out is Accumulated CO2
[p3A] = stdplot(Xp,AIsodM,AIsodS,'y','y','y',0.15);

% Respiration & Biodegradation in 2 parts
[Risod1,Bisod1,Misod1,Sisod1] = repbiodeg(Blank(1:504,:),Isod(1:504,:),Tiso);           % Without filter
[Risod1,BisodSm1,MisodSm1,SisodSm1] = repbiodegsm(Blank(1:504,:),Isod(1:504,:),Tiso);   % With filter

[Risod2,Bisod2,Misod2,Sisod2] = repbiodeg(Blank(505:1072,:),Isod(505:1072,:),Tiso2);     % Without filter
[Risod2,BisodSm2,MisodSm2,SisodSm2] = repbiodegsm(Blank(505:1072,:),Isod(505:1072,:),Tiso2);   % With filter

Risod = [Risod1 ; Risod2];

Bisod = [Bisod1 ; Bisod2];
Misod = [Misod1 ; Misod2];
Sisod = [Sisod1 ; Sisod2];

BisodSm = [BisodSm1 ; BisodSm2];
MisodSm = [MisodSm1 ; MisodSm2];
SisodSm = [SisodSm1 ; SisodSm2];

figure(9)
[p3] = stdplot(Xp,Misod,Sisod,'y','y','y',0.3);

% figure(10)
% % [p3sm] = stdplot(Xp,MisodSm,SisodSm,'y','y','y',0.3);
%% ISO then PISOX 
% Extract ISO then PISOX
Isop = Acc(:,16:18);

% Accumulated CO2 mean & std
[AIsopM,AIsopS] = accmeanstd(Isop);
figure(8) % Figure 8 from here out is Accumulated CO2
[p4A] = stdplot(Xp,AIsopM,AIsopS,'b','b','b',0.15);

% Respiration & Biodegradation in 2 parts
[Risop1,Bisop1,Misop1,Sisop1] = repbiodeg(Blank(1:504,:),Isop(1:504,:),Tiso);           % Without filter
[Risop1,BisopSm1,MisopSm1,SisopSm1] = repbiodegsm(Blank(1:504,:),Isop(1:504,:),Tiso);   % With filter

[Risop2,Bisop2,Misop2,Sisop2] = repbiodeg(Blank(505:1072,:),Isop(505:1072,:),Tpis2);     % Without filter
[Risop2,BisopSm2,MisopSm2,SisopSm2] = repbiodegsm(Blank(505:1072,:),Isop(505:1072,:),Tpis2);   % With filter

Risop = [Risop1 ; Risop2];

Bisop = [Bisop1 ; Bisop2];
Misop = [Misop1 ; Misop2];
Sisop = [Sisop1 ; Sisop2];

BisopSm = [BisopSm1 ; BisopSm2];
MisopSm = [MisopSm1 ; MisopSm2];
SisopSm = [SisopSm1 ; SisopSm2];
 
figure(9)
[p4] = stdplot(Xp,Misop,Sisop,'b','b','b',0.3);

% figure(10)
% % [p4sm] = stdplot(Xp,MisopSm,SisopSm,'b','b','b',0.3);
%% Oxalic Low
% Extract Oxalic Low
Oxa = Acc(:,[19 21]); % Exclude vessels 89, 90 (column 22, 23) due to shorter time and col20, v67, due to leakage
ToxaEx = Toxa(:,1:2); % Adjust T for exclusion

% Accumulated CO2 mean & std
% [AOxaM,AOxaS] = accmeanstd(Oxa);
% figure(8) % Figure 8 from here out is Accumulated CO2
% [p5A] = stdplot(Xp,AOxaM,AOxaS,[0.5,0,0.5],[0.5,0,0.5],[0.5,0,0.5],0.15);

% Respiration & Biodegradation
[Roxa,Boxa,Moxa,Soxa] = repbiodeg(Blank,Oxa,ToxaEx);           % Without filter
[Roxa,BoxaSm,MoxaSm,SoxaSm] = repbiodegsm(Blank,Oxa,ToxaEx);   % With filter

figure(9)
[p5] = stdplot(Xp,Moxa,Soxa,[0.5,0,0.5],[0.5,0,0.5],[0.5,0,0.5],0.3);

figure(10)
[p5sm] = stdplot(Xp,MoxaSm,SoxaSm,[0.5,0,0.5],[0.5,0,0.5],[0.5,0,0.5],0.3);
%% PET 
% Extract PET
Pet = Acc(:,24:26);

% Accumulated CO2 mean & std
% [APetM,APetS] = accmeanstd(Pet);
% figure(8) % Figure 8 from here out is Accumulated CO2
% [p6A] = stdplot(Xp,APetM,APetS,[0.5 0.5 0],[0.5 0.5 0],[0.5 0.5 0],0.15);

% Respiration & Biodegradation
[Rpet,Bpet,Mpet,Spet] = repbiodeg(Blank,Pet,Tpet);           % Without filter
[Rpet,BpetSm,MpetSm,SpetSm] = repbiodegsm(Blank,Pet,Tpet);   % With filter

figure(9)
[p6] = stdplot(Xp,Mpet,Spet,'k','k','k',0.3);

figure(11)
[p6sm] = stdplot(Xp,MpetSm,SpetSm,'k','k','k',0.3);
%% Monomers mixture 
% Etract Monomers mixture
Mix = Acc(:,27:29);

% Accumulated CO2 mean & std
% [AMixM,AMixS] = accmeanstd(Mix);
% figure(8) % Figure 8 from here out is Accumulated CO2
% [p7A] = stdplot(Xp,AMixM,AMixS,'c','c','c',0.15);

% Respiration & Biodegradation
[Rmix,Bmix,Mmix,Smix] = repbiodeg(Blank,Mix,Tmix);           % Without filter
[Rmix,BmixSm,MmixSm,SmixSm] = repbiodegsm(Blank,Mix,Tmix);   % With filter

figure(9)
[p7] = stdplot(Xp,Mmix,Smix,'c','c','c',0.3);

figure(11)
[p7sm] = stdplot(Xp,MmixSm,SmixSm,'c','c','c',0.3);

%% FDCA
Xp = 1:834; % Adjust horizontal points for shorter time

% Extract FDCA
Fdca = Acc(1:834,36:38); % Exclude ending NAN values to secure correct std fill
BlankSh = Blank(1:834); % Adjust Blank to correct time/length

% Accumulated CO2 mean & std
% [AFdcaM,AFdcaS] = accmeanstd(Fdca);
% figure(8) % Figure 8 from here out is Accumulated CO2
% [p9A] = stdplot(Xp,AFdcaM,AFdcaS,'b','b','b',0.3);

% Respiration & Biodegradation
[Rfdca,Bfdca,Mfdca,Sfdca] = repbiodeg(BlankSh,Fdca,Tfdca);           % Without filter
[Rfdca,BfdcaSm,MfdcaSm,SfdcaSm] = repbiodegsm(BlankSh,Fdca,Tfdca);   % With filter

figure(9)
[p8] = stdplot(Xp,Mfdca,Sfdca,[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5],0.3);

% figure(10)
% [p8sm] = stdplot(Xp,MfdcaSm,SfdcaSm,[0.5 0.5 0.5],[0.5 0.5 0.5],[0.5 0.5 0.5],0.3);      
%% NAOX
% Respiration
Naox = Acc(1:834,39:40); % Exclude ending NAN values to secure correct std fill

% Accumulated CO2 mean & std
% [ANaoxM,ANaoxS] = accmeanstd(Naox);
figure(8) % Figure 8 from here out is Accumulated CO2
% [p10A] = stdplot(Xp,ANaoxM,ANaoxS,'k','k','k',0.3);
title('Accumulated CO2')
ylabel('Accumulated CO2 (mg)')
ylim([0 120])
yticks([0:20:120])
legend([pBA p3A p4A],{'Blank','Iso-Isosorbide','Iso-PISOX'}) 
grid minor
ax = gca;
set(gca,'FontSize',30)

% Respiration & Biodegradation
[Rnaox,Bnaox,Mnaox,Snaox] = repbiodeg(BlankSh,Naox,Tnaox);           % Without filter
[Rnaox,BnaoxSm,MnaoxSm,SnaoxSm] = repbiodegsm(BlankSh,Naox,Tnaox);   % With filter

figure(9)
[p9] = stdplot(Xp,Mnaox,Snaox,[0.5 0.5 0],[0.5 0.5 0],[0.5 0.5 0],0.3);

figure(10)
[p9sm] = stdplot(Xp,MnaoxSm,SnaoxSm,[0.5 0.5 0],[0.5 0.5 0],[0.5 0.5 0],0.3);  

%% PISOX 
Xp = 1:1072; % Adjust length to Pisox again
% Extract PISOX
Pis = Acc(:,30:35);

dlmwrite('PisoxAcc.txt',Pis,'delimiter','\t','newline','pc')

% Accumulated CO2 mean & std
% [APisM,APisS] = accmeanstd(Pis);
% figure(8) % Figure 8 from here out is Accumulated CO2
% [p8A] = stdplot(Xp,APisM,APisS,'m','m','m',0.15);

% Respiration & Biodegradation
[Rpis,Bpis,Mpis,Spis] = repbiodeg(Blank,Pis,Tpis);           % Without filter
[Rpis,BpisSm,MpisSm,SpisSm] = repbiodegsm(Blank,Pis,Tpis);   % With filter

figure(9)
[p10] = stdplot(Xp,Mpis,Spis,'m','m','m',0.3);
legend([p1 p2 p3 p4 p5 p6 p7 p8 p9 p10],{'Cellulose','Hexanediol low','Iso then Iso','Iso then Pisox'...
    'Oxalic low','PET','Monomers mixture','FDCA','NAOX','Pisox'})
grid minor

figure(10)
[p10sm10] = stdplot(Xp,MpisSm,SpisSm,'m','m','m',0.3);
legend([p1sm10 p2sm p5sm p9sm p10sm10],{ 'Cellulose','Hexanediol Low','Oxalic Acid Low'...
    'NAOX','PISOX'},'FontSize',20) % change legend and p's accordingly, and scale
grid minor
ax = gca;
set(gca,'FontSize',30)
title('Biodegradation in Soil') 
yticks(-10:10:100);
ylim([-10 100]);

figure(11)
[p10sm11] = stdplot(Xp,MpisSm,SpisSm,'m','m','m',0.3);
legend([p1sm11 p6sm p7sm p10sm11],{ 'Cellulose','PET','Monomers Mix'...
    ,'PISOX'},'FontSize',20)
grid minor
ax = gca;
set(gca,'FontSize',30)
title('Biodegradation in Soil') 
yticks(-10:10:100);
ylim([-10 100]);

%% Overview
%close all % change this line to see all graphs
figure(12)
p = plot(BcelSm,'r');
hold on
p1 = plot(BfdcaSm,'b');
hold on
p2 = plot(BhexSm,'g');
hold on
p3 = plot(BisodSm,'y');
hold on
p4 = plot(BisopSm,'y--');
hold on
p5 = plot(BmixSm,'c');
hold on
p6 = plot(BnaoxSm,'k');
hold on
p7 = plot(BoxaSm,'Color',[0.5,0,0.5]);
hold on
p8 = plot(BpisSm,'m');

legend([p(1) p1(1) p2(1) p3(1) p4(1) p5(1) p6(1) p7(1) p8(1)],{'Cellulose','FDCA','Hexanediol Low','ISO then ISO','ISO then PISOX'...
    ,'Monomers mixture','NAOX','Oxalic Low','PISOX'})
title('Biodegradation per group')
xlabel('Time (h)')
ylabel('Biodegradation %')
yticks(-10:10:120)
xlim([0 1100])
ylim([-10 120])
grid on
grid minor


