%% add path

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')

%% set directory

myFolder = '/snel/share/share/data/kastner/pulvinar/lfp/M20170127/';
%Specify where to save the reorganized data
saveDir = '/snel/share/share/derived/kastner/data_processed/pulvinar/lfp/170127/';


%% load data
fileName = 'M20170127-g3-g4-g5-allFP-evokedLfps-v11.mat';
fullFileName = fullfile(myFolder, fileName);
data = load(fullFileName);
%% set directory

myFolder = '/snel/share/share/data/kastner/pulvinar/lfp/M20170127/';
%Specify where to save the reorganized data
saveDir = '/snel/share/share/derived/kastner/data_processed/pulvinar/lfp/170127/';


%% load data
fileName = 'M20170127-g3-g4-g5-allFP-evokedLfps-v11.mat';
fullFileName = fullfile(myFolder, fileName);
data = load(fullFileName);

%%

CO_tensor = data.cueOnsetLfp.lfp;
AO_tensor = data.arrayOnsetLfp.lfp;
TD_tensor = data.targetDimLfp.lfp;

%%

CO_struct(size(CO_tensor,2)).lfp = 1;
AO_struct(size(AO_tensor,2)).lfp = 1;
TD_struct(size(TD_tensor,2)).lfp = 1;

for trial = 1:size(CO_tensor,2)
    CO_struct(trial).lfp = squeeze(CO_tensor(:,trial,:));
    CO_struct(trial).cueLoc = data.UE.cueLoc(trial);
    CO_struct(trial).isHold = data.UE.isHoldTrial(trial);
    CO_struct(trial).alignType = 1;
end

for trial = 1:size(AO_tensor,2)
    AO_struct(trial).lfp = squeeze(AO_tensor(:,trial,:));
    AO_struct(trial).cueLoc = data.UE.cueLoc(trial);
    AO_struct(trial).isHold = data.UE.isHoldTrial(trial);
    AO_struct(trial).alignType = 2;
end

cueLocTD = data.UE.cueLoc(data.UE.isHoldTrial);
for trial = 1:size(TD_tensor,2)
    TD_struct(trial).lfp = squeeze(TD_tensor(:,trial,:));
    TD_struct(trial).cueLoc = cueLocTD(trial);
    TD_struct(trial).isHold = 1;
    TD_struct(trial).alignType = 3;
end

R = [CO_struct,AO_struct,TD_struct];
    
%% select the channels that correspond to the good neurons

nIndices = [1 3 4 6 7 8 11 13 14 15 16 17 19 20 21 22 23 24 25 26 29 30 31 32];
% 170127

for i = 1:length(R)
    R(i).lfp = R(i).lfp(nIndices,:);
    %R(i).rfloc = R(i).rfloc(nIndices,:);
end

%% save
cd(saveDir);
save('170127_cueOnArrayOnTargetDim_HoldRel_lfp.mat', 'R');