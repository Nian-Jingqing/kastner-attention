%% add paths here
addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/kastner_analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools')

addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/myTools/jPCA_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/old_pulvinarRepo/Kastner_Attention/codes/postAnalysisCodes/Manoj')

%%
buildRuns_20190226

%%
loadForPostAnalysis

%%

events = {'barOn', 'cueOn', 'targetOn'};
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};
window = round([-400  400] / binsize_rescaled);
timePoints = window(1):window(2);

%%
nCue = 1;
nBar = 1;
cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
trialIndicesPerCond(1).cueOn.(cue_cond_field) = (UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue);

%%
keepTrials_struct = alf{1}(trialIndicesPerCond(1).cueOn.(cue_cond_field));

for itrial = 1 : numel(keepTrials_struct)
    neuronData(itrial).forSpectrum = keepTrials_struct(itrial).spikes(22,:);
end

%% do pwelch
f2 = figure
cla
hold on
for itrial=1:20
    fs = 1000;
    %itrial = 30;
    tmp = neuronData(itrial).forSpectrum;
    tmp = tmp(~isnan(tmp));

    tmp = tmp -mean(tmp);
    [pxx,f] = pwelch(tmp,2^10,2^9,[],fs);

    %

    plot(f,10*log10(pxx))
    xlabel('Frequency (Hz)')
    ylabel('PSD (dB/Hz)')
end

%%
savedir = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/pulvinar/180208/postAnalysis/withExternalInput_20190226/spectrums/';
if ~isdir(savedir)
    mkdir(savedir)
end
cd(savedir);
name = 'MU-22-smoothedSpikes5ms_zoomedMore';
title(name)
print(f2, name, '-dpng');
