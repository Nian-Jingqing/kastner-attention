%% add path

addpath('/snel/home/fzhu23/bin/LFADS/lfads-run-manager/src')
addpath('/snel/home/fzhu23/bin/analysis_tools')
addpath('/snel/home/fzhu23/Projects/Pulvinar/codes')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools/eeglab14_1_1b/')
addpath('/snel/home/fzhu23/Projects/Pulvinar/myTools/')

%% set directory

myFolder = '/snel/share/share/derived/kastner/data_processed/pulvinar/lfp/170127/';


%% load data
fileName = '170127_cueOnArrayOnTargetDim_HoldRel_lfp.mat';
fullFileName = fullfile(myFolder, fileName);
data = load(fullFileName);
R = data.R;

%% filter the LFP data
filtHighCutoff = 15;
filtLowCutoff = 2;
Fs = 1000;

R = bandpassFilter( R, 'lfp', 'lfp_theta', filtHighCutoff, filtLowCutoff, Fs );

%eegplot(d, 'srate', 1000, 'winlength', 1.5)

%% set names and chop times
chopStart = 301;
chopEnd = 1100;
savedirBase = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withGoodNeurons_PBTRun_20180403/lfp/'];
condNames = {'CueOnsetCueLoc1', 'CueOnsetCueLoc3', 'ArrayOnsetCueLoc1_hold', 'ArrayOnsetCueLoc3_hold', 'ArrayOnsetCueLoc1_rel', 'ArrayOnsetCueLoc3_rel', 'TargetDimCueLoc1', 'TargetDimCueLoc3' };
channelLabel = cell(1, size(R(1).lfp, 1));
channelVector = size(R(1).lfp, 1) : -1 : 1;
for c = 1:size(R(1).lfp, 1)
    channelLabel{c} = int2str(channelVector(c));
end


 %% select condition type and get trials in that condition
cond = 1;
R_selected = R([data.R.cueLoc] == 1 & [data.R.alignType] == 1);

%% loop over first 40 trials and plots
for t = 1:40
    lfpThisTrial = R_selected(t).lfp_theta;
%     f1 = figure;
    chopped_lfp = lfpThisTrial(:,chopStart:chopEnd)';
    dis = 0.5*(max(chopped_lfp(:)) - min(chopped_lfp(:)));
    base = size(chopped_lfp, 2):-1:1;
    riseBy = base*dis;
    scrolled_lfp = chopped_lfp + riseBy;
    plot(scrolled_lfp, 'b');
    set(gca,'YTick', fliplr(riseBy));
    set(gca,'YTickLabels',channelLabel);
    
    
    
    eegplot(lfpThisTrial(:,chopStart:chopEnd), 'srate', 1000, 'winlength', 1, 'title', condNames{cond});
%     title(condNames{cond});
    savedirOne = fullfile(savedirBase, condNames{cond});
    cd(savedirOne)
    print(gca, ['Trial ' int2str(t)], '-dpng');
    close;
end
    
%% clear set
for condType = 1:length(condNames)
    trialsForThisCond = r_lfadsWhole.r(condIx == condType);
    trialsOfRealForThisCond = rbinned(condIx == condType);
    nTrialsThisCond = length(trialsForThisCond);
    for t = 1:40
        lfadsRatesThisTrial = trialsForThisCond(t).rates;
        inputThisTrial = trialsForThisCond(t).controller_outputs;
        realSpikingThisTrial = trialsOfRealForThisCond(t).spikeCounts;
        f1 = figure;
%         for i = 1:nInputs
%             sp(i) = subplot(7, 2, i*2 - 1);
%             plot(inputThisTrial(i, chopStart:chopEnd), 'b');
%             set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
%             set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
%             set(sp(i), 'FontSize', 7);
%             xlim([0 div*nTimesLFADS]);
%             title(sp(i), ['Input ' int2str(i)]);
%         end
%         sp(7) = subplot(7,2,13);
%         sp(7) = subplot('position', [0.1300 0.0339 0.3347 0.150]);
        sp(1) = subplot(2, 1, 1)
        imagesc( lfadsRatesThisTrial(:, chopStart:chopEnd) );
        set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
        set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
        title(sp(1), 'LFADS rates');
        ylabel('units');
%         set(sp(7), 'FontSize', 7);
%         set(sp(7), 'position', [0.1300 0.1039 0.7750 0.1])
%         sp(8) = subplot(7,2,14);
%         sp(8) = subplot('position', [0.5703 0.0339 0.3347 0.150]);
        sp(2) = subplot(2, 1, 2)
        imagesc( realSpikingThisTrial(:, chopStart:chopEnd) );
        set(gca,'XTick',[1 div*0.5*nTimesLFADS div*nTimesLFADS]);
        set(gca,'XTickLabels',{'-400','AlignedTime','+400'});
        title(sp(2), 'Real Spiking');
        ylabel('units');
%         set(sp(8), 'FontSize', 7);
        
        suptitle(condNames{condType});
%         set(f1, 'Position', [279 53 648 913]);
%         set(f1, 'Position', [400 38 1219 928]);
        savedirOne = fullfile(savedirBase, condNames{condType});
        cd(savedirOne)
        print(f1, ['Trial ' int2str(t)], '-dpng');
        close;
    end
end



%%
savedir = ['/snel/share/share/derived/kastner/LFADS_runs/pulvinar/Multi-day/' ...
    'multiDay_CO_AO_TD_HoldRel_JanToApr/postAnalysis/withGoodNeurons_PBTRun_20180403/lfp/cueOnsetCueLoc1/'];

for i = 1:length(R_selected)
    f1 = eegplot(R_selected(i).lfp(:,chopStart:chopEnd), 'srate', 1000, 'winlength', 1);
    title(['Trial ' int2str(i)]);
    print(f1, ['Trial ' int2str(i)], '-dpng');
end
%%


