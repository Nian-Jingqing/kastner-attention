%% datasets
%datasets(1).shortName = '0218';
%datasets(1).midName = '021819';
%datasets(1).longName = '02182019';
%datasets(2).shortName = '0306';
%datasets(2).midName = '030619';
%datasets(2).longName = '03062019';
%datasets(3).shortName = '0311';
%datasets(3).midName = '031119';
%datasets(3).longName = '03112019';
%datasets(4).shortName = '0314';
%datasets(4).midName = '031419';
%datasets(4).longName = '03142019';
%datasets(5).shortName = '0406';
%datasets(5).midName = '040619';
%datasets(5).longName = '04062019';
%datasets(6).shortName = '0425';
%datasets(6).midName = '042519';
%datasets(6).longName = '04252019';
%datasets(7).shortName = '0502';
%datasets(7).midName = '050219';
%datasets(7).longName = '05022019';

dataset(1).date = '02182019';
dataset(2).date = '03062019';
dataset(3).date = '03112019';
dataset(4).date = '03142019';
%dataset(5).date = '03272019'; % previously didn't work
dataset(5).date = '04062019';
dataset(6).date = '04252019';
dataset(7).date = '05022019';
dataset(8).date = '02082019';
dataset(9).date = '02132019';
dataset(10).date = '02142019';
dataset(11).date = '02152019';
dataset(12).date = '02162019';
dataset(13).date = '02262019';
dataset(14).date = '02282019';
dataset(15).date = '03012019';
dataset(16).date = '03022019';
dataset(17).date = '03032019';
dataset(18).date = '03312019';
dataset(19).date = '04012019';


%% good channel info
good_channels{1} = [8, 22, 23, 24, 25, 26, 27, 28, 29, 31, 32];
good_channels{2} = [22, 23, 24, 26, 27, 28, 30, 31];
good_channels{3} = [14, 24, 25, 26, 30, 31];
good_channels{4} = [8, 21, 24, 25, 26, 27, 31, 32];
good_channels{5} = [17, 18, 19, 21, 23];
good_channels{6} = [9, 10, 21, 30, 31, 32];
good_channels{7} = [17, 24];
good_channels{8} = [11, 14, 15, 18, 19, 20, 21, 22, 23, 24, 25, 26];
good_channels{9} = [4, 5, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 26, 27, 28, 29, 30, 32];
good_channels{10} = [6, 15, 17, 18, 19, 22, 23, 24, 28, 29, 30, 31, 32];
good_channels{11} = [13, 16, 18, 21, 23, 25, 27, 29, 32];
good_channels{12} = [9, 11, 12, 15, 16, 19, 21, 23, 24, 29];
good_channels{13} = [11, 13, 14, 15, 16, 17, 18, 19, 20, 21, 25, 26, 27, 28, 29];
good_channels{14} = [13, 14, 15, 16, 23, 24, 25, 30, 32];
good_channels{15} = [10, 20, 30, 31];
good_channels{16} = [21, 24, 27, 28, 30, 32];
good_channels{17} = [18, 25, 27, 29, 30, 31];
good_channels{18} = [15 16 21 22 23 29 32]; % decided to get rid of 22
good_channels{19} = [6, 7, 13, 14, 15, 19]; % decided to get rid of 14



%% load data
%loadpath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/withoutDataMasking/higherThreshold/';
loadpath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/spiking_data/withExtInp_lowerThresh/';

% iterate over days, load each day and add it to a 'olapChopped' cell array
clear tc_data
%for nday = 1:numel( datasets )
%for nday = 1:numel(dataset)
for nday = 19
    disp( sprintf( 'loading tc day %g / %g', nday, numel( dataset ) ) );
    fileName = ['LIP_spiking_r_', dataset(nday).date, '.mat'];
    loadFile = fullfile(loadpath, fileName);
    %fname = sprintf( '%s%s_v1.mat', loadpath, dataset( nday ).date );
    %tmp = load( fname );
    %tc_data{ nday } = tmp.combinedData;
    tc_data{ nday } = load(loadFile);
end

%% get rid of bad channels first (to speed up the analysis below)
%for nday = 1:numel(dataset)
for nday = 19
    for itrial = 1 : numel(tc_data{nday}.r.r)
        tc_data{nday}.r.r(itrial).spikes_selected = tc_data{nday}.r.r(itrial).spikes(good_channels{nday}, :);
    end
end



%%
%for nday = 1:numel(dataset)
for nday = 19
    allSpikes{nday} = [tc_data{nday}.r.r.spikes_selected];
    sequence_length = size(allSpikes{nday}, 2);
    %corr_matrix
    sX = full(allSpikes{nday}(:, 0.2*sequence_length:0.5*sequence_length));
    corr_r = zeros(1,size(sX,1)^2);

    tic;
    x = size(sX,1)^2;
    t = size(sX,1);
    idx_par = zeros(2, x);

    parfor i=1:x
        idx1 = ceil(i/t);
        idx2 = mod(i,t)+1;
        if idx1 < idx2
            corr_r(i+1) = xcorr(sX(idx1,:), sX(idx2,:), 0, 'coeff');
        else
            corr_r(i+1) = nan;
        end
        idx_par(:, i + 1) = [idx1;idx2];
        disp( sprintf( 'running %g / %g', i, x ) );
    end
    toc;
    %corr_r = corr_r(2:end);
    %idx_par = idx_par(:, 2:end);
    %corr_r = reshape(corr_r, [size(sX,1) size(sX,1)]);
    %for j = 1:size(sX,1)
    %    tmp2end = corr_r(1:end-1, j);
    %    corr_r(1,j) = corr_r(end,j);
    %    corr_r(2:end, j) = tmp2end;
    %end
    corr_days{nday} = corr_r;
    corrIdx_days{nday} = idx_par;
    %correlationMatrices{nday} = corrcoef(allSpikes{nday});

    % formatting correlation matrix
    corr_r = corr_r(2:end);
    idx_par = idx_par(:, 2:end);
    corr_r = reshape(corr_r, [size(sX,1) size(sX,1)]);
    for j = 1:size(sX,1)
        tmp2end = corr_r(1:end-1, j);
        corr_r(1,j) = corr_r(end,j);
        corr_r(2:end, j) = tmp2end;
    end
    corr_matrix{nday} = corr_r;
    corrIdx{nday} = idx_par;
    
end


%%
%saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/highCorrCheck/higherThreshold/'
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/highCorrCheck/newSP/onlyGoodChannels/';
if ~isdir(saveDir)
    mkdir(saveDir)
end

cd(saveDir)
%%
for nday = 19
    f1 = figure
    imagesc(corr_matrix{nday})
    colorbar
    title(['MU-correlations-day ' int2str(nday)])
    print(f1, ['MU_correlations_day ' int2str(nday)], '-dpng');
end

%%
f1 = figure
isub = 1;
%for nday = [1, 2, 4, 8, 9, 10, 11, 13, 14, 17]
for nday = 19
    %subplot(2, 5, isub);
    histogram(corr_days{nday});
    title(['MU-correlations-day ' int2str(nday)])
    %isub = isub + 1;
end
set(gcf, 'Position', [86, 53, 1612, 857])
%print(f1, 'Distribution_allDays', '-dpng');
print(f1, 'Distribution_day03312019', '-dpng');

%% find high corr channels:
day = 19;
high_corr_idx = find(corr_days{day} > 0.05);
which_pair = corrIdx_days{day}(:, high_corr_idx);

%
% 02182019
%hc_channels{1} = [24, 25, 26, 27];
%hc_channels{2} = [23, 24, 25, 27, 28];
%hc_channels{3} = [8, 9];
%hc_channels{4} = [];
%hc_channels{5} = [18, 19];
%hc_channels{6} = [];
%hc_channels{7} = [31, 32];

% find which channels have higher firing rate
meanFR = mean(allSpikes{day}(:, 1:round(0.5*end)), 2);

%%
% 02182019
rm_channels{1} = [24, 26];
rm_channels{2} = [24, 25, 27];
rm_channels{3} = [9];
rm_channels{4} = [];
rm_channels{5} = [18];
rm_channels{6} = [];
rm_channels{7} = [32];





%% clear up data and save
savepath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
            'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/180614data_rm_highCorr/'];
cd(savepath)
%%
for itr = 1:numel(olapChopped{2}.R)
    olapChopped{2}.R(itr).spikeCounts(16,:) = [];
    olapChopped{2}.R(itr).rfloc(16,:) = [];
end

for itr = 1:numel(olapChopped{2}.r.r)
    olapChopped{2}.r.r(itr).spikes(16,:) = [];
    olapChopped{2}.r.r(itr).spikes_smoothed(16,:) = [];
    olapChopped{2}.r.r(itr).rfloc(16,:) = [];
end

for itr = 1:numel(olapChopped{5}.R)
    olapChopped{5}.R(itr).spikeCounts([12, 14, 27],:) = [];
    olapChopped{5}.R(itr).rfloc([12, 14, 27],:) = [];
end

for itr = 1:numel(olapChopped{5}.r.r)
    olapChopped{5}.r.r(itr).spikes([12, 14, 27],:) = [];
    olapChopped{5}.r.r(itr).spikes_smoothed([12, 14, 27],:) = [];
    olapChopped{5}.r.r(itr).rfloc([12, 14, 27],:) = [];
end

for itr = 1:numel(olapChopped{6}.R)
    olapChopped{6}.R(itr).spikeCounts([18, 24],:) = [];
    olapChopped{6}.R(itr).rfloc([18, 24],:) = [];
end

for itr = 1:numel(olapChopped{6}.r.r)
    olapChopped{6}.r.r(itr).spikes([18, 24],:) = [];
    olapChopped{6}.r.r(itr).spikes_smoothed([18, 24],:) = [];
    olapChopped{6}.r.r(itr).rfloc([18, 24],:) = [];
end

%%
for nday = [2 5 6]
    clear combinedData
    combinedData = olapChopped{nday};
    saveName = [datasets(nday).shortName, '_cueOnArrayOnTargetDim_HoldRel'];
    save(saveName, 'combinedData');
end