%% datasets
datasets(1).shortName = '170127';
datasets(1).longName = '20170127';
datasets(2).shortName = '170130';
datasets(2).longName = '20170130';
datasets(3).shortName = '170201';
datasets(3).longName = '20170201';
datasets(4).shortName = '170211';
datasets(4).longName = '20170211';
datasets(5).shortName = '170308';
datasets(5).longName = '20170308';
datasets(6).shortName = '170311';
datasets(6).longName = '20170311';

%% load data
loadpath = ['/snel/share/share/derived/kastner/data_processed/pulvinar/' ...
            'multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/'];

% iterate over days, load each day and add it to a 'olapChopped' cell array
clear olapChopped
for nday = 1:numel( datasets )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_cueOnArrayOnTargetDim_HoldRel.mat', savepath, datasets( nday ).shortName );

    tmp = load( fname );
    olapChopped{ nday } = tmp.combinedData;
end

%%
%for nday = 1:numel(olapChopped)
for nday = 6
    allSpikes{nday} = [olapChopped{nday}.r.r.spikes];
    sequence_length = size(allSpikes{nday}, 2);
    %corr_matrix
    sX = full(allSpikes{nday}(:, 1:0.5*sequence_length));
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
end

%%
saveDir = '/snel/share/share/derived/kastner/data_processed/pulvinar/multi-unit/continuousOverlapChop/multiDay_JanToMar/withExternalInput_withLag/';
cd(saveDir)
for nday = 1:numel(olapChopped)
    f1 = figure
    imagesc(corr_days{nday})
    colorbar
    print(f1, ['MU_correlations_day ' int2str(nday)], '-dpng');
end

%%
f1 = figure
for nday = 1:numel(olapChopped)
    subplot(2, 3, nday);
    histogram(corr_days{nday});
    title(['MU-correlations-day ' int2str(nday)])
end
set(gcf, 'Position', [86, 53, 1612, 857])
print(f1, 'Distribution_allDays', '-dpng');

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