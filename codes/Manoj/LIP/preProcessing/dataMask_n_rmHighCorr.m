%% list the dates for all needed datasets
datasets(1).shortName = '0218';
datasets(1).midName = '021819';
datasets(1).longName = '02182019';
datasets(2).shortName = '0306';
datasets(2).midName = '030619';
datasets(2).longName = '03062019';
datasets(3).shortName = '0311';
datasets(3).midName = '031119';
datasets(3).longName = '03112019';
datasets(4).shortName = '0314';
datasets(4).midName = '031419';
datasets(4).longName = '03142019';
datasets(5).shortName = '0406';
datasets(5).midName = '040619';
datasets(5).longName = '04062019';
datasets(6).shortName = '0425';
datasets(6).midName = '042519';
datasets(6).longName = '04252019';
datasets(7).shortName = '0502';
datasets(7).midName = '050219';
datasets(7).longName = '05022019';

%% get drifting periods and bad channels due to drifting
drifty_info


%% bad channels due to high correlation
hc_channels{1} = [24, 26];
hc_channels{2} = [24, 25, 27];
hc_channels{3} = [9];
hc_channels{4} = [];
hc_channels{5} = [18];
hc_channels{6} = [];
hc_channels{7} = [32];

%% remove channels
clear rm_channels
allChannels = 1:32;
for i = 1:numel(datasets)
    mu_idx{i} = 1:32;
    rm_channels{i} = unique([hc_channels{i}, drift_channels{i}]);
    mu_idx{i}(ismember(allChannels, rm_channels{i})) = [];
end


%% % load all the chopped and recombined data (post-LFADS)
loadpath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/withoutDataMasking/selected/';

%%
% iterate over days, load each day and add it to a 'olapChopped' cell array
clear tc_origin
for nday = 1:numel( datasets )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_v1.mat', loadpath, datasets( nday ).midName );
    tmp = load( fname );
    tc_origin{ nday } = tmp.combinedData;
    %UE{ nday } = tmp.combinedData.UE;
    %barOn{ nday } = UE{ nday }.barOn - UE{ nday }.fixOn;
end

%%
%saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/data_masked_highCorr_rm/selected/';
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/highCorr_rm_noMasking/selected/';
if ~isdir(saveDir)
    mkdir(saveDir)
end


for nday = 1:7
    %    tmp = tc_origin{nday}.r.copy();
    clear combinedData
    %for nday = 1
    for itr = 1:numel(tc_origin{nday}.r.r)
        tc_origin{nday}.r.r(itr).spikes = tc_origin{nday}.r.r(itr).spikes(mu_idx{nday},:);
    end

    % data masking
    %    for n = 1:size(tc_origin{nday}.r.r(1).spikes,1)
    %    for itrial = 1:length(d{nday}{mu_idx{nday}(n)})
    %        tc_origin{nday}.r.r(d{nday}{mu_idx{nday}(n)}(itrial)).spikes(n,:) = NaN;
    %    end
    %end
    combinedData.r = tc_origin{nday}.r;
    combinedData.UE = tc_origin{nday}.UE;
    cd(saveDir);
    saveName = [datasets(nday).shortName '_v2.mat'];
    save(saveName, 'combinedData', '-v7.3')
end

    


