%% list the dates for all needed datasets
datasets(1).shortName = '020819';
datasets(1).longName = '20190208';
datasets(2).shortName = '021819';
datasets(2).longName = '20190218';
datasets(3).shortName = '022619';
datasets(3).longName = '20190226';
datasets(4).shortName = '022719';
datasets(4).longName = '20190227';
datasets(5).shortName = '030819';
datasets(5).longName = '20190308';
datasets(6).shortName = '031019';
datasets(6).longName = '20190310';
datasets(7).shortName = '031119';
datasets(7).longName = '20190311';

%% good MU indices
mu_idx{1} = [5:20, 22:32] - 4;
mu_idx{2} = [5, 7:32] - 4;
mu_idx{3} = [6:22, 24:32] - 4;
mu_idx{4} = [5:28, 30:32] - 4;
mu_idx{5} = [5:32] - 4;
mu_idx{6} = [5:32] - 4;
mu_idx{7} = [5:32] - 4;

%% 
goodTrialIdx

%% % load all the chopped and recombined data (post-LFADS)
loadpath = ['/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/pulvinar/multi-unit/continuous/withExternalInput_withLag/'];

%%
% iterate over days, load each day and add it to a 'olapChopped' cell array
clear olapChopped
for nday = 1:numel( datasets )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_v1.mat', loadpath, datasets( nday ).shortName );
    tmp = load( fname );
    olapChopped{ nday } = tmp.combinedData;
    %UE{ nday } = tmp.combinedData.UE;
    %barOn{ nday } = UE{ nday }.barOn - UE{ nday }.fixOn;
end

%%
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/pulvinar/multi-unit/continuous/withExternalInput_withLag/nan_replaced/';
for nday = 1:numel(datasets)
    clear combinedData
%for nday = 1
    for itr = 1:numel(olapChopped{nday}.r.r)
        olapChopped{nday}.r.r(itr).spikes = olapChopped{nday}.r.r(itr).spikes(mu_idx{nday},:);
    end
    for n = 1:size(olapChopped{nday}.r.r(1).spikes,1)
        for itrial = 1:length(d{nday}{n})
            olapChopped{nday}.r.r(d{nday}{n}(itrial)).spikes(n,:) = NaN;
        end
    end
    combinedData.r = olapChopped{nday}.r;
    combinedData.UE = olapChopped{nday}.UE;
    cd(saveDir);
    saveName = [datasets(nday).shortName '_v3.mat'];
    save(saveName, 'combinedData')
end

    


