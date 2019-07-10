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
    UE{ nday } = tmp.combinedData.UE;
    barOn{ nday } = UE{ nday }.barOn - UE{ nday }.fixOn;
end

%%
baseline = {};
for nday = 1:numel( datasets )
    dayBaseline = NaN(size(olapChopped{nday}.r.r(1).spikes, 1), numel(olapChopped{nday}.r.r));
    for itr = 1:numel(olapChopped{nday}.r.r)
        dayBaseline(:, itr) = 1000*mean(olapChopped{nday}.r.r(itr).spikes(:, (barOn{nday}(itr) - 500):barOn{nday}(itr)),2);
    end
    baseline{nday} = dayBaseline;
    clear dayBaseline;
end

%%
saveDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/pulvinar/multi-unit/continuous/withExternalInput_withLag/driftingTest';

for nday = 1:numel( datasets )
    daySaveDir = fullfile(saveDir, datasets(nday).shortName);
    if ~isdir(daySaveDir)
        mkdir(daySaveDir);
    end
    cd(daySaveDir)

    for n = 1: size(baseline{nday}, 1)
        x = 1: size(baseline{nday}, 2);
        f1 = figure
        scatter(x, baseline{nday}(n,:));
        xlabel('Trial')
        ylabel('Firing rate (spikes/s)')
        title(['MU ' int2str(n+4)])
        print(f1, ['MU ' int2str(n+4)], '-dpng');
        close;
    end
end


