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
end


baseWindow = [1 50];
leftB = 400:10:1200;

