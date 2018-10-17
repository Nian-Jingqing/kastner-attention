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

%% load LFP data for 2 locations
lfpBasePath = '/snel/share/share/derived/kastner/data_processed/pulvinar/lfp/continuousOverlapChop/';
for nday = 1:numel( datasets )
    disp( sprintf( 'loading LFPs day %g / %g', nday, numel( datasets ) ) );
    lfpPath = fullfile( lfpBasePath, datasets( nday ).shortName );
    fname = sprintf( '%s/%s_cueOnArrayOnTargetDim_HoldRel_lfp.mat', lfpPath, datasets( nday ).shortName );
    tmp = load( fname );
    lfp{ nday } = tmp.r.r;
end

%%
for nday = 1:numel( datasets )
    for ntr = 1:numel(lfp{ nday })
        lfp{ nday }(ntr).lfps = lfp{ nday }(ntr).lfps(1:end - 1, :);
    end
end

%%
%ind = 0;
%for nday = 1:numel( datasets )
%    outlierOrNot = [lfp{ nday }.r.r.isTrialOutlier];
%    nonOutlierInds = find(~outlierOrNot);
%    for ntr = 1 : numel( nonOutlierInds )
%        ind = nonOutlierInds(ntr);
%        lfp{ nday }.r.r( ind ).lfps