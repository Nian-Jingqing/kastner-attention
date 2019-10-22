%% load UEs
UE_baseDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/';
for nday = 1:numel(dataset)
    UE_file = ['UE_' dataset(nday).date '.mat'];
    load_UE_path = fullfile(UE_baseDir, UE_file);
    clear tmp
    tmp = load(load_UE_path);
    UE{ nday } = tmp.UE;
end

%% load and smooth/rebin spikes
baseDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/notch_filtering/notchFilterPlusBandPass/spiking_data/'; 
alf = {};
for nday = 1:numel(dataset)
    fileName = ['LIP_spiking_r_', dataset(nday).date, '.mat'];
    loadFile = fullfile(baseDir, fileName);
    disp( sprintf( 'loading day %g / %g', nday, numel( datasets ) ) );
    load(loadFile);
    %keyboard
    tc_r = R.Rstruct(r.r);
    % smooth and rebin
    sigma = 100;
    binsize_rescaled = 10;
    %
    tc_r.smoothFieldInR('spikes', 'spike_smoothed', sigma, 1);
    tc_rebinned = tc_r.binData({'spike_smoothed'}, [binsize_rescaled]);

    barStart = UE{nday}.barOn - UE{nday}.fixOn;
    cueStart = UE{nday}.cueOn - UE{nday}.fixOn;
    targetStart = UE{nday}.targetOn - UE{nday}.fixOn;
    %tc = [];
    for ntr = 1:numel(tc_rebinned)
        alf{nday}(ntr).spikes = tc_rebinned(ntr).spike_smoothed(good_channels{nday}, :);
        alf{nday}(ntr).barOnset = round( barStart( ntr ) /binsize_rescaled);
        alf{nday}(ntr).cueOnset = round( cueStart( ntr ) / binsize_rescaled);
        alf{nday}(ntr).targetStart = round( targetStart( ntr ) / binsize_rescaled);
    end
end