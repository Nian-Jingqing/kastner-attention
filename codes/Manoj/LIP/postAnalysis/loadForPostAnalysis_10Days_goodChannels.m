%% list the dates for all needed datasets
dataset(1).date = '02082019';
dataset(2).date = '02132019';
dataset(3).date = '02142019';
dataset(4).date = '02152019';
dataset(5).date = '02182019';
dataset(6).date = '02262019';
dataset(7).date = '02282019';
dataset(8).date = '03032019';
dataset(9).date = '03062019';
dataset(10).date = '03142019';

%% % load all the chopped and recombined data (post-LFADS)
%loadpath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/newSignalProcessing/withPreAligned_noExtInp/';
loadpath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/newSignalProcessing/withPreAligned_withExtInp/lowerThresh/';

%%
% iterate over days, load each day and add it to a 'olapChopped' cell array
clear olapChopped
for nday = 1:numel( dataset )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( dataset ) ) );
    fname = sprintf( '%s%s_v2.mat', loadpath, dataset( nday ).date );
    %fname = sprintf( '%s%s_v1.mat', loadpath, dataset( nday ).date );
    tmp = load( fname );
    olapChopped{ nday } = tmp.combinedData;
    UE{ nday } = tmp.combinedData.UE;
end

%% get the LFADS factors and rates for each day
% load LFADS factors and rates from saved processed data or from posterior mean directly
loadFromLFADS = 1; % need to change this when you want to change datafile for loading
rebinOrNot = 0; % need to change to 1 if you want to rebin your LFADS factors and rates
rfactor = 1;
binsize = par.spikeBinMs;
binsize_rescaled = binsize * rfactor;
loadFromSavedDir = 'lfads'; % need to specify
r_id = 1;
run = rc2.runs(r_id);

use_rates = false; % need to figure out what this means
rawSampleRate = 1000;
trial_time_ms = 500;
trial_olap_ms = 100;

if loadFromLFADS
    for nday = 1: numel( dataset )
        day_id = nday;
        disp( sprintf( 'getting LFADS firing rates day %g / %g', nday, numel( dataset ) ) );
        r_lfads2 = olapChopped{ nday }.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms, 'rates');
        assembled_lfads2 = R.Rstruct(r_lfads2);
        alf{ nday } = assembled_lfads2.r;
        
        disp( sprintf( 'getting factors day %g / %g', nday, numel( dataset ) ) );
        r_lfads3 = olapChopped{ nday }.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms, 'factors');
        assembled_lfads3 = R.Rstruct(r_lfads3);
        %    alf{ nday } = assembled_lfads3.r;
      
        barStart = UE{ nday }.barOn - UE{ nday }.fixOn;
        cueStart = UE{ nday }.cueOn - UE{ nday }.fixOn;
        targetStart = UE{ nday }.targetOn - UE{ nday }.fixOn;
        
        for ntr = 1:numel( alf{ nday } )
            if use_rates
                alf{ nday }( ntr ).rates = log( alf{ nday }( ntr ).rates );
            end
            if rebinOrNot
                rdim = size( alf{ nday }(ntr).rates );
                tmp1 = alf{ nday }(ntr).rates(:, 1:floor( rdim( 2 ) / rfactor ) * rfactor );
                tmp1 = reshape( tmp1, rdim(1), rfactor, floor( rdim( 2 ) / rfactor ) );
                tmp1 = squeeze( mean( tmp1, 2 ) );
                alf{ nday }(ntr).rates = tmp1;
                
                rdim2 = size( assembled_lfads3.r( ntr ).rates );
                tmp2 = assembled_lfads3.r( ntr ).rates(:, 1:floor( rdim2( 2 ) / rfactor ) * rfactor );
                tmp2 = reshape( tmp2, rdim2(1), rfactor, floor( rdim2( 2 ) / rfactor ) );
                tmp2 = squeeze( mean( tmp2, 2 ) );
                alf{ nday }(ntr).factors = tmp2;
            else
                alf{ nday }( ntr ).factors = assembled_lfads3.r( ntr ).rates;
            end
            alf{ nday }( ntr ).barOnset = round( barStart( ntr ) / binsize_rescaled);
            alf{ nday }( ntr ).cueOnset = round( cueStart( ntr ) / binsize_rescaled);
            alf{ nday }( ntr ).targetStart = round( targetStart( ntr ) / binsize_rescaled);
        end
    end
else
    load(loadFromSavedDir);
end

%% get real spikes for each day, and smooth / rebin
smoothOrNot = 1;
sigma = 20;
for nday = 1: numel( dataset )
    % get all real data
    tc_r = R.Rstruct(olapChopped{ nday }.r.r);
    tc_r.smoothFieldInR('spikes', 'spike_smoothed', sigma, 1);
    tc_rebinned = tc_r.binData({'spike_smoothed'}, [binsize_rescaled]);
    for ntr = 1:numel(tc_rebinned)
        alf{nday}(ntr).spikes = tc_rebinned(ntr).spike_smoothed;
    end
end

%% fix any weirdness with zeros in the ALF
for nd = 1:numel(alf)
    for ntr = 1:numel(alf{nd})
        alf{nd}(ntr).rates(alf{nd}(ntr).rates==0) = nan;
    end
end

%% get info for good channels
for nday = 1:numel(dataset)
    alf{nday}(1).channel_info = olapChopped{nday}.channel_info;
end

