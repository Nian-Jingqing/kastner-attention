%% list the dates for all needed datasets
datasets(1).shortName = '021819';
datasets(1).longName = '20190218';
datasets(2).shortName = '030619';
datasets(2).longName = '20190306';
datasets(3).shortName = '031119';
datasets(3).longName = '20190311';
datasets(4).shortName = '031419';
datasets(4).longName = '20190314';
datasets(5).shortName = '040619';
datasets(5).longName = '20190406';
datasets(6).shortName = '042519';
datasets(6).longName = '20190425';
datasets(7).shortName = '050219';
datasets(7).longName = '20190502';

%% % load all the chopped and recombined data (post-LFADS)
loadpath = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/thresholdCrossings/highCorr_rm_noMasking/preAlignedData_added/selected/';

%%
% iterate over days, load each day and add it to a 'olapChopped' cell array
clear olapChopped
for nday = 1:numel( datasets )
    disp( sprintf( 'loading chopped day %g / %g', nday, numel( datasets ) ) );
    fname = sprintf( '%s%s_v3.mat', loadpath, datasets( nday ).shortName );
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
    for nday = 1: numel( datasets )
        day_id = nday;
        disp( sprintf( 'getting LFADS firing rates day %g / %g', nday, numel( datasets ) ) );
        r_lfads2 = olapChopped{ nday }.r.get_output_from_lfads(run, day_id, trial_time_ms, trial_olap_ms, 'rates');
        assembled_lfads2 = R.Rstruct(r_lfads2);
        alf{ nday } = assembled_lfads2.r;
        
        disp( sprintf( 'getting factors day %g / %g', nday, numel( datasets ) ) );
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
overWriteRealOrNot = 1;
if overWriteRealOrNot
    for nday = 1: numel( datasets )
        % get all real data
        trueSpikes = [olapChopped{ nday }.r.r.spikes];
        nSpikes = size(trueSpikes, 2);
        out = olapChopped{ nday }.r.generate_overlap_chop_lfads_data( trial_time_ms, trial_olap_ms );
        tot_spikes = out.counts;
        % nPieces = size(tot_spikes, 1);
        %spikes = zeros(size(trueSpikes, 1), nSpikes);
        spikes(length(olapChopped{ nday }.r.r)).spikes = 1;
        %% assembly
        nChops = 0;
        for r_trial = 1: length(olapChopped{ nday }.r.r)
            spikes(r_trial).spikes = zeros(size(olapChopped{ nday }.r.r(r_trial).spikes));
            nPieces = ceil( (size( spikes(r_trial).spikes,2 ) - trial_time_ms ) / ( trial_time_ms - trial_olap_ms ) );
            for i = 1:nPieces
                spikesThisPiece = squeeze( tot_spikes( i+nChops, :, : ) );
                trialStart = 1 + (i-1)*trial_time_ms - (i-1)*trial_olap_ms;
                trialEnd = trialStart + trial_time_ms - 1;
                %         if i~= size( tot_spikes, 1 )
                %             trialStart = 1 + (i-1)*trial_time_ms - (i-1)*trial_olap_ms;
                %             trialEnd = trialStart + trial_time_ms - 1;
                %         else
                %             trialStart = nSpikes - trial_time_ms + 1;
                %             trialEnd = nSpikes;
                %         end
                spikes(r_trial).spikes( :, trialStart:trialEnd ) = spikesThisPiece;
            end
            nChops = nChops + nPieces;
        end
        assembled_real2 = R.Rstruct(spikes);
        if smoothOrNot
            assembled_real2.smoothFieldInR('spikes', 'spike_smoothed', sigma, 1);
            rbinned_real = assembled_real2.binData({'spike_smoothed'}, [binsize_rescaled]);
        else
            rbinned_real = assembled_real2.binData({'spikes'}, [binsize_rescaled]);
        end

        for ntr = 1:numel( alf{ nday } )
            if smoothOrNot
                alf{ nday }( ntr ).spikes = rbinned_real( ntr ).spike_smoothed;
            else            
                alf{ nday }( ntr ).spikes = rbinned_real( ntr ).spikes;
            end
        end
    end
end        

%% fix any weirdness with zeros in the ALF
for nd = 1:numel(alf)
    for ntr = 1:numel(alf{nd})
        alf{nd}(ntr).rates(alf{nd}(ntr).rates==0) = nan;
    end
end