classdef Dataset < LFADS.Dataset
% A single-day of raw data processed in a common way
    properties(SetAccess=protected)
        data
    end
    
    methods
        function ds = Dataset(collection, relPath, monkey)
            ds = ds@LFADS.Dataset(collection, relPath);

            % extract the dataset date
            putativeDatestring = relPath( 4:13 );
            
            % define the datestring based on filenames
            match = regexp(putativeDatestring, ['(?<datestr>\d{4}-\' ...
                                'd{2}-\d{2})'], 'names');
            % override the default name with the date
            ds.name = [collection.monkey '_' match.datestr ];
            ds.datenum  = datenum(match.datestr);
        end
        

        function data = loadData(ds) % overload the "loadData"
                                     % function defined in LFADS.Dataset
            in = load(ds.path);

            originalR=in.R;

            %% postprocess originalR to remove unhittable trials
            unusable_lidxs = [originalR.unhittable];    % Exclude trials based on the trial type, e.g. trials that can't be moved to.
            unusable_lidxs = unusable_lidxs | [originalR.photoBoxError];   % Exclude trials if the there was a photobox error.
            unusable_lidxs = unusable_lidxs | [originalR.possibleRTproblem];   % Exclude trials if there was a possible reaction
                                                                       % time issue.

            % note: original analysis (lfads_1) was run on 2435 trials.
            % subselection excluded: [unhittable | photoBoxError | possibleRTproblem]
            % DID NOT exclude [trialType<=0 | ~isConsistent] or unsuccessful trials
            % more recent analyses might exclude these trials...
            unusable_lidxs = unusable_lidxs | [originalR.trialType]<=0; % exclude trialtype==0
            unusable_lidxs = unusable_lidxs | [originalR.isConsistent]~=1; % exclude inconsistent trials
            
            s=[originalR.SUCCESS_SUMMARY]; % remove unsuccessful trials
            unusable_lidxs = unusable_lidxs | [s.errorNumber]~=0;

            originalR(unusable_lidxs)=[];
        
            data = [];
            fields2pass = {'blockID', 'trialID', 'runID', 'saveTag', ...
                          'offlineRT', 'offlineMoveOnsetTime', ...
                           'possibleRTproblem', ...
                           'actualFlyAppears', 'actualLandingTime'};

            nTrialsToProcess = numel( originalR );
            pbar = LFADS.Utils.ProgressBar(nTrialsToProcess, ...
                                           'Processing raw Rstruct');

            % for each trial
            for ntr = 1:nTrialsToProcess
                trialEnd = originalR(ntr).trialEndsTime;
                if isnan(trialEnd)
                    disp('oops');
                    keyboard
                end
                % create a spikeraster and store down
                % the hand position
                nUnits = numel(originalR(ntr).unit);
                spikes = sparse(zeros(nUnits, trialEnd));
                for nu = 1:nUnits
                    spikeTimes = originalR(ntr).unit(nu).spikeTimes;
                    %spikesMS = histc(spikeTimes, msedges);
                    spikes(nu, floor( spikeTimes( spikeTimes < trialEnd )) ...
                             + 1) = 1;
                end
                spikes = sparse(spikes);
                data(ntr).y = spikes;

                % save hand position
                data(ntr).X=[originalR(ntr).HAND.X(:,1)';...
                         originalR(ntr).HAND.Y(:,1)'];

                data(ntr).T = trialEnd;
                data(ntr).params.dtMS = 1;

                % calculate path length
                pdx = diff(originalR(ntr).HAND.X);
                pdy = diff(originalR(ntr).HAND.Y);
                pathdiff = [pdx(:)';
                            pdy(:)'];
                pd = sum( sum(pathdiff.^2) );
                data(ntr).pathLength = pd;

                % save the target position
                activeFly = 1;
                if originalR(ntr).numFlies > 1
                    activeFly = originalR(ntr).activeFly;
                end
                data(ntr).posTarget = [originalR(ntr).PARAMS.flyX(activeFly); originalR(ntr).PARAMS.flyY(activeFly)];
                data(ntr).allTargets = [originalR(ntr).PARAMS.flyX; originalR(ntr).PARAMS.flyY];
                data(ntr).activeFly = activeFly;
                data(ntr).whichFly = originalR(ntr).whichFly;

                % calculate and save the endpoint angle
                endpoint = data( ntr ).posTarget;
                data( ntr ).endpointAngle = angle( endpoint(1,:) + sqrt(-1)*endpoint(2,:) );


                
                % things to save down
                for nf = 1:numel(fields2pass)
                    if isfield(originalR(ntr), fields2pass{nf})
                        data(ntr).(fields2pass{nf}) = ...
                            originalR(ntr).(fields2pass{nf});
                    end
                end
                pbar.update(ntr);
            end
            disp(' ');

            % create 'conditionCode' - unique identifyer for each
            % trial type
            conditionCodes = 1000*[ originalR.trialType ] + [ originalR.trialVersion ];
            tmp2 = num2cell(conditionCodes);
            [data.conditionCode] = deal(tmp2{:});

            % single items to save down
            data(1).unitRatings = originalR(1).unitRatings;
            data(1).SaveInfo = originalR(1).SaveInfo;
            data(1).nUnits = nUnits;
            data(1).nTrials = numel(data);



            % now use the Rstruct class to store the data
            ds.data = R.Rstruct( data );
            data = ds.data;
        end

        function seq = makeColormap(ds, angleField)
            %% make a colormap
            if ~exist('angleField','var')
                angleField = 'endpointAngle';
            end
            c = hsv;
            for nt = 1:numel( ds.data.r )
                ang = ds.data.r(nt).( angleField );
                angind = (ang - -pi) / (2*pi);
                clrind = floor(angind * 64)+1;
                ds.data.r(nt).color = c(clrind,:);
                ds.data.r(nt).colorInd = clrind;
            end
        end

        function makeMeanTrajectories( ds, prePostWindow )
        % calculate the average trajectory for each condition
            if ~exist( 'prePostWindow', 'var')
                prePostWindow = [450 650];
            end
            avgStruct = ds.data.alignAndAverage( 'conditionCode', 'X', prePostWindow, 'offlineMoveOnsetTime' );
            
            for nc = 1:numel( avgStruct )
                theseTrials = find( [ ds.data.r.conditionCode ] == avgStruct( nc ).conditionCode );
                for it = 1:numel( theseTrials )
                    nt = theseTrials( it );
                    ds.data.r( nt ).meanTrajectory = avgStruct( nc ).X;
                end
            end
        end

        function plotTrajectories( ds, prePostWindow, showMarkers )
        % plot all the trajectories
            if ~exist( 'prePostWindow', 'var')
                prePostWindow = [450 650];
            end
            % check if colormap is defined
            if ~isfield( ds.data.r, 'color' )
                ds.makeColormap();
            end
            

            avgStruct = ds.data.alignAndAverage( 'conditionCode', 'X', prePostWindow, 'offlineMoveOnsetTime' );
            
            for nc = 1:numel( avgStruct )
                traj = avgStruct( nc ).X;

                cc = avgStruct( nc ).conditionCode;
                exampleTrial = find( [ds.data.r.conditionCode] == cc, 1 );
                clr = ds.data.r( exampleTrial ).color;
                posTarget = ds.data.r( exampleTrial ).posTarget;
                l = Plot.patchline(traj(1, :), traj(2, :), ...
                              'edgealpha',0.5);
                set(l, 'edgecolor', clr);
                set(l, 'linewidth', 2);
                hold on;
            end            
            
        end


        function [ excludedTrials, keptTrials ] = createAngleExclusion(ds, angleRange, data)
            assert( numel( angleRange ) == 2, 'must provide a high and low angle');
            assert( angleRange(2) > angleRange(1), 'must provide a high and low angle');

            if ~exist('data', 'var')
                data = ds.loadData();
            end


            endpoint = [ data.r.posTarget ];
            endang = angle( endpoint(1,:) + sqrt(-1)*endpoint(2,:) );
            
            % this is not needed angle(..) is always 0< < 2*pi
            %endang = mod(endang, 2 * pi - eps);

            keepTrials = true( size( data.r ) );

            trialsThisRange = find(endang >= angleRange(1) & endang < angleRange(2) );
            excludedTrials = trialsThisRange;

            keepTrials(trialsThisRange) = false;
            keptTrials = find( keepTrials );
        end

        function [ targetCurvatureId ] = groupByTargetAndCurvature( ds, data )
            if ~exist('data', 'var')
                data = ds.loadData();
            end
            %% a note on the maze data
            % the most complex dataset has 108 "conditions"
            % there are actually 36 different targets (36 * 3 = 108)
            % for each target, there are three conditions - one straight reach to the target, and two curved
            %  the only difference between the curved reaches is that there may or may not be a distractor targ
            %  this has virtually 0 effect on kinematics
            % therefore, we want to come out with 36 *2 = 72 conditions
            % to match the variety of kinematics present
            % so first group by target, then group by curved v. straight

            ccodes = [data.conditionCode];
            targetInd = floor(ccodes / 1000);
            % there should be 36 unique targets. if not, something's wrong
            assert( numel( unique(targetInd) ) == 36, 'error');

            isStraight = mod(ccodes, 1000) == 0;

            tcId = targetInd + sqrt(-1) * isStraight;

            [~, ~, targetCurvatureId] = unique( tcId );

        end



        function loadInfoFromData(ds, data)
            ds.subject = data(1).SaveInfo.monkey;
            assert( ~isempty(ds.datenum), ['this should already ' ...
                                'have been defined'] );
            % find save tag
            ds.saveTags = unique([data.saveTag]);

            ds.nChannels = data(1).nUnits;
            ds.nTrials = data(1).nTrials;
        end
    end
end
