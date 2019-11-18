%%
% basic PSTH save root
saveRoot = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/postAnalysis/tc_newSP_PBT_191107/PSTH/';
if ~isdir(saveRoot)
    mkdir(saveRoot);
end

events = {'barOn', 'cueOn', 'targetOn'};
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};
window = round([-300  500] / binsize_rescaled);
timePoints = window(1):window(2);
%%
for nday = 1 : numel(alf)
    for nBar = 1:numel(conditions.barOn)
        trialIndicesPerCond(nday).barOn.(conditions.barOn{nBar}) = UE{nday}.barType == nBar;
        % compute real psth
        [psth, raster_tensor, stderr] = preparePSTH_Manoj(alf{nday}, 'spikes', trialIndicesPerCond(nday).barOn.(conditions.barOn{nBar}), [alf{nday}.barOnset], timePoints, binsize_rescaled);
        PSTH(nday).barOn.(conditions.barOn{nBar}).real.raster_tensor = raster_tensor;
        PSTH(nday).barOn.(conditions.barOn{nBar}).real.psth = psth;
        PSTH(nday).barOn.(conditions.barOn{nBar}).real.stderr = stderr;

        % compute lfads psth
        [psth, raster_tensor, stderr] = preparePSTH_Manoj(alf{nday}, 'rates', trialIndicesPerCond(nday).barOn.(conditions.barOn{nBar}), [alf{nday}.barOnset], timePoints, binsize_rescaled);
        PSTH(nday).barOn.(conditions.barOn{nBar}).lfads.raster_tensor = raster_tensor;
        PSTH(nday).barOn.(conditions.barOn{nBar}).lfads.psth = psth;
        PSTH(nday).barOn.(conditions.barOn{nBar}).lfads.stderr = stderr;
        
        for nCue = 1 : numel(conditions.cueOn)
            cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
            trialIndicesPerCond(nday).cueOn.(cue_cond_field) = (UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue);

            % compute real psth
            [psth, raster_tensor, stderr] = preparePSTH_Manoj(alf{nday}, 'spikes', trialIndicesPerCond(nday).cueOn.(cue_cond_field), [alf{nday}.cueOnset], timePoints, binsize_rescaled);
            PSTH(nday).cueOn.(cue_cond_field).real.raster_tensor = raster_tensor;
            PSTH(nday).cueOn.(cue_cond_field).real.psth = psth;
            PSTH(nday).cueOn.(cue_cond_field).real.stderr = stderr;

            % compute lfads psth
            [psth, raster_tensor, stderr] = preparePSTH_Manoj(alf{nday}, 'rates', trialIndicesPerCond(nday).cueOn.(cue_cond_field), [alf{nday}.cueOnset], timePoints, binsize_rescaled);
            PSTH(nday).cueOn.(cue_cond_field).lfads.raster_tensor = raster_tensor;
            PSTH(nday).cueOn.(cue_cond_field).lfads.psth = psth;
            PSTH(nday).cueOn.(cue_cond_field).lfads.stderr = stderr;
        end
    end
end

%% calculate max value for each plot series
% need to do so


%% first let's plot barOn
%for nday = 1:numel(datasets)


%% plotting PSTHs
pre_post_times = [300, 500] / binsize_rescaled;
t = pre_post_times;
t = [-t(1):-1 0:t(2)] * binsize_rescaled; %t(end)=[];

for nday = 1:10
    %for nday = 1:numel(dataset)
    saveRootDay = fullfile(saveRoot, dataset(nday).date);
    if ~isdir(saveRootDay)
        mkdir(saveRootDay);
    end
    cd(saveRootDay)
    nNeurons = size(PSTH(nday).barOn.Vert.lfads.psth, 1);
    % let's plot cueOn

    allIndicators(1).cond = {'Exo_TL_Vert', 'Exo_BL_Vert', 'Exo_TR_Vert', 'Exo_BR_Vert'};
    allIndicators(2).cond = {'Exo_TL_Hori', 'Exo_BL_Hori', 'Exo_TR_Hori', 'Exo_BR_Hori'};
    cueCondType{1} = 'Vert';
    cueCondType{2} = 'Hori';
    %allIndicators(3).cond = {'Endo_TL_Vert', 'Endo_BL_Vert', 'Endo_TR_Vert', 'Endo_BR_Vert'};
    %allIndicators(4).cond = {'Endo_TL_Hori', 'Endo_BL_Hori', 'Endo_TR_Hori', 'Endo_BR_Hori'};
    for nIndicator = 1:2
        plotIndicator = allIndicators(nIndicator).cond;
        %    saveRootDay = fullfile(saveRoot, datasets(nday).shortName);
        cmap = lines();
        for n = 1:nNeurons
            f2 = figure;
            y_max1 = max([PSTH(nday).cueOn.(plotIndicator{1}).lfads.psth(n, :), PSTH(nday).cueOn.(plotIndicator{2}).lfads.psth(n, :), PSTH(nday).cueOn.(plotIndicator{3}).lfads.psth(n, :), PSTH(nday).cueOn.(plotIndicator{4}).lfads.psth(n, :)]);
            y_max2 = max([PSTH(nday).cueOn.(plotIndicator{1}).real.psth(n, :), PSTH(nday).cueOn.(plotIndicator{2}).real.psth(n, :), PSTH(nday).cueOn.(plotIndicator{3}).real.psth(n, :), PSTH(nday).cueOn.(plotIndicator{4}).real.psth(n, :)]);
            y_max = max([y_max1, y_max2]);
            s1 = subplot(5,2,1);
            for j = 1:numel(plotIndicator)
                condStr = plotIndicator{j};
                whereUnderline = strfind(condStr, '_');
                condStr(whereUnderline) = '-';
                s_tmp = PSTH(nday).cueOn.(plotIndicator{j}).real.psth(n, :);
                s_tmp_upper = s_tmp + PSTH(nday).cueOn.(plotIndicator{j}).real.stderr(n, :);
                s_tmp_lower = s_tmp - PSTH(nday).cueOn.(plotIndicator{j}).real.stderr(n, :);             
                plot(s_tmp, 'Color', cmap(j,:), 'LineWidth', 2, 'DisplayName', condStr);
                hold on
                plot(s_tmp_upper, 'Color', cmap(j,:), 'LineWidth', 1);
                hold on
                plot(s_tmp_lower, 'Color', cmap(j,:), 'LineWidth', 1);
                hold on
            end
            set(gca,'XTick',[1 (3/8)*length(timePoints) length(timePoints)]);
            set(gca,'XTickLabels',{'-300ms', 'cueOnset' ,'+500ms'});
            title(s1, 'Real', 'FontSize', 16);
            ylabel('Firing Rate');
            axes(s1)
            ylim([0 y_max*(6/5)])
            xlim([0 81]) %need to change there
                         %xlim([0 41]) %need to change there
            set(gca, 'fontsize', 12);    

            s2 = subplot(5,2,2);
            for j = 1:numel(plotIndicator)
                condStr = plotIndicator{j};
                whereUnderline = strfind(condStr, '_');
                condStr(whereUnderline) = '-';
                r_tmp = PSTH(nday).cueOn.(plotIndicator{j}).lfads.psth(n, :);
                r_tmp_upper = r_tmp + PSTH(nday).cueOn.(plotIndicator{j}).lfads.stderr(n, :);
                r_tmp_lower = r_tmp - PSTH(nday).cueOn.(plotIndicator{j}).lfads.stderr(n, :);            
                h(j) = plot(r_tmp, 'Color', cmap(j,:), 'LineWidth', 2, 'DisplayName', condStr);
                hold on
                plot(r_tmp_upper, 'Color', cmap(j,:), 'LineWidth', 1);
                hold on
                plot(r_tmp_lower, 'Color', cmap(j,:), 'LineWidth', 1);
                hold on       
            end
            set(gca,'XTick',[1 (3/8)*length(timePoints) length(timePoints)]);
            set(gca,'XTickLabels',{'-300ms', 'cueOnset' ,'+500ms'});
            title(s2, 'LFADS', 'FontSize', 16);
            ylabel('Firing Rate');
            axes(s2)
            ylim([0 y_max*(6/5)])
            xlim([0 81]) % need to change here to not hard coded
            %xlim([0 41]) % need to change here to not hard coded
            set(gca, 'fontsize', 12);
            legend(h(1:4))
            set(legend, 'Position', [0.8254 0.8557 0.1048 0.1014]);

            for j = 1:numel(plotIndicator)

                % first, let's set up the title string
                condStr = plotIndicator{j};
                whereUnderline = strfind(condStr, '_');
                condStr(whereUnderline) = '-';
                
                subplot(5,2,2*j+1);
                a = squeeze(PSTH(nday).cueOn.(plotIndicator{j}).real.raster_tensor(n,:,:));
                imagesc(a)
                set(gca,'XTick',[1 (3/8)*length(timePoints) length(timePoints)]);
                set(gca,'XTickLabels',{'-300ms', 'cueOnset' ,'+500ms'});
                title(condStr);
                ylabel('Trials');
                set(gca, 'fontsize', 10);

                subplot(5,2,2*j+2);
                a = squeeze(PSTH(nday).cueOn.(plotIndicator{j}).lfads.raster_tensor(n,:,:));
                imagesc(a)
                set(gca,'XTick',[1 (3/8)*length(timePoints) length(timePoints)]);
                set(gca,'XTickLabels',{'-300ms', 'cueOnset' ,'+500ms'});
                title(condStr);
                ylabel('Trials');
                set(gca, 'fontsize', 10);    
            end

            suptitle(['Multi-unit ' int2str(alf{nday}(1).channel_info(n))]);
            set(f2, 'Position', [428 4 1193 962]);
            print(f2,['Multi-unit ' int2str(alf{nday}(1).channel_info(n)) ' - ' cueCondType{nIndicator}], '-dpng');

            %printpdf(f2,int2str(nIndices(n)) )
            close;
        end
    end
end
