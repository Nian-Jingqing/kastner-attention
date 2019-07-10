
%%
buildRuns_20190411

%%
loadForPostAnalysis_20190411

%%
% basic PSTH save root
saveRoot = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/pulvinar/180208/postAnalysis/withExternalInput_20190411/PSTH/';
if ~isdir(saveRoot)
    mkdir(saveRoot);
end

events = {'barOn', 'cueOn', 'targetOn'};
conditions.barOn = {'Vert', 'Hori'};
conditions.cueOn = {'Exo_TL', 'Exo_BL', 'Exo_TR', 'Exo_BR', 'Endo_TL', 'Endo_BL', 'Endo_TR', 'Endo_BR'};
window = round([-400  400] / binsize_rescaled);
timePoints = window(1):window(2);
%%
for nday = 1 : numel(alf)
    for nBar = 1:numel(conditions.barOn)
        trialIndicesPerCond(nday).barOn.(conditions.barOn{nBar}) = UE{nday}.barType == nBar;
        % compute real psth
        [psth, raster_tensor] = preparePSTH_Manoj(alf{nday}, 'spikes', trialIndicesPerCond(nday).barOn.(conditions.barOn{nBar}), [alf{nday}.barOnset], timePoints, binsize_rescaled);
        PSTH(nday).barOn.(conditions.barOn{nBar}).real.raster_tensor = raster_tensor;
        PSTH(nday).barOn.(conditions.barOn{nBar}).real.psth = psth;

        % compute lfads psth
        [psth, raster_tensor] = preparePSTH_Manoj(alf{nday}, 'rates', trialIndicesPerCond(nday).barOn.(conditions.barOn{nBar}), [alf{nday}.barOnset], timePoints, binsize_rescaled);
        PSTH(nday).barOn.(conditions.barOn{nBar}).lfads.raster_tensor = raster_tensor;
        PSTH(nday).barOn.(conditions.barOn{nBar}).lfads.psth = psth;
        
        for nCue = 1 : numel(conditions.cueOn)
            cue_cond_field = [conditions.cueOn{nCue},'_', conditions.barOn{nBar}];
            trialIndicesPerCond(nday).cueOn.(cue_cond_field) = (UE{nday}.barType == nBar) & (UE{nday}.cueType == nCue);

            % compute real psth
            [psth, raster_tensor] = preparePSTH_Manoj(alf{nday}, 'spikes', trialIndicesPerCond(nday).cueOn.(cue_cond_field), [alf{nday}.cueOnset], timePoints, binsize_rescaled);
            PSTH(nday).cueOn.(cue_cond_field).real.raster_tensor = raster_tensor;
            PSTH(nday).cueOn.(cue_cond_field).real.psth = psth;

            % compute lfads psth
            [psth, raster_tensor] = preparePSTH_Manoj(alf{nday}, 'rates', trialIndicesPerCond(nday).cueOn.(cue_cond_field), [alf{nday}.cueOnset], timePoints, binsize_rescaled);
            PSTH(nday).cueOn.(cue_cond_field).lfads.raster_tensor = raster_tensor;
            PSTH(nday).cueOn.(cue_cond_field).lfads.psth = psth;            
        end
    end
end

%% calculate max value for each plot series
% need to do so


%% first let's plot barOn
for nday = 2:numel(datasets)
%for nday = 1
    saveRootDay = fullfile(saveRoot, datasets(nday).shortName);
    barSaveDir = fullfile(saveRootDay, 'barOn');
    if ~isdir(barSaveDir)
        mkdir(barSaveDir);
    end
    cd(barSaveDir)
    %saveDirTmp = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/pulvinar/180208/postAnalysis/withExternalInput_20190226/PSTH/bar_tmp';
    %cd(saveDirTmp)
    nNeurons = size(PSTH(nday).barOn.Vert.lfads.psth, 1);
  
    for n = 1:nNeurons
        f1 = figure;
        y_max = max([PSTH(nday).barOn.Vert.lfads.psth(n, :), PSTH(nday).barOn.Vert.real.psth(n, :), PSTH(nday).barOn.Hori.lfads.psth(n, :), PSTH(nday).barOn.Hori.real.psth(n, :)]);
        s1 = subplot(3,2,1);
        % make the first subplot for plotting avg Reak rates for different cue locations
        plot(PSTH(nday).barOn.Vert.real.psth(n,:), 'r', 'LineWidth', 2, 'DisplayName','Vertical Bars');
        hold on
        plot(PSTH(nday).barOn.Hori.real.psth(n,:), 'b', 'LineWidth', 2, 'DisplayName','Horizontal Bars');
        hold on 
        set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
        set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
        title(s1, 'Real', 'FontSize', 16);
        ylabel('Firing Rate');
        axes(s1)
        ylim([0 y_max])
        xlim([0 81]) %need to change there
        set(gca, 'fontsize', 12);    

        s2 = subplot(3,2,2);
        % make the first subplot for plotting avg Reak rates for different cue locations
        plot(PSTH(nday).barOn.Vert.lfads.psth(n,:), 'r', 'LineWidth', 2, 'DisplayName','Vertical Bars');
        hold on
        plot(PSTH(nday).barOn.Hori.lfads.psth(n,:), 'b', 'LineWidth', 2, 'DisplayName','Horizontal Bars');
        hold on 
        set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
        set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
        title(s2, 'LFADS', 'FontSize', 16);
        ylabel('Firing Rate');
        axes(s2)
        ylim([0 y_max])
        xlim([0 81]) % need to change here to not hard coded
        set(gca, 'fontsize', 12);    
        legend('show')
        set(legend, 'Location', 'southeast');

        s3 = subplot(3,2,3);
        a = squeeze(PSTH(nday).barOn.Vert.real.raster_tensor(n,:,:));
        imagesc(a)
        %        caxis([0 150])
        set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
        set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
        set(gca, 'fontsize', 12);    
        title(s3, 'Vert');
        ylabel('Trials');

        s4 = subplot(3,2,4);
        a = squeeze(PSTH(nday).barOn.Vert.lfads.raster_tensor(n,:,:));
        imagesc(a)
        %        caxis([0 150])        
        set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
        set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
        set(gca, 'fontsize', 12);    
        title(s4, 'Vert');
        ylabel('Trials');

        s5 = subplot(3,2,5);
        a = squeeze(PSTH(nday).barOn.Hori.real.raster_tensor(n,:,:));
        imagesc(a)
        % caxis([0 150])        
        set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
        set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
        set(gca, 'fontsize', 12);    
        title(s5, 'Hori');
        ylabel('Trials');

        s6 = subplot(3,2,6);
        a = squeeze(PSTH(nday).barOn.Hori.lfads.raster_tensor(n,:,:));
        imagesc(a)
        %caxis([0 150])        
        set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
        set(gca,'XTickLabels',{'-400ms', 'barOnset' ,'+400ms'});
        set(gca, 'fontsize', 12);    
        title(s6, 'Hori');
        ylabel('Trials');
        
        suptitle(['Multi-unit ' int2str(n+4)]);
        set(f1, 'Position', [428 200 1215 811]);
        print(f1,['Multi-unit ' int2str(n+4)], '-dpng');
        %printpdf(f1,int2str(nIndices(n)) )
        close;
    end

    % let's plot cueOn

    allIndicators(1).cond = {'Exo_TL_Vert', 'Exo_BL_Vert', 'Exo_TR_Vert', 'Exo_BR_Vert'};
    allIndicators(2).cond = {'Exo_TL_Hori', 'Exo_BL_Hori', 'Exo_TR_Hori', 'Exo_BR_Hori'};
    allIndicators(3).cond = {'Endo_TL_Vert', 'Endo_BL_Vert', 'Endo_TR_Vert', 'Endo_BR_Vert'};
    allIndicators(4).cond = {'Endo_TL_Hori', 'Endo_BL_Hori', 'Endo_TR_Hori', 'Endo_BR_Hori'};
    for nIndicator = 1:4
        plotIndicator = allIndicators(nIndicator).cond;
        %    saveRootDay = fullfile(saveRoot, datasets(nday).shortName);
        cueSaveDir = fullfile(saveRootDay, num2str(nIndicator));
        if ~isdir(cueSaveDir)
            mkdir(cueSaveDir);
        end
        cd(cueSaveDir)

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
                plot(PSTH(nday).cueOn.(plotIndicator{j}).real.psth(n,:), 'Color', cmap(j,:), 'DisplayName', condStr);
                hold on
            end
            set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
            set(gca,'XTickLabels',{'-400ms', 'cueOnset' ,'+400ms'});
            title(s1, 'Real', 'FontSize', 16);
            ylabel('Firing Rate');
            axes(s1)
            ylim([0 y_max])
            xlim([0 81]) %need to change there
            set(gca, 'fontsize', 12);    

            s2 = subplot(5,2,2);
            for j = 1:numel(plotIndicator)
                condStr = plotIndicator{j};
                whereUnderline = strfind(condStr, '_');
                condStr(whereUnderline) = '-';
                plot(PSTH(nday).cueOn.(plotIndicator{j}).lfads.psth(n,:), 'Color', cmap(j,:), 'DisplayName', condStr);
                hold on
            end
            set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
            set(gca,'XTickLabels',{'-400ms', 'cueOnset' ,'+400ms'});
            title(s2, 'LFADS', 'FontSize', 16);
            ylabel('Firing Rate');
            axes(s2)
            ylim([0 y_max])
            xlim([0 81]) % need to change here to not hard coded
            set(gca, 'fontsize', 12);    
            legend('show')
            set(legend, 'Position', [0.8254 0.8557 0.1048 0.1014]);

            for j = 1:numel(plotIndicator)

                % first, let's set up the title string
                condStr = plotIndicator{j};
                whereUnderline = strfind(condStr, '_');
                condStr(whereUnderline) = '-';
                
                subplot(5,2,2*j+1);
                a = squeeze(PSTH(nday).cueOn.(plotIndicator{j}).real.raster_tensor(n,:,:));
                imagesc(a)
                set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
                set(gca,'XTickLabels',{'-400ms', 'cueOnset' ,'+400ms'});
                title(condStr);
                ylabel('Trials');
                set(gca, 'fontsize', 10);

                subplot(5,2,2*j+2);
                a = squeeze(PSTH(nday).cueOn.(plotIndicator{j}).lfads.raster_tensor(n,:,:));
                imagesc(a)
                set(gca,'XTick',[1 0.5*length(timePoints) length(timePoints)]);
                set(gca,'XTickLabels',{'-400ms', 'cueOnset' ,'+400ms'});
                title(condStr);
                ylabel('Trials');
                set(gca, 'fontsize', 10);    
            end

            suptitle(['Multi-unit ' int2str(n+4)]);
            set(f2, 'Position', [428 4 1193 962]);
            print(f2,['Multi-unit ' int2str(n+4)], '-dpng');
            %printpdf(f2,int2str(nIndices(n)) )
            close;
        end
    end
end

%%