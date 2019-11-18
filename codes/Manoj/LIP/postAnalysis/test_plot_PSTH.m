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
            %xlim([0 81]) %need to change there
            xlim([0 41]) %need to change there
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
            set(gca,'XTickLabels',{'-400ms', 'cueOnset' ,'+400ms'});
            title(s2, 'LFADS', 'FontSize', 16);
            ylabel('Firing Rate');
            axes(s2)
            ylim([0 y_max*(6/5)])
            %xlim([0 81]) % need to change here to not hard coded
            xlim([0 41]) % need to change here to not hard coded
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
                set(gca,'XTickLabels',{'-400ms', 'cueOnset' ,'+400ms'});
                title(condStr);
                ylabel('Trials');
                set(gca, 'fontsize', 10);

                subplot(5,2,2*j+2);
                a = squeeze(PSTH(nday).cueOn.(plotIndicator{j}).lfads.raster_tensor(n,:,:));
                imagesc(a)
                set(gca,'XTick',[1 (3/8)*length(timePoints) length(timePoints)]);
                set(gca,'XTickLabels',{'-400ms', 'cueOnset' ,'+400ms'});
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