function make_pbt_run_plots( pbt_dir, results_save_dir, image_format, plot_hps, benchmarks_fitlog_dir, benchmark_type )
% pbt_dir, directory where PBT workers are saved
%
% results_save_dir, directory to save the generated plots by this function
%
% image_format, '-dpng', 'pdf', '-djpg'
%
% benchmarks_fitlog_dir, cell array containing the directories containing
% fitlog.csv for benchmark runs (to plot on top of PBT plots). Leave empty
% if not used
%
% benchmark_type, type of benchmark run: 'lfadslite' (default) or 'LFADS'

if ~exist('plot_hps')
    plot_hps = true;
end

if ~exist('benchmarks_fitlog_dir')
    benchmarks_fitlog_dir = {};
end
if ~exist('benchmark_type')
    benchmark_type = 'lfadslite';
end
if ~exist('image_format')
    image_format = '-dpng';
end


% save directory
%
param_file = fullfile(pbt_dir, 'pbt_params.json');
pbt_params = jsondecode(fileread(param_file));
run_name = pbt_params.name;

results_save_dir = fullfile(results_save_dir, run_name);
if ~exist(results_save_dir, 'file')
    mkdir(results_save_dir)
end


[runs, epoch_per_gen] = PBT_analysis.load_pbt_results( pbt_dir );

epoch_per_gen = pbt_params.steps_to_ready;


%% plot some
r = runs(min([size(runs,1),3]):end,:); %ignore first 3 generation for ylim
idx = cellfun(@(x) length(x), {r.valid});
idx = idx == epoch_per_gen;

v = [r(idx).valid]; v = v(:);
try
    sv = [r(idx).valid_samp]; sv = sv(:);
    if any(isnan(sv))
        sv = inf;
    end
catch
    sv = inf;
end
tr = [r(idx).train]; tr = tr(:);
minY = min([min(v), min(sv), min(tr)]);
if sv==inf, sv=-inf; end
maxY = max([max(v), max(sv), max(tr)]);

%ylims = [mean(tr) - 0.2*std(tr), mean(v) + 0.3*std(v)];
ylims = [minY, maxY + 0.3*max([std(v), std(sv), std(tr)])];

opacity = 0.9;

%
fh = figure(1); clf;
set(gcf, 'color', [1 1 1]);
ah(1) = subplot(3,1,1);
hold on
ah(2) = subplot(3,1,2);
hold on
ah(3) = subplot(3,1,3);
hold on

%% PBT results
axes(ah(1));
PBT_analysis.plot_pbt_results( runs, 'train', opacity );
ylabel('Train');
title('Costs')

axes(ah(2));
PBT_analysis.plot_pbt_results( runs, 'valid', opacity );
ylabel('Valid-Trial');

% plot valid samp
try
axes(ah(3));
PBT_analysis.plot_pbt_results( runs, 'valid_samp', opacity );
ylabel('Valid-Samp');
xlabel('Generation');
end

linkaxes(ah)

%set(ah, 'fontsize', 10);
ylim( ylims )
xlim([1 size(runs,1)])

% removes spaces around plots
Plot.tightfig(fh);

%%
if ~isempty(benchmarks_fitlog_dir)
   
    num_benchmarks = numel(benchmarks_fitlog_dir);
    colmap = lines(num_benchmarks);

    for i = 1:num_benchmarks 
        % overlay benchmark runs
        benchmark = sprintf('%s/fitlog.csv', benchmarks_fitlog_dir{i});
        log = PBT_analysis.read_fitlog( benchmark );
        train = cellfun(@(x) str2num(x), log(:, 9) );
        epoch = 1:length(train);
        if strcmp(benchmark_type, 'LFADS')    
            % this is for LFADS log file
            valid = cellfun(@(x) str2num(x), log(:, 10) );   
        else
            % this is for lfadslite log file
            valid = cellfun(@(x) str2num(x), log(:, 11) );
        end
        axes(ah(1))
        h = plot((epoch+epoch_per_gen)/epoch_per_gen, train, 'Color', colmap(i, :), 'LineWidth', 0.5);
        axes(ah(2))
        h=plot((epoch+epoch_per_gen)/epoch_per_gen, valid, 'Color', colmap(i, :), 'LineWidth', 0.5);
        alpha(opacity)
    end
end

% save plots
fname = fullfile(results_save_dir, 'costs');
if strcmp(image_format, 'pdf')
    export_fig([fname '.pdf'])
else
    print(image_format, '-r300', fname);
end

if plot_hps
    % Get the searched HPs
    all_hps = fields(runs(1,1).hps)';
    hp_list = {};
    for hp = all_hps
        hp = hp{1};
        if isnumeric(runs(1,1).hps.(hp)) && runs(1,1).hps.(hp) ~= runs(end,1).hps.(hp)
            hp_list = {hp_list{:} hp};
        end
    end
    %
    % plots HPs and save them
    opacity = 0.90;
    %hp_list = {'keep_prob', 'learning_rate_init', 'l2_gen_scale', 'l2_ic_enc_scale', 'l2_ci_enc_scale', 'l2_con_scale', 'kl_co_weight', 'kl_ic_weight' };
    for hp = hp_list
        fh = figure; clf;
        set(gcf, 'color', [1 1 1]);
        PBT_analysis.plot_pbt_results( runs, hp{1}, opacity );
        set(gca, 'yscale', 'log')
        axis('tight')
        ylabel(hp{1}, 'Interpreter', 'none');
        xlabel('Generation');
        %set(gca, 'fontsize', 12);
        xlim([1 size(runs,1)])
        %title(run_name, 'Interpreter', 'none', 'FontSize', 9)
        title(hp{1}, 'Interpreter', 'none')
        Plot.tightfig(fh);
        if strcmp(image_format, 'pdf')
            export_fig(fullfile(results_save_dir, [hp{1} '.pdf']))
        else
            print(image_format, '-r300', fullfile(results_save_dir, hp{1}));
        end
    end
end

end

