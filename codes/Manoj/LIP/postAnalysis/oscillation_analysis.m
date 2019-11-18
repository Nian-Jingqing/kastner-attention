%% loop through each day, each channel, and calculate a trial-avg Pxx for each channel
Fs = 100;
preTime = 1000;
whichfield = 'targetStart';
window = round([-1*preTime 0] / binsize_rescaled);
timePoints = window(1):window(2);

%%
day_Pxx = {};
day_stderr = {};
for nday = 1:numel(alf)
%for nday = 2
    longDelayTrialInd = find([alf{nday}.targetStart] - [alf{nday}.cueOnset] > round(preTime/binsize_rescaled));
    avg_Pxx = [];
    stderr = [];
    for nNeuron = 1:size(alf{nday}(1).rates, 1)
        Pxx = [];
        for itr = 1:numel(longDelayTrialInd)
            ntr = longDelayTrialInd(itr);
            selected_data = normalize(alf{nday}(ntr).rates(nNeuron, alf{nday}(ntr).(whichfield) + timePoints), 'centered');
            [Pxx(:, itr), w] = pwelch(selected_data, 50);
        end
        avg_Pxx(nNeuron, :) = mean(10*log10(Pxx), 2);
        tmp_std = std(10*log10(Pxx), 0, 2);
        stderr(nNeuron, :) = tmp_std / sqrt(size(Pxx, 2));
    end
    day_Pxx{nday} = avg_Pxx;
    day_stderr{nday} = stderr;
end

%%
saveRoot = '/snel/share/share/derived/kastner/LFADS_runs/Manoj/LIP/postAnalysis/tc_newSP_PBT_191107/powerAnalysis/';
if ~isdir(saveRoot)
    mkdir(saveRoot);
end
cd(saveRoot)

%%
for nday = 1:numel(alf)
%for nday = 1
    f1 = figure
    for nNeuron = 1:size(day_Pxx{nday}, 1)
        if mod(size(day_Pxx{nday}, 1), 2)
            num_sub = floor(size(day_Pxx{nday}, 1) / 2) + 1;
        else
            num_sub = size(day_Pxx{nday}, 1) / 2;
        end        
        subplot(2, num_sub, nNeuron)
        tmp_signal = day_Pxx{nday}(nNeuron,:);
        up = tmp_signal + day_stderr{nday}(nNeuron, :);
        down = tmp_signal - day_stderr{nday}(nNeuron, :);
        plot(w*Fs/(2*pi), tmp_signal,'b', 'lineWidth', 2);
        hold on
        plot(w*Fs/(2*pi), up,'b', 'lineWidth', 1);
        hold on
        plot(w*Fs/(2*pi), down,'b', 'lineWidth', 1);
        hold on
        xlim([0 50])
    end
    set(gcf, 'Position', [28 269 1877 663])
    title(['Day ' dataset(nday).date])
    print(f1,['Day ' dataset(nday).date], '-dpng');
end


%%
day_data = zeros(numel(alf), numel(timePoints));
for nday = 1:numel(alf)
    trial_data = zeros(numel(alf{nday}), numel(timePoints));
    for itr = 1:numel(alf{nday})
        selected_data = normalize(alf{nday}(itr).rates(:, alf{nday}(itr).(whichfield) + timePoints), 'centered');
        trial_data(itr, :) = mean(selected_data);
    end
    day_data(nday, :) = mean(trial_data);
end
%%
overall = mean(day_data);

imagesc(day_data)