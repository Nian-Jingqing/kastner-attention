%%
dataset(1).date = '02182019';
dataset(2).date = '03062019';
dataset(3).date = '03112019';
dataset(4).date = '03142019';
dataset(5).date = '03272019';
dataset(6).date = '04062019';
dataset(7).date = '04252019';
dataset(8).date = '05022019';

dataset(9).date = '02082019';
dataset(10).date = '02132019';
dataset(11).date = '02142019';
dataset(12).date = '02152019';
dataset(13).date = '02162019';
dataset(14).date = '02262019';
dataset(15).date = '02282019';
dataset(16).date = '03012019';
dataset(17).date = '03022019';
dataset(18).date = '03032019';

dataset(19).date = '02272019';
dataset(20).date = '03042019';
dataset(21).date = '03072019';
dataset(22).date = '03092019';
dataset(23).date = '03102019';
dataset(24).date = '03122019';
dataset(25).date = '03132019';
dataset(26).date = '03152019';

dataset(27).date = '03162019';
dataset(28).date = '03182019';
dataset(29).date = '03292019';
dataset(30).date = '03312019';
dataset(31).date = '04012019';
dataset(32).date = '04032019';
dataset(33).date = '04052019';
dataset(34).date = '04242019';
dataset(35).date = '04262019';
dataset(36).date = '04292019';

%% load UEs
UE_baseDir = '/snel/share/share/derived/kastner/data_processed/ManojData/singleArea/LIP/UEs/';
for nday = 1:numel(dataset)
    UE_file = ['UE_' dataset(nday).date '.mat'];
    load_UE_path = fullfile(UE_baseDir, UE_file);
    clear tmp
    tmp = load(load_UE_path);
    UE{ nday } = tmp.UE;
end

%%
trial_num = zeros(1, numel(dataset));
successful_num = zeros(1, numel(dataset));
hitRate = zeros(1, numel(dataset));
hitRate_validTarget = zeros(1, numel(dataset));
hitRate_sameObjTarget = zeros(1, numel(dataset));
hitRate_diffObjTarget = zeros(1, numel(dataset));
for nday = 1:numel(dataset)
    trial_num(nday) = numel(UE{nday}.fixOn);
    successful_num(nday) = trial_num(nday) - sum(UE{nday}.isErrorTrial);
    hitRate(nday) = successful_num(nday)/trial_num(nday);
    hitRate_validTarget(nday) = 1 - sum(UE{nday}.isErrorTrial(UE{nday}.isValidTarget)) / sum(UE{nday}.isValidTarget);
    hitRate_sameObjTarget(nday) = 1 - sum(UE{nday}.isErrorTrial(UE{nday}.isSameObjTarget)) / sum(UE{nday}.isSameObjTarget);
    hitRate_diffObjTarget(nday) = 1 - sum(UE{nday}.isErrorTrial(UE{nday}.isDiffObjTarget)) / sum(UE{nday}.isDiffObjTarget);
end

%% sort by trial num?
sortByTrialNum = 1;
if sortByTrialNum
    [sorted_trial_num, I] = sort(trial_num, 'descend');
    trial_num = sorted_trial_num;
    successful_num = successful_num(I);
    hitRate = hitRate(I);
    hitRate_validTarget = hitRate_validTarget(I);
    hitRate_sameObjTarget = hitRate_sameObjTarget(I);
    hitRate_diffObjTarget = hitRate_diffObjTarget(I);
end
    

%%
x = 1:numel(dataset);
sz = 50;
figure
subplot(5,1,1)
scatter(x, trial_num, sz, 'MarkerEdgeColor','b','MarkerFaceColor',[0.5,0.5,0.5])

hold on

scatter(x, successful_num, sz, 'MarkerEdgeColor','g','MarkerFaceColor',[0.5,0.5,0.5])

legend({'Num of total trials', 'Num of hit trials'})
ylim([0 1100])
xlabel('Session')
ylabel('Num of Trials')

subplot(5,1,2)
scatter(x, hitRate, sz, 'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5])
legend({'Hit rate'})
ylim([0 1.1])
ylabel('Hit rate')
xlabel('Session')

subplot(5,1,3)
scatter(x, hitRate_validTarget, sz, 'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5])
title('Cued-Loc target hit rate')
ylim([0 1.1])
ylabel('Hit rate')
xlabel('Session')

subplot(5,1,4)
scatter(x, hitRate_sameObjTarget, sz, 'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5])
title('Same-Obj target hit rate')
ylim([0 1.1])
ylabel('Hit rate')
xlabel('Session')

subplot(5,1,5)
scatter(x, hitRate_diffObjTarget, sz, 'MarkerEdgeColor','k','MarkerFaceColor',[0.5,0.5,0.5])
title('Diff-Obj target hit rate')
ylim([0 1.1])
ylabel('Hit rate')
xlabel('Session')


