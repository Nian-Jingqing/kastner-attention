% ****** These are indices after the bad MUs are removed
%day 1 - 02182019
% all trials
all_trials_MU = [1, 3:10, 12, 14, 17:22, 24:25, 27, 29, 31:32];
for i = 1:length(all_trials_MU)
    d{1}{all_trials_MU(i)} = 1:666;
end
% individual trials
d{1}{2} = NaN;
d{1}{11} = 100:666;
d{1}{13} = [1:614, 628:666];
d{1}{15} = 1:500;
d{1}{16} = NaN;
d{1}{23} = 1:480;
d{1}{26} = [1:583, 585:666];
d{1}{28} = 1:380;
d{1}{30} = 1:520;

drift_channels{1} = [2, 16];

%%
%day 2 - 03062019
% all
% all trials
all_trials_MU = [1:2, 4, 7, 22:23, 25, 27:32];
for i = 1:length(all_trials_MU)
    d{2}{all_trials_MU(i)} = 1:623;
end
%%
% individual trials
d{2}{3} = 200:623;
d{2}{5} = 200:623;
d{2}{6} = 200:623;
d{2}{8} = [1:250, 360:623];
d{2}{9} = [1:170, 280:623];
d{2}{10} = 200:623;
d{2}{11} = 300:623;
d{2}{12} = 250:623;
d{2}{13} = 200:623;
d{2}{14} = [1:80 250:623];
d{2}{15} = 300:623;
d{2}{16} = 1:550;
d{2}{17} = 1:450;
d{2}{18} = [1:370, 570:623];
d{2}{19} = 200:623;
d{2}{20} = 200:623;
d{2}{21} = 180:623;
d{2}{24} = 90:600;
d{2}{26} = 130:623;

drift_channels{2} = [];

%%
%day 3 - 03112019
% all trials
all_trials_MU = [2:4, 6:8, 11, 13:14, 16, 19:21, 23:24, 26:29, 31:32];
for i = 1:length(all_trials_MU)
    d{3}{all_trials_MU(i)} = 1:886;
end

%%
% individual trials
d{3}{1} = [200:380, 440:886];
d{3}{5} = 1:480;
d{3}{9} = NaN;
d{3}{10} = 1:800;
d{3}{12} = 250:886;
d{3}{15} = 300:886;
d{3}{17} = 140:886;
d{3}{18} = 1:800;
d{3}{22} = 270:886;
d{3}{25} = 200:886;
d{3}{30} = 120:800;

drift_channels{3} = 9;

%%
%day 4 - 03142019
% all trials
all_trials_MU = [2:5, 8:11, 13, 15:24, 26:28, 31:32];
for i = 1:length(all_trials_MU)
    d{4}{all_trials_MU(i)} = 1:639;
end

%%
% individual trials
d{4}{1} = 1:560;
d{4}{6} = 1:560;
d{4}{7} = 120:639;
d{4}{12} = 100:560;
d{4}{14} = 250:639;
d{4}{25} = 130:639;
d{4}{29} = 1:560;
d{4}{30} = 1:560;

drift_channels{4} = [];

%%
%day 5 - 04062019
% all trials
all_trials_MU = [3:4, 6:12, 29:31];
for i = 1:length(all_trials_MU)
    d{5}{all_trials_MU(i)} = 1:841;
end

%%
% individual trials
d{5}{1} = NaN;
d{5}{2} = 350:841;
d{5}{5} = 250:760;
d{5}{13} = 1:730;
d{5}{14} = 350:841;
d{5}{15} = NaN;
d{5}{16} = 1:450;
d{5}{17} = 130:841;
d{5}{18} = 150:841;
d{5}{19} = 150:760;
d{5}{20} = 70:700;
d{5}{21} = 280:841;
d{5}{22} = 1:720;
d{5}{23} = 200:841;
d{5}{24} = 150:841;
d{5}{25} = 1:780;
d{5}{26} = [1:470, 560:841];
d{5}{27} = [140:320, 450:841];
d{5}{28} = [1:180, 240:780];
d{5}{32} = 220:841;

drift_channels{5} = [1, 15];


%%
%day 6 - 04252019
% all trials
all_trials_MU = [13, 16, 19, 21:24, 26, 28, 30, 32];
for i = 1:length(all_trials_MU)
    d{6}{all_trials_MU(i)} = 1:618;
end

%%
% individual trials
d{6}{1} = [1:60, 150:560];
d{6}{2} = [1:130, 220:618];
d{6}{3} = [1:360, 480:618];
d{6}{4} = [1:300, 450:570];
d{6}{5} = 1:450;
d{6}{6} = NaN;
d{6}{7} = 1:500;
d{6}{8} = 200:618;
d{6}{9} = 180:618;
d{6}{10} = 1:540;
d{6}{11} = 250:618;
d{6}{12} = [1:250, 400:550];
d{6}{14} = 100:550;
d{6}{15} = 220:618;
d{6}{17} = 180:618;
d{6}{18} = 1:550;
d{6}{20} = [1:50, 120:550];
d{6}{25} = 200:618;
d{6}{27} = 230:618;
d{6}{29} = 120:618;
d{6}{31} = [1:410, 490:618];

drift_channels{6} = 6;

%%
%day 7 - 05022019
% all trials
all_trials_MU = [1:8, 10:11, 16, 21, 23:25, 27:28, 30, 32];
for i = 1:length(all_trials_MU)
    d{7}{all_trials_MU(i)} = 1:791;
end

%%
% individual trials
d{7}{9} = [1:300, 360:650];
d{7}{12} = [1:360, 520:750];
d{7}{13} = 1:600;
d{7}{14} = 60:791;
d{7}{15} = 290:650;
d{7}{17} = 1:670;
d{7}{18} = 1:650;
d{7}{19} = 200:791;
d{7}{20} = [1:200, 400:791];
d{7}{22} = 250:791;
d{7}{26} = 180:791;
d{7}{29} = 1:600;
d{7}{31} = 200:791;

drift_channels{7} = [];