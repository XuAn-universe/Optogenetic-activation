%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Oct 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
[filename, pathname] = uigetfile('*.mat', 'Pick all of the individual data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    return;
end
N = numel(filename);
dist1_start = nan(1, N);
dist1_stop = nan(1, N);
dist2_start = nan(1, N);
dist2_stop = nan(1, N);
for i = 1:N
    temp = load([pathname filename{i}]);
    result = temp.result;
    dist1_trials_start = [];
    dist1_trials_stop = [];
    dist2_trials_start = [];
    dist2_trials_stop = [];
    for j = 1:size(result.nose_start1, 1)
        mouseclick = [result.nose_start1(j, 3) result.leftpaw_start1(j, 3) result.nose_stop1(j, 3) result.leftpaw_stop1(j, 3)];
        if ~any(mouseclick == 3)
            dist1_trials_start = [dist1_trials_start vecnorm(result.nose_start1(j, 1:2)-result.leftpaw_start1(j, 1:2))];
            dist1_trials_stop = [dist1_trials_stop vecnorm(result.nose_stop1(j, 1:2)-result.leftpaw_stop1(j, 1:2))];
        end
        mouseclick = [result.nose_start2(j, 3) result.leftpaw_start2(j, 3) result.nose_stop2(j, 3) result.leftpaw_stop2(j, 3)];
        if ~any(mouseclick == 3)
            dist2_trials_start = [dist2_trials_start vecnorm(result.nose_start2(j, 1:2)-result.leftpaw_start2(j, 1:2))];
            dist2_trials_stop = [dist2_trials_stop vecnorm(result.nose_stop2(j, 1:2)-result.leftpaw_stop2(j, 1:2))];
        end
    end
    dist1_start(i) = mean(dist1_trials_start);
    dist1_stop(i) = mean(dist1_trials_stop);
    dist2_start(i) = mean(dist2_trials_start);
    dist2_stop(i) = mean(dist2_trials_stop);
end
figure;
paired_plot(dist1_start, dist1_stop, 'hand to nose distance (pixel)', filename);
set(gca, 'XTick', [1 2], 'XTickLabel', {'before', 'after'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);
figure;
paired_plot(dist2_start, dist2_stop, 'hand to nose distance (pixel)', filename);
set(gca, 'XTick', [1 2], 'XTickLabel', {'before', 'after'}, 'XLim', [0.4 2.6], 'TickLength', [0 0], 'FontSize', 12);