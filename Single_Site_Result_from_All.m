%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
f10n = 6;
f50n = 9;
p10n = 6;
p50n = 6;
spot = 102;

%% Fezf2 10 Hz
x_all = zeros(size(dataf10_1.trajectory_all, 2), f10n)/0;
y_all = zeros(size(dataf10_1.trajectory_all, 2), f10n)/0;
z1_all = zeros(size(dataf10_1.trajectory_all, 2), f10n)/0;
z2_all = zeros(size(dataf10_1.trajectory_all, 2), f10n)/0;
for i = 1:f10n
    eval(['x = dataf10_' num2str(i) '.trajectory_all(:, :, spot, 1);']);
    x(isnan(x(:, 1)), :) = [];
    x = mean(x);
    x_all(:, i) = x;
    eval(['y = dataf10_' num2str(i) '.trajectory_all(:, :, spot, 2);']);
    y(isnan(y(:, 1)), :) = [];
    y = mean(y);
    y_all(:, i) = y;
    eval(['z1 = dataf10_' num2str(i) '.trajectory_all(:, :, spot, 3);']);
    z1(isnan(z1(:, 1)), :) = [];
    z1 = mean(z1);
    z1_all(:, i) = z1;
    eval(['z2 = dataf10_' num2str(i) '.trajectory_all(:, :, spot, 4);']);
    z2(isnan(z2(:, 1)), :) = [];
    z2 = mean(z2);
    z2_all(:, i) = z2;
end

%% Fezf2 50 Hz
x_all = zeros(size(dataf50_1.trajectory_all, 2), f50n)/0;
y_all = zeros(size(dataf50_1.trajectory_all, 2), f50n)/0;
z1_all = zeros(size(dataf50_1.trajectory_all, 2), f50n)/0;
z2_all = zeros(size(dataf50_1.trajectory_all, 2), f50n)/0;
for i = 1:f50n
    eval(['x = dataf50_' num2str(i) '.trajectory_all(:, :, spot, 1);']);
    x(isnan(x(:, 1)), :) = [];
    x = mean(x);
    x_all(:, i) = x;
    eval(['y = dataf50_' num2str(i) '.trajectory_all(:, :, spot, 2);']);
    y(isnan(y(:, 1)), :) = [];
    y = mean(y);
    y_all(:, i) = y;
    eval(['z1 = dataf50_' num2str(i) '.trajectory_all(:, :, spot, 3);']);
    z1(isnan(z1(:, 1)), :) = [];
    z1 = mean(z1);
    z1_all(:, i) = z1;
    eval(['z2 = dataf50_' num2str(i) '.trajectory_all(:, :, spot, 4);']);
    z2(isnan(z2(:, 1)), :) = [];
    z2 = mean(z2);
    z2_all(:, i) = z2;
end

%% PlexinD1 10 Hz
x_all = zeros(size(datap10_1.trajectory_all, 2), p10n)/0;
y_all = zeros(size(datap10_1.trajectory_all, 2), p10n)/0;
z1_all = zeros(size(datap10_1.trajectory_all, 2), p10n)/0;
z2_all = zeros(size(datap10_1.trajectory_all, 2), p10n)/0;
for i = 1:p10n
    eval(['x = datap10_' num2str(i) '.trajectory_all(:, :, spot, 1);']);
    x(isnan(x(:, 1)), :) = [];
    x = mean(x);
    x_all(:, i) = x;
    eval(['y = datap10_' num2str(i) '.trajectory_all(:, :, spot, 2);']);
    y(isnan(y(:, 1)), :) = [];
    y = mean(y);
    y_all(:, i) = y;
    eval(['z1 = datap10_' num2str(i) '.trajectory_all(:, :, spot, 3);']);
    z1(isnan(z1(:, 1)), :) = [];
    z1 = mean(z1);
    z1_all(:, i) = z1;
    eval(['z2 = datap10_' num2str(i) '.trajectory_all(:, :, spot, 4);']);
    z2(isnan(z2(:, 1)), :) = [];
    z2 = mean(z2);
    z2_all(:, i) = z2;
end

%% PlexinD1 50 Hz
x_all = zeros(size(datap50_1.trajectory_all, 2), p50n)/0;
y_all = zeros(size(datap50_1.trajectory_all, 2), p50n)/0;
z1_all = zeros(size(datap50_1.trajectory_all, 2), p50n)/0;
z2_all = zeros(size(datap50_1.trajectory_all, 2), p50n)/0;
for i = 1:p50n
    eval(['x = datap50_' num2str(i) '.trajectory_all(:, :, spot, 1);']);
    x(isnan(x(:, 1)), :) = [];
    x = mean(x);
    x_all(:, i) = x;
    eval(['y = datap50_' num2str(i) '.trajectory_all(:, :, spot, 2);']);
    y(isnan(y(:, 1)), :) = [];
    y = mean(y);
    y_all(:, i) = y;
    eval(['z1 = datap50_' num2str(i) '.trajectory_all(:, :, spot, 3);']);
    z1(isnan(z1(:, 1)), :) = [];
    z1 = mean(z1);
    z1_all(:, i) = z1;
    eval(['z2 = datap50_' num2str(i) '.trajectory_all(:, :, spot, 4);']);
    z2(isnan(z2(:, 1)), :) = [];
    z2 = mean(z2);
    z2_all(:, i) = z2;
end

%% 2D plot
linecolors = jet(256);
markersize = 8;
figure;
subplot(2, 2, 1);
p2y = zeros(1, size(y_all, 2));
for j = 1:size(y_all, 2)
    p2y(j) = plot(y_all(:, j)-y_all(1, j), z1_all(:, j)-z1_all(1, j), '-', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'LineWidth', 2);
    hold on;
    plot(y_all(end, j)-y_all(1, j), z1_all(end, j)-z1_all(1, j), 's', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerSize', markersize);
end
if size(y_all, 2) > 1
    plot(mean(y_all, 2)-mean(y_all(1, :)), mean(z1_all, 2)-mean(z1_all(1, :)), '-k', 'LineWidth', 2);
    plot(mean(y_all(end, :)-y_all(1, :)), mean(z1_all(end, :)-z1_all(1, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
end
box off;
xlabel('dY position (mm)');
ylabel('dZ position (mm)');
title(['N = ' num2str(size(y_all, 2)) ' | Side Camera']);
subplot(2, 2, 2);
sp2y = zeros(1, size(y_all, 2));
for j = 1:size(y_all, 2)
    sp2y(j) = plot(y_all(:, j), z1_all(:, j), '-', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'LineWidth', 2);
    hold on;
    plot(y_all(1, j), z1_all(1, j), 'o', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerSize', markersize);
    plot(y_all(end, j), z1_all(end, j), 's', 'Color', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(y_all, 2)*256), :), 'MarkerSize', markersize);
end
if size(y_all, 2) > 1
    plot(mean(y_all, 2), mean(z1_all, 2), '-k', 'LineWidth', 2);
    plot(mean(y_all(1, :)), mean(z1_all(1, :)), 'ok', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    plot(mean(y_all(end, :)), mean(z1_all(end, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
end
box off;
xlabel('Y position (mm)');
ylabel('Z position (mm)');
title(['N = ' num2str(size(y_all, 2)) ' | Side Camera']);
subplot(2, 2, 3);
p2x = zeros(1, size(x_all, 2));
for j = 1:size(x_all, 2)
    p2x(j) = plot(x_all(:, j)-x_all(1, j), z2_all(:, j)-z2_all(1, j), '-', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'LineWidth', 2);
    hold on;
    plot(x_all(end, j)-x_all(1, j), z2_all(end, j)-z2_all(1, j), 's', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerSize', markersize);
end
if size(x_all, 2) > 1
    plot(mean(x_all, 2)-mean(x_all(1, :)), mean(z2_all, 2)-mean(z2_all(1, :)), '-k', 'LineWidth', 2);
    plot(mean(x_all(end, :)-x_all(1, :)), mean(z2_all(end, :)-z2_all(1, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
end
box off;
xlabel('dX position (mm)');
ylabel('dZ position (mm)');
title(['N = ' num2str(size(x_all, 2)) ' | Front Camera']);

subplot(2, 2, 4);
sp2x = zeros(1, size(x_all, 2));
for j = 1:size(x_all, 2)
    sp2x(j) = plot(x_all(:, j), z2_all(:, j), '-', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'LineWidth', 2);
    hold on;
    plot(x_all(1, j), z2_all(1, j), 'o', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerSize', markersize);
    plot(x_all(end, j), z2_all(end, j), 's', 'Color', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerFaceColor', linecolors(round(j/size(x_all, 2)*256), :), 'MarkerSize', markersize);
end
if size(x_all, 2) > 1
    plot(mean(x_all, 2), mean(z2_all, 2), '-k', 'LineWidth', 2);
    plot(mean(x_all(1, :)), mean(z2_all(1, :)), 'ok', 'MarkerFaceColor','k', 'MarkerSize', markersize);
    plot(mean(x_all(end, :)), mean(z2_all(end, :)), 'sk', 'MarkerFaceColor','k', 'MarkerSize', markersize);
end
box off;
xlabel('X position (mm)');
ylabel('Z position (mm)');
title(['N = ' num2str(size(x_all, 2)) ' | Front Camera']);

%% End2Ref
group = [];
f50_end2ref = zeros(1, f50n)/0;
for i = 1:f50n
    eval(['end2ref = dataf50_' num2str(i) '.end2ref_all(spot, :);'])
    end2ref(isnan(end2ref)) = [];
    f50_end2ref(i) = mean(end2ref);
    group{i} = 'F 50';
end
p50_end2ref = zeros(1, p50n)/0;
for i = 1:p50n
    eval(['end2ref = datap50_' num2str(i) '.end2ref_all(spot, :);'])
    end2ref(isnan(end2ref)) = [];
    p50_end2ref(i) = mean(end2ref);
    group{numel(group)+1} = 'P 50';
end
f10_end2ref = zeros(1, f10n)/0;
for i = 1:f10n
    eval(['end2ref = dataf10_' num2str(i) '.end2ref_all(spot, :);'])
    end2ref(isnan(end2ref)) = [];
    f10_end2ref(i) = mean(end2ref);
    group{numel(group)+1} = 'F 10';
end
p10_end2ref = zeros(1, p10n)/0;
for i = 1:p10n
    eval(['end2ref = datap10_' num2str(i) '.end2ref_all(spot, :);'])
    end2ref(isnan(end2ref)) = [];
    p10_end2ref(i) = mean(end2ref);
    group{numel(group)+1} = 'P 10';
end

figure;
boxplot([f50_end2ref p50_end2ref f10_end2ref p10_end2ref], group);
ylabel('End2Ref (mm)');

figure;
plot([ones(1, f50n) ones(1, p50n)*2 ones(1, f10n)*3 ones(1, p10n)*4]+(rand(1, f50n+p50n+f10n+p10n)-0.5)/2,...
    [f50_end2ref p50_end2ref f10_end2ref p10_end2ref], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(f50_end2ref) mean(f50_end2ref)], '-r', 'LineWidth', 2);
plot([2-0.25 2+0.25], [mean(p50_end2ref) mean(p50_end2ref)], '-r', 'LineWidth', 2);
plot([3-0.25 3+0.25], [mean(f10_end2ref) mean(f10_end2ref)], '-r', 'LineWidth', 2);
plot([4-0.25 4+0.25], [mean(p10_end2ref) mean(p10_end2ref)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'F 50' 'P 50' 'F 10' 'P 10'}); 
ylabel('End2Ref (mm)');

%% Total Distance
group = [];
f50_distance = zeros(1, f50n)/0;
for i = 1:f50n
    eval(['distance = dataf50_' num2str(i) '.distance_all(spot, :);'])
    distance(isnan(distance)) = [];
    f50_distance(i) = mean(distance);
    group{i} = 'F 50';
end
p50_distance = zeros(1, p50n)/0;
for i = 1:p50n
    eval(['distance = datap50_' num2str(i) '.distance_all(spot, :);'])
    distance(isnan(distance)) = [];
    p50_distance(i) = mean(distance);
    group{numel(group)+1} = 'P 50';
end
f10_distance = zeros(1, f10n)/0;
for i = 1:f10n
    eval(['distance = dataf10_' num2str(i) '.distance_all(spot, :);'])
    distance(isnan(distance)) = [];
    f10_distance(i) = mean(distance);
    group{numel(group)+1} = 'F 10';
end
p10_distance = zeros(1, p10n)/0;
for i = 1:p10n
    eval(['distance = datap10_' num2str(i) '.distance_all(spot, :);'])
    distance(isnan(distance)) = [];
    p10_distance(i) = mean(distance);
    group{numel(group)+1} = 'P 10';
end

figure;
boxplot([f50_distance p50_distance f10_distance p10_distance], group);
ylabel('Total Distance (mm)');

figure;
plot([ones(1, f50n) ones(1, p50n)*2 ones(1, f10n)*3 ones(1, p10n)*4]+(rand(1, f50n+p50n+f10n+p10n)-0.5)/2,...
    [f50_distance p50_distance f10_distance p10_distance], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(f50_distance) mean(f50_distance)], '-r', 'LineWidth', 2);
plot([2-0.25 2+0.25], [mean(p50_distance) mean(p50_distance)], '-r', 'LineWidth', 2);
plot([3-0.25 3+0.25], [mean(f10_distance) mean(f10_distance)], '-r', 'LineWidth', 2);
plot([4-0.25 4+0.25], [mean(p10_distance) mean(p10_distance)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'F 50' 'P 50' 'F 10' 'P 10'}); 
ylabel('Total Distance (mm)');

%% Linear Distance
group = [];
f50_distanceshort = zeros(1, f50n)/0;
for i = 1:f50n
    eval(['distanceshort = dataf50_' num2str(i) '.distanceshort_all(spot, :);'])
    distanceshort(isnan(distanceshort)) = [];
    f50_distanceshort(i) = mean(distanceshort);
    group{i} = 'F 50';
end
p50_distanceshort = zeros(1, p50n)/0;
for i = 1:p50n
    eval(['distanceshort = datap50_' num2str(i) '.distanceshort_all(spot, :);'])
    distanceshort(isnan(distanceshort)) = [];
    p50_distanceshort(i) = mean(distanceshort);
    group{numel(group)+1} = 'P 50';
end
f10_distanceshort = zeros(1, f10n)/0;
for i = 1:f10n
    eval(['distanceshort = dataf10_' num2str(i) '.distanceshort_all(spot, :);'])
    distanceshort(isnan(distanceshort)) = [];
    f10_distanceshort(i) = mean(distanceshort);
    group{numel(group)+1} = 'F 10';
end
p10_distanceshort = zeros(1, p10n)/0;
for i = 1:p10n
    eval(['distanceshort = datap10_' num2str(i) '.distanceshort_all(spot, :);'])
    distanceshort(isnan(distanceshort)) = [];
    p10_distanceshort(i) = mean(distanceshort);
    group{numel(group)+1} = 'P 10';
end

figure;
boxplot([f50_distanceshort p50_distanceshort f10_distanceshort p10_distanceshort], group);
ylabel('Linear Distance (mm)');

figure;
plot([ones(1, f50n) ones(1, p50n)*2 ones(1, f10n)*3 ones(1, p10n)*4]+(rand(1, f50n+p50n+f10n+p10n)-0.5)/2,...
    [f50_distanceshort p50_distanceshort f10_distanceshort p10_distanceshort], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(f50_distanceshort) mean(f50_distanceshort)], '-r', 'LineWidth', 2);
plot([2-0.25 2+0.25], [mean(p50_distanceshort) mean(p50_distanceshort)], '-r', 'LineWidth', 2);
plot([3-0.25 3+0.25], [mean(f10_distanceshort) mean(f10_distanceshort)], '-r', 'LineWidth', 2);
plot([4-0.25 4+0.25], [mean(p10_distanceshort) mean(p10_distanceshort)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'F 50' 'P 50' 'F 10' 'P 10'}); 
ylabel('Linear Distance (mm)');

%% Peak Speed
group = [];
f50_peakspeed = zeros(1, f50n)/0;
for i = 1:f50n
    eval(['peakspeed = dataf50_' num2str(i) '.peakspeed(spot);'])
    f50_peakspeed(i) = peakspeed;
    group{i} = 'F 50';
end
p50_peakspeed = zeros(1, p50n)/0;
for i = 1:p50n
    eval(['peakspeed = datap50_' num2str(i) '.peakspeed(spot);'])
    p50_peakspeed(i) = peakspeed;
    group{numel(group)+1} = 'P 50';
end
f10_peakspeed = zeros(1, f10n)/0;
for i = 1:f10n
    eval(['peakspeed = dataf10_' num2str(i) '.peakspeed(spot);'])
    f10_peakspeed(i) = peakspeed;
    group{numel(group)+1} = 'F 10';
end
p10_peakspeed = zeros(1, p10n)/0;
for i = 1:p10n
    eval(['peakspeed = datap10_' num2str(i) '.peakspeed(spot);'])
    p10_peakspeed(i) = peakspeed;
    group{numel(group)+1} = 'P 10';
end

figure;
boxplot([f50_peakspeed p50_peakspeed f10_peakspeed p10_peakspeed], group);
ylabel('Peak Speed (mm/s)');

figure;
plot([ones(1, f50n) ones(1, p50n)*2 ones(1, f10n)*3 ones(1, p10n)*4]+(rand(1, f50n+p50n+f10n+p10n)-0.5)/2,...
    [f50_peakspeed p50_peakspeed f10_peakspeed p10_peakspeed], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(f50_peakspeed) mean(f50_peakspeed)], '-r', 'LineWidth', 2);
plot([2-0.25 2+0.25], [mean(p50_peakspeed) mean(p50_peakspeed)], '-r', 'LineWidth', 2);
plot([3-0.25 3+0.25], [mean(f10_peakspeed) mean(f10_peakspeed)], '-r', 'LineWidth', 2);
plot([4-0.25 4+0.25], [mean(p10_peakspeed) mean(p10_peakspeed)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'F 50' 'P 50' 'F 10' 'P 10'}); 
ylabel('Peak Speed (mm/s)');

%% Delay
group = [];
f50_delay = zeros(1, f50n)/0;
for i = 1:f50n
    eval(['delay = dataf50_' num2str(i) '.delay(spot);'])
    f50_delay(i) = delay;
    group{i} = 'F 50';
end
p50_delay = zeros(1, p50n)/0;
for i = 1:p50n
    eval(['delay = datap50_' num2str(i) '.delay(spot);'])
    p50_delay(i) = delay;
    group{numel(group)+1} = 'P 50';
end
f10_delay = zeros(1, f10n)/0;
for i = 1:f10n
    eval(['delay = dataf10_' num2str(i) '.delay(spot);'])
    f10_delay(i) = delay;
    group{numel(group)+1} = 'F 10';
end
p10_delay = zeros(1, p10n)/0;
for i = 1:p10n
    eval(['delay = datap10_' num2str(i) '.delay(spot);'])
    p10_delay(i) = delay;
    group{numel(group)+1} = 'P 10';
end

figure;
boxplot([f50_delay p50_delay f10_delay p10_delay], group);
ylabel('Delay (ms)');

figure;
plot([ones(1, f50n) ones(1, p50n)*2 ones(1, f10n)*3 ones(1, p10n)*4]+(rand(1, f50n+p50n+f10n+p10n)-0.5)/2,...
    [f50_delay p50_delay f10_delay p10_delay], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(f50_delay(~isnan(f50_delay))) mean(f50_delay(~isnan(f50_delay)))], '-r', 'LineWidth', 2);
plot([2-0.25 2+0.25], [mean(p50_delay(~isnan(p50_delay))) mean(p50_delay(~isnan(p50_delay)))], '-r', 'LineWidth', 2);
plot([3-0.25 3+0.25], [mean(f10_delay(~isnan(f10_delay))) mean(f10_delay(~isnan(f10_delay)))], '-r', 'LineWidth', 2);
plot([4-0.25 4+0.25], [mean(p10_delay(~isnan(p10_delay))) mean(p10_delay(~isnan(p10_delay)))], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'F 50' 'P 50' 'F 10' 'P 10'}); 
ylabel('Delay (ms)');

%% Frequency
group = [];
f50_frequency = zeros(1, f50n)/0;
for i = 1:f50n
    eval(['frequency = dataf50_' num2str(i) '.frequency_all(spot, :);'])
    frequency(isnan(frequency)) = [];
    f50_frequency(i) = mean(frequency);
    group{i} = 'F 50';
end
p50_frequency = zeros(1, p50n)/0;
for i = 1:p50n
    eval(['frequency = datap50_' num2str(i) '.frequency_all(spot, :);'])
    frequency(isnan(frequency)) = [];
    p50_frequency(i) = mean(frequency);
    group{numel(group)+1} = 'P 50';
end
f10_frequency = zeros(1, f10n)/0;
for i = 1:f10n
    eval(['frequency = dataf10_' num2str(i) '.frequency_all(spot, :);'])
    frequency(isnan(frequency)) = [];
    f10_frequency(i) = mean(frequency);
    group{numel(group)+1} = 'F 10';
end
p10_frequency = zeros(1, p10n)/0;
for i = 1:p10n
    eval(['frequency = datap10_' num2str(i) '.frequency_all(spot, :);'])
    frequency(isnan(frequency)) = [];
    p10_frequency(i) = mean(frequency);
    group{numel(group)+1} = 'P 10';
end

figure;
boxplot([f50_frequency p50_frequency f10_frequency p10_frequency], group);
ylabel('Frequency (Hz)');

figure;
plot([ones(1, f50n) ones(1, p50n)*2 ones(1, f10n)*3 ones(1, p10n)*4]+(rand(1, f50n+p50n+f10n+p10n)-0.5)/2,...
    [f50_frequency p50_frequency f10_frequency p10_frequency], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(f50_frequency) mean(f50_frequency)], '-r', 'LineWidth', 2);
plot([2-0.25 2+0.25], [mean(p50_frequency) mean(p50_frequency)], '-r', 'LineWidth', 2);
plot([3-0.25 3+0.25], [mean(f10_frequency) mean(f10_frequency)], '-r', 'LineWidth', 2);
plot([4-0.25 4+0.25], [mean(p10_frequency) mean(p10_frequency)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'F 50' 'P 50' 'F 10' 'P 10'}); 
ylabel('Frequency (Hz)');

%% Power
group = [];
f50_power = zeros(1, f50n)/0;
for i = 1:f50n
    eval(['power = dataf50_' num2str(i) '.power_all(spot, :);'])
    power(isnan(power)) = [];
    f50_power(i) = mean(power);
    group{i} = 'F 50';
end
p50_power = zeros(1, p50n)/0;
for i = 1:p50n
    eval(['power = datap50_' num2str(i) '.power_all(spot, :);'])
    power(isnan(power)) = [];
    p50_power(i) = mean(power);
    group{numel(group)+1} = 'P 50';
end
f10_power = zeros(1, f10n)/0;
for i = 1:f10n
    eval(['power = dataf10_' num2str(i) '.power_all(spot, :);'])
    power(isnan(power)) = [];
    f10_power(i) = mean(power);
    group{numel(group)+1} = 'F 10';
end
p10_power = zeros(1, p10n)/0;
for i = 1:p10n
    eval(['power = datap10_' num2str(i) '.power_all(spot, :);'])
    power(isnan(power)) = [];
    p10_power(i) = mean(power);
    group{numel(group)+1} = 'P 10';
end

figure;
boxplot([f50_power p50_power f10_power p10_power], group);
ylabel('Power (mm^{2})');

figure;
plot([ones(1, f50n) ones(1, p50n)*2 ones(1, f10n)*3 ones(1, p10n)*4]+(rand(1, f50n+p50n+f10n+p10n)-0.5)/2,...
    [f50_power p50_power f10_power p10_power], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(f50_power) mean(f50_power)], '-r', 'LineWidth', 2);
plot([2-0.25 2+0.25], [mean(p50_power) mean(p50_power)], '-r', 'LineWidth', 2);
plot([3-0.25 3+0.25], [mean(f10_power) mean(f10_power)], '-r', 'LineWidth', 2);
plot([4-0.25 4+0.25], [mean(p10_power) mean(p10_power)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', [1 2 3 4], 'XTickLabel', {'F 50' 'P 50' 'F 10' 'P 10'}); 
ylabel('Power (mm^{2})');
