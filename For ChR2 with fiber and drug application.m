%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%% 2D plot
n = 6;
for i = 1:n
    eval(['temp = result' num2str(i) ';']);
    x_all(:, i) = temp.meanx;
    y_all(:, i) = temp.meany;
    z1_all(:, i) = temp.meanz;
    z2_all(:, i) = temp.meanz;
end
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
axis equal;
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
axis equal;
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
axis equal;
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
axis equal;
box off;
xlabel('X position (mm)');
ylabel('Z position (mm)');
title(['N = ' num2str(size(x_all, 2)) ' | Front Camera']);

%%
for i = 1:n
    eval(['resultb' num2str(i) ' = result' num2str(i) ';']);
end

%%
n = 6;
FrameRate = 100;
linecolors = jet(256);
figure;
hold on;
for i = 1:n
    eval(['temp = resultb' num2str(i) ';']);
    t = (-1:numel(temp.meanx)-2)/FrameRate*1000;
    plot(t, temp.meanx-temp.meanx(1), '-', 'Color', linecolors(round(i/n*256), :), 'LineWidth', 2);
    meandx_all(:, i) = temp.meanx-temp.meanx(1);
end
dxbefore = meandx_all;
meandx = mean(meandx_all, 2);
plot(t, meandx, '-k', 'LineWidth', 4);
title(['N = ' num2str(n) ' | before: dx']);
box off;
xlabel('Time after light on (ms)');
ylabel('dx (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

meandxsem = std(meandx_all, 0, 2)/sqrt(n);;
figure;
errorbar(t, meandx, meandxsem, 'k-', 'LineWidth', 1);
ax = gca;

figure;
hold on;
for i = 1:n
    eval(['temp = resulta' num2str(i) ';']);
    t = (-1:numel(temp.meanx)-2)/FrameRate*1000;
    plot(t, temp.meanx-temp.meanx(1), '-', 'Color', linecolors(round(i/n*256), :), 'LineWidth', 2);
    meandx_all(:, i) = temp.meanx-temp.meanx(1);
end
dxafter = meandx_all;
meandx = mean(meandx_all, 2);
plot(t, meandx, '-k', 'LineWidth', 4);
title(['N = ' num2str(n) ' | after: dx']);
box off;
xlabel('Time after light on (ms)');
ylabel('dx (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

meandxsem = std(meandx_all, 0, 2)/sqrt(n);
axes(ax);
hold(ax, 'on');
errorbar(ax, t, meandx, meandxsem, 'r-', 'LineWidth', 1);
title(['N = ' num2str(n) ' | dx']);
box off;
xlabel('Time after light on (ms)');
ylabel('dx (mm)');
xlim([t(1) t(end)+1000/FrameRate]);
legend(ax, {'before' 'after'});
legend boxoff;

figure;
dxchange = dxafter-dxbefore;
dxchangemean = mean(dxchange, 2);
dxchangesem = std(dxchange, 0, 2)/sqrt(n);
errorbar(t, dxchangemean, dxchangesem, 'k-', 'LineWidth', 1);
title(['N = ' num2str(n) ' | dx(after)-dx(before)']);
box off;
xlabel('Time after light on (ms)');
ylabel('dx(after)-dx(before) (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

%%
n = 6;
FrameRate = 100;
linecolors = jet(256);
figure;
hold on;
for i = 1:n
    eval(['temp = resultb' num2str(i) ';']);
    t = (-1:numel(temp.meany)-2)/FrameRate*1000;
    plot(t, temp.meany-temp.meany(1), '-', 'Color', linecolors(round(i/n*256), :), 'LineWidth', 2);
    meandy_all(:, i) = temp.meany-temp.meany(1);
end
dybefore = meandy_all;
meandy = mean(meandy_all, 2);
plot(t, meandy, '-k', 'LineWidth', 4);
title(['N = ' num2str(n) ' | before: dy']);
box off;
xlabel('Time after light on (ms)');
ylabel('dy (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

meandysem = std(meandy_all, 0, 2)/sqrt(n);;
figure;
errorbar(t, meandy, meandysem, 'k-', 'LineWidth', 1);
ax = gca;

figure;
hold on;
for i = 1:n
    eval(['temp = resulta' num2str(i) ';']);
    t = (-1:numel(temp.meany)-2)/FrameRate*1000;
    plot(t, temp.meany-temp.meany(1), '-', 'Color', linecolors(round(i/n*256), :), 'LineWidth', 2);
    meandy_all(:, i) = temp.meany-temp.meany(1);
end
dyafter = meandy_all;
meandy = mean(meandy_all, 2);
plot(t, meandy, '-k', 'LineWidth', 4);
title(['N = ' num2str(n) ' | after: dy']);
box off;
xlabel('Time after light on (ms)');
ylabel('dy (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

meandysem = std(meandy_all, 0, 2)/sqrt(n);
axes(ax);
hold(ax, 'on');
errorbar(ax, t, meandy, meandysem, 'r-', 'LineWidth', 1);
title(['N = ' num2str(n) ' | dy']);
box off;
xlabel('Time after light on (ms)');
ylabel('dy (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

legend(ax, {'before' 'after'});
legend boxoff;

figure;
dychange = dyafter-dybefore;
dychangemean = mean(dychange, 2);
dychangesem = std(dychange, 0, 2)/sqrt(n);
errorbar(t, dychangemean, dychangesem, 'k-', 'LineWidth', 1);
title(['N = ' num2str(n) ' | dy(after)-dy(before)']);
box off;
xlabel('Time after light on (ms)');
ylabel('dy(after)-dy(before) (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

%%
n = 6;
FrameRate = 100;
linecolors = jet(256);
figure;
hold on;
for i = 1:n
    eval(['temp = resultb' num2str(i) ';']);
    t = (-1:numel(temp.meanz)-2)/FrameRate*1000;
    plot(t, temp.meanz-temp.meanz(1), '-', 'Color', linecolors(round(i/n*256), :), 'LineWidth', 2);
    meandz_all(:, i) = temp.meanz-temp.meanz(1);
end
dzbefore = meandz_all;
meandz = mean(meandz_all, 2);
plot(t, meandz, '-k', 'LineWidth', 4);
title(['N = ' num2str(n) ' | before: dz']);
box off;
xlabel('Time after light on (ms)');
ylabel('dz (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

meandzsem = std(meandz_all, 0, 2)/sqrt(n);;
figure;
errorbar(t, meandz, meandzsem, 'k-', 'LineWidth', 1);
ax = gca;

figure;
hold on;
for i = 1:n
    eval(['temp = resulta' num2str(i) ';']);
    t = (-1:numel(temp.meanz)-2)/FrameRate*1000;
    plot(t, temp.meanz-temp.meanz(1), '-', 'Color', linecolors(round(i/n*256), :), 'LineWidth', 2);
    meandz_all(:, i) = temp.meanz-temp.meanz(1);
end
dzafter = meandz_all;
meandz = mean(meandz_all, 2);
plot(t, meandz, '-k', 'LineWidth', 4);
title(['N = ' num2str(n) ' | after: dz']);
box off;
xlabel('Time after light on (ms)');
ylabel('dz (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

meandzsem = std(meandz_all, 0, 2)/sqrt(n);
axes(ax);
hold(ax, 'on');
errorbar(ax, t, meandz, meandzsem, 'r-', 'LineWidth', 1);
title(['N = ' num2str(n) ' | dz']);
box off;
xlabel('Time after light on (ms)');
ylabel('dz (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

legend(ax, {'before' 'after'});
legend boxoff;

figure;
dzchange = dzafter-dzbefore;
dzchangemean = mean(dzchange, 2);
dzchangesem = std(dzchange, 0, 2)/sqrt(n);
errorbar(t, dzchangemean, dzchangesem, 'k-', 'LineWidth', 1);
title(['N = ' num2str(n) ' | dz(after)-dz(before)']);
box off;
xlabel('Time after light on (ms)');
ylabel('dz(after)-dz(before) (mm)');
xlim([t(1) t(end)+1000/FrameRate]);

%%
n = 6;
FrameRate = 100;
linecolors = jet(256);
figure;
hold on;
for i = 1:n
    eval(['temp = resultb' num2str(i) ';']);
    t = (-1:numel(temp.speedmean)-2)/FrameRate*1000;
    plot(t, temp.speedmean, '-', 'Color', linecolors(round(i/n*256), :), 'LineWidth', 2);
    speedmean_all(:, i) = temp.speedmean;
end
speedbefore = speedmean_all;
speedmean = mean(speedmean_all, 2);
plot(t, speedmean, '-k', 'LineWidth', 4);
title(['N = ' num2str(n) ' | before: speed']);
box off;
xlabel('Time after light on (ms)');
ylabel('Speed (mm/s)');
xlim([t(1) t(end)+1000/FrameRate]);

speedmeansem = std(speedmean_all, 0, 2)/sqrt(n);;
figure;
errorbar(t, speedmean, speedmeansem, 'k-', 'LineWidth', 1);
ax = gca;

figure;
hold on;
for i = 1:n
    eval(['temp = resulta' num2str(i) ';']);
    t = (-1:numel(temp.speedmean)-2)/FrameRate*1000;
    plot(t, temp.speedmean, '-', 'Color', linecolors(round(i/n*256), :), 'LineWidth', 2);
    speedmean_all(:, i) = temp.speedmean;
end
speedafter = speedmean_all;
speedmean = mean(speedmean_all, 2);
plot(t, speedmean, '-k', 'LineWidth', 4);
title(['N = ' num2str(n) ' | after: speed']);
box off;
xlabel('Time after light on (ms)');
ylabel('Speed (mm/s)');
xlim([t(1) t(end)+1000/FrameRate]);

speedmeansem = std(speedmean_all, 0, 2)/sqrt(n);
axes(ax);
hold(ax, 'on');
errorbar(ax, t, speedmean, speedmeansem, 'r-', 'LineWidth', 1);
title(['N = ' num2str(n) ' | Speed']);
box off;
xlabel('Time after light on (ms)');
ylabel('Speed (mm/s)');
xlim([t(1) t(end)+1000/FrameRate]);

legend(ax, {'before' 'after'});
legend boxoff;

figure;
speedchange = speedafter-speedbefore;
speedchangemean = mean(speedchange, 2);
speedchangesem = std(speedchange, 0, 2)/sqrt(n);
errorbar(t, speedchangemean, speedchangesem, 'k-', 'LineWidth', 1);
title(['N = ' num2str(n) ' | speed(after)-speed(before)']);
box off;
xlabel('Time after light on (ms)');
ylabel('speed(after)-speed(before) (mm/s)');
xlim([t(1) t(end)+1000/FrameRate]);

%%
n = 5;
for i = 1:n
    eval(['temp = resultb' num2str(i) ';']);
    try
        end2ref_before(i) = temp.end2ref;
    end
    distance_before(i) = temp.distance;
    distanceshort_before(i) = temp.distanceshort;
    straightnessindex_before(i) = temp.straightnessindex;
    peakspeed_before(i) = temp.peakspeed;
    delay_before(i) = temp.delay;
    frequency_before(i) = temp.frequency;
    power_before(i) = temp.power;
    
    eval(['temp = resulta' num2str(i) ';']);
    try
        end2ref_after(i) = temp.end2ref;
    end
    distance_after(i) = temp.distance;
    distanceshort_after(i) = temp.distanceshort;
    straightnessindex_after(i) = temp.straightnessindex;
    peakspeed_after(i) = temp.peakspeed;
    delay_after(i) = temp.delay;
    frequency_after(i) = temp.frequency;
    power_after(i) = temp.power;
end
alpha = 0.05;
try
    figure;
    hold on;
    line(repmat([1;2], 1, n), [end2ref_before; end2ref_after], 'Color', 'k');
    line(repmat([3;4], 1, n), [distance_before; distance_after], 'Color', 'k');
    line(repmat([5;6], 1, n), [distanceshort_before; distanceshort_after], 'Color', 'k');
    line([(1:2:6)-0.35; (1:2:6)+0.35], repmat([mean(end2ref_before) mean(distance_before) mean(distanceshort_before)], 2, 1), 'Color', 'k', 'LineWidth', 3);
    line([(2:2:6)-0.35; (2:2:6)+0.35], repmat([mean(end2ref_after) mean(distance_after) mean(distanceshort_after)], 2, 1), 'Color', 'r', 'LineWidth', 3);
    plot(repmat(1:2:6, n, 1), [end2ref_before' distance_before' distanceshort_before'], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot(repmat(2:2:6, n, 1), [end2ref_after' distance_after' distanceshort_after'], 'or', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    % bar([mean(end2ref_after-end2ref_before) mean(distance_after-distance_before) mean(distanceshort_after-distanceshort_before)], 'Width', 0.5)
    box off;
    set(gca, 'XTick', 1:6, 'XTickLabel', {'End2Ref before' 'End2Ref after' 'Rpath before' 'Rpath after' 'Ro before' 'Ro after'}, 'XTickLabelRotation', 45);
    xlim([0.5 6.5]);
    ylabel('mm');
    
    [H, P] = ttest(end2ref_before, end2ref_after, 'alpha', 0.05);
    disp(['end2ref before vs after: H = ' num2str(H) '; P = ' num2str(P)]);
    [H, P] = ttest(distance_before, distance_after, 'alpha', 0.05);
    disp(['Rpath before vs after: H = ' num2str(H) '; P = ' num2str(P)]);
    [H, P] = ttest(distanceshort_before, distanceshort_after, 'alpha', 0.05);
    disp(['Ro before vs after: H = ' num2str(H) '; P = ' num2str(P)]);
    
    figure;
    plot([ones(1, n) ones(1, n)*2 ones(1, n)*3]+(rand(1, 3*n)-0.5)/2,...
        [end2ref_after-end2ref_before distance_after-distance_before distanceshort_after-distanceshort_before], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
    hold on;
    plot([1-0.25 1+0.25], [mean(end2ref_after-end2ref_before) mean(end2ref_after-end2ref_before)], '-r', 'LineWidth', 2);
    plot([2-0.25 2+0.25], [mean(distance_after-distance_before) mean(distance_after-distance_before)], '-r', 'LineWidth', 2);
    plot([3-0.25 3+0.25], [mean(distanceshort_after-distanceshort_before) mean(distanceshort_after-distanceshort_before)], '-r', 'LineWidth', 2);
    % bar([mean(end2ref_after-end2ref_before) mean(distance_after-distance_before) mean(distanceshort_after-distanceshort_before)], 'Width', 0.5)
    box off;
    set(gca, 'XTick', [1 2 3], 'XTickLabel', {'\DeltaEnd2Ref' '\DeltaRpath' '\DeltaRo'});
    ylabel('mm');
catch
    close
    figure;
    hold on;
    line(repmat([1;2], 1, n), [distance_before; distance_after], 'Color', 'k');
    line(repmat([3;4], 1, n), [distanceshort_before; distanceshort_after], 'Color', 'k');
    line([(1:2:4)-0.35; (1:2:4)+0.35], repmat([mean(distance_before) mean(distanceshort_before)], 2, 1), 'Color', 'k', 'LineWidth', 3);
    line([(2:2:4)-0.35; (2:2:4)+0.35], repmat([mean(distance_after) mean(distanceshort_after)], 2, 1), 'Color', 'r', 'LineWidth', 3);
    plot(repmat(1:2:4, n, 1), [distance_before' distanceshort_before'], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot(repmat(2:2:4, n, 1), [distance_after' distanceshort_after'], 'or', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
    % bar([mean(end2ref_after-end2ref_before) mean(distance_after-distance_before) mean(distanceshort_after-distanceshort_before)], 'Width', 0.5)
    box off;
    set(gca, 'XTick', 1:4, 'XTickLabel', {'Rpath before' 'Rpath after' 'Ro before' 'Ro after'}, 'XTickLabelRotation', 45);
    xlim([0.5 4.5]);
    ylabel('mm');
    
    [H, P] = ttest(distance_before, distance_after, 'alpha', 0.05);
    disp(['Rpath before vs after: H = ' num2str(H) '; P = ' num2str(P)]);
    [H, P] = ttest(distanceshort_before, distanceshort_after, 'alpha', 0.05);
    disp(['Ro before vs after: H = ' num2str(H) '; P = ' num2str(P)]);
    
    figure;
    plot([ones(1, n) ones(1, n)*2]+(rand(1, 2*n)-0.5)/2,...
        [distance_after-distance_before distanceshort_after-distanceshort_before], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
    hold on;
    plot([1-0.25 1+0.25], [mean(distance_after-distance_before) mean(distance_after-distance_before)], '-r', 'LineWidth', 2);
    plot([2-0.25 2+0.25], [mean(distanceshort_after-distanceshort_before) mean(distanceshort_after-distanceshort_before)], '-r', 'LineWidth', 2);
    % bar([mean(end2ref_after-end2ref_before) mean(distance_after-distance_before) mean(distanceshort_after-distanceshort_before)], 'Width', 0.5)
    box off;
    set(gca, 'XTick', [1 2], 'XTickLabel', {'\DeltaRpath' '\DeltaRo'});
    ylabel('mm');
end

figure;
hold on;
line(repmat([1;2], 1, n), [straightnessindex_before; straightnessindex_after], 'Color', 'k');
line([1-0.35 1+0.35], [mean(straightnessindex_before) mean(straightnessindex_before)], 'Color', 'k', 'LineWidth', 3);
line([2-0.35 2+0.35], [mean(straightnessindex_after) mean(straightnessindex_after)], 'Color', 'r', 'LineWidth', 3);
plot(ones(1, n), straightnessindex_before, 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(ones(1, n)*2, straightnessindex_after, 'or', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% bar([mean(end2ref_after-end2ref_before) mean(distance_after-distance_before) mean(distanceshort_after-distanceshort_before)], 'Width', 0.5)
box off;
set(gca, 'XTick', 1:2, 'XTickLabel', {'Straightness Index before' 'Straightness Index after'}, 'XTickLabelRotation', 45);
xlim([0.5 2.5]);
[H, P] = ttest(straightnessindex_before, straightnessindex_after, 'alpha', 0.05);
disp(['Straightness Index before vs after: H = ' num2str(H) '; P = ' num2str(P)]);
figure;
plot(ones(1, n)+(rand(1, n)-0.5)/2,...
    [straightnessindex_after-straightnessindex_before], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(straightnessindex_after-straightnessindex_before) mean(straightnessindex_after-straightnessindex_before)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', 1, 'XTickLabel', '\DeltaStraightness Index');

figure;
hold on;
line(repmat([1;2], 1, n), [delay_before; delay_after], 'Color', 'k');
line([1-0.35 1+0.35], [mean(delay_before) mean(delay_before)], 'Color', 'k', 'LineWidth', 3);
line([2-0.35 2+0.35], [mean(delay_after) mean(delay_after)], 'Color', 'r', 'LineWidth', 3);
plot(ones(1, n), delay_before, 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(ones(1, n)*2, delay_after, 'or', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% bar([mean(end2ref_after-end2ref_before) mean(distance_after-distance_before) mean(distanceshort_after-distanceshort_before)], 'Width', 0.5)
box off;
set(gca, 'XTick', 1:2, 'XTickLabel', {'Delay before' 'Delay after'}, 'XTickLabelRotation', 45);
xlim([0.5 2.5]);
ylabel('ms');
[H, P] = ttest(delay_before, delay_after, 'alpha', 0.05);
disp(['Delay before vs after: H = ' num2str(H) '; P = ' num2str(P)]);
figure;
plot(ones(1, n)+(rand(1, n)-0.5)/2,...
    [delay_after-delay_before], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(delay_after-delay_before) mean(delay_after-delay_before)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', 1, 'XTickLabel', '\DeltaDelay'); 
ylabel('ms');

figure;
hold on;
line(repmat([1;2], 1, n), [peakspeed_before; peakspeed_after], 'Color', 'k');
line([1-0.35 1+0.35], [mean(peakspeed_before) mean(peakspeed_before)], 'Color', 'k', 'LineWidth', 3);
line([2-0.35 2+0.35], [mean(peakspeed_after) mean(peakspeed_after)], 'Color', 'r', 'LineWidth', 3);
plot(ones(1, n), peakspeed_before, 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
plot(ones(1, n)*2, peakspeed_after, 'or', 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
% bar([mean(end2ref_after-end2ref_before) mean(distance_after-distance_before) mean(distanceshort_after-distanceshort_before)], 'Width', 0.5)
box off;
set(gca, 'XTick', 1:2, 'XTickLabel', {'Peak Speed before' 'Peak Speed after'}, 'XTickLabelRotation', 45);
xlim([0.5 2.5]);
ylabel('mm/s');
[H, P] = ttest(peakspeed_before, peakspeed_after, 'alpha', 0.05);
disp(['Peak Speed before vs after: H = ' num2str(H) '; P = ' num2str(P)]);
figure;
plot(ones(1, n)+(rand(1, n)-0.5)/2,...
    [peakspeed_after-peakspeed_before], 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k')
hold on;
plot([1-0.25 1+0.25], [mean(peakspeed_after-peakspeed_before) mean(peakspeed_after-peakspeed_before)], '-r', 'LineWidth', 2);
box off;
set(gca, 'XTick', 1, 'XTickLabel', '\DeltaPeak Speed'); 
ylabel('mm/s');
