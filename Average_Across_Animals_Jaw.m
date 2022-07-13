%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
%%
[filename, pathname] = uigetfile('*.mat', 'Pick all of the individual data set', 'MultiSelect', 'on');
if isequal(filename, 0)
    return;
end
N = numel(filename);
for i = 1:N
    temp = load([pathname filename{i}]);
    eval(['data' num2str(i) ' = temp.data;']);
end

%%
n = N;
nspots = size(data1.trajectory_all, 3);
scalefactor = 50;
data.rows = data1.rows;
data.columns = data1.columns;
%1D maps
% data.startpoints_all = [];
% data.endpoints_all = [];
% for i = 1:n
%     eval(['data.startpoints_all = [data.startpoints_all data' num2str(i) '.startpoints_all];']);
%     eval(['data.endpoints_all = [data.endpoints_all data' num2str(i) '.endpoints_all];']);
% end
% for i = 1:3
%     [d_mean, pvalue, end_mean, end_std] = takecareNaN(data.startpoints_all(:, :, i), data.endpoints_all(:, :, i));
%     d_absmax = max(abs(d_mean));
%     im(:, :, 1) = mat2gray(reshape(d_mean, data.rows, data.columns), [-d_absmax d_absmax])*0.67;
%     im(:, :, 3) = 1-mat2gray(reshape(pvalue, data.rows, data.columns), [0 1]);
%     figure;
%     subplot(1, 2, 1);
%     imshow(imresize(hsv2rgb(im), scalefactor, 'nearest'));
%     title(title_text{i*2-1});
%     cbaxish = subplot(1, 2, 2);
%     imshow(hsv2rgb(colorbarim))
%     axis on;
%     box off;
%     set(cbaxish, 'XTick', [0.5 1*scalefactor+0.5], 'XTickLabel', {'1' '0'})
%     set(cbaxish, 'YTick', [0.5 0.5+data.rows*scalefactor/2 data.rows*scalefactor+0.5], 'YTickLabel', {num2str(d_absmax, '%.1f') '0' num2str(-d_absmax, '%.1f')})
%     set(cbaxish, 'YAxisLocation', 'right')
%     set(cbaxish, 'Position', [0.31 cbaxish.Position(2:4)]);
%     
%     im(:, :, 1) = mat2gray(reshape(end_mean, data.rows, data.columns))*0.67;
%     im(:, :, 3) = 1-mat2gray(reshape(end_std, data.rows, data.columns));
%     figure;
%     subplot(1, 2, 1);
%     imshow(imresize(hsv2rgb(im), scalefactor, 'nearest'));
%     title(title_text{i*2});
%     cbaxish = subplot(1, 2, 2);
%     imshow(hsv2rgb(colorbarim))
%     axis on;
%     box off;
%     set(cbaxish, 'XTick', [0.5 1*scalefactor+0.5], 'XTickLabel', {num2str(max(end_std), '%.1f') num2str(min(end_std), '%.1f')})
%     set(cbaxish, 'YTick', [0.5 0.5+data.rows*scalefactor/2 data.rows*scalefactor+0.5], 'YTickLabel', {num2str(max(end_mean), '%.1f') num2str(max(end_mean)/2+min(end_mean)/2, '%.1f') num2str(min(end_mean), '%.1f')})
%     set(cbaxish, 'YAxisLocation', 'right')
%     set(cbaxish, 'Position', [0.31 cbaxish.Position(2:4)]);
% end
colorbarim = ones(data.rows*scalefactor, 1*scalefactor, 3);
colorbarim(:, :, 1) = repmat(linspace(0.67, 0, data.rows*scalefactor).', 1, 1*scalefactor);
colorbarim(:, :, 3) = repmat(linspace(0, 1, 1*scalefactor), data.rows*scalefactor, 1);
im = zeros(data.rows, data.columns, 3)/0;
im(:, :, 2) = 1;
title_text = {'dx map' 'dz map'};
for i = 1:2
    d_mean_all = zeros(nspots, n)/0;
    pvalue_all = zeros(nspots, n)/0;
    for j = 1:n
        eval(['[d_mean_all(:, j), pvalue_all(:, j), ~, ~] = takecareNaN(data' num2str(j) '.startpoints_all(:, :, i^2), data' num2str(j) '.endpoints_all(:, :, i^2));']);
%         d_mean_all(:, j) = d_mean_all(:, j)/max(abs(d_mean_all(:, j)));
    end
    d_mean = mean(d_mean_all, 2);
    d_absmax = max(abs(d_mean));
    im(:, :, 1) = mat2gray(reshape(d_mean, data.rows, data.columns), [-d_absmax d_absmax])*0.67;
    pvalue = mean(pvalue_all, 2);
    im(:, :, 3) = 1-reshape(pvalue, data.rows, data.columns);
    figure;
    subplot(1, 2, 1);
    imshow(imresize(hsv2rgb(im), scalefactor, 'nearest'));
    title(title_text{i});
    cbaxish = subplot(1, 2, 2);
    imshow(hsv2rgb(colorbarim))
    axis on;
    box off;
    set(cbaxish, 'XTick', [0.5 1*scalefactor+0.5], 'XTickLabel', {'1' '0'})
    set(cbaxish, 'YTick', [0.5 0.5+data.rows*scalefactor/2 data.rows*scalefactor+0.5], 'YTickLabel', {num2str(d_absmax, '%.1f') '0' num2str(-d_absmax, '%.1f')})
    set(cbaxish, 'YAxisLocation', 'right')
    set(cbaxish, 'Position', [0.31 cbaxish.Position(2:4)]);
end

%%
%2D maps
d_mean = zeros(nspots, 4)/0;
d_mean_all = zeros(nspots, n, 4)/0;
interval = 1;
[X,Y] = meshgrid(interval/2:interval:(data.columns-1)*interval+interval/2, (data.rows-1)*interval+interval/2:-interval:interval/2);
for i = 1:4
    for j = 1:n
        eval(['[d_mean_all(:, j, i), ~, ~, ~] = takecareNaN(data' num2str(j) '.startpoints_all(:, :, i), data' num2str(j) '.endpoints_all(:, :, i));']);
    end
end
for i = 1:4
    d_mean(:, i) = mean(d_mean_all(:, :, i), 2);
end
x = d_mean(:, 1);
y = d_mean(:, 2);
z1 = d_mean(:, 3);
z2 = d_mean(:, 4);
x = x/max(abs([x; z2]));
z2 = z2/max(abs([x; z2]));
y = y/max(abs([y; z1]));
z1 = z1/max(abs([y; z1]));
figure;
% subplot(1, 2, 1);
hold on;
for i = 1:length(x)
    quiver(X(i), Y(i), x(i), z2(i), 'Color', 'black', 'LineWidth', 1, 'MaxHeadSize', 1);
end
axis image;
box on;
set(gca, 'XTick', [], 'YTick', [], 'XLim', [0 interval*data.columns], 'YLim', [0 interval*data.rows]);
xlabel('X');
ylabel('Z');
% subplot(1, 2, 2);
% hold on;
% for i = 1:length(y)
%     quiver(X(i), Y(i), y(i), z1(i), 'Color', 'black', 'LineWidth', 1, 'MaxHeadSize', 1);
% end
% axis image;
% box on;
% set(gca, 'XTick', [], 'YTick', [], 'XLim', [0 interval*data.columns], 'YLim', [0 interval*data.rows]);
% xlabel('Y');
% ylabel('Z');

%%
%Scatter map
% scatter_all = zeros(nspots, n);
% endpoints3_all = zeros(nspots, n, 3);
% startpoints3_all = zeros(nspots, n, 3);
% for j = 1:n
%     eval(['tempe = data' num2str(j) '.endpoints3_all;']);
%     eval(['temps = data' num2str(j) '.startpoints3_all;']);
%     for i = 1:size(tempe, 1)
%         x_row = temps(i, :, 1);
%         x_row(isnan(x_row)) = [];
%         y_row = temps(i, :, 2);
%         y_row(isnan(y_row)) = [];
%         z_row = temps(i, :, 3);
%         z_row(isnan(z_row)) = [];
%         startpoints3_all(i, j, 1) = mean(x_row);
%         startpoints3_all(i, j, 2) = mean(y_row);
%         startpoints3_all(i, j, 3) = mean(z_row);
%         
%         x_row = tempe(i, :, 1);
%         x_row(isnan(x_row)) = [];
%         y_row = tempe(i, :, 2);
%         y_row(isnan(y_row)) = [];
%         z_row = tempe(i, :, 3);
%         z_row(isnan(z_row)) = [];
%         scatter_all(i, j) = mean(sqrt((x_row-mean(x_row)).^2+(y_row-mean(y_row)).^2+(z_row-mean(z_row)).^2));
%         endpoints3_all(i, j, 1) = mean(x_row);
%         endpoints3_all(i, j, 2) = mean(y_row);
%         endpoints3_all(i, j, 3) = mean(z_row);
%     end
% end
% figure;
% imshow(imresize(reshape(mean(scatter_all, 2), data.rows, data.columns), scalefactor, 'nearest'), []);
% title('Scatter map');
% colormap parula;
% colorbar;

%%
%3D plot
% temp = repmat(linspace(0, 1, data.columns), data.rows, 1);
% facecolor_all(:, 1) = temp(:);
% temp = repmat(linspace(0, 1, data.rows)', 1, data.columns);
% facecolor_all(:, 2) = temp(:);
% facecolor_all(:, 3) = 0;
% figure;
% subplot(1, 3, 1);
% imshow(imresize(reshape(facecolor_all, data.rows, data.columns, 3), scalefactor, 'nearest'), []);
% p3h1 = subplot(1, 3, 2);
% for i = 1:nspots
%     plot3(mean(endpoints3_all(i, :, 1)-startpoints3_all(i, :, 1), 2), mean(endpoints3_all(i, :, 2)-startpoints3_all(i, :, 2), 2), mean(endpoints3_all(i, :, 3)-startpoints3_all(i, :, 3), 2),...
%         'o', 'MarkerFaceColor',facecolor_all(i, :), 'MarkerEdgeColor', facecolor_all(i, :));
%     hold on;
% end
% slocationh = plot3(0, 0, 0, 'o', 'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor', [0 0 1], 'MarkerSize', 12);
% axis equal;
% xlabel('dX position (mm)');
% ylabel('dY position (mm)');
% zlabel('dZ position (mm)');
% title('Distribution of end points in space');
% legend(slocationh, 'Rest Location');
% legend boxoff;
% p3h2 = subplot(1, 3, 3);
% for i = 1:nspots
%     plot3(mean(endpoints3_all(i, :, 1), 2), mean(endpoints3_all(i, :, 2), 2), mean(endpoints3_all(i, :, 3), 2), 'o', 'MarkerFaceColor',facecolor_all(i, :), 'MarkerEdgeColor', facecolor_all(i, :));
%     hold on;
% end
% slocationh = plot3(mean2(startpoints3_all(:, :, 1)), mean2(startpoints3_all(:, :, 2)), mean2(startpoints3_all(:, :, 3)), 'o', 'MarkerFaceColor',[0 0 1], 'MarkerEdgeColor', [0 0 1], 'MarkerSize', 12);
% axis equal;
% xlabel('X position (mm)');
% ylabel('Y position (mm)');
% zlabel('Z position (mm)');
% title('Distribution of end points in space');
% legend(slocationh, 'Rest Location');
% legend boxoff;
% % colormap(gray)
% % colorbar('YTick', [0 0.5 1], 'YTickLabel', {num2str(min(endpoints3(:, 4)), '%.1f'), num2str((min(endpoints3(:, 4))+max(endpoints3(:, 4))/2), '%.1f'), num2str(max(endpoints3(:, 4)), '%.1f')});
% set(p3h1, 'ButtonDownFcn', @(hObject, eventdata)view_adjust(hObject, eventdata, p3h2));
% set(p3h2, 'ButtonDownFcn', @(hObject, eventdata)view_adjust(hObject, eventdata, p3h1));

%%
%distance maps
distance_all = zeros(nspots, n)/0;
distanceshort_all = zeros(nspots, n)/0;
straightness_index = zeros(nspots, n)/0;
for i = 1:n
    for j = 1:nspots
        eval(['distance_row = data' num2str(i) '.distance_all(j, :);']);
        distance_row(isnan(distance_row)) = [];
        distance_all(j, i) = mean(distance_row);
        eval(['distanceshort_row = data' num2str(i) '.distanceshort_all(j, :);']);
        distanceshort_row(isnan(distanceshort_row)) = [];
        distanceshort_all(j, i) = mean(distanceshort_row);
    end
    straightness_index(:, i) = distanceshort_all(:, i)./distance_all(:, i);
end
distance = mean(distance_all, 2);
distanceshort = mean(distanceshort_all, 2);
figure;
subplot(1, 4, 1);
imshow(imresize(reshape(distanceshort, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Linear distance map');
colormap parula;
colorbar;
subplot(1, 4, 2);
imshow(imresize(reshape(distance, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Travel distance map');
colormap parula;
colorbar;
subplot(1, 4, 3);
imshow(imresize(reshape(mean(straightness_index, 2), data.rows, data.columns), scalefactor, 'nearest'), []);
title('Straightness index map');
colormap parula;
colorbar;
subplot(1, 4, 4);
imshow(imresize(reshape(distanceshort./distance, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Straightness index map');
colormap parula;
colorbar;

%%
%end2ref map
% end2ref_all = zeros(nspots, n)/0;
% for j = 1:n
%     for i = 1:nspots
%         eval(['end2ref_row = data' num2str(j) '.end2ref_all(i, :);']);
%         end2ref_row(isnan(end2ref_row)) = [];
%         end2ref_all(i, j) = mean(end2ref_row);
%     end
% end
% figure;
% end2ref = mean(end2ref_all, 2);
% imshow(imresize(reshape(end2ref, data.rows, data.columns), scalefactor, 'nearest'), []);
% title('End2Ref map');
% colormap parula;
% colorbar;

%%
%fft maps
pw_all = zeros(nspots, n)/0;
frequency_all = zeros(nspots, n)/0;
for j = 1:n
    for i = 1:nspots
        eval(['pw_row = data' num2str(j) '.power_all(i, :);']);
        pw_row(isnan(pw_row)) = [];
        pw_all(i, j) = mean(pw_row);
        eval(['frequency_row = data' num2str(j) '.frequency_all(i, :);']);
        frequency_row(isnan(frequency_row)) = [];
        frequency_all(i, j) = mean(frequency_row);
    end
end
figure;
frequency = mean(frequency_all, 2);
subplot(1, 2, 1);
imshow(imresize(reshape(frequency, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Frequency map');
colormap parula;
colorbar;
pw = mean(pw_all, 2);
subplot(1, 2, 2);
imshow(imresize(reshape(pw, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Rhythmicity map');
colormap parula;
colorbar;

%%
%rotation angle maps
% pangle_all = zeros(nspots, n)/0;
% nangle_all = zeros(nspots, n)/0;
% for j = 1:n
%     for i = 1:nspots
%         eval(['pangle_row = data' num2str(j) '.pangle_all(i, :);']);
%         pangle_row(isnan(pangle_row)) = [];
%         pangle_all(i, j) = mean(pangle_row);
%         eval(['nangle_row = data' num2str(j) '.nangle_all(i, :);']);
%         nangle_row(isnan(nangle_row)) = [];
%         nangle_all(i, j) = mean(nangle_row);
%     end
% end
% pangle = mean(pangle_all, 2);
% nangle = mean(nangle_all, 2);
% figure;
% subplot(1, 2, 1);
% imshow(imresize(reshape(pangle, data.rows, data.columns), scalefactor, 'nearest'), []);
% title('Forward rotation angle map');
% colormap(gca, 'parula');
% colorbar;
% subplot(1, 2, 2);
% imshow(imresize(reshape(nangle, data.rows, data.columns), scalefactor, 'nearest'), []);
% title('Backward rotation angle map');
% colormap(gca, 'parula');
% colorbar;

%%
%peak speed map
peakspeed_all = zeros(nspots, n)/0;
for j = 1:n
    eval(['peakspeed_all(:, j) = data' num2str(j) '.peakspeed;']);
end
figure;
peakspeed = mean(peakspeed_all, 2);
imshow(imresize(reshape(peakspeed, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Peak Speed map');
colormap parula;
colorbar;

%%
%delay map
delay_all = zeros(nspots, n)/0;
for j = 1:n
    eval(['delay_all(:, j) = data' num2str(j) '.delay;']);
end
delay = zeros(1, nspots)/0;
for i = 1:nspots
    temp = delay_all(i, :);
    temp(isnan(temp)) = [];
    if ~isempty(temp) && numel(temp) > n/2
        delay(i) = mean(temp);
    end
end
figure;
alpha = imresize(reshape(delay, data.rows, data.columns), scalefactor, 'nearest');
alpha(~isnan(alpha)) = 1;
alpha(isnan(alpha)) = 0;
imshow(zeros([size(alpha) 3]), []);
hold on;
h = imshow(imresize(reshape(delay, data.rows, data.columns), scalefactor, 'nearest'), []);
title('Delay map');
colormap(gca, 'parula');
colorbar;
set(h, 'AlphaData', alpha);

%%
%discreteness index map
im = get(get(gca, 'Children'), 'CData');

%%
n = 3;
rsquare = zeros(size(im1));
for i = 1:n
    eval(['rsquare = rsquare+im' num2str(i) ';']);
end
rsquare = rsquare/n;
figure;
imshow(rsquare, []);
title('Discreteness index map');
colormap(gca, 'parula');
colorbar;