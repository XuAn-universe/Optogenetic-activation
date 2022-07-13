%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Sep 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
front_check = 0;
thr = 0.5;
if thr ~= 1
    vidFrame = rgb2hsv(RGB);
    intensity = vidFrame(:, :, 3);
    intensity = mat2gray(intensity, [0 thr]);
    vidFrame(:, :, 3) = intensity;
    vidFrame = hsv2rgb(vidFrame);
end
figure;
imshow(vidFrame);

markersize = 5;
leftcolor = [0.6 0.6 0];
lavgcolor = [1 1 0];
rightcolor = [0 0.6 0];
ravgcolor = [0 1 0];
linecolors = jet(256);

% for i = 1:size(result_left.x_pixel, 2)
%     hold on;
%     lc = linecolors(round(i/size(result_left.x_pixel, 2)*256), :);
%     plot(result_left.x_pixel(:, i), result_left.z2_pixel(:, i), '-', 'LineWidth', 1, 'Color', lc);
%     plot(result_left.x_pixel(1, i), result_left.z2_pixel(1, i), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', lc, 'Color', lc);
%     plot(result_left.x_pixel(end, i), result_left.z2_pixel(end, i), 's', 'MarkerSize', markersize, 'MarkerFaceColor', lc, 'Color', lc);
% end
% for i = 1:size(result_right.x_pixel, 2)
%     hold on;
%     rc = linecolors(round(i/size(result_right.x_pixel, 2)*256), :);
%     plot(result_right.x_pixel(:, i), result_right.z2_pixel(:, i), '-', 'LineWidth', 1, 'Color', rc);
%     plot(result_right.x_pixel(1, i), result_right.z2_pixel(1, i), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', rc, 'Color', rc);
%     plot(result_right.x_pixel(end, i), result_right.z2_pixel(end, i), 's', 'MarkerSize', markersize, 'MarkerFaceColor', rc, 'Color', rc);
% end

if front_check
    for i = 1:size(result_left.x_pixel, 2)
        hold on;
        plot(result_left.x_pixel(:, i), result_left.z2_pixel(:, i), '-', 'LineWidth', 1, 'Color', leftcolor);
        plot(result_left.x_pixel(1, i), result_left.z2_pixel(1, i), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', leftcolor, 'Color', leftcolor);
        plot(result_left.x_pixel(end, i), result_left.z2_pixel(end, i), 's', 'MarkerSize', markersize, 'MarkerFaceColor', leftcolor, 'Color', leftcolor);
    end
    try
        for i = 1:size(result_right.x_pixel, 2)
            hold on;
            plot(result_right.x_pixel(:, i), result_right.z2_pixel(:, i), '-', 'LineWidth', 1, 'Color', rightcolor);
            plot(result_right.x_pixel(1, i), result_right.z2_pixel(1, i), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', rightcolor, 'Color', rightcolor);
            plot(result_right.x_pixel(end, i), result_right.z2_pixel(end, i), 's', 'MarkerSize', markersize, 'MarkerFaceColor', rightcolor, 'Color', rightcolor);
        end
    end
    plot(mean(result_left.x_pixel, 2), mean(result_left.z2_pixel, 2), '-', 'LineWidth', 3, 'Color', lavgcolor);
    plot(mean(result_left.x_pixel(1, :), 2), mean(result_left.z2_pixel(1, :), 2), 'o', 'MarkerSize', markersize*2, 'MarkerFaceColor', lavgcolor, 'Color', lavgcolor);
    plot(mean(result_left.x_pixel(end, :), 2), mean(result_left.z2_pixel(end, :), 2), 's', 'MarkerSize', markersize*2, 'MarkerFaceColor', lavgcolor, 'Color', lavgcolor);
    pixelpermm = (mean(result_left.x_pixel(end, :), 2)-mean(result_left.x_pixel(1, :), 2))/(result_left.meanx(end)-result_left.meanx(1));
    plot([20 20+pixelpermm*5], [size(RGB, 1)-20 size(RGB, 1)-20], '-w', 'LineWidth', 3); % 5mm scale bar
    try
        plot(mean(result_right.x_pixel, 2), mean(result_right.z2_pixel, 2), '-', 'LineWidth', 3, 'Color', ravgcolor);
        plot(mean(result_right.x_pixel(1, :), 2), mean(result_right.z2_pixel(1, :), 2), 'o', 'MarkerSize', markersize*2, 'MarkerFaceColor', ravgcolor, 'Color', ravgcolor);
        plot(mean(result_right.x_pixel(end, :), 2), mean(result_right.z2_pixel(end, :), 2), 's', 'MarkerSize', markersize*2, 'MarkerFaceColor', ravgcolor, 'Color', ravgcolor);
    end
else
    for i = 1:size(result_left.y_pixel, 2)
        hold on;
        plot(result_left.y_pixel(:, i), result_left.z1_pixel(:, i), '-', 'LineWidth', 1, 'Color', leftcolor);
        plot(result_left.y_pixel(1, i), result_left.z1_pixel(1, i), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', leftcolor, 'Color', leftcolor);
        plot(result_left.y_pixel(end, i), result_left.z1_pixel(end, i), 's', 'MarkerSize', markersize, 'MarkerFaceColor', leftcolor, 'Color', leftcolor);
    end
    try
        for i = 1:size(result_right.y_pixel, 2)
            hold on;
            plot(result_right.y_pixel(:, i), result_right.z1_pixel(:, i), '-', 'LineWidth', 1, 'Color', rightcolor);
            plot(result_right.y_pixel(1, i), result_right.z1_pixel(1, i), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', rightcolor, 'Color', rightcolor);
            plot(result_right.y_pixel(end, i), result_right.z1_pixel(end, i), 's', 'MarkerSize', markersize, 'MarkerFaceColor', rightcolor, 'Color', rightcolor);
        end
    end
    plot(mean(result_left.y_pixel, 2), mean(result_left.z1_pixel, 2), '-', 'LineWidth', 3, 'Color', lavgcolor);
    plot(mean(result_left.y_pixel(1, :), 2), mean(result_left.z1_pixel(1, :), 2), 'o', 'MarkerSize', markersize*2, 'MarkerFaceColor', lavgcolor, 'Color', lavgcolor);
    plot(mean(result_left.y_pixel(end, :), 2), mean(result_left.z1_pixel(end, :), 2), 's', 'MarkerSize', markersize*2, 'MarkerFaceColor', lavgcolor, 'Color', lavgcolor);
    pixelpermm = (mean(result_left.y_pixel(end, :), 2)-mean(result_left.y_pixel(1, :), 2))/(result_left.meany(end)-result_left.meany(1));
    plot([20 20+pixelpermm*5], [size(RGB, 1)-20 size(RGB, 1)-20], '-w', 'LineWidth', 3); % 5mm scale bar
    try
        plot(mean(result_right.y_pixel, 2), mean(result_right.z1_pixel, 2), '-', 'LineWidth', 3, 'Color', ravgcolor);
        plot(mean(result_right.y_pixel(1, :), 2), mean(result_right.z1_pixel(1, :), 2), 'o', 'MarkerSize', markersize*2, 'MarkerFaceColor', ravgcolor, 'Color', ravgcolor);
        plot(mean(result_right.y_pixel(end, :), 2), mean(result_right.z1_pixel(end, :), 2), 's', 'MarkerSize', markersize*2, 'MarkerFaceColor', ravgcolor, 'Color', ravgcolor);
    end
end