%*---------------------------------------------------------------------*
% Huang Lab
% Duke University
% Author : Xu An, Sep 2021
% xu.an@duke.edu
% Version: 1.0
%*---------------------------------------------------------------------*
%%
markersize = 8;
thr = 0.6;
frames = find(~isnan(leftpaw(:, 1)));
nframe = numel(frames);
if thr ~= 1
    vidFrame = rgb2hsv(RGB);
    intensity = vidFrame(:, :, 3);
    intensity = mat2gray(intensity, [0 thr]);
    vidFrame(:, :, 3) = intensity;
    vidFrame = hsv2rgb(vidFrame);
end
figure;
imshow(rgb2gray(vidFrame));
hold on;
plot(leftpaw(frames, 1), leftpaw(frames, 2), '-', 'Color', [1 1 0], 'LineWidth', 1);
plot(leftpaw(frames(1), 1), leftpaw(frames(1), 2), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', [1 1 0], 'Color', [1 1 0]);
plot(leftpaw(frames(end), 1), leftpaw(frames(end), 2), 's', 'MarkerSize', markersize, 'MarkerFaceColor', [1 1 0], 'Color', [1 1 0]);
plot(rightpaw(frames, 1), rightpaw(frames, 2), '-', 'Color', [0 1 0], 'LineWidth', 1);
plot(rightpaw(frames(1), 1), rightpaw(frames(1), 2), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', [0 1 0], 'Color', [0 1 0]);
plot(rightpaw(frames(end), 1), rightpaw(frames(end), 2), 's', 'MarkerSize', markersize, 'MarkerFaceColor', [0 1 0], 'Color', [0 1 0]);
plot(nose(frames, 1), nose(frames, 2), '-', 'Color', [1 0 0], 'LineWidth', 1);
plot(nose(frames(1), 1), nose(frames(1), 2), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', [1 0 0], 'Color', [1 0 0]);
plot(nose(frames(end), 1), nose(frames(end), 2), 's', 'MarkerSize', markersize, 'MarkerFaceColor', [1 0 0], 'Color', [1 0 0]);
plot(righthindpaw(frames, 1), righthindpaw(frames, 2), '-', 'Color', [0 0 1], 'LineWidth', 1);
plot(righthindpaw(frames(1), 1), righthindpaw(frames(1), 2), 'o', 'MarkerSize', markersize, 'MarkerFaceColor', [0 0 1], 'Color', [0 0 1]);
plot(righthindpaw(frames(end), 1), righthindpaw(frames(end), 2), 's', 'MarkerSize', markersize, 'MarkerFaceColor', [0 0 1], 'Color', [0 0 1]);