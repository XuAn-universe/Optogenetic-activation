%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
n = size(data, 1);
width = 0.5;
figure;
hold on;
% bar(1:size(data, 2), mean(data), width, 'FaceColor', 'w');
for i = 1:size(data, 2)
    plot(ones(1, n)*i+(rand(1, n)-0.5)*width, data(:, i)*10^3, 'ok', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
    plot([i-width/2 i+width/2], [mean(data(:, i)) mean(data(:, i))]*10^3, 'r', 'LineWidth', 2);
end
set(gca, 'XTick', [1:size(data, 2)], 'XTickLabel', {'Plexin AP' 'Fezf AP' 'Plexin ML' 'Fezf ML'}, 'XTickLabelRotation', 45)
ylabel('\mum')
title('Sensory Forelimb Location')