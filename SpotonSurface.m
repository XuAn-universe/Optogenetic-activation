%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
xs = 0;
xe = 800;
ys = 0;
ye = 1600;
nheight = 16;
nwidth = 8;
index = 1;
figure;
hold on;
for i = 1:nheight+1
    plot([xs xe], [ys+(i-1)*(ye-ys)/nheight ys+(i-1)*(ye-ys)/nheight], 'Color', [0 0 0], 'LineWidth', 1);
end
for i = 1:nwidth+1
    plot([xs+(i-1)*(xe-xs)/nwidth xs+(i-1)*(xe-xs)/nwidth], [ys ye], 'Color', [0 0 0], 'LineWidth', 1);
end
if index ~= 0
    row = mod(index, nheight);
    if row == 0
        row = nheight;
    end
    column = ceil(index/nheight);
    hp = patch([xs+(column-1)*(xe-xs)/nwidth xs+(column-1)*(xe-xs)/nwidth xs+column*(xe-xs)/nwidth xs+column*(xe-xs)/nwidth],...
        [ye-(row-1)*(ye-ys)/nheight ye-row*(ye-ys)/nheight ye-row*(ye-ys)/nheight ye-(row-1)*(ye-ys)/nheight], [0 0 1]);
    set(hp, 'EdgeColor', 'none', 'FaceAlpha', 0.5);
end
axis image;