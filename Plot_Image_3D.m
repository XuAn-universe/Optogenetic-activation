%%
ax = get(gca, 'Children');
t = get(ax, 'CData');
%%
[X,Y] = meshgrid(0:1:399, 799:-1:0);
figure;
hold on;
surf(X, Y, ones(800,400)*0, t, 'LineStyle', 'none');
surf(X, Y, ones(800,400)*250, f, 'LineStyle', 'none');
surf(X, Y, ones(800,400)*500, p, 'LineStyle', 'none');
colormap jet;
colorbar;
axis equal;
axis tight;
axis off;
line([0 0 0 0 0 0; 399 399 399 399 399 399], [-1 -1 -1 800 800 800; -1 -1 -1 800 800 800], [0 250 500 0 250 500; 0 250 500 0 250 500], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);
line([-1 -1 -1 400 400 400; -1 -1 -1 400 400 400], [-2 -2 -2 -2 -2 -2; 801 801 801 801 801 801], [0 250 500 0 250 500; 0 250 500 0 250 500], 'Color', [0.5 0.5 0.5], 'LineWidth', 2);