%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
%%
rows = 16;
columns = 8;
mag = 50;
c = [1 2];
c = setdiff([1:3], c);
sites = [66:69 82:86 98:103 114:120];
im = ones(rows, columns, 3);
for i = 1:numel(sites)
    row = mod(sites(i), rows);
    row(row == 0) = rows;
    column = ceil(sites(i)/rows);
    im(row, column, c) = 0;
end
figure;
imshow(imresize(im, mag, 'nearest'), []);
hold on;
x = [0.5; rows*mag+0.5];
x = repmat(x, 1, rows+1);
y = 0.5:mag:rows*mag+0.5;
y = [y; y];
line(x, y, 'Color', [0 0 0], 'LineWidth', 2)
line(y, x, 'Color', [0 0 0], 'LineWidth', 2)

%%
rows = 16;
columns = 8;
mag = 50;
sites = [82:86 98:103 114:120];
im = zeros(rows, columns, 3);
% for i = 1:numel(sites)
%     row = mod(sites(i), rows);
%     row(row == 0) = rows;
%     column = ceil(sites(i)/rows);
%     im(row, column, :) = 1;
% end
figure;
imshow(imresize(im, mag, 'nearest'), []);
hold on;
x = [0.5; rows*mag+0.5];
x = repmat(x, 1, rows+1);
y = 0.5:mag:rows*mag+0.5;
y = [y; y];
line(x, y, 'Color', [1 1 1], 'LineWidth', 2)
line(y, x, 'Color', [1 1 1], 'LineWidth', 2)