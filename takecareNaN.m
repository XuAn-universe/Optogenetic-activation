function [d_mean, pvalue, end_mean, end_std] = takecareNaN(start_all, end_all)
%*---------------------------------------------------------------------*
% Huang Lab
% Cold Spring Harbor Laboratory
% Author : Xu An, March 2017
% xan@cshl.edu
%*---------------------------------------------------------------------*
d_mean = zeros(size(start_all, 1), 1)/0;
pvalue = zeros(size(start_all, 1), 1)/0;
end_mean = zeros(size(start_all, 1), 1)/0;
end_std = zeros(size(start_all, 1), 1)/0;
for i = 1:size(start_all, 1)
    start_row = start_all(i, :);
    start_row(isnan(start_row)) = [];
    end_row = end_all(i, :);
    end_row(isnan(end_row)) = [];
    d_mean(i) = mean(end_row-start_row);
    [~, pvalue(i)] = ttest(start_row, end_row);
    end_mean(i) = mean(end_row);
    end_std(i) = std(end_row);
end